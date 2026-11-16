"""Check that the browser app and the Python package give identical predictions.

    python -m probe.verify --scope broad

Checks, for one model scope:
1. The model bundled in the app has the same parameters as the model in models/.
2. Python reproduces the frozen predictions recorded before the outcomes were read.
3. app/probe.js (run in Node, with no outcome files present) matches Python on full,
   partial and single-missing input panels, including error bands.
"""
import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from .evaluate import ROOT, read_bundle
from .model import Model

NODE = '''const fs = require('fs'); const P = require('./probe.js');
const model = JSON.parse(fs.readFileSync('model.json')), support = JSON.parse(fs.readFileSync('support.json'));
const cases = JSON.parse(fs.readFileSync('cases.json'));
fs.writeFileSync('out.json', JSON.stringify(cases.map(x => P.infer(model, x, true, support))));'''


def verify(scope, root=ROOT):
    model_path = root/'models'/scope/'model.json'
    model = Model.load(model_path)
    data = read_bundle(root/'app'/f'{scope}-data.js')
    assert data['model'] == model.metadata, 'App model differs from models/'
    assert read_bundle(root/'app/support.js')[model.version] == json.loads(model_path.with_name('support.json').read_text())

    reference = json.loads((root/'tests/fixtures'/f'{scope}-reference.json').read_text())
    inputs = np.asarray(reference['inputs'], float)
    assert [c['id'] for c in data['cases']] == reference['model_ids']
    np.testing.assert_array_equal([c['values'] for c in data['cases']], inputs)
    frozen = np.asarray(reference['predictions'][f'targeted_{len(model.probes)}'])
    np.testing.assert_allclose(model.panel.predict(inputs), frozen, rtol=0, atol=1e-12)

    n = len(model.probes)
    cases = []
    for k in [n, n // 2, 4, 0]:
        cases += [np.where(np.arange(n) < k, row, np.nan) for row in inputs]
    for row in inputs[:5]:
        cases += [np.where(np.arange(n) == i, np.nan, row) for i in range(n)]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        shutil.copy(root/'app/probe.js', tmp)
        shutil.copy(model_path, tmp)
        shutil.copy(model_path.with_name('support.json'), tmp)
        (tmp/'cases.json').write_text(json.dumps([[None if np.isnan(v) else float(v) for v in c] for c in cases]))
        (tmp/'run.js').write_text(NODE)
        subprocess.run(['node', 'run.js'], cwd=tmp, check=True)
        browser = json.loads((tmp/'out.json').read_text())

    worst = 0.0
    for x, js in zip(cases, browser):
        py = model.predict(x)
        diff = np.abs(np.array(js['scores']) - [r['predicted_effect'] for r in py])
        worst = max(worst, float(diff.max()))
        assert js['status'] == py[0]['calibration_status']
        assert (js['bands'] is not None) == (py[0]['lower'] is not None)
        if js['bands']:
            np.testing.assert_allclose([b['lower'] for b in js['bands']], [r['lower'] for r in py], rtol=0, atol=1e-12)
            np.testing.assert_allclose([b['upper'] for b in js['bands']], [r['upper'] for r in py], rtol=0, atol=1e-12)
            assert [b['dependency'] for b in js['bands']] == [r['dependency_call'] for r in py]
    assert worst <= 1e-12, worst
    return dict(scope=scope, inputs=n, outputs=len(model.targets), held_out_profiles=len(inputs),
                panels_checked=len(cases), max_difference=worst)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--scope', choices=['broad', 'compact'], default='broad')
    r = verify(parser.parse_args().scope)
    print(f"Passed: {r['scope']} model ({r['inputs']} inputs → {r['outputs']} outputs). "
          f"Python and browser agree on {r['panels_checked']} input panels (max difference {r['max_difference']:.1e}).")


if __name__ == '__main__':
    main()
