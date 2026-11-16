import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from probe import Model, read_measurements
from probe.evaluate import ROOT, evaluate, read_bundle
from probe.panel import fit, select
from probe.verify import verify

SCOPES = ['broad', 'compact']
load = lambda scope: Model.load(ROOT/'models'/scope/'model.json')


@pytest.mark.skipif(shutil.which('node') is None, reason='needs Node.js')
@pytest.mark.parametrize('scope', SCOPES)
def test_browser_and_python_agree(scope):
    assert verify(scope)['max_difference'] <= 1e-12


@pytest.mark.parametrize('scope', SCOPES)
def test_recomputed_results_match_recorded_test(scope):
    r = evaluate(scope)
    recorded = {(row['policy'], int(row['budget'])): float(row['normalized_mse'])
                for row in csv.DictReader(open(ROOT/'results'/scope/'comparison.csv'))}
    n = len(load(scope).probes)
    errors = list(r['errors'].values())
    assert errors[0] == pytest.approx(recorded['targeted', n], abs=1e-12)
    assert errors[1] == pytest.approx(recorded['d_optimal', 2 * n if scope == 'broad' else 32], abs=1e-12)
    evidence = read_bundle(ROOT/'app'/f'{scope}-data.js')['evidence']
    for rule, key in [('conservative', 'conservative'), ('wider', 'point')]:
        assert r[rule]['flags'] == evidence[key]['calls'] and r[rule]['correct'] == evidence[key]['correct']
    assert r['dependencies'] == evidence['point']['total_dependencies']


def test_headline_numbers():
    r = evaluate('broad')
    assert (r['profiles'], r['donors'], r['observed'], r['dependencies']) == (124, 120, 24177, 5881)
    assert (r['conservative'], r['wider']) == (dict(flags=788, correct=733), dict(flags=4409, correct=3081))
    assert round(r['coverage'], 3) == 0.876


def test_bands_only_for_complete_measured_in_range_panel():
    model = load('broad')
    x = np.asarray(read_bundle(ROOT/'app/broad-data.js')['cases'][1]['values'], float)
    assert model.predict(x)[0]['calibration_status'] == 'complete_panel'
    assert model.predict(x, measured=False)[0]['lower'] is None
    partial = x.copy(); partial[0] = np.nan
    assert model.predict(partial)[0]['calibration_status'] == 'missing_inputs'
    far = x.copy(); far[0] = model.upper[0] + 1
    assert model.predict(far)[0]['calibration_status'] == 'outside_training_range'
    assert np.allclose([r['predicted_effect'] for r in model.predict(np.full(len(x), np.nan))], model.panel.output_mean)


def test_read_measurements():
    probes = load('compact').probes
    blank = (ROOT/'models/compact/blank_inputs.csv').read_text()
    assert np.isnan(read_measurements(blank, probes)).all()
    values = read_measurements((ROOT/'models/compact/example_inputs.csv').read_text(), probes)
    assert np.isfinite(values).all()
    bad = ['gene_id,gene_effect\n', blank.replace(probes[0]['id'], '999999999', 1),
           blank.replace(f"{probes[0]['id']},{probes[0]['symbol']},", f"{probes[0]['id']},{probes[0]['symbol']},abc")]
    for text in bad:
        with pytest.raises(ValueError):
            read_measurements(text, probes)


def test_select_finds_informative_inputs():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(300, 20))
    y = x[:, [7, 12]] @ rng.normal(size=(2, 6)) + 0.1 * rng.normal(size=(300, 6))
    assert set(select(x, y, 2, 0.1)) == {7, 12}
    panel = fit(x[:, [7, 12]], y, 0.01)
    assert np.corrcoef(panel.predict(x[:, [7, 12]]).ravel(), y.ravel())[0, 1] > 0.99


def test_cli(tmp_path):
    out = tmp_path/'predictions.csv'
    subprocess.run([sys.executable, '-m', 'probe', '--model', 'models/broad/model.json',
                    '--input', 'models/broad/example_inputs.csv', '--output', str(out)], cwd=ROOT, check=True)
    rows = list(csv.DictReader(open(out)))
    assert len(rows) == 195 and float(rows[0]['predicted_effect']) <= float(rows[-1]['predicted_effect'])


def test_local_links_resolve():
    for page in [ROOT/'README.md', ROOT/'app/index.html', ROOT/'app/app.js', ROOT/'index.html']:
        for link in re.findall(r'(?:\]\(|href="|src=")(?!https?:)([^)"#]+)', page.read_text()):
            for scope in SCOPES:
                assert (page.parent/link.replace('${scope}', scope)).exists(), f'{page.name}: {link}'


def test_missing_and_invalid_training_ranges():
    model = load('broad')
    x = np.zeros(len(model.probes))
    unguarded = Model(model.metadata)
    assert unguarded.predict(x)[0]['calibration_status'] == 'unverified_range'
    assert unguarded.predict(x)[0]['lower'] is None
    support = json.loads((ROOT/'models/broad/support.json').read_text())
    for change in [{'model_version': 'wrong'}, {'lower': []}, {'upper': [float('nan')] * len(x)}]:
        with pytest.raises(ValueError):
            Model(model.metadata, support | change)


def test_invalid_model_does_not_silently_drop_predictions():
    payload = json.loads((ROOT/'models/broad/model.json').read_text())
    payload['calibration']['absolute_error_radii'] = [0.1]
    with pytest.raises(ValueError, match='radius'):
        Model(payload)
    payload = json.loads((ROOT/'models/broad/model.json').read_text())
    payload['model']['input_scale'][0] = 0
    with pytest.raises(ValueError, match='scales'):
        Model(payload)


@pytest.mark.parametrize('bad', ['missing_column', 'one_row', 'infinite', 'constant_output', 'nan_noise'])
def test_invalid_training_data_is_rejected(bad):
    x = np.arange(12, dtype=float).reshape(4, 3)
    y = x[:, :1].copy()
    noise = 0.1
    if bad == 'missing_column': x[:, 0] = np.nan
    if bad == 'one_row': x, y = x[:1], y[:1]
    if bad == 'infinite': x[0, 0] = np.inf
    if bad == 'constant_output': y[:] = 1
    if bad == 'nan_noise': noise = np.nan
    with pytest.raises(ValueError):
        select(x, y, 1, noise)


@pytest.mark.skipif(shutil.which('node') is None, reason='needs Node.js')
@pytest.mark.parametrize('scope', SCOPES)
def test_browser_scoring_matches_python(scope):
    script = '''
const assert = require('node:assert/strict');
const fs = require('node:fs');
const P = require('./app/probe.js'), E = require('./app/evaluate.js');
const scope = process.argv[1];
require(`./app/${scope}-data.js`); require(`./app/${scope}-outcomes.js`);
require('./app/validation.js'); require('./app/support.js');
const data = globalThis[scope.toUpperCase() + '_DATA'];
const outcomes = globalThis[scope.toUpperCase() + '_OUTCOMES'];
const support = globalThis.PROBE_SUPPORT[data.model.version];
const x = data.cases[1].values;
assert.equal(P.infer(data.model, x).status, 'unverified_range');
assert.equal(P.infer(data.model, x).bands, null);
assert.throws(() => P.infer(data.model, Array(x.length)));
assert.throws(() => P.infer(data.model, x, true, {...support, model_version:'wrong'}));
assert.throws(() => P.infer({...data.model, calibration:{absolute_error_radii:[0.1]}}, x));
const csv = fs.readFileSync(`models/${scope}/example_inputs.csv`, 'utf8');
assert.equal(P.readMeasurements(csv, data.model.probes).length, x.length);
assert.throws(() => P.readMeasurements(csv.replace(data.model.probes[0].symbol, 'wrong'), data.model.probes));
const predictions = data.cases.map(c => P.infer(data.model, c.values, true, support).scores);
const validation = globalThis.PROBE_VALIDATION[scope];
assert.throws(() => E.score(data.model, [...data.cases].reverse(), predictions, outcomes, validation));
console.log(JSON.stringify(E.score(data.model, data.cases, predictions, outcomes, validation)));
'''
    js = json.loads(subprocess.check_output(['node', '-e', script, scope], cwd=ROOT, text=True))
    py = evaluate(scope)
    assert [r['error'] for r in js['losses']] == pytest.approx(list(py['errors'].values()), abs=1e-12)
    assert js['coverage'] == pytest.approx(py['coverage'], abs=1e-12)
    assert (js['profiles'], js['donors'], js['values']) == (py['profiles'], py['donors'], py['observed'])
    for rule in ['conservative', 'wider']:
        assert (js[rule], js[rule + 'Correct']) == (py[rule]['flags'], py[rule]['correct'])


def test_python_scoring_rejects_reordered_profiles(monkeypatch):
    def reordered(path):
        bundle = read_bundle(path)
        if Path(path).name == 'validation.js':
            bundle['broad']['ids'].reverse()
        return bundle
    monkeypatch.setattr('probe.evaluate.read_bundle', reordered)
    with pytest.raises(ValueError, match='profiles'):
        evaluate('broad')
