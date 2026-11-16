"""Recompute the held-out test from the frozen model, the held-out inputs and the recorded outcomes.

    python -m probe.evaluate --scope broad
"""
import argparse
import json
from pathlib import Path

import numpy as np

from .model import DEPENDENCY, Model

ROOT = Path(__file__).resolve().parent.parent


def read_bundle(path):
    """Read an app data file of the form `globalThis.NAME = {...};`."""
    return json.loads(Path(path).read_text().split(' = ', 1)[1].rstrip().removesuffix(';'))


def normalized_error(predictions, outcomes, variance, donors):
    """Squared error / training variance, averaged over genes, then profiles per donor, then donors."""
    per_profile = np.nanmean((predictions - outcomes) ** 2 / variance, axis=1)
    return float(np.mean([per_profile[donors == d].mean() for d in np.unique(donors)]))


def evaluate(scope, root=ROOT):
    model = Model.load(root/'models'/scope/'model.json')
    cases = read_bundle(root/'app'/f'{scope}-data.js')['cases']
    outcome_map = read_bundle(root/'app'/f'{scope}-outcomes.js')
    validation = read_bundle(root/'app/validation.js')[scope]
    if [c['id'] for c in cases] != validation['ids'] or len(validation['donors']) != len(cases):
        raise ValueError('Validation profiles do not match the model data')
    inputs = np.array([[np.nan if v is None else v for v in c['values']] for c in cases])
    predictions = model.panel.predict(inputs)
    outcomes = np.array([[np.nan if v is None else v for v in outcome_map[c['id']]] for c in cases])
    variance, donors = np.asarray(validation['variance']), np.asarray(validation['donors'])

    observed = np.isfinite(outcomes)
    truth = outcomes < DEPENDENCY
    conservative = (predictions + model.radius < DEPENDENCY) & observed
    wider = (predictions < DEPENDENCY) & observed
    errors = {f'PROBE ({len(model.probes)} inputs)': normalized_error(predictions, outcomes, variance, donors)}
    for control in validation['controls']:
        errors[control['label']] = normalized_error(np.asarray(control['predictions']), outcomes, variance, donors)
    return dict(
        scope=scope, profiles=len(cases), donors=len(np.unique(donors)), observed=int(observed.sum()),
        dependencies=int((truth & observed).sum()), errors=errors,
        coverage=float((np.abs(predictions - outcomes) <= model.radius)[observed].mean()),
        conservative=dict(flags=int(conservative.sum()), correct=int((conservative & truth).sum())),
        wider=dict(flags=int(wider.sum()), correct=int((wider & truth).sum())),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--scope', choices=['broad', 'compact'], default='broad')
    r = evaluate(parser.parse_args().scope)
    print(f"{r['profiles']} held-out profiles, {r['donors']} donors, {r['observed']:,} observed responses\n")
    print('Normalized error (lower is better)')
    for name, value in r['errors'].items():
        print(f'  {name:<40} {value:.3f}')
    print(f"\n90% band coverage   {100 * r['coverage']:.1f}%")
    for rule in ['conservative', 'wider']:
        c = r[rule]
        print(f"{rule.capitalize():<12} {c['correct']} of {c['flags']} flags correct "
              f"({100 * c['correct'] / c['flags']:.1f}% precision, {100 * c['correct'] / r['dependencies']:.1f}% recall)")


if __name__ == '__main__':
    main()
