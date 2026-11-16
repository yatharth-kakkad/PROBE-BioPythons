"""Load a frozen PROBE model and predict with calibrated 90% error bands."""
import csv
import io
import json
import re
from pathlib import Path

import numpy as np

from .panel import Panel

NUMBER = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?\Z', re.ASCII)
DEPENDENCY = -0.5  # gene effect below this counts as a dependency


class Model:
    def __init__(self, payload, support=None):
        if payload.get('schema') != 'probe.functional.v1':
            raise ValueError('Unsupported model file')
        self.metadata = payload
        self.version = payload['version']
        self.probes, self.targets = payload['probes'], payload['targets']
        m = payload['model']
        self.panel = Panel(**{k: np.asarray(m[k], float) for k in
                              ['input_mean', 'input_scale', 'output_mean', 'covariance', 'cross_covariance', 'noise']})
        n, t = len(self.probes), len(self.targets)
        if self.panel.covariance.shape != (n, n) or self.panel.cross_covariance.shape != (n, t):
            raise ValueError('Model dimensions do not match its gene lists')
        self.radius = np.asarray(payload['calibration']['absolute_error_radii'], float)
        if self.radius.shape != (t,) or not np.isfinite(self.radius).all() or (self.radius < 0).any():
            raise ValueError('Expected one finite, nonnegative error radius per response gene')
        self.lower = self.upper = None
        if support is not None:
            if (support.get('schema') != 'probe.input-support.v1' or
                    support.get('model_version') != self.version or support.get('probe_ids') != [g['id'] for g in self.probes]):
                raise ValueError('Training-range file belongs to a different model')
            self.lower, self.upper = np.asarray(support['lower'], float), np.asarray(support['upper'], float)
            if (self.lower.shape != (n,) or self.upper.shape != (n,) or
                    not np.isfinite([self.lower, self.upper]).all() or (self.lower >= self.upper).any()):
                raise ValueError('Expected finite lower and upper training ranges for every input')

    @classmethod
    def load(cls, path):
        path = Path(path)
        support = path.with_name('support.json')
        return cls(json.loads(path.read_text()), json.loads(support.read_text()) if support.exists() else None)

    def status(self, x, measured=True):
        if not measured:
            return 'edited_inputs'
        if not np.isfinite(x).all():
            return 'missing_inputs'
        if self.lower is None:
            return 'unverified_range'
        if ((x < self.lower) | (x > self.upper)).any():
            return 'outside_training_range'
        return 'complete_panel'

    def predict(self, values, measured=True):
        """One profile (NaN = not measured). Bands only for a complete, measured, in-range panel."""
        x = np.asarray(values, float)
        with np.errstate(over='ignore', invalid='ignore'):
            prediction = self.panel.predict(x[None])[0]
        if not np.isfinite(prediction).all():
            raise ValueError('Measurements exceed the numerical range of this model')
        status = self.status(x, measured)
        banded = status == 'complete_panel'
        observed = int(np.isfinite(x).sum())
        return [dict(gene_id=g['id'], gene=g['symbol'], predicted_effect=float(p),
                     lower=float(p - r) if banded else None, upper=float(p + r) if banded else None,
                     dependency_call=bool(p + r < DEPENDENCY) if banded else None,
                     measured_inputs=observed, calibration_status=status)
                for g, p, r in zip(self.targets, prediction, self.radius)]


def read_measurements(text, probes):
    """CSV with gene_id, gene_effect and optional gene. Every input gene exactly once; blank = missing."""
    rows = [r for r in csv.reader(io.StringIO(text.lstrip('﻿')), strict=True) if any(c.strip() for c in r)]
    if not rows:
        raise ValueError('Empty CSV')
    header, *rows = rows
    if len(set(header)) != len(header) or not {'gene_id', 'gene_effect'} <= set(header) <= {'gene_id', 'gene', 'gene_effect'}:
        raise ValueError('CSV requires gene_id and gene_effect; optional gene column')
    expected = {g['id']: g['symbol'] for g in probes}
    found = {}
    for row in rows:
        if len(row) != len(header):
            raise ValueError('Misaligned CSV row')
        item = dict(zip(header, row))
        gid = item['gene_id']
        if gid not in expected or gid in found:
            raise ValueError('Unknown or duplicate gene ID: ' + gid)
        if item.get('gene', expected[gid]) != expected[gid]:
            raise ValueError('Gene name does not match ID: ' + gid)
        raw = item['gene_effect'].strip()
        if raw and (not NUMBER.fullmatch(raw) or not np.isfinite(float(raw))):
            raise ValueError('Use decimal gene effects, or blanks for missing values')
        found[gid] = float(raw) if raw else np.nan
    if found.keys() != expected.keys():
        raise ValueError(f'Include all {len(probes)} genes exactly once; leave missing effects blank')
    return np.array([found[g['id']] for g in probes])
