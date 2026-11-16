// Browser inference. Mirrors probe/model.py exactly; probe/verify.py checks the two agree.
(function (root) {
  'use strict';
  const NUMBER = /^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$/;

  function cholesky(a) {
    const n = a.length, l = Array.from({length: n}, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
      let v = a[i][j];
      for (let k = 0; k < j; k++) v -= l[i][k] * l[j][k];
      if (i === j && !(v > 0)) throw Error('Model covariance is not positive definite');
      l[i][j] = i === j ? Math.sqrt(v) : v / l[j][j];
    }
    return l;
  }

  function solve(l, b) {
    const n = l.length, y = Array(n).fill(0), x = Array(n).fill(0);
    for (let i = 0; i < n; i++) { y[i] = b[i]; for (let j = 0; j < i; j++) y[i] -= l[i][j] * y[j]; y[i] /= l[i][i]; }
    for (let i = n - 1; i >= 0; i--) { x[i] = y[i]; for (let j = i + 1; j < n; j++) x[i] -= l[j][i] * x[j]; x[i] /= l[i][i]; }
    return x;
  }

  // values: one gene effect per input gene, or null when not measured.
  // measured=false marks hand-edited inputs, which turns error bands off.
  function infer(model, values, measured = true, support = null) {
    const n = model.probes.length, m = model.model;
    if (!Array.isArray(values) || values.length !== n || Array.from(values).some(v => v !== null && !Number.isFinite(v)))
      throw Error(`Supply ${n} finite gene effects, or blanks for missing measurements`);
    const radii = model.calibration.absolute_error_radii;
    if (radii.length !== model.targets.length || radii.some(r => !Number.isFinite(r) || r < 0))
      throw Error('Expected one finite, nonnegative error radius per response gene');
    if (support && (support.schema !== 'probe.input-support.v1' || support.model_version !== model.version ||
        support.probe_ids?.length !== n || support.probe_ids.some((id, i) => id !== model.probes[i].id) ||
        support.lower?.length !== n || support.upper?.length !== n ||
        support.lower.some((v, i) => !Number.isFinite(v) || !Number.isFinite(support.upper[i]) || v >= support.upper[i])))
      throw Error('Training-range file does not match the model or contains invalid bounds');
    const observed = values.map((v, i) => v === null ? -1 : i).filter(i => i >= 0);
    const l = observed.length ? cholesky(observed.map(i => observed.map(j => m.covariance[i][j] + (i === j ? m.noise : 0)))) : [];
    // weights[t][i]: change in prediction t per unit change in input i.
    const weights = model.targets.map((_, t) => {
      const w = Array(n).fill(0);
      if (observed.length) solve(l, observed.map(i => m.cross_covariance[i][t]))
        .forEach((v, k) => { w[observed[k]] = v / m.input_scale[observed[k]]; });
      return w;
    });
    const scores = m.output_mean.map((mean, t) => mean + observed.reduce((s, i) => s + weights[t][i] * (values[i] - m.input_mean[i]), 0));
    if (scores.some(v => !Number.isFinite(v))) throw Error('Measurements exceed the numerical range of this model');
    const outside = support ? observed.filter(i => values[i] < support.lower[i] || values[i] > support.upper[i]) : [];
    const status = !measured ? 'edited_inputs' : observed.length !== n ? 'missing_inputs' : !support ? 'unverified_range' : outside.length ? 'outside_training_range' : 'complete_panel';
    const bands = status === 'complete_panel'
      ? scores.map((p, t) => ({lower: p - radii[t], upper: p + radii[t], dependency: p + radii[t] < -.5})) : null;
    return {scores, weights, bands, observed: observed.length, calibrated: !!bands, status, outside};
  }

  function parseCSV(text) {
    const rows = []; let row = [], cell = '', state = 'start';
    const field = () => { row.push(cell); cell = ''; state = 'start'; };
    const line = () => { field(); if (row.some(v => v.trim())) rows.push(row); row = []; };
    for (let i = 0; i < text.length; i++) {
      const c = text[i];
      if (state === 'quoted') { if (c !== '"') cell += c; else if (text[i + 1] === '"') { cell += '"'; i++; } else state = 'closed'; }
      else if (c === ',') field();
      else if (c === '\n' || c === '\r') { if (c === '\r' && text[i + 1] === '\n') i++; line(); }
      else if (state === 'closed') throw Error('Unexpected characters after a quoted CSV field');
      else if (c === '"') { if (state !== 'start') throw Error('Quote inside an unquoted CSV field'); state = 'quoted'; }
      else { cell += c; state = 'plain'; }
    }
    if (state === 'quoted') throw Error('Unclosed CSV quotation');
    line(); return rows;
  }

  // CSV with columns gene_id, gene_effect and optional gene; every input gene exactly once.
  function readMeasurements(text, probes) {
    const rows = parseCSV(text.replace(/^﻿/, '')), header = rows.shift() || [];
    if (new Set(header).size !== header.length || !header.includes('gene_id') || !header.includes('gene_effect') ||
        header.some(h => !['gene_id', 'gene', 'gene_effect'].includes(h))) throw Error('CSV requires gene_id and gene_effect; optional gene column');
    const expected = new Map(probes.map(g => [g.id, g.symbol])), found = new Map();
    for (const row of rows) {
      if (row.length !== header.length) throw Error('Misaligned CSV row');
      const item = Object.fromEntries(header.map((h, i) => [h, row[i]])), id = item.gene_id;
      if (!expected.has(id) || found.has(id)) throw Error('Unknown or duplicate gene ID: ' + id);
      if ('gene' in item && item.gene !== expected.get(id)) throw Error('Gene name does not match ID: ' + id);
      const raw = item.gene_effect.trim();
      if (raw && (!NUMBER.test(raw) || !Number.isFinite(Number(raw)))) throw Error('Use decimal gene effects, or blanks for missing values');
      found.set(id, raw ? Number(raw) : null);
    }
    if (found.size !== probes.length) throw Error(`Include all ${probes.length} genes exactly once; leave missing effects blank`);
    return probes.map(g => found.get(g.id));
  }

  const toCSV = rows => rows.map(r => r.map(v => '"' + String(v ?? '').replaceAll('"', '""') + '"').join(',')).join('\n') + '\n';

  const api = {infer, readMeasurements, toCSV};
  root.Probe = api;
  if (typeof module !== 'undefined') module.exports = api;
})(globalThis);
