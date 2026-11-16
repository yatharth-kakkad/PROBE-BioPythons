// Scores predictions on the held-out confirmation profiles. Mirrors probe/evaluate.py.
(function (root) {
  'use strict';
  // Normalized error: squared error / training variance per gene, averaged over genes,
  // then over profiles within a donor, then over donors.
  function normalizedError(cases, predictions, outcomes, validation) {
    const donors = new Map();
    cases.forEach((c, i) => {
      let loss = 0, count = 0;
      outcomes[c.id].forEach((a, j) => { if (a !== null) { loss += (predictions[i][j] - a) ** 2 / validation.variance[j]; count++; } });
      const d = validation.donors[i];
      donors.set(d, [...(donors.get(d) || []), loss / count]);
    });
    const means = [...donors.values()].map(v => v.reduce((a, b) => a + b, 0) / v.length);
    return means.reduce((a, b) => a + b, 0) / means.length;
  }

  function score(model, cases, predictions, outcomes, validation) {
    if (cases.length !== validation.ids.length || cases.length !== validation.donors.length ||
        cases.some((c, i) => c.id !== validation.ids[i])) throw Error('Validation profiles do not match the model data');
    const radii = model.calibration.absolute_error_radii;
    const t = {values: 0, dependencies: 0, covered: 0, conservative: 0, conservativeCorrect: 0, wider: 0, widerCorrect: 0};
    cases.forEach((c, i) => outcomes[c.id].forEach((a, j) => {
      if (a === null) return;
      const p = predictions[i][j], dep = a < -.5;
      t.values++; t.dependencies += dep;
      t.covered += Math.abs(p - a) <= radii[j];
      t.wider += p < -.5; t.widerCorrect += p < -.5 && dep;
      t.conservative += p + radii[j] < -.5; t.conservativeCorrect += p + radii[j] < -.5 && dep;
    }));
    const methods = [{label: `PROBE · ${model.probes.length} inputs`, predictions}, ...validation.controls];
    return {
      profiles: cases.length, donors: new Set(validation.donors).size, ...t,
      losses: methods.map(m => ({label: m.label, error: normalizedError(cases, m.predictions, outcomes, validation)})),
      coverage: t.covered / t.values,
    };
  }

  root.ProbeEvaluate = {score};
  if (typeof module !== 'undefined') module.exports = {score};
})(globalThis);
