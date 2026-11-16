(function () {
  'use strict';
  const P = window.Probe, E = window.ProbeEvaluate;
  const $ = id => document.getElementById(id);
  const fmt = v => v == null ? 'N/A' : v.toFixed(3), pct = v => v == null ? 'N/A' : (100 * v).toFixed(1) + '%';
  const esc = v => String(v).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const DATA = {broad: window.BROAD_DATA, compact: window.COMPACT_DATA};
  const OUTCOME_FILES = {broad: ['broad-outcomes.js', 'BROAD_OUTCOMES'], compact: ['compact-outcomes.js', 'COMPACT_OUTCOMES']};

  let scope, D, M, profile, values, original, imported, edited, revealed, selected, result;
  let token = 0; // increments whenever the profile changes, so stale async work is ignored

  const support = () => window.PROBE_SUPPORT[M.version];
  const outcomes = () => window[OUTCOME_FILES[scope][1]];
  const effective = () => values.map((v, i) => i < Number($('budget').value) ? v : null);

  const scripts = new Map();
  function loadScript(src) {
    if (!scripts.has(src)) scripts.set(src, new Promise((resolve, reject) => {
      const s = document.createElement('script');
      s.src = src; s.onload = resolve;
      s.onerror = () => { scripts.delete(src); s.remove(); reject(Error('Could not load ' + src)); };
      document.head.appendChild(s);
    }));
    return scripts.get(src);
  }

  function download(name, rows) {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([P.toCSV(rows)], {type: 'text/csv'}));
    a.download = name; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  function showError(message) {
    token++; result = null;
    $('notice').textContent = message; $('notice').className = 'notice error';
    ['ranking', 'detail', 'atlas', 'distribution', 'response-overview'].forEach(id => $(id).replaceChildren());
    $('export').disabled = $('reveal').disabled = true;
    $('available').textContent = 'INPUT ERROR'; $('decision-evidence').textContent = '';
    $('count').textContent = 'Correct the input to calculate predictions'; $('calls').textContent = '';
  }

  function changeScope() {
    scope = $('scope').value; D = DATA[scope]; M = D.model; token++;
    const e = D.evidence;
    const panel = e.contrasts.find(c => c.comparator === 'd_optimal');
    const molecular = e.contrasts.filter(c => c.comparator.startsWith('molecular'))
      .sort((a, b) => a.relative_error_reduction - b.relative_error_reduction)[0];
    const comparison = [
      [`PROBE (${e.probes} inputs)`, 1],
      [`Strongest panel (${panel.comparator_budget} inputs)`, 1 / (1 - panel.relative_error_reduction)],
      ['Best RNA + mutation model', 1 / (1 - molecular.relative_error_reduction)],
    ];
    const maxError = Math.max(...comparison.map(([, error]) => error));
    $('benchmark-chart').innerHTML = '<p class="muted">Relative prediction error. Lower is better; PROBE = 1.</p>' +
      comparison.map(([label, error]) => `<div class="validation-bar"><div>${label}</div><div><span style="width:${100 * error / maxError}%"></span></div><code>${error.toFixed(2)}×</code></div>`).join('');
    $('headline-evidence').innerHTML =
      `<div><strong>${pct(panel.relative_error_reduction)}</strong><p>Lower error than the strongest ${panel.comparator_budget}-input panel</p></div>` +
      `<div><strong>${pct(molecular.relative_error_reduction)}</strong><p>Lower error than the best RNA + mutation model</p></div>` +
      `<div><strong>${e.donors}</strong><p>Held-out donors, ${e.targets} response genes</p></div>`;
    $('budget').innerHTML = [M.probes.length, M.probes.length / 2, 4, 0].map((n, i) =>
      `<option value="${n}">${i === 0 ? `All ${n} measurements` : n ? `First ${n} measurements` : 'No measurements · training average'}</option>`).join('');
    const tissues = [...new Set(D.cases.map(c => c.tissue))].sort();
    $('tissue').innerHTML = '<option value="">All lineages</option>' + tissues.map(t => `<option value="${esc(t)}">${esc(t.replaceAll('_', ' '))}</option>`).join('');
    $('evidence-description').textContent = `Trained on 692 cancer cell models from 665 donors, then frozen and tested on ${e.profiles} models from ${e.donors} other donors, predicting ${e.targets} responses each.`;
    $('evidence-links').innerHTML = `<a href="../models/${scope}/model.json">Frozen model</a><a href="../results/${scope}/comparison.csv">Every comparator</a><a href="../results/${scope}/per_gene.csv">Every gene</a>`;
    $('validation-result').replaceChildren(); $('validate').disabled = false;
    populateProfiles();
  }

  function populateProfiles() {
    const tissue = $('tissue').value, cases = D.cases.filter(c => !tissue || c.tissue === tissue);
    $('cell').innerHTML = cases.map(c => `<option value="${esc(c.id)}">${esc(c.name)} · ${esc(c.id)}</option>`).join('');
    // Default to the first profile whose inputs fall inside the training range (never chosen by outcomes).
    const first = cases.find(c => P.infer(M, c.values, true, support()).calibrated);
    if (first) $('cell').value = first.id;
    loadProfile();
  }

  function loadProfile() {
    token++;
    profile = D.cases.find(c => c.id === $('cell').value);
    values = profile.values.slice(); original = values.slice();
    imported = edited = revealed = false; selected = 0;
    $('budget').value = String(M.probes.length); $('search').value = ''; $('decision').value = 'all';
    $('profile-description').textContent = `${D.cases.length} held-out profiles available. Calculations run on this device.`;
    $('cell').querySelector('option[value=imported]')?.remove();
    renderInputs(); calculate();
  }

  function renderInputs() {
    const budget = Number($('budget').value);
    $('queries').innerHTML = M.probes.map((g, i) =>
      `<div class="query-row"><label for="q-${i}" title="Input ${i + 1}, Entrez ${esc(g.id)}">${esc(g.symbol)}</label>` +
      `<input id="q-${i}" data-index="${i}" type="number" step="any" value="${values[i] ?? ''}" placeholder="Missing" ${i >= budget ? 'disabled' : ''} aria-label="${esc(g.symbol)} gene effect"></div>`).join('');
    $('queries').querySelectorAll('input').forEach(input => input.addEventListener('input', () => {
      if (input.validity.badInput) return showError('Use a number or leave the field empty.');
      token++;
      values[Number(input.dataset.index)] = input.value === '' ? null : Number(input.value);
      edited = values.some((v, i) => v !== original[i]);
      calculate();
    }));
  }

  function calculate() {
    try {
      result = P.infer(M, effective(), !edited, support());
      const status = {
        edited_inputs: 'What-if scenario. Error bands are off.',
        missing_inputs: `Partial panel. Error bands need all ${M.probes.length} measurements.`,
        outside_training_range: `${result.outside.map(i => M.probes[i].symbol).join(', ')} outside the training range. Point predictions only.`,
        complete_panel: 'Complete panel. 90% error bands are on.',
      }[result.status];
      $('notice').className = 'notice';
      $('notice').textContent = `${status}${revealed && edited ? ' Observed outcomes belong to the original profile.' : ''}`;
      $('available').textContent = `${result.observed} / ${M.probes.length}`;
      $('count').textContent = `${result.observed} measured inputs → ${M.targets.length} predictions`;
      $('calls').textContent = result.bands ? `${result.bands.filter(b => b.dependency).length} conservative candidates` : 'Point predictions · bands off';
      $('response-overview').innerHTML = [
        [result.observed, 'Measured knockouts'], [M.targets.length, 'Predicted responses'],
        [result.bands ? result.bands.filter(b => b.dependency).length : 'N/A', 'Conservative candidates'],
      ].map(([value, label]) => `<div><strong>${value}</strong><span>${label}</span></div>`).join('');
      $('export').disabled = false; $('reveal').disabled = imported;
      render();
    } catch (e) { showError(e.message); }
  }

  const ranked = () => M.targets.map((_, i) => i).sort((a, b) => result.scores[a] - result.scores[b] || a - b);
  const matches = (i, q) => M.targets[i].symbol.includes(q) || M.targets[i].id.includes(q);
  const actual = i => revealed && !imported ? outcomes()[profile.id][i] : null;

  function shown() {
    const q = $('search').value.trim().toUpperCase(), rule = $('decision').value;
    return ranked().filter(i => matches(i, q) &&
      (rule === 'all' || rule === 'wider' && result.scores[i] < -.5 || rule === 'conservative' && result.bands?.[i].dependency));
  }

  function render() { renderTable(); renderAtlas(); renderDetail(); }
  function choose(i) {
    const container = document.activeElement.closest('#ranking, #atlas')?.id;
    selected = i; render();
    if (container) $(container).querySelector(`[data-gene="${i}"]`)?.focus({preventScroll: true});
  }

  function renderAtlas() {
    const visible = new Set(shown());
    renderDistribution([...visible].map(i => result.scores[i]));
    if (!visible.size) { $('atlas').innerHTML = '<p class="muted">No responses match this view.</p>'; return; }
    $('atlas').innerHTML = M.targets.map((g, i) => {
      const strength = Math.min(1, Math.max(0, -result.scores[i] / 1.1));
      const colors = ['#f0e4d7', '#ddbc9f', '#c28c63', '#9a5a32', '#63381f'];
      const level = Math.min(4, Math.floor(strength * 5));
      const title = `${g.symbol}: predicted ${fmt(result.scores[i])}${revealed ? ', observed ' + fmt(actual(i)) : ''}`;
      return `<button data-gene="${i}" aria-pressed="${selected === i}" class="${result.bands?.[i].dependency ? 'is-candidate' : ''}" title="${esc(title)}" ` +
        `style="background:${colors[level]};color:${level >= 3 ? '#fff8f0' : '#25211d'}" ${!visible.has(i) ? 'hidden' : ''}>${esc(g.symbol)}</button>`;
    }).join('');
    $('atlas').querySelectorAll('[data-gene]').forEach(b => b.addEventListener('click', () => choose(Number(b.dataset.gene))));
  }

  function renderDistribution(scores) {
    if (!scores.length) { $('distribution').innerHTML = '<p class="empty-state">No responses match this view.</p>'; return; }
    const min = Math.floor(Math.min(...scores, -.5) * 4) / 4;
    const max = Math.max(min + .25, Math.ceil(Math.max(...scores, -.5) * 4) / 4);
    const bins = Array(12).fill(0), width = (max - min) / bins.length;
    scores.forEach(v => bins[Math.min(11, Math.floor((v - min) / width))]++);
    const peak = Math.max(...bins), dependent = scores.filter(v => v < -.5).length;
    $('distribution').innerHTML = `<div class="distribution-title"><h3>Response distribution</h3><span class="muted">${dependent} of ${scores.length} below -0.5</span></div>` +
      `<div class="histogram" role="img" aria-label="Distribution of ${scores.length} predicted effects; ${dependent} below the dependency threshold of -0.5">` +
      bins.map((count, i) => `<div class="histogram-column" title="${fmt(min + i * width)} to ${fmt(min + (i + 1) * width)}: ${count} responses"><i style="height:${100 * count / peak}%"><span>${count || ''}</span></i></div>`).join('') +
      `<div class="threshold" style="left:${100 * (-.5 - min) / (max - min)}%" title="Dependency threshold: -0.5"></div></div>` +
      `<div class="histogram-axis"><span>${fmt(min)}</span><span>Predicted effect</span><span>${fmt(max)}</span></div>`;
  }

  function renderTable() {
    const rule = $('decision').value, rows = shown(), stats = D.evidence[rule === 'conservative' ? 'conservative' : 'point'];
    $('outcome-heading').textContent = revealed ? 'Observed outcome' : 'Outcome hidden';
    $('reveal').textContent = revealed ? 'Hide outcomes' : 'Reveal outcomes';
    $('decision-evidence').textContent =
      rule === 'all' ? `Showing ${rows.length} of ${M.targets.length} responses. Choose a follow-up view to filter and export a shortlist.`
      : !result.calibrated ? 'Shortlist accuracy below applies only to complete, unedited, in-range panels.'
      : `On held-out data this rule was ${pct(stats.precision)} precise (${stats.correct} of ${stats.calls} flags correct) and caught ${pct(stats.recall)} of ${stats.total_dependencies} true dependencies.`;
    $('ranking').innerHTML = rows.map((i, rank) => {
      const b = result.bands?.[i], g = M.targets[i];
      return `<tr class="${i === selected ? 'selected' : ''}"><td>${rank + 1}</td><td><button class="row-button" data-gene="${i}" aria-pressed="${i === selected}" title="Entrez ${esc(g.id)}">${esc(g.symbol)}</button>` +
        `${b?.dependency ? '<span class="candidate" title="Upper error band below -0.5">Candidate</span>' : ''}</td><td>${fmt(result.scores[i])}</td>` +
        `<td>${b ? fmt(b.lower) + ' to ' + fmt(b.upper) : 'Off'}</td><td>${revealed ? fmt(actual(i)) : 'Hidden'}</td></tr>`;
    }).join('') || '<tr><td colspan="5">No responses match this view.</td></tr>';
    $('ranking').querySelectorAll('[data-gene]').forEach(b => b.addEventListener('click', () => choose(Number(b.dataset.gene))));
  }

  function renderDetail() {
    const x = effective(), g = M.targets[selected], e = D.evidence.per_gene[selected];
    const parts = result.weights[selected].map((w, i) => ({i, v: x[i] === null ? 0 : w * (x[i] - M.model.input_mean[i])}));
    const top = parts.slice().sort((a, b) => Math.abs(b.v) - Math.abs(a.v)).slice(0, 4);
    const largest = Math.max(...top.map(p => Math.abs(p.v)), 1e-12);
    $('detail').innerHTML =
      `<div class="detail-heading"><h3>${esc(g.symbol)}</h3><span>Entrez ${esc(g.id)}</span></div>` +
      `<dl class="detail-metrics"><div><dt>Predicted effect</dt><dd>${fmt(result.scores[selected])}</dd></div>` +
      `<div><dt>Training average</dt><dd>${fmt(M.model.output_mean[selected])}</dd></div>` +
      `<div><dt>Input contribution</dt><dd>${fmt(parts.reduce((s, p) => s + p.v, 0))}</dd></div>` +
      `${revealed ? `<div><dt>Observed outcome</dt><dd>${fmt(actual(selected))}</dd></div>` : ''}</dl>` +
      '<p class="muted">Largest input contributions</p><div class="contribution-chart"><div class="chart-axis"><span>Stronger dependency</span><span>Weaker</span></div>' +
      top.map(p => `<div class="contribution-row"><span>${esc(M.probes[p.i].symbol)}</span><div class="contribution-track" aria-hidden="true"><span class="${p.v < 0 ? 'negative' : 'positive'}" style="width:${100 * Math.abs(p.v) / largest}%"></span></div><code>${p.v >= 0 ? '+' : ''}${fmt(p.v)}</code></div>`).join('') + '</div>' +
      `<p class="muted">Held-out error ${pct(Math.abs(e.reduction_vs_mean))} ${e.reduction_vs_mean < 0 ? 'higher' : 'lower'} than the training average. Correlation ${fmt(e.correlation)}.</p>`;
  }

  const inputRows = x => [['gene_id', 'gene', 'gene_effect'], ...M.probes.map((g, i) => [g.id, g.symbol, x[i]])];

  $('import').addEventListener('click', () => $('upload').click());
  $('scope').addEventListener('change', changeScope);
  $('tissue').addEventListener('change', populateProfiles);
  $('cell').addEventListener('change', loadProfile);
  $('reset').addEventListener('click', () => { if (imported) $('cell').value = profile.id; loadProfile(); });
  $('budget').addEventListener('change', () => { token++; renderInputs(); calculate(); });
  $('search').addEventListener('input', () => { if (result) { renderTable(); renderAtlas(); } });
  $('decision').addEventListener('change', () => { if (result) { renderTable(); renderAtlas(); } });
  $('inputs-export').addEventListener('click', () => download(`probe-${scope}-inputs.csv`, inputRows(effective())));
  $('blank-template').addEventListener('click', () => download(`probe-${scope}-template.csv`, inputRows(M.probes.map(() => null))));

  $('upload').addEventListener('change', async event => {
    const file = event.target.files[0], request = ++token;
    if (!file) return;
    try {
      if (file.size > 1e6) throw Error('The CSV is larger than 1 MB');
      const text = await file.text();
      if (request !== token) return;
      values = P.readMeasurements(text, M.probes); original = values.slice();
      imported = true; edited = revealed = false;
      $('cell').querySelector('option[value=imported]')?.remove();
      $('cell').insertAdjacentHTML('beforeend', '<option value="imported" disabled>Imported measurements</option>');
      $('cell').value = 'imported'; $('budget').value = String(M.probes.length);
      $('profile-description').textContent = `Imported ${file.name}. Your data stays on this device.`;
      renderInputs(); calculate();
    } catch (e) { if (request === token) showError(e.message); }
    finally { $('upload').value = ''; }
  });

  $('export').addEventListener('click', () => download(`probe-${scope}-${$('decision').value}-predictions.csv`, [
    ['gene_id', 'gene', 'predicted_effect', 'lower', 'upper', 'dependency_call', 'measured_inputs', 'calibration_status', 'model_version', 'profile'],
    ...shown().map(i => [M.targets[i].id, M.targets[i].symbol, result.scores[i], result.bands?.[i].lower, result.bands?.[i].upper,
      result.bands?.[i].dependency, result.observed, result.status, M.version, imported ? 'imported' : profile.id])]));

  $('reveal').addEventListener('click', async () => {
    if (imported) return;
    if (revealed) { revealed = false; return calculate(); }
    const request = token;
    try {
      await loadScript(OUTCOME_FILES[scope][0]);
      if (request === token) { revealed = true; calculate(); }
    } catch (e) { if (request === token) showError(e.message); }
  });

  $('validate').addEventListener('click', async () => {
    const key = scope, data = D, out = $('validation-result');
    $('validate').disabled = true;
    out.innerHTML = '<p>Loading held-out outcomes and comparator predictions…</p>';
    try {
      await Promise.all([loadScript(OUTCOME_FILES[key][0]), loadScript('validation.js')]);
      const predictions = [];
      for (let i = 0; i < data.cases.length; i++) {
        if (scope !== key) return;
        predictions.push(P.infer(data.model, data.cases[i].values).scores);
        if (i % 16 === 0) {
          out.innerHTML = `<p>Predicting profile ${i + 1} of ${data.cases.length}…</p><progress value="${i + 1}" max="${data.cases.length}"></progress>`;
          await new Promise(r => setTimeout(r, 0));
        }
      }
      const s = E.score(data.model, data.cases, predictions, window[OUTCOME_FILES[key][1]], window.PROBE_VALIDATION[key]);
      const max = Math.max(...s.losses.map(l => l.error));
      out.innerHTML = `<p><strong>Recomputed on ${s.profiles} profiles, ${s.donors} donors, ${s.values.toLocaleString()} observed responses.</strong> Lower error is better.</p>` +
        `<div class="validation-bars">${s.losses.map(l => `<div class="validation-bar"><div>${esc(l.label)}</div><div><span style="width:${100 * l.error / max}%"></span></div><code>${fmt(l.error)}</code></div>`).join('')}</div>` +
        `<div class="validation-summary"><div><strong>${pct(s.coverage)}</strong><p>of outcomes inside the 90% error bands</p></div>` +
        `<div><strong>${pct(s.conservativeCorrect / s.conservative)} / ${pct(s.conservativeCorrect / s.dependencies)}</strong><p>Conservative precision / recall<br>${s.conservativeCorrect} of ${s.conservative} flags correct</p></div>` +
        `<div><strong>${pct(s.widerCorrect / s.wider)} / ${pct(s.widerCorrect / s.dependencies)}</strong><p>Wider shortlist precision / recall<br>${s.widerCorrect} of ${s.wider} flags correct</p></div></div>`;
    } catch (e) { if (scope === key) out.textContent = e.message; }
    finally { if (scope === key) $('validate').disabled = false; }
  });

  changeScope();
})();
