# PROBE

**Measure 32 gene knockouts in a cancer cell model. Predict 195 more.**

## [Open PROBE — live website](https://yatharth-kakkad.github.io/PROBE-BioPythons/)

Use the app directly in your browser. No download, installation, or sign-in required.

**Team:** Bio Pythons · Andy Pham & Yatharth Kakkad · Pine View School
**Track:** AI + Health · Veritas AI International High School Hackathon 2026

> This project is a prototype and is not medical advice.

## The problem

A cancer researcher cannot knock out every gene in every cell model. They have to choose which few experiments to run first, and which genes to follow up as possible drug targets. Chosen badly, the first experiments say little about the rest.

## What we built

PROBE is a browser tool that turns a small panel of measured knockout effects into predictions for knockouts that were not run.

- **Input:** knockout effects for 32 chosen genes in one cell model (or 16, in the compact model).
- **Output:** predicted effects for 195 other genes (46 in the compact model). Complete, measured panels within the training ranges also receive nominal 90% error bands.
- **Decision:** two shortlists of likely dependencies (genes the cell needs to survive) to test next.

The user can import their own CSV, edit inputs and watch predictions update, see which inputs drive each prediction, reveal the real outcome for held-out cell models, and rerun the full held-out test in the browser. It runs offline.

**User:** a researcher deciding which knockout experiments to run next on a cancer cell model.

## What we added

Predicting many knockout effects from a few has prior work (Zhao et al. 2022; Hossain et al. 2025, goloco). We add: response genes qualified as predictable on training data only, a measurement panel chosen for those genes, a frozen test on held-out donors against larger panels and RNA/mutation models, and a tool that shows every prediction's inputs and recomputes its own test.

## How it works

1. **Data.** DepMap 24Q2 CRISPR knockout effects (Chronos scores) for cancer cell lines. Cell models are split by donor, so no donor appears in more than one split.
2. **Choose outputs.** On training data, keep genes that vary across cell lines and that a single other knockout already predicts well: 195 of 808 candidates.
3. **Choose inputs.** Greedily pick the 32 genes that most reduce the remaining uncertainty across those 195 outputs (A-optimal design).
4. **Predict.** A regularized linear model (conditional Gaussian) predicts each output from the inputs. If some inputs are missing, it uses the ones measured.
5. **Error bands.** Sized on a separate development set (71 models, 70 donors).
6. **Test.** The model, gene lists and test plan were committed before the held-out outcomes were read.

Training used 692 cell models from 665 donors.

## Results

Held-out test: **124 cell models from 120 new donors, 24,177 observed responses.** Normalized error is squared error divided by each gene's training variance; lower is better.

| Predictor | Measured knockouts | Normalized error |
|---|---:|---:|
| Training average | 0 | 1.211 |
| Best RNA + mutation model | 0 | 1.030 |
| Random panel (5 seeds) | 64 | 1.053 |
| Strongest other panel (D-optimal) | 64 | 0.796 |
| **PROBE** | **32** | **0.725** |

- **8.9% lower error** than the strongest panel with twice the measurements (95% familywise interval 6.0 to 11.8%).
- **29.6% lower error** than the best RNA + mutation model (23.7 to 35.5%).
- All 9 planned comparisons passed.
- 90% error bands contained 87.6% of outcomes.
- **Conservative shortlist:** 733 of 788 flags correct (93.0% precision).
- **Wider shortlist:** 3,081 of 4,409 flags correct (69.9% precision).

**Compact model** (16 inputs → 46 outputs, 85 held-out donors): error 0.425 versus 0.565 for the strongest 32-input panel (24.8% lower) and 0.851 for the best RNA + mutation model (50.0% lower). Conservative shortlist: 156 of 161 flags correct.

### Where it fails

- **It misses most dependencies.** Of 5,881 true dependencies, the conservative shortlist catches 12.5% and the wider one 52.4%.
- **Some genes remain poorly predicted.** Eight of the 195 genes are predicted worse than their training average. Broad and compact errors use different target sets and confirmation cohorts, so 0.725 versus 0.425 is not a controlled comparison of coverage.
- **Outputs were chosen for being predictable.** The result says nothing about the 613 candidates that did not qualify.
- **Many input–output pairs are known partners**, such as TSC1 and TSC2, or genes in the same mitochondrial complex. PROBE uses this shared structure; it does not discover new biology.
- **No lab saving is shown.** The inputs are archived values from genome-wide pooled screens, not a cheaper assay.
- **Cancer cell lines only.** The held-out donors are new, but still from DepMap. No patient data or clinical benefit.

## What we would do next

1. Measure the 32-gene panel in a small arrayed knockout experiment and test whether those values work as inputs.
2. Repeat the test on cell models from a source outside DepMap.

## Optional Python tools

**Python** (3.10+). Run from the repository root. Install Node.js to run browser-parity checks; pytest skips those checks when Node.js is absent:

```sh
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -e ".[test]"
python -m probe --model models/broad/model.json --input models/broad/example_inputs.csv --output predictions.csv
python -m probe.evaluate --scope broad   # recompute the held-out results above
python -m probe.verify --scope broad     # check the app and Python give identical predictions
python -m probe.verify --scope compact
python -m pytest -q
```

To make a blank CSV, run `python -m probe --model models/broad/model.json --template --output inputs.csv`. Keep every gene ID exactly once; blank effects mean missing measurements. Genetic inputs are Chronos effects, not RNA levels or drug viability values.

Keep each `support.json` beside its `model.json`. Without the training ranges, inference still returns point estimates but reports `unverified_range` and omits bands. For edited scenarios through Python, call `Model.predict(values, measured=False)`; the CLI treats CSV values as measured inputs.

| Path | Contents |
|---|---|
| `app/` | The browser tool |
| `probe/` | Python model, panel selection, evaluation and verification |
| `models/broad/`, `models/compact/` | Frozen models, training input ranges, example input CSVs |
| `results/broad/`, `results/compact/` | Recorded held-out test: every comparator, every gene, exposure logs |

## Credits

**Data**
- DepMap 24Q2 Public, Broad Institute (2024). [doi:10.25452/figshare.plus.25880521.v1](https://doi.org/10.25452/figshare.plus.25880521.v1). CC BY 4.0. Derived files in `app/`, `models/` and `results/` are subsets and model outputs of this release.

**Scientific references**
- Dempster et al. (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biology*. [doi:10.1186/s13059-021-02540-7](https://doi.org/10.1186/s13059-021-02540-7)
- Zhao et al. (2022). A pan-CRISPR analysis of mammalian cell specificity identifies ultra-compact sgRNA subsets for genome-scale experiments. *Nature Communications*. [doi:10.1038/s41467-022-28045-w](https://doi.org/10.1038/s41467-022-28045-w)
- Hossain et al. (2025). goloco: a web application to create genome scale information from surprisingly small experiments. *BMC Bioinformatics*. [doi:10.1186/s12859-025-06070-y](https://doi.org/10.1186/s12859-025-06070-y)

**AI tools**
- OpenAI Codex: research, code, review.
- Anthropic Claude: research, code, review.

**Software**
- This repository: Python and NumPy; pytest for tests.
- Training also used SciPy, pandas and scikit-learn.
- JavaScript (no libraries) for the app; Node.js for verification.
