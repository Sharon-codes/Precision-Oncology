# TCR benchmark audit — full author reproducibility archive

The main manuscript is `paper/manuscript.pdf`; the editable source is `paper/manuscript.tex`. The work is a computational audit of synthetic-negative construction, not validation of an unavailable neural binding predictor. Start with `paper/author_only/review_report.pdf` for the adversarial review and evidence ledger.

This is an **author archive**, not an anonymized submission supplement. It preserves identifying source licenses and provenance. The manuscript itself omits author identities, and the separate title-page template identifies declarations that only the authors can confirm. Provide this complete source-evidence bundle to the editor for designated reviewer access with identifying copyright and provenance intact. No external submission or correspondence has occurred.

## Contents

| Folder | Evidence |
|---|---|
| `paper/` | LaTeX, compiled main and supplementary PDFs, generated tables, five main figures and one supplementary figure, author report |
| `audit_work/` | Executed original-sampler tests, synthetic experiment, exhaustive matching tests, rejection reference and independent numerical verifier |
| `public_audit/` | Cohort construction, source-order/batch/confidence experiments, distance partitions, annotation sensitivity, aggregate outputs and source-row manifests |
| `crossfit_audit/` | Three independently executed receptor-held-out metadata protocols, complete metrics, frozen manifests and compressed out-of-fold probabilities |
| `research/` | Untouched public data archive and separate geometric counterexample used to audit the original attention claim |
| `reviews/` | Reviewer-role reports and debates; original manuscript audit; verified bibliography; annotation evidence and validators |
| `precision_repo/` | Unmodified reference source files and textual outputs from the original declared Git commit, with MIT license |

## Data and source provenance

- VDJdb release: https://github.com/antigenomics/vdjdb-db/releases/tag/2023-06-01
- Original ZIP SHA-256: `7c4994ce32eba8e0ff6c23d1ec53b85a982be69651ed1f467ffe64c86c66c1d1`
- Extracted `vdjdb.txt` SHA-256: `2823281ece02f37b15cdd9187be374e99ee7b46e041572fc2cde4b295a9bf52e`
- Source sampler: https://github.com/Sharon-codes/Precision-Oncology
- Exact commit: `0ba2d538065e2d3a8b4642547eb77447ecc89a5b`
- Frozen official HLA nomenclature: https://github.com/ANHIG/IMGTHLA/tree/38398a75e9762ff070d8e9bd714d074332646cd7 (file headers version 3.52.0, 17 April 2023).

The original VDJdb archive is included unchanged with its license. Modified sequence CSVs are excluded from distribution and regenerated locally. Released derived constructions use source-row identifiers. Official nomenclature files remain unchanged. Original source-code copyright notices are preserved. New analysis scripts are supplied as source for inspection and reproduction; no new blanket license over third-party materials is asserted.

## Reproduce

Python 3.12.14 was used. Install the exact dependencies in an isolated environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements_reproducibility.txt
python run_all.py --verify-only
```

`--verify-only` reads the archived numerical outputs and recomputes all pooled and fold-level metrics from stored probabilities. It does not retrain models or certify biology. The delivered execution passed 19,726 assertions across 2,473,986 prediction rows, 1,860 pooled and 9,300 fold-level metric sets. The accompanying JSON states the exact checks.

For a full regeneration of the deterministic scientific outputs from the untouched data archive:

```bash
python run_all.py
```

This regenerates local sequence materializations, every reported main experiment and sensitivity, model controls, source-row manifests, figures, and independent checks. It runs sequentially on CPU; duration depends on hardware. It overwrites generated analysis outputs in this extracted copy, so preserve the delivered ZIP for comparison. Scientific values and fixed manifests are reproducible under the pinned environment; gzip headers, PDF timestamps and provenance metadata need not be byte-identical. `SHA256SUMS.txt` fingerprints the delivered bytes, not every possible later rebuild. No GPU, PyTorch, ESM model, checkpoint, or network download is needed once dependencies are installed. Intermediate sequence files are local computational products and should not be redistributed as a modified copy of the source dataset.

Add `--compile-pdfs` when a suitable TeX Live and `latexmk` installation is available. Alternatively compile only the documents as described in `paper/README.md`.

## Interpretation constraints

1. Identity-oracle AUROCs are same-dataset category-frequency diagnostics, distinct from held-out classifier metrics.
2. Trained models use peptide/HLA metadata and held-out exact receptor groups within the antigen collection. They are not unseen-peptide binding models.
3. Negative-sampler seed variation is not a biological confidence interval. No test-label tuning or score inversion was used.
4. Known-positive exclusion means consistency with a stated beta/peptide/source-HLA registry. It does not prove nonbinding or validate a full alpha-beta receptor.
5. Closed matching preserves marginals only on a selected retained population, with substantial HLA/peptide coverage loss and possible duplicate negative triples.
6. Raw source-label syntax is distinguished from nomenclature validity. The two flagged labels were excluded in separate fully rerun robustness analyses.
7. The initial globally constructed cross-validation protocol, including below-chance results, is preserved as a sensitivity rather than suppressed.

`SHA256SUMS.txt` fingerprints the delivered files. The execution and source provenance reports distinguish newly executed results from unverified claims in the supplied original manuscript.

## Reconstruction identifiers and assignment witnesses

Antigen donor IDs in reconstruction manifests are the first source row carrying a required peptide/HLA value. They reconstruct tuples but need not identify the internally selected donor and do not certify closed cycles. `public_audit/assignment_witnesses/` supplies the actual solver-selected receiver/donor maps: 58 witnesses, 33,839 assignments, including all metadata folds and distance sensitivities. Each is independently checked for closure, registry exclusion and exact agreement with reconstructed tuples. The original audited sampler remains unmodified.

The synthetic fitting positive control is executable as `python crossfit_audit/validate_synthetic_control.py` and is included in the full runner. It archives 400 prediction rows from two receptor-fold seeds: peptide and peptide/HLA AUROC 1.0; HLA and receptor-only AUROC 0.5.

## Large files in the GitHub copy

The downloadable complete research ZIP contains ordinary complete compressed prediction files. The GitHub copy stores these files as lossless numbered parts to fit the upload transport limit. `python run_all.py --verify-only` automatically restores them before checking results; alternatively run `python restore_large_artifacts.py`. Whole-file and per-part SHA-256 values must match before a reconstructed file is installed. This changes storage only, not a single prediction byte. After restoration, the original `SHA256SUMS.txt` applies. `GITHUB_SHA256SUMS.txt` describes the distributed GitHub files before restoration.

Restoration preserves already-existing outputs, which the numerical verifier then checks. Use `python restore_large_artifacts.py --verify-existing` only when checking the exact delivered GitHub snapshot; later regenerated gzip headers may differ without changing predictions.
