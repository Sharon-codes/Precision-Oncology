# Release changes after the archived scientific runs

The final release distinguishes the source used to execute the stored experiments from later changes to command-line handling, packaging, and additional verification. Existing execution metadata has not been rewritten to imply that all models were retrained after these changes.

## Annotation-wrapper output default

`crossfit_audit/run_annotation_filtered.py` originally inherited the standalone primary experiment's `split_first_results` output directory if invoked without `--out`. The full reproduction pipeline had always supplied `--out crossfit_audit/annotation_filtered_results`, so this default did not affect the stored annotation sensitivity or any reported scientific value.

The release wrapper now resolves its own `--out` argument first, defaults to `annotation_filtered_results`, and delegates the resolved destination to the unchanged shared experiment. Explicit `--out PATH`, `--out=PATH`, and unambiguous argparse abbreviations remain supported. No eligibility, matching, fold construction, model, seed, or score computation changed.

| Record | SHA-256 |
|---|---|
| Historical wrapper recorded when annotation models executed | `0e9a8273431e0bbf8873b153db4f407eb43f9495ee5d83029746a9e9bbc5ed10` |
| Wrapper after safe-default fix | `31f65137f4acebb2eafdd32febf10ee0d1fb1720cc90c446abf29fb91670d833` |
| Unchanged archived `annotation_filtered_results/run_metadata.json` | `b84d592da099c54692e07476a3f6bcccc0617d5169b6717c8eeb49af9d74e4dc` |
| Unchanged shared `run_crossfit.py` | `5bdbc4c805f570a3879f9d226bba6be2b64ad8b83ad5a41b6a51462311f84d89` |
| Unchanged shared `run_split_first.py` | `852d9004409a5c6366fda02776adf681681f01bec19884c4b2ea40149cb3df51` |

The historical wrapper can be reconstructed exactly from the release wrapper by removing `import argparse`, removing the eight-line argument-resolution block preceding `experiment.prepare = annotation_filtered_prepare`, and replacing that removed output-variable assignment with these two original lines immediately after `experiment.main()`:

```python
    args = sys.argv
    output = Path(args[args.index("--out")+1]) if "--out" in args else experiment.HERE/"split_first_results"
```

The comments preceding that block are part of the added block and must also be removed to obtain the historical byte hash. The historical metadata fingerprint intentionally differs from the current wrapper fingerprint. A future full regeneration will correctly record the current wrapper fingerprint in its newly produced metadata.

An isolated argument-dispatch check used a stub shared experiment and temporary directories, avoiding all archived outputs. Default, explicit-space, explicit-equals, and abbreviated output arguments all selected the intended destination; the primary default directory was never created. The archived annotation metadata hash remained unchanged. This verifies CLI handling, not a fresh full annotation-model fit.

## Newly executed verification artifacts

The executable synthetic positive-control fixture was added under `crossfit_audit/validate_synthetic_control.py`; its entirely artificial input, predictions, metrics, and validation record are in `crossfit_audit/synthetic_control_results/`. It exercises the actual shared model routine, confirms AUROC 1.0 for perfect peptide signal and 0.5 for constant-HLA/receptor controls, and adds no biological observations.

The assignment-witness exporter was added under `audit_work/export_assignment_witnesses.py`. It passively captures actual donor occurrences from the unchanged assignment function and compares their reconstructed decoy values with archived representative-antigen manifests. It does not replace the original scientific results. The fresh reproducibility reviewer independently verified all 58 witness files and 33,839 assignments directly from frozen source rows.

`run_all.py` now includes `metadata_positive_control` and `assignment_witnesses` stages. The witness stage follows the source, confidence, distance, primary split-first, and annotation experiments whose manifests it checks. An isolated entrypoint-dispatch test confirmed both stage names and their ordering; both target scripts had already been executed successfully on their real inputs. No claim is made that the full grid was rerun merely to test stage registration.

The fresh-stage packaging fix prevents obsolete files or locally regenerated sequence materializations from surviving into rebuilt archives. Its adversarial sentinel test is recorded separately in `third_round_packaging_check.json`. Scientific output reproducibility and delivered-file byte fingerprints remain distinct: later timestamps and runtime metadata can change bytes without changing scientific values.

## Lossless GitHub transport

The GitHub-specific numbered parts are a storage representation of the three original compressed score archives. The final `restore_large_artifacts.py` helper was tested against all three complete-file hashes: 119,181,206 bytes reconstructed exactly. A corrupted part is rejected before installation and leaves no partial output. Existing outputs are preserved by default so a subsequent numerical verification can inspect a locally regenerated run; `--verify-existing` additionally enforces the exact delivered snapshot. The final helper SHA-256 is `f6363e6d4abaacbb1f84880529078badee85b16dceb68df00ce6be0f4e25d7ed`.

The final runner calls restoration before either archived-result verification or full rebuilding. The complete downloadable ZIP contains the ordinary whole files and does not need a part manifest. `third_round_large_artifact_transport.json` records final-helper tests and manifest/runner fingerprints. Whole reconstructed files were removed from the GitHub transport staging folder after the tests; the thirteen parts and manifest remain intact.
