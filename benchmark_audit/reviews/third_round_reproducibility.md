# Fresh third-round reproducibility review

This review was performed after completion of the revised manuscript. I read the main and supplementary LaTeX sources, data and source provenance, `run_all.py`, the pinned requirements, package documentation, model/construction scripts, and figure/table generators. I did not rely on the previous 19,726-assertion verdict to conclude that fitting or source execution reproduced. My fresh checks and their limits are distinguished below. No external publication or repository push was performed by this reviewer.

## Verdict

The numerical construction audit is reproducible in the tested environment. I found no discrepancy between the independently regenerated subset and its archived construction rows, fold diagnostics, fitted model probabilities, or metrics. All generated table sources and PNG figure assets also reproduce exactly. The full experimental grid was **not** retrained in this review; the statistical referee independently checked all archived score-based metrics.

The review identified practical packaging, witness-provenance, executable-positive-control, and output-default gaps. These have been fixed and independently checked; the final documentation supplies an explicit editor-mediated source-access route and distinguishes representative antigen IDs from actual donor witnesses. **No unresolved scientific or executable-code blocker remains within the stated construction-audit scope.** Final archive checksums and external repository publication must reflect the actual final files. This verdict does not certify acceptance or claim a full-grid retraining in this round.

## Independently executed evidence

The repeatable checker is `reviews/third_round_reproducibility_check.py`; its results are in `reviews/third_round_reproducibility_checks.json`. It imports neither `verify_archived_results.py` nor the assignment-witness exporter.

| Check | Independently observed result |
|---|---|
| Frozen input ZIP and its TSV member | Exact declared SHA-256 values match |
| Reference checkout | Correct commit `0ba2d538065e2d3a8b4642547eb77447ecc89a5b`; checked source files equal their committed Git blobs |
| Archive source without `.git` | A fresh four-file source copy plus `SOURCE_PROVENANCE.json` runs successfully |
| Source tamper detection | A harmless appended comment in `dataset.py` is rejected by the archive fingerprint validator |
| Fresh synthetic source audit | Seeds 0 and 1 reproduce all four saved rows exactly; source counterexamples execute successfully |
| Fresh split-first construction and model fitting | Fold seed 17 and sampling seed 0, all four samplers and all four feature sets |
| Exact replay of that fitted subset | 21,336 probability rows, 16 pooled metric sets, 80 fold diagnostic rows, 20 construction metric rows, five retention rows |
| Construction manifests for fitted subset | Four files byte-identical to the archive |
| Figure/table regeneration | Five LaTeX tables and six PNG figures byte-identical after generation in a fresh temporary directory |
| Actual assignment witnesses | 58 constructions and 33,839 assignments independently reconstructed from raw TSV rows; closure, no fixed points, same HLA, registry exclusion, and exact equality with archived negative antigen values all pass |

The saved predictions were read with `float_precision="round_trip"`, preserving serialized floating-point values. The statistical referee separately observed that default CSV parsing can change tied-score AUROCs at about the sixth decimal place. That parsing effect must not be mistaken for a model discrepancy.

`reviews/third_round_archive_portability.json` records a separate execution of the **actual** `run_all.py --verify-only` entrypoint in a copied archive with spaces in its directory name, while the current directory was outside the archive. It completed without an extracted TSV, Git history, local pair CSVs, or downloads. This confirms entrypoint portability; its 19,726 assertions remain an archived-output consistency check, not evidence of fresh training of all models.

`reviews/third_round_packaging_check.json` records an adversarial packaging test. I planted stale local-pair, extracted-TSV, and obsolete-code sentinels in an isolated prior stage. The corrected packager removed all three by creating a fresh staging tree. All 694 manifest-listed files in that test matched their SHA-256 values. This count is a snapshot, not a required final archive count.

## Findings, required fixes, and their disposition

1. **Stale staging files could enter the delivered ZIP — fixed and tested.** The original packager copied into an existing tree and then archived every remaining file. A user could reproduce there and subsequently redistribute locally generated sequence tables or deleted obsolete files. The replacement stages into a new temporary tree and replaces only its task-owned destination. The sentinel experiment confirms that stale excluded content is removed.

2. **Representative antigen rows did not certify actual donor-row closure — fixed computationally and documented.** Existing `negative_antigen_source_row_0based` fields identify a representative row having the selected peptide/HLA value. They can refer outside a receptor fold even when the antigen value was sampled entirely within the fold. Therefore these fields reconstruct values but do not alone prove an assignment permutation on retained rows. The new `assignment_witnesses` files explicitly preserve `actual_donor_source_row_0based`, with cycle annotations. I independently reconstructed all 58 witnesses and verified their values against the existing manifests. The final supplementary methods and reproduction README explicitly distinguish these representative IDs from actual donor occurrences.

3. **A supplementary synthetic positive-control claim lacked an executable archived fixture — fixed, executed, and integrated.** I added `crossfit_audit/validate_synthetic_control.py`, which passes entirely artificial category strings through the actual `run_crossfit.run_dataset` function. Its 100 receptor categories supply one positive and one negative each; peptide category perfectly identifies the label and HLA is constant. Across two fold seeds and five folds, peptide and peptide-plus-HLA AUROCs are 1.0; HLA and receptor AUROCs and probabilities are exactly 0.5. The synthetic input, 400 score rows, eight pooled metric sets, 40 fold sets, and validation JSON are saved under `crossfit_audit/synthetic_control_results/`. Both this fixture and the assignment-witness exporter are registered in `run_all.py`; an isolated dispatch check confirmed the stage ordering. These data test an implementation and are not additional biological evidence.

4. **Determinism language needed to distinguish values from file bytes — fixed in the main reproduction README.** Scientific values and fixed manifests can be replayed under the declared software/data order. Gzip headers, PDF dates, and runtime provenance need not have identical bytes on a later build. `SHA256SUMS.txt` verifies the delivered snapshot. The revised README states this explicitly; do not characterize it as a universal bit-for-bit environment guarantee.

5. **Executable source access must survive submission handling — route now specified.** A fingerprint in a separate inaccessible author archive is not enough for reviewers to reproduce an exact-AST audit. The complete source bundle with preserved original copyright/license notices is acceptable evidence. The final main manuscript, supplement, and reproduction README specify an editor-mediated source/evidence route and identify the archive as nonanonymous. No copyright notices were removed. Actual editor access occurs at submission, and the public repository URL/revision should be recorded after publication succeeds; no such external action is assumed by this review.

6. **Annotation-wrapper standalone default could overwrite primary results — fixed and tested.** `run_all.py` already supplied a safe explicit annotation output directory. The standalone wrapper now defaults to `annotation_filtered_results` as well. Isolated argument-dispatch checks passed for default, explicit-space, explicit-equals, and abbreviated output arguments; the primary directory was not created. No scientific logic changed and historical execution metadata was preserved byte-for-byte. `reviews/release_provenance.md` records historical/current hashes and explains why they legitimately differ.

7. **GitHub transport must preserve the full probability archive — tested losslessly.** Thirteen numbered parts reconstruct the three compressed prediction files, totaling 119,181,206 bytes. I tested the final helper against all original complete-file SHA-256 values and strict existing-snapshot verification; all passed. A deliberately corrupted part was rejected without leaving a partial target. The default preserves already-regenerated outputs, while `--verify-existing` rejects bytes that differ from the delivered snapshot. This resolves the initially identified conflict between new gzip timestamps and a second verification run. Test results and the final helper hash are in `reviews/third_round_large_artifact_transport.json`. Reconstructed whole files were removed from the transport staging folder afterward, leaving all parts and the manifest intact.

## Leakage and interpretation audit

The actual model routine fits both `OneHotEncoder` and logistic regression only on `folds != held_out_fold`, transforms held-out categories using that fitted encoder, and checks exact receptor-group separation. The fresh replay reproduces the archived probabilities exactly, providing stronger support than a metrics-only check. No test-label parameter tuning, threshold search, or score inversion was found.

The *construction* uses the known-positive registry and positive records within each held-out fold to select a feasible retained population and generate decoys. It is therefore a retrospective construction diagnostic conditional on known records. This is compatible with the explicitly bounded manuscript question. It cannot be reinterpreted as prospective unseen-peptide prediction or independent biological nonbinding validation. An antigen lookup ID outside a fold is not itself leakage when the referenced antigen value equals an in-fold candidate; the true row witnesses now make the distinction inspectable.

Exact receptor-string grouping prevents exact-string overlap; it does not establish separation of related receptors, donors, studies, or peptide families. The manuscript does not claim those stronger generalizations. The globally constructed protocol couples complementary fold imbalances by design; preserving and explaining its below-chance values is more informative than treating them as biological reverse prediction.

## Limits of this review

This is a fresh technical reviewer-role audit, not journal peer review or an independent biological experiment. I did not recreate the missing original neural study, certify source annotations, retrain the complete experiment grid, or test another OS/BLAS/TeX installation. Five scientific Python packages are pinned; this is not a complete transitive dependency or operating-system lockfile. Final PDF rendering and the external GitHub publication are separate release tasks. The definitive tested script/manuscript hashes are recorded in the result JSON, so later edits can be identified without pretending that an earlier check verified unexamined bytes.
