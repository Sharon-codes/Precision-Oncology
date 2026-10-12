# Final independent execution review

Reviewed `paper/manuscript.tex` against the executable sampler, archived public/synthetic/reference outputs, matching validation, cross-fitting code, and all archived out-of-fold probabilities. I did not edit the manuscript.

**Execution verdict: PASS for the stated construction-audit claims.** No blocking numerical mismatch or misleading claim of neural-model execution was found. This is not a guarantee of scientific novelty, journal acceptance, biological label validity, or correctness of the former manuscript's unavailable learned-model results.

## Verification performed

The standalone script `audit_work/verify_archived_results.py` completed **19,726 numeric and certificate assertions**. Output: `audit_work/archive_verification_results.json`.

It independently:

- Reconstructed every primary and annotation-filtered public aggregate from 722 per-seed rows (468 aggregate records).
- Checked synthetic summaries, all 30 rejection-reference runs, and the 166-fixture/924,960-permutation matching-validation totals.
- Recomputed AUROC, average precision and Brier score from every archived OOF score vector, without retraining.
- Recomputed all 9,300 fold-level metric sets and all 1,860 pooled metric sets across three protocol/cohort analyses.
- Checked 2,473,986 stored score rows: valid probability ranges, unique example identifiers, complete fold coverage, disjoint receptor hashes across folds, equal positive/negative multiplicity for every receptor, and identical positive populations across the compared samplers and sampler seeds within each fold assignment.
- Recomputed fold-averaged seed summaries and manuscript rounding for the main trained controls.

The verifier reads floating-point predictions with round-trip precision to preserve ties. It needs NumPy, pandas and scikit-learn but no raw biological sequences, source repository import, networking or retraining. Preserve the `audit_work/`, `public_audit/` and `crossfit_audit/` relative directory layout, then run:

```
python audit_work/verify_archived_results.py --root .
```

This establishes internal consistency of archived predictions and declared fold certificates. Stored scores alone cannot prove how training happened; the separately inspected `run_crossfit.py` fits both `OneHotEncoder` and logistic regression strictly on `train = folds != fold`, predicts on its complement, treats convergence warnings as errors and leaves scores uninverted. `run_split_first.py` constructs negatives separately inside the preassigned receptor folds before passing those assignments to that training routine.

## Targeted main numeric checks

| Claim | Independently recovered value |
|---|---:|
| Primary release-order contradictions | 2.98868636% ± 0.10220160 percentage points |
| Primary peptide-sorted contradictions | 0.88204484% ± 0.07605975 percentage points |
| Primary shuffled contradictions | 0.22836790% ± 0.06950899 percentage points |
| Positive retention | 2,691 / 3,182 |
| HLA strata retained | 22 / 41 |
| Original-to-retained HLA TV | 0.1274123933 |
| Split-first peptide/HLA AUROC, sorted | 0.5558443585 |
| Split-first peptide/HLA AUROC, release order | 0.5168048952 |
| Split-first peptide/HLA AUROC, shuffled / matched | 0.500000 / 0.500000 |
| Global-first peptide/HLA AUROC, sorted | 0.7146393700 |
| Global-first peptide/HLA AUROC, shuffled | 0.3297992593 |
| Global-first peptide/HLA AUROC, matched | 0.3261593726 |

The main reference-rejection, confidence-reference, exact-expected-TV, matching-diversity and coverage-cost numbers agree with their saved evidence. The global-first below-chance results are retained, and the manuscript appropriately distinguishes an explanatory coupling mechanism from a controlled single-factor causal comparison.

## Source mechanism check

The exact `dataset.py:TCRCollate` class:

1. Stores the positive registry but never reads it.
2. For batch size greater than one, repeatedly shuffles indices until it obtains an index derangement.
3. Sequentially attempts swaps that avoid equal-peptide assignments and positional fixed points. If no valid local swap is found, it marks that row for fallback.
4. For a fallback row, makes at most 20 independent uniform choices from the supplied global pool, accepting the first peptide different from its original peptide. The public audit passes distinct peptide–HLA values to this pool, so this is category-weighted fallback rather than row-weighted reassignment.
5. Uses peptide `AAAAA` and HLA `A*02:01` if no accepted fallback is found.
6. Duplicates receptor embeddings for positive and negative blocks; it does not consult a learned model when choosing decoys.

The manuscript's source-execution explanation is faithful. Identity tags and the NumPy shim test selection and assembly, not ESM or Mamba inference; the draft makes that distinction explicitly. The fallback mechanism changes donor marginals, while positive-registry nonuse is a separate defect. Neither the public audit nor the revised draft claims the earlier study's missing learned-model outputs have been reproduced.

## Minor editorial precision

Two small clarifications identified in the earlier draft are now resolved in the final manuscript:

- The methods now explicitly restrict index derangement to batches larger than one and describe singleton fallback.
- The common-population table caption now identifies the entries as means, with standard deviations shown selectively and complete summaries archived.

No additional experiment is required to resolve an execution inconsistency uncovered in this review.

## Annotation-sensitivity extension

The final verifier also recomputes every prediction metric in `crossfit_audit/annotation_filtered_results`, every aggregate in `public_audit/annotation_sensitivity`, and checks the independent construction-only review JSON. The annotation-filtered population contains 3,108 primary positives; closed assignment retains 2,665 and 21/39 source strata. The corrected manifest exactly equals the original corrected manifest with the flagged source strata removed. Original raw labels are excluded without recoding.

- Annotation-filtered primary source contradiction percentages: 2.8882453882, 0.9330759331 and 0.2445302445 for release, sorted and shuffled loading.
- Corresponding peptide identity-oracle AUROCs: 0.6422945249, 0.7935729511 and 0.5.
- Independently recomputed split-first peptide/HLA classifier means: sorted 0.5618565011, release order 0.5187752368, shuffled 0.5, matching 0.5.
- The added OOF archive contains 814,804 prediction rows, 155 complete five-fold runs, 620 pooled metric sets and 3,100 fold-level sets; all checks passed.

The original `reviews/annotation_sensitivity.json` accurately describes its earlier construction-only task. The later public-sampler and metadata reruns have separate result directories and are included in the final verification; that earlier JSON's `not_rerun` field should be understood in its local task scope.

The verifier records SHA-256 of the actual manuscript and all three OOF score archives. Final verified manuscript hash at this run: `f104024c5059949e89adbb5c6f0bdab0068d09d4a27f058ad119c0dce3eef1f0`.
