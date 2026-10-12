# Fresh statistical referee review

Date: 12 October 2026. Reviewed the revised manuscript and Supplementary Methods, their rendered PDFs, the original sampler class, the source-audit shim, public-audit construction and rejection code, both cross-validation protocols, annotation filtering, distance partitions, figure-generation code, and underlying score/manifests. Earlier reviewer verdicts were not read before this assessment. This is a review of an exploratory sampler-construction audit, not a judgment about neural binding prediction or journal acceptance.

## Verdict

**PASS within the manuscript's stated construction-audit scope.** I found no fatal inference error, mathematical optimizer error, undisclosed test-feature fitting, or headline numerical failure. All four conditions identified below are now closed: the final PDFs explain the two donor-identifier types, use the correct rounding, show Figure 5b without overlapping text, and the synthetic positive-control claim has an executed actual-pipeline fixture. This is a scoped statistical and reproducibility verdict, not certification of biological nonbinding, neural prediction superiority, or journal acceptance.

## Concrete issues found and disposition

1. **Actual donor rows versus antigen-value lookup rows, resolved.** The original `negative_antigen_source_row_0based` manifest field selects the first eligible source record having the chosen peptide/HLA value. It is therefore an antigen-value reconstruction pointer, not necessarily the row occurrence selected by the assignment. A pointer outside a receptor fold does not establish candidate-value leakage, but those original identifiers cannot by themselves certify donor-row closure. The parent requested and I implemented a stronger correction: `audit_work/export_assignment_witnesses.py` passively captures the final arrays from the unmodified assignment implementation. `public_audit/assignment_witnesses/` now contains **58 actual donor-cycle witnesses, covering 33,839 retained assignments**, plus checksums and verification results. All actual receiver/donor multisets close; HLA agrees; every negative is outside its declared registry; all declared marginals agree; every negative tuple exactly matches its existing representative-value manifest. The final supplement's rendered page 2 explicitly distinguishes both manifest types and documents all 58 witnesses.

2. **Rounding mismatch, resolved.** The complete-primary-cohort sorted peptide-oracle mean is **0.7916499788051978**, which rounds directly to **0.7916**, not 0.7917. Supplementary Table S4 was already correct. I verified that the final rebuilt main PDF now reports 0.7916 as well.

3. **Figure 5b text overlap, resolved.** In the initial rendered main-page 11, the `0.715` label overlapped the second line of the protocol-placement annotation. The final figure on rendered main-page 12 separates the annotation above the label. I visually inspected that page and confirmed the overlap is gone.

4. **Positive-control fixture, resolved during review.** The supplement claims a separately executed perfect-peptide-signal check with AUROC 1.0. The reproduction referee could not initially find its executable fixture or saved result. They subsequently executed `crossfit_audit/validate_synthetic_control.py` through the actual `run_dataset` path and saved its input, 400 OOF rows, eight pooled metric sets, 40 fold metric sets, and validation report under `crossfit_audit/synthetic_control_results`. I inspected this fixture and report: artificial peptide/joint features yield AUROC 1.0 in both complete runs and all ten folds; constant HLA and receptor features yield exactly 0.5. The supplementary sentence is now supported. This fixture tests implementation, not biology.

## Independent numerical checks

- Read all three released compressed OOF score archives: original global construction, primary split-first, and annotation-filtered split-first. These contain **2,473,986 score rows in 465 complete runs**, yielding **1,860 feature-level evaluations**. Recomputed AUROC, average precision and Brier score for every evaluation. All three metrics agree **exactly** with every saved result when read with `pandas.read_csv(..., float_precision='round_trip')`.
- An initial default-parser read produced tie-sensitive differences no larger than 1.26e-6. These vanished completely with round-trip floating-point parsing. This is a CSV parsing/near-tie issue, not a model discrepancy. Reproduction instructions should use round-trip parsing for exact comparisons.
- Independently checked AUROC with the rank-sum U-statistic, including tied ranks. Checked unique example identifiers within each complete run, one fold per receptor hash, and equal positive/negative multiplicity for each receptor. All passed.
- Primary split-first peptide/HLA AUROC: sorted **0.5558443585**, release **0.5168048952**, shuffled **0.5**, assignment **0.5**. Sorted all-run range: **0.5434817258–0.5633408591**. Its sampler-seed SD after averaging fold-seed runs is **0.001104**, consistent with the manuscript's 0.0011.
- Original global-construction peptide/HLA AUROC: sorted **0.71463937**, shuffled **0.32979926**, assignment **0.32615937**. The below-chance outputs were not silently reversed.
- Annotation-filtered split-first peptide/HLA AUROC: sorted **0.5618565011**, release **0.5187752368**, shuffled and assignment **0.5**. These agree with the supplement.
- Independently restated the raw-source eligibility masks and confirmed **35,237 eligible source rows**, **3,182 primary tuples**, **2,691 retained rows**, **22 retained HLA strata**, **236 retained peptides**, **2,624 unique negative tuples**, and original-versus-retained HLA TV **0.1274123933**.
- Reran the source sampler for primary seed 0 under all three orders. Registry contradiction counts and peptide TV exactly match the stored results: release 93 and 0.1285355123; sorted 29 and 0.4352608422; shuffled 9 and 0. Source fallback counts are 428, 1,498 and 0, respectively.
- Passively captured and independently checked the actual solver-selected row occurrences in all 58 exported constructions, then reopened every saved witness, verified its checksum, and traversed its recorded cycles. All passed.
- Independently tested the assignment objective on **84 new random directed-graph fixtures**, sizes 1–7, with **70,956 exhaustive permutations**. Every minimum-cost solution attained the exhaustive maximum retained closed subset and avoided forbidden edges. This is additional evidence independent of the manuscript's already reported 166 fixtures; it does not replace the proof or establish biological validity.

Supporting outputs from these checks are in `reviews/tmp_statistics/`: `oof_roundtrip_recalculation.json`, `primary_reconstruction.json`, and `independent_assignment_proof_check.json`. The round-trip result is authoritative for exact metric replay; the earlier default-parser file is retained as a diagnostic.

## Why the inference is defensible within scope

The optimizer proof is correct: the identity permutation costs n; every forbidden edge alone costs n+1; every feasible optimum consequently contains only admissible off-diagonals and diagonal abstentions. Its cost equals the number discarded. The nonfixed indices form complete cycles, proving identical retained receiver/donor row sets. This optimizes the expressly declared closed-row assignment contract, not all possible weighted, replacement-based, or biologically valid negative designs. The manuscript makes this distinction.

The peptide identity oracle is explicitly labeled a same-dataset resubstitution diagnostic. Assignment AUROC 0.5 is openly treated as a construction identity, and shuffled sampling's identical result is reported. The trained sorted condition supplies an independent, modest label-learnability result; the paper does not turn structural chance performance into a claimed prediction victory.

The comparison with independent within-HLA rejection has identical retained positives and the same registry. Its nonuniform accepted donor probabilities are derived correctly. TV of the expected negative distribution is distinguished from mean realized TV. Source-versus-assignment comparisons change several constraints; the manuscript acknowledges that they do not isolate coordination alone.

Cross-validation fits encoders and logistic regressions only on training folds, keeps exact receptor identities disjoint, retains unseen antigen categories, and does not invert scores or tune on test results. The global registry is explicitly a frozen retrospective sanitation reference, rather than a deployable predictive input. The global-construction coupling identity is correct. Because protocol placement also changes retained populations and batches, the manuscript appropriately avoids a single-factor causal interpretation of the total AUROC difference.

Uncertainty is narrowly and accurately described as conditional seed variation on frozen data. Fold-seed runs are collapsed before sampler-seed SDs; full-run ranges expose additional partition variation. These are not biological confidence intervals or independent biological replicates. More sampler seeds would not resolve database ascertainment, omitted receptor context, or population-selection limits, and the paper says so.

The original-versus-retained population change is substantial but disclosed quantitatively rather than hidden by common-cohort comparison. Registry absence is never certified as nonbinding. These restrictions must remain intact in the abstract, discussion, conclusions, title-page material and any public description of the study.

## Final artifact verification

All review conditions were closed after reading the rebuilt PDF text and inspecting the final Figure 5 page and supplementary witness explanation. The reviewed final PDFs have SHA-256 checksums:

- `paper/manuscript.pdf`: `e08f1ac6272460037759ee7e2f413b05b4932b3f027ee8f3e11ffb6f82b5aed2`
- `paper/supplementary.pdf`: `9c3472a40c3b6c56315daacda7ec8ef64760d4f392440eb567e72cde54b267c7`

No statistical blocker remains for the stated exploratory construction-audit conclusions.
