# Final methods/statistics review of revised manuscript

Reviewed `paper/manuscript.tex`, the executed public-audit code/outputs, the rejection-reference results, and both trained-control pipelines. This review concerns the replacement construction-audit paper, not the original unsupported biological-superiority manuscript.

## Verdict

**The revised scientific story is internally defensible after the specific corrections below.** I found no remaining blocking mathematical or statistical defect in its bounded claims. The manuscript now distinguishes construction invariants, empirical source behavior, trained metadata controls and biological prediction. It preserves negative/less favorable results and does not present simulations as clinical evidence.

The work remains a focused audit of one research-code sampler and one database release. Journal acceptance and field-level impact are uncertain; the supportable contribution is the executable construction contract and measured feasibility/coverage tradeoff. I would not restore the original predictor-superiority, rank-collapse or clinical claims.

## Required corrections

1. **HLA validation language (important).** Biology review identified raw-release HLA-A*08:01 and the deleted HLA-A*24:01 label. The eligibility test is a regular expression, not a nomenclature-dictionary check. Replace “unambiguous two-field HLA-A/B/C strings” with “source labels matching two-field HLA-A/B/C syntax.” Explicitly distinguish syntactic eligibility from verified allele identity. The sensitivity excluding both flagged labels removes 74 primary rows, leaving 3,108 positives and 39 HLA source-label strata; global closed assignment retains 2,665. Do not silently replace original main results. Deleting rows before receptor-fold assignment changes the foldwise analyses, so I reran the complete five-fold-map/ten-decoy-seed study rather than assuming invariance; results below retain the same qualitative conclusions.
2. **Multiple changes in source-versus-correction comparison.** The source sampler moves peptide and HLA together and can change the receiver's HLA; both rejection and assignment hold HLA fixed. Add one explicit sentence in Methods: “Source-to-reference comparisons change the HLA and exclusion constraints as well as donor coordination; the rejection-versus-assignment comparison shares HLA restriction and registry exclusion.” This prevents interpreting the full source-to-assignment difference as the isolated effect of coordinated assignment.
3. **Table common-population caption.** It currently promises means plus SDs but several source TV/oracle cells show means alone. Either say SDs are supplied where shown and all SDs are in the archive, or add the missing values. Verified source values: release peptide-TV SD=0.000931, oracle SD=0.000865; sorted TV SD=0.002691, oracle SD=0.001901; shuffled TV SD=0.000068, oracle SD=0.000066. Rounding the tiny shuffled mean to zero would lose a real distinction; keep the current precision or scientific notation.
4. **Unknown biological status.** Prefer “generated pairs absent from K have unknown binding status” over “remain unverified nonbinders.” The latter still sounds like nonbinding is the default biological conclusion. The manuscript otherwise handles this correctly.

## Specific contracts independently checked

### Receptor–HLA preservation is valid

Within each HLA stratum, a retained receiver keeps its original beta string and HLA; its negative receives a donor peptide from that same HLA. Each receiver appears once in both classes. Therefore the joint beta–HLA marginal is preserved, not merely the beta marginal. I independently checked the actual common-cohort positive/negative counters and reconstructed all five split-first matched manifests from untouched source rows. Each matched negative retained receiver HLA, and each fold had exactly equal joint peptide–HLA class marginals. Results are recorded in `crossfit_audit/final_contract_validation.json`.

### Assignment maximum is correctly scoped

The identity assignment costs n, any forbidden edge alone costs n+1, and an admissible nonfixed edge costs zero. An optimal assignment therefore avoids forbidden edges and minimizes fixed points. Nonfixed cycles have equal source/donor row sets. The manuscript correctly limits “maximum” to the declared closed within-HLA row-assignment contract and distinguishes ordinary bipartite matching as an upper bound. Its mathematical proof is sound. The reported 166-fixture/924,960-permutation counts match the independent audit output; my separate 120 random-graph enumeration also passed.

### Chance AUROC statement is valid in the primary protocol

If class multisets of a feature are equal inside each test fold, any fixed per-fold deterministic score of that feature produces identical positive and negative score distributions in that fold. Because each fold contains equal class counts, pooling those fold-specific scores still yields equal class-conditional score distributions and AUROC 0.5. The joint peptide–HLA matching contract meets this requirement. The primary logistic controls also empirically returned constant 0.5 scores.

Global-only marginal equality is different: training/test category imbalances are complementary, and category scores trained separately can invert held-out imbalances. The manuscript correctly states this distinction, retains the 0.3262 and 0.3298 outcomes without reversal, and does not explain all differences between placement protocols through a single mechanism. No further correction is needed here.

### Full and retained populations are separated correctly

The 3,182-row full-cohort source results, the 2,691-row common-population comparison, and the 2,652–2,670 foldwise retained populations are clearly labeled and not interchanged. The maximum score-3 peptide frequency, 139-row matched retention and the 22-of-41 HLA coverage finding are supported by outputs. Preserve the caveat above concerning source-label nomenclature.

## Numerical/statistical checks

- Full primary contradiction rates, peptide TVs and identity-oracle AUROCs match `public_audit/aggregate_results.csv`.
- Table common-population values match the retained-cohort rows in that same file.
- Rejection TV 0.1348569305, SD 0.0058775595, oracle AUROC 0.5990066834, SD 0.0037711414 and expected-distribution TV 0.1019994442 match `reference_rejection_results.json`.
- Duplicate counts (2,624 matched unique negatives; 2,656.7 rejection mean) match the saved results.
- Primary metadata AUROC 0.555844, all-run range 0.543482–0.563341, source-order 0.516805, shuffled and matched 0.5 match the independently executed crossfit output.
- The reported SDs average fold-map results within each decoy seed before summarizing ten decoy-seed averages. They are appropriately described as conditional variability, not biological confidence intervals. Matching has no sampling-seed SD.
- The exclusion registry includes all eligible scores, whereas the main positive population uses score>=1. The manuscript states this clearly and avoids equating a registry contradiction with a true false-negative binding assay result.
- No label-dependent learner tuning, threshold selection or sign reversal was used. Group separation and encoder fitting are correct. Antigen overlap across receptor folds is disclosed; this is not mislabeled as an unseen-antigen benchmark.

## Agreement with other reviewers

The editorial reviewer and I agree the paper should explicitly identify a single custom research-code implementation rather than imply systematic prevalence across the literature. The biology reviewer and I agree on the HLA syntax/identity distinction, unknown-negative wording and the multifactor source/reference comparison. These changes strengthen the bounded contribution without adding unsupported claims.

## Final confidence

High confidence in the executed construction and metadata-control results, assignment proof, and stated diagnostic scope. Remaining uncertainty concerns generalization to other releases/samplers, biological utility of the selected negatives, and editorial novelty—not an unresolved calculation being concealed. No further simulation is required to support the manuscript's present limited claims after the listed edits.

## Executed annotation sensitivity added after nomenclature review

The same split-first control was rerun using `crossfit_audit/run_annotation_filtered.py`, excluding HLA-A*08:01 and HLA-A*24:01 from every positive cohort and the known-positive registry, without recoding either label. The all-score reference declined from 29,306 to 29,146 triples, and the primary cohort from 3,182 to 3,108. The score-3 cohort was unchanged. Receptor folds, within-fold feasible cohorts and all fitted controls were regenerated. Original outputs remain separate.

| Peptide+HLA control | Mean AUROC | Range over all 50 runs | SD across ten fold-map-averaged decoy seeds |
|---|---:|---:|---:|
| Release order | 0.518775 | 0.513858–0.523180 | 0.000457 |
| Peptide sorted | 0.561857 | 0.553202–0.573754 | 0.001069 |
| Row shuffled | 0.500000 | 0.500000–0.500000 | 0.000000 |
| Closed assignment | 0.500000 | 0.500000–0.500000 | Not applicable |

The five fold-map seeds retained 2,618, 2,625, 2,636, 2,642 and 2,621 positives, respectively. This sensitivity supports the qualitative trained-control conclusion while explicitly acknowledging annotation imperfections. It is not a comprehensive HLA nomenclature validation and does not justify silently calling all remaining source labels verified alleles. Exact score files, manifests, metadata, raw metrics and summary are in `crossfit_audit/annotation_filtered_results/`.
