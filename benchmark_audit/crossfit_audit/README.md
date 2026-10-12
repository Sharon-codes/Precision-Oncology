# Independently executed trained negative controls

This experiment asks whether peptide/HLA identity can predict the **constructed** positive/decoy label without receptor sequence input. It is a within-antigen-collection construction diagnostic with held-out receptor groups. It is **not** a peptide-OOD or biological binding prediction benchmark.

## Recommended primary result: split before negative construction

Use `split_first_results/summary.csv` and `split_first_results/oof_metrics.csv`.

We assigned exact CDR3-beta sequence groups to five folds, then constructed negatives separately within each fold. For each fold, the cycle-cover procedure selected the largest closed subset admitting a within-HLA, catalogued-positive-excluding antigen permutation. Every comparator sampler then used **that same retained positive cohort** in that fold. Training and test receptor groups were disjoint; antigen identities could recur. Five group-assignment seeds were 17, 41, 73, 107 and 211. There were ten negative-sampling seeds for each source-collator regime; correction was deterministic for a fixed fold/cohort.

For every test fold, a one-hot encoder was fitted only on the other four folds. Logistic regression used L2 regularization, C=1, `liblinear`, tolerance 1e-8, max_iter=2,000 and random_state=1729; no tuning used test labels. Unknown test categories were ignored. Separate controls used peptide+HLA, peptide only, HLA only and receptor identity only. All positive/negative pairs for a receptor stayed together. The receptor-only control had balanced examples for each training receptor and unknown test receptors, and returned 0.5 throughout.

**Peptide+HLA control:**

| Construction on each common fold population | Mean out-of-fold AUROC | SD across 10 decoy-seed means | Range over all fold-map/decoy runs |
|---|---:|---:|---:|
| Source order | 0.516805 | 0.000582 | 0.511448–0.522349 |
| Peptide-sorted rows | 0.555844 | 0.001104 | 0.543482–0.563341 |
| Row-shuffled | 0.500000 | 0.000000 | 0.500000–0.500000 |
| Closed within-HLA matching | 0.500000 | Not applicable | 0.500000–0.500000 |

The SD first averages the five fold-map runs for each decoy seed, then measures variability between those ten averages. It is **not** a confidence interval over biological populations. The all-run range also includes fold assignment variation. Prefer means with the all-run ranges if one simple summary is needed.

Retained positives for the five fold seeds were 2,667, 2,659, 2,652, 2,670 and 2,658 out of 3,182. These are 5,334, 5,318, 5,304, 5,340 and 5,316 total examples after 1:1 decoy construction. All methods share the same positive set within a fold seed; different fold seeds change the feasible retained cohort. Mean unseen test peptide-category fractions were 4.08% (matched), 4.88% (sorted), 4.34% (source order) and 4.08% (shuffled); these are reported, not silently discarded.

**Interpretation:** modest but reproducible peptide/HLA-only discrimination survives receptor-group holdout under peptide-sorted construction. Shuffling removed this marginal signal in the present matched population, as did exact matching. The result does not show that every shuffled sampler is unbiased, that the correction improves a biological model, or that real binding has been learned. Source-collator known-positive collisions can persist despite AUROC=0.5 for the metadata control.

## Mandatory sensitivity: construct globally before cross-fitting

The initial experiment `results/` froze negatives globally on the common 2,691-positive population, then made receptor folds. It was preserved unchanged after it revealed a protocol-placement effect.

| Construction | Mean peptide+HLA out-of-fold AUROC | Range over all runs |
|---|---:|---:|
| Source order | 0.510860 | 0.504216–0.518748 |
| Peptide-sorted | 0.714639 | 0.710338–0.718214 |
| Row-shuffled | 0.329799 | 0.321197–0.341083 |
| Closed matching | 0.326159 | 0.321143–0.335068 |

The below-chance results were neither hidden nor score-flipped. They arise because exact global class marginals couple the training and held-out class counts. If category c has global positive and negative counts equal, then its held-out count difference is the negative of its training difference: Delta_test(c) = -Delta_train(c). A fitted category-frequency model can therefore learn the opposite of the held-out imbalance. This finite-dataset dependence is a reason to construct negative examples independently within evaluation partitions. It is not biological reverse prediction and does not demonstrate model failure under ordinary IID sampling.

Differences between the two panels are descriptive protocol sensitivities: the foldwise construction changes feasible retention, batch composition and the placement of constraints. Do not attribute the entire AUROC change to one cause. `metadata_crossfit_protocol.pdf` displays all results.

## Registry and population restrictions

The catalogued-positive exclusion registry is the fixed collection of eligible VDJdb records at all confidence scores, including score zero. Main positive examples require score>=1. The registry is used for constructing and auditing presumed-negative triples, not fitted as a predictive feature. The original collator receives the same registry but does not consult it. The correction enforces exclusion only at the exact beta/peptide/HLA level; unknown binding and omitted biological context remain. Registry exclusion does not certify a negative assay result.

This is a retrospective public-database construction audit, not a prospective deployment simulation. The correction's scope includes HLA matching and explicit omission of infeasible records; those exclusions can change the target population. It is not a universally unbiased sampler. Parent public-audit outputs report HLA retention and marginal population shifts.

## Independent checks

- `independent_matching_validation.json`: 120 random directed graphs of sizes 2–7 checked against exhaustive permutations; maximum retained closed subset agreed in every case.
- Independent checks on the public corrected CSV confirmed exact beta, peptide, HLA and joint peptide–HLA marginals and zero catalogued-positive triple collisions.
- Crossfit functions were tested on a synthetic dataset with perfect peptide-only signal: antigen controls gave AUROC=1, null HLA/receptor controls 0.5, with receptor separation asserted.
- Every fit asserts group disjointness and complete score coverage. All samplers in a cohort are checked for identical positive multisets. Convergence warnings are fatal.

## Files and reproduction

1. Run `../public_audit/run_public_audit.py` on the pinned VDJdb release to reconstruct common local pair CSVs.
2. Run `python run_crossfit.py --manifest input_manifest.json --output results` for the global-construction sensitivity.
3. Run `python run_split_first.py` for the recommended split-before-construction control.
4. Run `python plot_crossfit.py` to regenerate the two-panel figure.

Scores, labels and receptor hashes are saved in `oof_scores.csv.gz` in each results folder. Fold assignment, category coverage, means, raw metrics, fixed source-index negative manifests and algorithm settings are included. Raw source sequence materializations remain in local public-audit scratch; package source indices and reproducible scripts with the unchanged licensed source obtained separately. The figure is descriptive; no statistical significance test or biological confidence interval was manufactured.

## HLA-annotation sensitivity

Biological review subsequently flagged raw-source HLA-A*08:01 and HLA-A*24:01 labels. `run_annotation_filtered.py --out annotation_filtered_results` excludes both labels from all cohorts and the exclusion registry, without recoding, and reruns the complete split-first experiment with freshly assigned receptor folds. This is a two-label sensitivity, not comprehensive allele nomenclature validation.

The primary cohort becomes 3,108 positives with 39 HLA labels; the all-score registry becomes 29,146 triples. Retained positives across the five fold seeds are 2,618, 2,625, 2,636, 2,642 and 2,621. Peptide+HLA AUROC is 0.561857 for sorted construction (all-run range 0.553202–0.573754), 0.518775 for release order (0.513858–0.523180), and 0.500000 for both shuffled construction and closed matching. Results and complete OOF scores are saved separately in `annotation_filtered_results`; all original raw-annotation outputs remain available. Use `primary_peptide_hla_summary.csv` for the mean, all-run range and decoy-seed SD.
