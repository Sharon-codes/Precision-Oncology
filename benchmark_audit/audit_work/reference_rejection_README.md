# Within-HLA rejection reference on the common retained cohort

Run `python reference_rejection.py` from this directory after the public audit. It reconstructs exactly the same 2,691 retained positives from the corrected source-row manifest and the untouched VDJdb release. The same all-score registry supplies exclusion keys.

For **each receiver independently**, draw uniformly **with replacement over eligible donor rows** within its HLA. An eligible donor produces a beta–peptide–HLA triple absent from the known-positive registry. This is a newly defined transparent reference sampler, not an official published biological method. Every receiver is feasible because this cohort already admits a closed matching.

Thirty NumPy `default_rng` seeds (0–29) are used. Source-row manifests are saved separately for every seed. The identity oracle is a marginal resubstitution diagnostic, not an independently tested predictor.

| Diagnostic | Reference mean ± seed SD | Fixed closed matching |
|---|---:|---:|
| Peptide total variation | 0.134857 ± 0.005878 | 0 |
| Peptide identity-oracle AUROC | 0.599007 ± 0.003771 | 0.5 |
| Peptide–HLA total variation | 0.135947 ± 0.005443 | 0 |
| Peptide–HLA identity-oracle AUROC | 0.599817 ± 0.003653 | 0.5 |
| Known-positive registry collisions | 0 for every seed | 0 |
| HLA / beta total variation | 0 / 0 for every seed | 0 / 0 |
| Unique negative triples | 2656.7 ± 3.33 | 2624 |
| Duplicate negative rows | 34.3 ± 3.33 | 67 |

The negative dataset always contains 2,691 rows. Repeated receptor or antigen identities can produce duplicate output triplets even under a donor-row permutation. Thus matching does not guarantee triplet uniqueness, and this particular fixed matching has more duplicates than the reference.

If E_i is receiver i's eligible donor-row set, donor row j's expected sampling probability is

`q_j = (1/N) sum_i I(j in E_i) / |E_i|`.

These incoming probabilities need not equal 1/N. Exact expectation over this reference sampler gives peptide TV **0.1019994442**, demonstrating a systematic imbalance beyond finite-draw noise. For example, a single-HLA cohort with peptide frequencies (1/2,1/4,1/4), unique receptors and only own-peptide known positives yields rejection frequencies (1/3,1/3,1/3), although it admits a full marginal-preserving permutation.

Zero known-positive collisions do not prove biological nonbinding. This comparison controls the retained positive cohort; it does not remove or conceal the selection cost that produced it.
