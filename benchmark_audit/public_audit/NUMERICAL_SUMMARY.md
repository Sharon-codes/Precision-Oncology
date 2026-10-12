# Numerical summary — public-release audit

All source-collator estimates below summarize30 RNG seeds. SD describes construction randomness, not biological uncertainty. Identity-oracle AUROC is a same-dataset marginal diagnostic, not held-out model performance.

## Primary cohort:3,182 positive triples, B64

| Order | Known-positive collision % (mean ±SD) | Fallback % | Peptide TV | Identity-oracle AUROC |
|---|---:|---:|---:|---:|
| release_order | 2.9887 ± 0.1022 | 13.4507 | 0.129415 | 0.604779 |
| peptide_sorted | 0.8820 ± 0.0761 | 47.0773 | 0.434192 | 0.791650 |
| row_shuffled | 0.2284 ± 0.0695 | 0.0000 | 0.000000 | 0.500000 |

## Common retained cohort:2,691 positive triples, B64

| Order | Known-positive collision % (mean ±SD) | Fallback % | Peptide TV | Identity-oracle AUROC |
|---|---:|---:|---:|---:|
| release_order | 2.5319 ± 0.1086 | 12.5975 | 0.121120 | 0.594188 |
| peptide_sorted | 0.6838 ± 0.0694 | 39.7250 | 0.377939 | 0.761059 |
| row_shuffled | 0.2564 ± 0.0886 | 0.0012 | 0.000012 | 0.500012 |
| Corrected closed assignment |0|0|0|0.500000|

## Confidence/reference sensitivities

- All-score cohort29,306: shuffled B64 collision5.1545%; peptide TV0.008535; oracle0.506612.
- Score3 cohort320: majority peptide186/320(58.125%); shuffled B64 collision0.09375%; TV0.167083; oracle0.606925. Corrected within-HLA assignment retains139(43.44%).
- Main correction using all-score reference retains2691; using only score>=1 reference retains2696 but its assignment contains8 collisions against the broader reference.

| Main cohort B64 order | Broad-reference collisions % | Score>=1-reference collisions % |
|---|---:|---:|
|peptide_sorted|0.8820|0.8003|
|release_order|2.9887|2.6639|
|row_shuffled|0.2284|0.1446|

## Coverage

Primary exact correction retains84.57% of rows,236/253 peptides,22/41 HLA alleles and2584 distinct beta sequences. Missing HLA strata and dominant-peptide restrictions are an explicit population-selection limitation. The closed assignment is optimal for retained row count under the stated within-HLA permutation and reference-exclusion constraints; it does not optimize biological representativeness.

## Strict peptide-component sensitivity

253 peptides yield221 edit-distance<=2 connected components. Five component splits (132/44/45train/validation/test) have cross-partition minimum edit distance3. No predictive model is trained in this experiment.

| Split seed | Train rows | Validation rows | Test rows | Retained test rows | Retention % |
|---|---:|---:|---:|---:|---:|
|0|1629|531|1022|509|49.80|
|1|1771|534|877|546|62.26|
|2|1808|731|643|380|59.10|
|3|1703|1107|372|271|72.85|
|4|2153|436|593|190|32.04|

## Interpretation

Shuffling removes peptide/HLA marginal bias in the primary cohort but leaves catalogued-positive contradictions. Under a stringent confidence threshold, a majority peptide makes exact distinct-peptide marginal-preserving permutation impossible on the full cohort, and shuffled loading alone retains substantial marginal separability. Corrected construction removes recorded-positive contradictions and the stated marginal disparities on its retained subset, with meaningful and explicitly quantified coverage loss. Remaining decoys are unverified nonbinders. These results establish construction properties; they do not reproduce original neural-model results or establish clinical utility.

## Annotation sensitivity

The primary population is a raw-annotation audit. Its HLA syntax filter does
not establish current nomenclature validity. The independent nomenclature review
identified unresolved HLA-A*08:01 and deleted legacy HLA-A*24:01. The secondary
script run_annotation_sensitivity.py excludes both exact labels from each cohort
AND the broad known-positive reference, with no recoding or imputation. It reruns
B64 with three loader orders and30 seeds on the clean full3108-row cohort and the
same2665-row retained cohort. Correction retains85.75% with0 reference-positive
collisions and exact marginal equality. Clean full-cohort release/sorted/shuffled
collision rates remain2.8882/0.9331/0.2445%; peptide TV is0.182937/0.439168/0.
Nomenclature status is not used to infer binding biology. Removing records
changes batch boundaries, so result changes are not interpreted biologically.
Evidence URLs and row exclusions are in annotation_sensitivity/provenance.json.
