# Public VDJdb negative-construction audit

This is a new audit of a pinned public release. It does **not** reproduce the
original manuscript's training cohorts or run a protein language model, Mamba,
or a trained binding predictor. The source collator is executed unchanged from
its class AST using numerical identity tags in place of embeddings. Consequently
the selected sequences and labels are audited exactly, while embedding quality
and model predictions remain outside scope.

## Reproduce

Expected surrounding paths in the reproduction bundle:

```
audit_work/run_source_audit.py
precision_repo/dataset.py
precision_repo/preprocess.py
precision_repo/train.py
precision_repo/run_official_sota_replications.py
research/data/vdjdb.txt
public_audit/run_public_audit.py
```

Python packages: numpy, pandas, scipy, scikit-learn. Exact execution versions and
the original repository commit/source SHA256 hashes are in `provenance.json`.

```
python public_audit/run_public_audit.py --seeds 30
python public_audit/run_split_sensitivity.py
python public_audit/run_reference_sensitivity.py
python public_audit/run_annotation_sensitivity.py
```

The main script also accepts `--source PATH --repo PATH --out PATH`. Extract the
unaltered VDJdb 2023-06-01 release before running. The source TSV SHA256 must match
the recorded provenance. Use the original source repository commit recorded in
provenance; the auxiliary AST wrapper verifies no neural weights.

## Population and references

The source contains 92,771 records. Sequential filters retain human sequences,
TRB, MHCI, syntactically exact two-field `HLA-[ABC]*DD:DD` strings (no imputed fields), canonical
amino-acid CDR3 beginning C and ending F/W, and canonical peptides of length8–15.
No CDR3 length restriction is silently imposed. Filtering precedes exact
beta/peptide/HLA deduplication. The main cohort is score>=1, with score>=0 and
score==3 sensitivities. All quality-eligible score>=0 positive triples constitute
the broad known-positive reference; the separate confidence-reference analysis
checks whether collisions remain using only score>=1 reference positives.

Source release ordering is retained. Peptide-sorted order uses peptide,HLA,then
source order; loader-shuffled order uses NumPy RNG100000+seed. These changes do
not alter the source collator's global fallback pool, which retains first-seen
unique peptide/HLA order. Collator Python RNG seed is separate. Main sampling
uses30 seeds0–29; all cohorts B64, main alsoB16 andB256. Each run creates one
negative record per positive. Monte Carlo variation describes construction
randomness, **not** uncertainty over biological populations.

## Metrics and limitations

- Known-positive collisions count generated decoys equal to a recorded positive
  beta/peptide/HLA triple. Beta/peptide matches ignoring HLA are reported separately.
  Their frequencies are not estimates of all biological false negatives.
- Peptide, HLA, joint peptide/HLA, and beta marginal total-variation distances
  compare equally weighted positive/negative classes.
- The identity-oracle AUROC assigns each category its positive fraction in the
  **same generated dataset**. It is a resubstitution diagnostic of marginal
  separability. It is not a trained model, held-out performance, or evidence of
  generalization to unseen peptides or receptors.
- Removing known-positive contradictions does not establish that remaining
  recombinations are experimentally verified nonbinders.

## Corrected finite-cohort construction

Within each HLA stratum, rows are both receptor sources and antigen donors.
An allowed off-diagonal assignment has cost0 iff the resulting beta/peptide/HLA
triple is absent from the broad positive reference. A diagonal self-assignment
costs1 and denotes abstention. Every forbidden assignment costsN+1, whereN is
stratum size. SciPy linear_sum_assignment solves the full assignment.

The identity permutation always costsN. Any forbidden edge alone costsN+1, so
none can be selected in an optimum. Every selected non-self edge belongs to a
cycle of the assignment permutation. Removing all diagonal self-loops therefore
retains exactly the same source and donor row set. Minimizing cost minimizes
the number of abstentions and maximizes retained row count among closed
within-HLA permutations satisfying the known-positive exclusion constraint.
This is an application of the standard assignment problem, not a claim of a new
optimization algorithm or a universally representative negative population.

We directly assert exact beta, peptide, HLA, and joint peptide/HLA class marginal
equality and zero known-positive collisions. The main solution retains2691/3182
rows,236/253 peptides,22/41 HLA types. Retention is biologically selective; rare or
single-peptide HLA strata can be entirely excluded. Original-vs-corrected
comparisons on the **same retained positive rows** accompany full-cohort source
diagnostics. A mere maximum partial bipartite matching would not suffice because
its matched source and donor subsets need not agree; the cycle-cover objective
explicitly enforces that agreement.

## Split sensitivity

All253 main-cohort peptides form221 connected components at unit-cost
Levenshtein distance<=2. For five prespecified RNG seeds,132/44/45 components go
to train/validation/test. These are approximately60/20/20 fractions of components,
not of records; actual record counts are reported. Exhaustive cross-partition
peptide comparison confirms minimum distance3. Source and corrected construction
are rerun exclusively inside each test partition, with10 sampler seeds per
ordering. This verifies construction behavior under strict peptide separation;
it remains a dataset audit, not an OOD model-performance experiment.

## Files and redistribution

- `seed_results.csv`, `aggregate_results.csv`:541 main runs (540 stochastic source
  runs plus one deterministic correction); means, SD, range, and empirical2.5/97.5
  percentiles. These percentile ranges are not population confidence intervals.
- `matching_retention_by_hla.csv`:full/retained stratum sizes and diagnostic
  maximum partial matching cardinalities.
- `manifests/`:source-row references reconstruct every frozen sample. Index0 is
  the first TSV data row after the header. Each negative uses the receptor from
  `receiver_source_row_0based` and peptide/HLA from
  `negative_antigen_source_row_0based`. Antigen row references may point to the
  first equivalent antigen pair in the eligible release, not the literal
  donor object selected inside an otherwise indistinguishable batch. A value-1
  denotes the original source's hardcoded `AAAAA`/`A*02:01` sentinel.
- `local_pairs/`:local sequence materializations for analysis; exclude these
  derived sequence tables from redistributed artifacts. Recreate them by running
  the scripts from the unchanged VDJdb release, whose license is CC BY-ND4.0.
- `manifest_validation.json`:independent reconstruction checks against the TSV.
- `split_sensitivity/`:five strict component-split diagnostics and manifests.
- `confidence_collision_sensitivity.csv`, `additional_sensitivity.json`:
  reference-confidence and score3 correction sensitivity.

The source data are heavily based on experimentally studied viral antigens.
This audit supports computational evaluation claims and does not validate
patient-specific neoantigen selection or clinical precision oncology utility.

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
