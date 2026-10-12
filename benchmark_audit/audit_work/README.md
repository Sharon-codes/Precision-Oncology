# Executable source audit and synthetic sampler controls

This directory audits source logic at repository commit `0ba2d538065e2d3a8b4642547eb77447ecc89a5b`.

## Scope

The audit executes the unmodified Python AST of `dataset.py:TCRCollate`, selected exact preprocessing methods, and the exact disease-holdout branch of `train.py`. PyTorch is not installed. A deliberately narrow NumPy shim provides only `stack`, `cat`, and `tensor`, and identity-tag embeddings preserve sequence/HLA identity. Thus this is an executable test of sequence selection and labels, **not neural-model inference, ESM reproduction, or verification of the paper's original metrics**. The actual source code is not modified.

The synthetic batch-order experiment uses 4,096 unique synthetic beta identifiers and nine distinct peptide categories: one comprises 50% of records and each of the other eight 6.25%. All records share one HLA. There are no underlying biological labels or simulated molecular binding mechanism. Labels are generated exclusively by the sampler. Batch size is 64. The prespecified score is one for the majority peptide and zero otherwise. No model is fitted and no labels are used to select the score. One hundred seeds (0–99) are run for each ordering; the shuffled order uses seed + 100,000. The sorted case consists of homogeneous batches by construction.

The fixed cyclic matching control shifts the globally sorted records by the largest category count. In this particular feasible design it preserves peptide and TCR marginals and avoids the known-positive pairs. This construction is a control, **not a proposed universally feasible matching algorithm**.

## Reproduce

Requires Python 3 with NumPy, pandas, scikit-learn, and (for figures) matplotlib:

```
python run_source_audit.py --repo ../precision_repo --seeds 100
python plot_source_audit.py
```

`source_audit_results.json` includes exact source hashes, version metadata, deterministic counterexamples, and Monte Carlo summaries. `source_audit_seed_results.csv` contains all 200 individual results.

## Findings and limits

1. The manuscript's described preliminary-filter/complement procedure admits a three-peptide counterexample at distance threshold 2. This does **not** establish that its original data actually contain this violation.
2. In contrast, the repository's actual connected-component preprocessing keeps that chain intact and passes all 72 cross-partition peptide comparisons of an additional 21-record synthetic fixture (minimum distance 9). The manuscript pseudocode and repository method must not be conflated.
3. The separate `leave_one_disease` branch filters disease records, then clusters only the remainder. The exact branch allows a constructed train/test peptide distance of 1. Disease holdout and peptide homology holdout are distinct requirements.
4. `TCRCollate` writes `self.pos_set` once and never reads it. A controlled two-record cross-positive fixture returns two known-positive triplets labeled as negatives. With the same random seed, providing the known-positive set or an empty set yields identical pair selections. This establishes a logical failure case, not the prevalence of false negatives in the original dataset.
5. A one-peptide pool exhausts 20 choices and uses the literal peptide `AAAAA`. Original embeddings might not contain it; the source audit deliberately uses identity tags to test selection independently of that separate cache failure.
6. The synthetic majority-peptide score achieves AUROC 0.718730 ± 0.001597 (mean ± SD over seeds) for peptide-sorted rows, 0.524587 ± 0.002293 for shuffled rows, and exactly 0.5 with marginal-preserving matching. The global-fallback fractions are 100% and 4.91748% respectively. These demonstrate a mechanism by which order changes the negative distribution. They are **not measured performance of the original TCR models or an estimate of the original study's bias**.

Quantiles in JSON summarize Monte Carlo variation across seeds. They are not population confidence intervals, and seeds are not independent biological cohorts.
