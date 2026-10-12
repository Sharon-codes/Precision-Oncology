# Independent methods and statistics review

Reviewer role: methods/statistics reviewer for *Interdisciplinary Sciences: Computational Life Sciences*.

Material reviewed: supplied `Precision_Oncology(2).pdf` and the exact declared public repository commit `0ba2d538065e2d3a8b4642547eb77447ecc89a5b`, cloned to `precision_repo`. Findings below distinguish manuscript inconsistencies, directly inspected implementation defects, archived outputs, and experiments still needed. Archived metrics were inspected, not reproduced from biological data.

## Recommendation

**Do not submit the current empirical superiority paper.** The problems are substantive rather than cosmetic: the released implementation contradicts the stated cohort, model, negative filtering, comparator provenance, and statistical evaluation. Headline results cannot be linked to released predictions, checkpoints, or cohort manifests. A careful executable audit and corrected evaluation protocol can be supported now; a strong biological prediction paper requires the missing artifacts and a clean rerun.

No inspected evidence warrants an allegation that predictions were fabricated. The more precise finding is that the released material does not substantiate several headline claims, while some archived results directly disagree with the manuscript. One historical comparator value is explicitly hardcoded in a report generator; this is not a substitute for regenerating that comparison.

## Blocking findings with direct source evidence

### M1. The empirical result lineage is broken

- Abstract/Table 3 report Direct Concatenation ROC-AUC `0.8753` and MCC `0.6145` on `7,416` test pairs.
- Figure 3 itself labels the ROC curve `0.6810`.
- `Evaluation/bootstrap_metrics_results.txt` reports ROC-AUC `0.681 [0.666, 0.696]` and PR-AUC `0.655 [0.635, 0.674]`.
- `Evaluation/unimodal_ablation_results.txt` reports intact-model AUC `0.6810`, with the same masked AUCs (`0.5000`, `0.5613`) that the paper attaches to an intact-model value of `0.8753`.
- `Evaluation/clinical_metrics_report.md` reports MCC `0.2977` at threshold `0.0766`. Its confusion counts are TN=1,363, FP=1,047, FN=655, TP=1,755, summing to **4,820**, not 7,416.
- `bmc_revision_metrics.txt` reports the same 0.6810 result and five-decoy-seed AUC `0.6892 ± 0.0043`.
- Repository-wide text search found no `0.875`, `0.891`, `0.912`, `0.827`, `0.836`, ERGO, EPIC, or DeLong implementation/results in released Python/text/Markdown/JSON files.
- `.gitignore` excludes `Dataset/`, `Processed/`, `Checkpoints/`, logs, attention arrays and occlusion arrays. These are absent from the checked-out release; neither the claimed split counts nor headline estimates can be regenerated.

**Required resolution:** frozen pair manifests, cohort provenance, label definitions, score files keyed by immutable pair ID, checkpoints, run configurations and logs for each table row. Regenerate every empirical table/figure from those exact artifacts. If absent, remove those results from the claimed new evidence and label historical outputs as unverified archived outputs.

### M2. The released model is not the reported model

- `model.py:21–38`: the model uses **30 beta + 30 alpha + 50 peptide–HLA tokens**, not a single-chain beta window of 95 and a peptide window of 15.
- The forward pass actively projects and concatenates alpha embeddings. It is therefore not valid to claim deliberate beta-only scope without evidence all alpha inputs are uniformly absent/zero in the actual cohort.
- `dataset.py:153–159` embeds peptide concatenated to an HLA pseudosequence.
- `modules.py:67–83`: head is flatten → Linear → **ReLU** → Dropout → Linear → ReLU → Dropout → Linear. It does not use the manuscript's GELU/LayerNorm combination. No mean pooling occurs.
- The manuscript also states mean pooling followed by a position-preserving 7,040-dimensional flatten, which is dimensionally inconsistent.

**Required resolution:** write the exact implemented architecture and modality scope; report chain availability and HLA mapping coverage. If a different beta-only model produced Table 3, release that exact implementation and run outputs rather than adapting the prose to an unrelated model.

### M3. Known-positive rejection is absent in the implemented decoy sampler

- `dataset.py:124–132` accepts/stores `positive_triplets_set` as `self.pos_set`.
- The subsequent collator does not read the set. It enforces index derangement and different peptide strings (`178–227`), but never tests whether a generated receptor/peptide/HLA combination is another known positive.
- `run_official_sota_replications.py:82–143` has the same unused positive-set issue.
- Fallback is a literal `AAAAA` peptide if no alternative is found (`dataset.py:218–224`; comparator script `138–140`). This contradicts the claim that decoys exclusively preserve natural sequence distributions and can create an artificial easy-negative signal or missing-embedding error.
- The only active sequence hash in `dataset.py` is MD5 for embedding cache filenames (`27–28`), not the manuscript's SHA256 known-positive trie.
- Even a correctly implemented known-positive exclusion can only prevent *known* false negatives. Unassayed shuffled pairs remain presumed nonbinders; cryptographic hashing does not prove biological incompatibility.

**Required resolution:** reject every known positive at the appropriate identity level; freeze evaluation pairs; quantify rejection and fallback rates; fail loudly if a valid decoy cannot be sampled; prohibit artificial fallback sequences. State decoy labels as assumptions. Include independent experimentally tested negatives where possible.

### M4. Dynamic evaluation negatives prevent a well-defined paired comparison

Every pass through `TCRCollate` draws new negatives. Different models, seeds, batch sizes and row orders can therefore see different test pairs. A common label vector is insufficient for paired testing: the same pair IDs and same labels must appear in the same order.

The first half of every batch is positive and the second half negative. Constant-score top-k results therefore depend on ordering and tie policy. A constant scorer does not deterministically achieve exactly 50 positives among its first 100 ranked entries; that is its expectation under random tie resolution in a balanced cohort.

**Required resolution:** materialize a single test-pair manifest before model scoring; check hash/order equality across predictors. Report uncertainty across alternative *frozen* negative manifests separately from training uncertainty. Make tie handling explicit for P@k.

### M5. Test-set threshold tuning is directly implemented

- `clinical_metrics.py:67–105` loads the test set, constructs test scores, chooses Youden's J threshold on that same test set, then reports MCC/F1/PPV/recall/specificity on the test set.
- `run_official_sota_replications.py:341–350` selects the threshold maximizing **test MCC**, whereas the manuscript says thresholds use Youden's J.

**Required resolution:** select calibration parameters and decision threshold exclusively on validation data; lock them before test evaluation. Historical threshold-optimized test MCC must be marked exploratory and optimistic, not held-out operating performance.

### M6. “Five-seed robustness” does not represent five independently trained models

`bmc_revision.py:344–369` instantiates one model, loads one checkpoint and sets eval mode *before* the seed loop. Seeds `[42,100,2024,7,999]` then resample evaluation negatives. This measures conditional variability under test-decoy sampling, not training reproducibility or population uncertainty.

**Required resolution:** distinguish split seeds, optimization seeds, decoy seeds and bootstrap resampling units. If claiming training robustness, independently train all seeds and preserve their manifests/checkpoints; otherwise rename the result to repeated decoy evaluation of one checkpoint.

### M7. The “official upstream” comparator claim is not supported by the released script

`run_official_sota_replications.py:183–231` and `236–287` define local PyTorch classes called NetTCR2 and TITAN. The file imports no official upstream implementation and contains no native published dataset evaluation. Its `main` uses only the project's Processed train/test files. No released native-split score artifact establishes the paper's 0.8274/0.8362 values.

Matching some architectural motifs is not sufficient to establish an official, parameter-matched replication. Native performance alone would still not rule out domain-shift-specific preprocessing or loading errors.

**Required resolution:** relabel these as local proxy implementations unless exact architecture/parity evidence is available. For official comparisons, pin upstream revisions, use their supported modality/data configurations, document any adaptation, validate native baselines, equalize tuning budgets, and evaluate identical frozen pairs.

### M8. The manuscript splitting algorithm is wrong, although released clustering code uses a valid construction

Paper Algorithm 1 checks preliminary test peptides only against preliminary training peptides, discards close test candidates into the final training set, and never rechecks retained test peptides against those moved records. A chain of distances 2, 2 and 4 creates a counterexample at threshold 2.

In contrast, `preprocess.py:513–620` builds the full threshold graph and allocates complete connected components to train/validation/test, which does guarantee no cross-split peptide edge at the chosen threshold if implemented inputs are consistent. This is a **documentation/implementation discrepancy**, not proof that the released split construction leaks.

The code targets approximately 80/10/10 positive-record ratios, whereas the paper describes a fixed 80/20 split after balanced pairing. Neither actual CSVs nor counts are present to resolve this.

**Required resolution:** replace pseudocode with the connected-component implementation; audit every final split pair; report achieved proportions, component sizes, peptide counts and train/test TCR identity/similarity. Peptide separation alone is not patient/study/TCR separation.

### M9. Cohort provenance and biological scope differ from the paper

- `preprocess.py:628–637` combines VDJdb, McPAS, **CEDAR, donor1 and T-PLL**, rather than the described three repositories alone.
- VDJdb processing (`155–212`) filters score 3 but does not filter host species. The PDF's mouse H-2 groups are consistent with this omitted restriction.
- `standardize_hla` maps coarse types such as A2/B27 to `:01` subtypes (`104–115`), which invents resolution that the observation does not provide.
- `dataset.py:79–87` substitutes locus consensus/one class-I pseudosequence when exact mapping is absent. HLA information therefore may be missing or misrepresented.

**Required resolution:** explicitly filter intended species, preserve uncertainty in HLA resolution, count exact versus fallback mappings, document data source/assay/donor and duplicates, and report chain availability. Do not use latent HLA clustering as evidence of biological discovery if HLA is an explicit input and source/type confounds remain.

### M10. Attention observations do not establish the stated mechanism

`quantify_oversmoothing.py:135–165` trains a comparator for five epochs and inspects only the first test batch. `compute_pairwise_cosine_similarities` (`96–110`) excludes self-pairs but includes padding/separator/missing-alpha tokens. The attention has residual addition and LayerNorm (`65`) and deliberately omits standard query-key scaling (`modules.py:54`).

High uncentered cosine similarity can reflect a shared offset or anisotropy without loss of discriminative residual information. Uniform attention alone does not force residual architectures to lose rank. The pure-attention theorem cited in the manuscript does not directly prove failure of this implementation. Hessian eigenspectrum bounds are neither the sole nor a necessary path to testing the causal mechanism.

**Required resolution:** masked, centered token metrics and effective rank; matched-capacity heads; residual/no-residual and positional/masking controls; full-test distributions; interventions that change homogenization while holding optimization and other factors fixed. Until then say “observed similarity concentration in one comparator,” not thermodynamic feature preservation or established causal rank degradation.

## Statistical and interpretation corrections

1. **Calibration:** a balanced-label prevalence-only predictor has Brier 0.25; the reported 0.3103 is worse. A monotonic reliability curve or low optimized threshold does not establish adequate calibration. Platt scaling must be fitted on validation and evaluated on test; report Brier skill and calibration slope/intercept as appropriate.
2. **Dependence:** paired rows share receptors and peptides. Row bootstrap (`bootstrap_metrics.py:90–103`) and ordinary DeLong inference treat observations more independently than the grouped biological design warrants. Use peptide-component clustered uncertainty, macro-peptide performance and, when data permit, crossed receptor/peptide resampling. Report the number of independent groups.
3. **Comparisons:** DeLong is not an “exact” finite-sample test. Extremely small p-values do not repair a changed test distribution, proxy comparator or unverified pipeline. Use effect sizes/paired intervals first; correct planned multiple comparisons.
4. **Leakage attribution:** `bmc_revision.py:334–335` hardcodes LODO AUC 0.6810, while random-split performance is freshly computed. Different test sets and train distributions confound difficulty and leakage. The 0.0825 contrast is descriptive and cannot by itself identify a causal leakage contribution.
5. **Constant logits:** OOD sequence distance does not mathematically force a linear model/CNN to output a constant. Inspect logit variance, saturation, features, gradients and checkpoint identity. Do not interpret constant scores as proof of a biological extrapolation ceiling.
6. **Below-chance MLP:** 0.2256 becomes 0.7744 under score reversal. Audit labels, score sign, class indexing and join/order mismatches; do not select a sign using the test set.
7. **Ablations:** inference masking creates a different input distribution and cannot alone prove modality necessity or interaction learning. Retrain TCR-only, peptide-only and HLA-only controls, including source/length/peptide-frequency predictors and within-HLA receptor permutations.
8. **Screening:** P@100=0.99 at 50% constructed prevalence is conditional benchmark discrimination. Its EF=1.98 is nearly the balanced benchmark ceiling of 2, not evidence of a comparable real clinical yield. Add low-prevalence sensitivity analysis and true negative data before clinical wording.
9. **Runtime:** CPU measurements cannot establish GPU throughput. Include model/embedding loading, frozen-feature preprocessing, batch size, thread count and repeated timings. Cross-attention complexity is O(Lq Lk d); flatten-MLP also depends on its input length and width.
10. **Biological language:** embeddings are not measurements of thermodynamic variance. Occlusion predicts sensitivity of a classifier, not physical contact. t-SNE/UMAP are descriptive, stochastic displays, not significance tests.

## Minimal decisive experiments, in priority order

### Executable now without original biological artifacts

1. Execute the manuscript split algorithm on a three-peptide edit-distance chain; verify the violation, and verify the released connected-component approach passes the same invariant.
2. Execute the exact released collator on two records whose swapped receptor/peptide pairs are both known positives. Demonstrate known-positive negative labels despite a supplied positive set.
3. Compare sorted versus shuffled input and vary batch size on a synthetic skewed peptide distribution. Quantify fallback, marginal peptide shift and peptide-only apparent discrimination; clearly label all results synthetic diagnostics.
4. AST/source audit: one checkpoint before seed loop, test-tuned threshold, local comparator definitions, ignored positive set, alpha modality and exact token lengths.
5. Exact mathematical controls: add a common vector to distinct token representations to increase cosine without destroying their informative residuals; show why cosine is insufficient for a collapse claim.

These experiments can establish software/protocol failure modes. They **cannot** validate the headline TCR predictor, recreate biological performance, establish a clinical benefit or repair missing provenance.

### Required when original artifacts become available

1. **Freeze lineage:** source row IDs → filtered positives → component split → fixed validation/test decoys → immutable pair IDs → model scores. Verify no known positives among decoys and no cross-split distance violations.
2. **Curation audit:** species, HLA exact mapping, chain availability, data source, donors, assays and deduplication. Compare intended beta-only and paired-chain tasks separately.
3. **Minimal credible benchmark:** independently retrained direct head, capacity-matched attention with masking/residual controls, peptide/HLA/length/frequency-only controls and at least one faithfully executed modern relevant comparator. Shared manifests, validation-only tuning and at least five training seeds.
4. **Operating metrics:** group-level intervals, macro-peptide AUROC/AP, fixed-threshold MCC, held-out calibration and bootstrap early retrieval. Report per-source/HLA coverage and positive-unlabeled caveats.
5. **External test:** independent assay/source with experimentally tested nonbinding pairs if feasible; keep test selection independent of model results.

## Dialogue with other reviewers

Biology reviewer independently confirmed omitted VDJdb species filtering, HLA subtype invention, model modality mismatch and missing positive-set use. They also identified multiple incorrect bibliographic entries and a missing large 2025/2026 benchmark. Editor reviewer agrees the executable-protocol audit is defensible if original artifacts cannot be obtained, but it must explicitly distinguish simulations from biological evidence and explain novelty relative to previous peptide-distance benchmarks.

## Confidence statement

High confidence in the source-code findings and internal document contradictions. High confidence that the current superiority/clinical claims are not ready for submission. Confidence in the eventual empirical hypothesis remains undetermined until the correct data and prediction artifacts are available. Rewriting alone cannot establish those claims.
