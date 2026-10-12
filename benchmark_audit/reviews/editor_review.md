# Independent editorial stress test

Reviewer role: journal editor / computational-life-sciences reviewer. This is an adversarial simulated review, not feedback from the journal. Date: 12 October 2026. Source reviewed: the complete supplied 21-page PDF, with particular attention to abstract, methods, results, figures, declarations, and references. Repository-level findings below are explicitly attributed to the team reviewers rather than independently re-executed by this reviewer.

## Editorial judgment

The topic fits Interdisciplinary Sciences: Computational Life Sciences, whose official scope includes original computational life science research. The current manuscript is not defensible as a submission. Its central claim of a superior, clinically useful, leakage-free predictor is undermined by internal contradictions and insufficient traceability. A substantial rewrite can produce an auditable computational-methodology study, but new simulations do not validate old biological results.

A strong replacement would investigate **negative-sampling integrity and peptide-marginal shortcuts in TCR–peptide evaluation**, use versioned public biological sequences for empirical demonstrations, and offer executable checks and corrected constructions. Its contribution would be the specific coupling between collision rejection, batch composition, sampler-induced marginals and evaluation outcomes. It must not claim novelty for unseen-epitope evaluation or distance splitting, which already have close precedents.

Recommended working title: **Auditing negative sampling and sequence separation in T-cell receptor–peptide benchmarks**. If results strongly establish batch dependence, a sharper alternative is **Batch-dependent decoy sampling creates peptide shortcuts in T-cell receptor–peptide benchmarks**. Use the latter only when real-sequence experiments demonstrate the effect across meaningful conditions.

## Submission-blocking findings in the supplied paper

1. **Experiment identity is unresolved.** Figure 3 identifies the direct-concatenation ROC curve as AUC 0.6810, whereas Table 3, the abstract and narrative identify the headline result as 0.8753. Table 1 separately contains 0.6810. The final paper needs one registry that maps each figure and table to data hash, split hash, sampler, model checkpoint, training seed and prediction file. Do not repair this by selecting the strongest number.
2. **LODO pseudocode fails its advertised invariant.** It tests candidates against the preliminary training set and moves rejected candidates into the final training set. A retained candidate can be near a rejected candidate newly added to training. Actual split manifests must be audited. The methods reviewer subsequently verified that the repository preprocessing uses connected components and can correctly enforce separation; the proven flaw is in the manuscript pseudocode, not evidence of observed repository split leakage. Correct mathematical separation should be proved and checked on all final train/test peptides, including generated negatives.
3. **Negative labels exceed available knowledge.** Permuting known positive sequences and rejecting recorded positive pairs does not establish biological non-binding. The phrase that cryptographic hashing eliminates false-negative biological collisions is false. Use presumed-negative decoys and quantify known-positive collision rates, including collisions after attempted repair.
4. **Architecture is internally contradictory.** Mean pooling discards the token dimension and cannot precede the stated 110-by-64 flattened matrix. Report the actual implementation and input fields. A two-chain input cannot be described as a single-beta-chain model. The biology reviewer reports repository code incorporating alpha and beta chains plus HLA despite the written single-beta scope; parent agent should confirm before correcting.
5. **Population claims conflict with the plotted data.** Figure 5 lists mouse H-2 types while the methods claim human-only filtering. Remove the figure unless its cohort can be proved, or regenerate from a declared correctly filtered cohort. A mixed-species exploratory plot cannot support human immunogenetic interpretation.
6. **Comparator failures are not explained by a valid theorem.** Constant outputs are not a deterministic property of unseen sequences or distance greater than two. Constant predictions in several unrelated baselines need executable input/output diagnostics, orientation checks, tokenizer checks, valid sequence lengths, trained weight verification and label balance checks. A native-split sanity check does not establish fidelity of the out-of-distribution data adapter.
7. **Clinical and mechanistic inferences are too strong.** Balanced synthetic-decoy discrimination does not establish repertoire-scale precision, true non-binding specificity, clinical utility or thermodynamic information preservation. P@100=0.99 on prevalence 0.5 is an early-ranking benchmark statistic, not 99% precision in deployment. A Brier score of 0.3103 is worse than the balanced-cohort constant-probability reference of 0.25. A low threshold does not repair calibration.
8. **No defensible statistical experimental unit is specified.** Repeated decoy draws do not constitute five independently trained model seeds. Pair-level bootstrap and DeLong inference can be overconfident with reused peptides, receptors or donors. Thresholds and calibration must be fitted on validation data. Distinguish training randomness, split randomness and negative-sampling randomness.
9. **Attention geometry is not a causal explanation.** High cosine similarity and uniform attention are diagnostics, not proof of information destruction or of cause for reduced generalization. Padding, sequence lengths, residual pathways and normalization must be controlled. Headwise Hessian bounds are not a necessary route to test this empirical causal claim.
10. **Runtime claims mix conditions.** The written CPU timing experiment cannot establish GPU throughput, and end-to-end embedding extraction must be separated from cached-feature head evaluation. Cross-attention cost should be described using receptor and peptide lengths rather than presenting every configuration as a generic square sequence-length cost.
11. **Bibliography is materially inaccurate.** Core methods are misattributed (see verified examples below), making novelty and comparator interpretation unreliable. Every retained entry must be verified against its DOI/publisher record. Delete filler and unverified papers.
12. **Declarations are incomplete and anonymity is violated.** The submitted PDF names authors and a person acknowledged; it gives an identifying repository URL. Funding, competing interests, author contributions, data/code availability and applicable ethics statements require author-approved factual wording, not invented assurances. Substantive AI-assisted analysis and drafting must be disclosed under the publisher's policy.

## What would make the replacement study significant enough?

The biology reviewer and I agree that simulations alone fixing errors in an unpublished manuscript would be a weak journal research contribution. A real public-sequence audit is essential. The parent agent is obtaining a versioned VDJdb release; that can support a stronger empirical computational study if the following minimum package is achieved:

- An explicit human TRA/TRB policy, peptide bounds, assay-quality threshold, MHC policy and reproducible deduplication, with counts at each step and retained source identifiers. Do not invent alleles to fill missing fields.
- Quantification of duplicates, known-positive collisions, peptide-frequency distributions, sequence-distance violations and rejected-pair fractions on the actual benchmark records.
- A controlled sampler comparison on the same underlying positive records: shuffled versus peptide-sorted batches; naive index derangement; known-positive rejection; corrected constrained matching that preserves declared marginals where feasible. Report sampling failures explicitly rather than silently falling back.
- A peptide-only frequency baseline, TCR-only baseline and interaction-aware baseline trained under identical, leakage-checked splits. A peptide-only shortcut score must use training information only; computing class likelihoods on the test labels produces an oracle diagnostic, not predictive performance.
- Multiple split/sampler seeds and peptide-grouped uncertainty. Dependence-aware summaries and transparent effect sizes matter more than many tiny p-values.
- A small analytic demonstration explaining when index derangement leaves the same peptide sequence unchanged, and when rejection distorts accepted negative marginals. Demonstrate the corrected postcondition by a reusable assertion.
- A concise scope statement: measured effects apply to these synthetic-decoy tasks and frozen public records, not a validated clinical binder classifier.

If only sampler audit and proofs are completed, the output is a coherent protocol-validation manuscript but I would describe its acceptance prospects as uncertain; I cannot confirm a dedicated research-note article category from the official guidance. If public-data model controls support reproducible sampler-induced inflation and correction, an original research article is plausible. Do not market the study as a new deep-learning predictor when its verified contribution is an audit.

## Reviewer debate and adjudication

Methods reviewer: rejects original inferential chain because the written split guarantee is broken, five decoy replicates may masquerade as training seeds, and statistical independence is unestablished. The reviewer confirms the repository connected-component implementation can be correct; this distinction prevents falsely alleging actual split leakage. I agree. The correct remedy is to withdraw those claims pending artifacts, not weaken the wording while retaining their numbers.

Biology reviewer: warns that a simulation-only rescue is unlikely substantial given Grazioli 2022, Drost 2025 and Lu 2026; asks for public-sequence evidence. I agree and recommend a biological-sequence audit with honest synthetic labels. There is no need to claim binding physics or treatment benefit for a computational benchmark-integrity contribution.

Constructive disagreement with an overly strict rejection: missing the old neural checkpoints does not prevent a new publishable study. A preregistered public-data comparison that isolates the sampler mechanism can stand independently. Its title, methods, results and conclusion must be rewritten around the newly executed study; it cannot inherit the old SOTA AUC table.

## Verified journal requirements (publisher source, checked 12 October 2026)

Source: https://link.springer.com/journal/12539/submission-guidelines#Instructions%20for%20Authors (retrieval ref turn17view2).

- Double-blind: anonymize manuscript and associated files; provide separate title page with author/contact details, acknowledgments, disclosures and funding.
- Editable sources required; mathematical manuscripts may use LaTeX; Springer Nature template recommended.
- Abstract: 150–250 words; no undefined abbreviations or unspecified references. Keywords: 4–6.
- Decimal headings: maximum three levels.
- Numeric square-bracket citations; consecutive reference numbering; include available DOI links; sn-basic.bst permitted.
- Tables: numbered, cited sequentially, captioned.
- Figure fonts embedded; EPS preferred for vectors, TIFF for halftones. Minimum resolutions: line 1200 dpi, halftone 300 dpi, combination 600 dpi. Online color is free.
- Relevant competing-interest statements required under Statements and Declarations.
- Substantive LLM use belongs in Methods; the copyediting exemption excludes autonomous content creation.

Supplemental official requirement PDF: https://media.springer.com/full/springer-instructions-for-authors-assets/pdf/1514467_Journal%20special%20requirments.pdf . Recommended order: Abstract; Keywords; Introduction; Materials and Methods; Results and Discussion; Conclusions; Acknowledgments; References. Corresponding author uses an asterisk; equal contribution uses pilcrow (first set) or ampersand (second set).

Official template landing page: https://www.springernature.com/gp/authors/campaigns/latex-author-support . Downloaded December 2024 package: https://cms-resources.apps.public.k8s.springernature.io/springer-cms/rest/v1/content/18782940/data/v12 . Local copy: reviews/springer_nature_template.zip.

Official scope: https://link.springer.com/journal/12539/aims-and-scope (turn13view1). Official ethics page confirms double-anonymous review: https://link.springer.com/journal/12539/ethics-and-disclosures (turn15view2).

Publication cost: official https://link.springer.com/journal/12539/how-to-publish-with-us (turn13view4) confirms a hybrid model and **no APC for subscription publication**. Optional open access is chargeable. This verifies the user's chosen no-APC route, not acceptance.

## Verified reference corrections and prior-art pointers

- TITAN: Weber A, Born J, Rodríguez Martínez M (2021), *TITAN: T-cell receptor specificity prediction with bimodal attention networks*. Bioinformatics 37(Suppl 1):i237–i244. https://doi.org/10.1093/bioinformatics/btab294 . Publisher: https://academic.oup.com/bioinformatics/article/37/Supplement_1/i237/6319659 (turn20search1).
- EPIC-TRACE: Korpela D et al (2023), *EPIC-TRACE: predicting TCR binding to unseen epitopes using attention and contextualized embeddings*. Bioinformatics 39:btad743. https://doi.org/10.1093/bioinformatics/btad743 . Publisher: https://academic.oup.com/bioinformatics/article/39/12/btad743/7467056 (turn20search0).
- Grazioli F et al (2022), *On TCR binding predictors failing to generalize to unseen peptides*. Front Immunol. https://doi.org/10.3389/fimmu.2022.1014256 . Publisher: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.1014256/full (turn20search23).
- Closely related distance-based evaluation: *Assessing the Generalization Capabilities of TCR Binding Predictors via Peptide Distance Analysis*, Castorina/Grazioli and colleagues: https://www.biorxiv.org/content/10.1101/2023.07.29.551100v3 (turn20search27). Have the literature reviewer verify final publication/version and include the relevant comparison. Novelty must exceed proposing peptide-distance separation.

## Decision threshold for the final deliverable

I can be confident in a carefully delimited, reproducible computational audit once the new scripts, datasets, tables and figures agree and the controls survive independent rerun. I cannot be confident in acceptance or in the original biomedical performance claims without their missing experimental artifacts. The final package should make this distinction explicit and include an author-facing evidence ledger describing withdrawn, corrected, retained and newly established claims.
