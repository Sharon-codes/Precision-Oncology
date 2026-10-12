# Final independent editorial review of the rewritten paper

**Reviewer role:** adversarial journal editor for Interdisciplinary Sciences: Computational Life Sciences (simulated, not journal feedback).

**Reviewed file:** `paper/manuscript.tex`, alongside the public-audit and trained-control README files and the verified negative-sampling prior-art addendum. This review assesses evidence/claim consistency and scientific framing; it does not claim an independent rerun of every analysis.

## Decision

The rewritten manuscript is a credible, carefully bounded computational benchmark-construction study. It is substantially stronger scientifically than the supplied manuscript because its central claims now correspond to executed, inspectable analyses. I found no remaining fatal biological or mathematical claim in the current framing. Complete the small corrections below before delivery, then I would support author review and consideration for submission as original research.

This is an incremental methodology contribution. Its journal acceptance prospects remain uncertain because it audits one research-code implementation and does not demonstrate improved biological prediction. Those are limitations of significance, not evidence that the current results are invalid. Recovering the unavailable original neural checkpoints is not required for this independently framed study and should not block delivery.

## Required changes

### 1. Remove the residual implication that absent registry pairs are nonbinders

Current sentence:

> Generated pairs absent from $K$ remain unverified nonbinders.

Replace with:

> Generated pairs absent from $K$ have unknown biological binding status.

This makes the empirical target precise. The biology reviewer independently agrees.

### 2. State the audited implementation's scope at first introduction

Current opening in the sampler subsection:

> The audited source sampler receives a positive registry and a pool of distinct peptide--HLA values.

Suggested replacement:

> We audited one research-code implementation of a batch sampler, without assuming that it represents unrelated published methods. The sampler receives a positive registry and a pool of distinct peptide--HLA values.

The limitation later in the Discussion is already good. Putting this qualifier in Methods prevents the reader initially interpreting the source audit as evidence about the prevalence of a field-wide problem. An anonymous copy must actually be available to reviewers; a fingerprint in an unavailable author-only archive would not provide reproducibility.

### 3. Identify the HLA-conditioning difference in the comparator interpretation

Insert after the independent-rejection description or at the beginning of the common-population results:

> The source sampler can reassign both peptide and HLA from a donor ligand, whereas rejection and closed assignment retain the receiver's HLA stratum. Consequently source-versus-assignment comparisons combine registry exclusion, HLA conditioning and coordinated marginal preservation. The within-HLA rejection-versus-assignment comparison uses the same retained positives, registry and admissible donor rule, more directly isolating the effect of coordinating donor assignments.

The current manuscript is not falsely attributing all source differences to one cause, but this explicit distinction closes an avoidable causal-comparison objection. The biology reviewer agrees this is the most important remaining interpretation clarification.

### 4. Align Table 2's uncertainty caption with its actual cells

The caption says entries are means plus/minus SD, but source peptide TV and oracle AUROC cells display means only. Either add their SD from the existing aggregate outputs or revise:

> Entries are means over 30 sampler seeds; selected standard deviations are shown and the complete per-seed summaries are supplied in Supplementary Tables. The assignment is deterministic.

Prefer adding SD if the table remains legible. This is a presentation correction, not a need for more experiments.

## Recommended improvements if convenient

1. Add one sentence in the Discussion stating the positive significance claim directly: **marginal balance, registry exclusion, cohort coverage and fold placement are distinct, testable construction properties; passing one does not imply passing the others**. This is the strongest transferable lesson.
2. If retaining the sentence attributing standard assignment/matching to Hall (1935), add a standard linear-assignment reference or the exact SciPy solver reference. Hall is a defensible source for matching feasibility, but it is not the algorithm used by the implementation. Do not add a citation without checking its bibliographic record.
3. Avoid overemphasizing SD across sampler-seed means in the main trained-control plot. The manuscript already explains it correctly and also gives all-run ranges. If one simple uncertainty summary is desired, the all-run range is more transparent about partition variation.
4. Keep the source-versus-assignment figures and the real-data coverage plot. Do not restore the old neural performance, attention geometry or clinical-utility figures; they are irrelevant to the validated new study.

## Novelty test against close prior art

Dens (2023) already establishes that negative-data choice can create apparent TCR specificity performance and that shuffled negatives are informative controls. Mi (2026) already proposes a constructive negative strategy and evaluates predictive effects. The rewritten introduction explicitly acknowledges both, so it does not claim that negative bias or negative construction is newly discovered.

The narrower contribution survives this comparison: a concrete fallback from donor rows to unique ligand categories is executed and traced; registry contradictions and marginal separation are measured separately; a closed within-HLA assignment provides an exact retained-population contract; biological coverage lost to that contract is quantified; and a partition-placement sensitivity exposes a complementary-fold dependence in metadata controls. The work is strongest as a transparent demonstration with reusable assertions, not as a universally preferred negative sampler.

There is no need to claim that an AUROC of 0.556 represents impressive predictive performance. It is informative precisely because a modest construction signal is detected without receptor features under the specified held-out-receptor protocol. The much larger same-data oracle AUROC is correctly segregated from generalization claims.

## Adversarial reviewer debate

The biology reviewer and I agree that the original paper's biological interpretation had to be withdrawn. We now agree that the replacement's limitations are specific and sufficient: beta-level registry semantics are stated, missing HLA is not invented, the population is filtered, synthetic decoys are not laboratory negatives, and the assignment's coverage loss is prominently reported.

I explicitly reject demanding unavailable checkpoint reruns as a condition for this rewritten paper: the manuscript no longer claims to replicate or validate those experiments. A reviewer may reasonably ask for another dataset or interaction-aware predictor as a significance extension, but the current construction properties do not require a neural binding model to be true.

The methods reviewer should provide the final code/figure numerical consistency check. My editorial conclusion is conditional on the delivered archive actually containing the declared source copy, manifests, raw aggregate outputs, fixed package versions and executable figure scripts.

## Submission status

The scientific manuscript can be delivered confidently within its explicitly bounded claims after these corrections. Author-specific funding, competing-interest and contribution statements still require factual completion and author approval; do not invent them. Verify all bibliography records, compile, inspect figures and ensure the review archive is anonymized. Acceptance cannot be guaranteed, and the file should not be described as already peer-reviewed by the journal.

## Late cross-review finding: source HLA annotation validation

The biology reviewer identified two primary-cohort rows annotated `HLA-A*08:01`, a label flagged as invalid in the VDJdb source's own patching/validation resources. Both rows are excluded by the assignment. The regular expression therefore validates two-field **syntax**, not all biological allele identities. Required correction: replace “unambiguous two-field HLA-A/B/C strings” with “source labels matching two-field HLA-A/B/C syntax,” describe HLA counts as source-annotation strata where appropriate, and explicitly disclose this unresolved annotation. A reported sensitivity dropping those two rows (3,180 positives, 40 label strata, 2,691 retained if rerun confirms) is preferable to silently recoding them. Verify the sensitivity before reporting it. This does not refute the construction theorem but prevents overstating biological curation.

The methods reviewer subsequently confirmed the headline numbers against the actual outputs and found no blocking mathematical/statistical claim. The reviewer independently agreed with the HLA-conditioning comparison clarification and Table 2 caption correction above.
