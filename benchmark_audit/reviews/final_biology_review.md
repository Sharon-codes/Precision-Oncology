# Final independent biological review

Reviewer: computational immunology and biological interpretation. Reviewed the rewritten `paper/manuscript.tex`, frozen VDJdb cohort, construction code and derived manifests. This review concerns the new benchmark-construction study, not the discarded claims of a superior neural model.

## Decision

The rewritten study has a defensible scientific question and substantially better evidential discipline. I can support it as a narrow, exploratory construction audit once the HLA-annotation issue below is explicitly resolved in text and the new sensitivity results are incorporated. I cannot certify journal acceptance or call its novelty major. The strongest plausible contribution is the combination of an executable failure mechanism, a precisely stated closed-assignment contract, its measurable population exclusions, and partition-placement controls on real public records.

The results do not demonstrate that the constructed decoys are biological nonbinders or that any interaction-aware model generalizes better. The draft now states this accurately. A request to restore affinity, structural-contact, clinical specificity, precision-oncology or state-of-the-art model claims should be rejected unless new direct evidence is supplied.

## Essential remaining correction: source HLA strings are annotations, not validated alleles

The eligibility regex checks the syntax of raw HLA-A/B/C strings. It does not establish that each label denotes a currently recognized allele or that the original assay was correctly annotated. The phrase “unambiguous two-field HLA-A/B/C strings” should therefore become “raw HLA-A/B/C annotations matching a two-field syntax.” Counts should refer to HLA annotation strata. A table label of “HLA strata” is acceptable after this scope is made explicit.

An independent comparison of all 41 primary labels to the official IPD-IMGT/HLA 3.52.0 WMDA nomenclature file, dated 17 April 2023, found two labels with no active allele sharing the recorded two-field prefix:

| Raw source annotation | Primary rows | Originally retained rows | Interpretation |
|---|---:|---:|---|
| HLA-A*08:01 | 2 | 0 | Unresolved source annotation absent from the official nomenclature comparison; no replacement inferred. |
| HLA-A*24:01 | 72 | 26 | Corresponds to legacy A*2401, which the official deletion list records as removed in March 1995 because its sequence contained errors. No replacement inferred. |

The official files were retrieved from commit `38398a75e9762ff070d8e9bd714d074332646cd7` in the IPD/IMGT repository. The git tag is `v3.52.0-alpha`; the files themselves identify release 3.52.0 and date 2023-04-17. This is a dated reference preceding the June 2023 VDJdb release, rather than an inference from later nomenclature changes. Presence of the remaining 39 prefixes does not independently verify the assay annotations.

Official sources:

- https://raw.githubusercontent.com/ANHIG/IMGTHLA/38398a75e9762ff070d8e9bd714d074332646cd7/wmda/hla_nom.txt
- https://raw.githubusercontent.com/ANHIG/IMGTHLA/38398a75e9762ff070d8e9bd714d074332646cd7/Deleted_alleles.txt
- Human-readable current deletion list: https://github.com/ANHIG/IMGTHLA/blob/Latest/Deleted_alleles.txt

Evidence, hashes, all 41 comparisons and independent sensitivity-manifest reconstruction checks are in `official_annotation_validation.json`; the checker is `validate_annotation_evidence.py`. The official files are saved beside them. The dated nomenclature SHA-256 is `8f77746876f13981d3f593e360e9cc35fd6b6b6b746db1f2f5e20d1a21c0cbd2`; the deletion-file SHA-256 is `de304e760edf70e60ee7524cf21572da96f99e132c52f62aa31b825895cfc1eb`.

### Independently executed annotation sensitivity

I removed each flagged label from both the primary positives and the broad known-positive registry, without recoding it, and reran the exact closed-assignment implementation:

| Exclusion | Primary rows | Registry triples | Retained rows | Retention | Retained/source HLA strata |
|---|---:|---:|---:|---:|---:|
| None, raw annotation audit | 3,182 | 29,306 | 2,691 | 84.57% | 22/41 |
| A*08:01 only | 3,180 | 29,241 | 2,691 | 84.62% | 22/40 |
| Both flagged annotations | 3,108 | 29,146 | 2,665 | 85.75% | 21/39 |

The A*08:01-only retained manifest is exactly unchanged. Removing both labels yields exactly the original retained source-row/negative-tuple manifest with the flagged strata removed. All declared marginals remain identical, and registry contradictions remain zero. After both exclusions, 235 of 251 peptides remain; original-to-retained HLA-distribution TV is 0.11884. The analysis code and machine-readable result are `run_annotation_sensitivity.py` and `annotation_sensitivity.json`.

The public-data reviewer separately reran full and common-retained source-sampler experiments for the two-label exclusion. I independently reconstructed all seven frozen seed-0 manifests from original source rows, checked positive and negative counts, verified absence of both labels, and verified exact equality of the corrected pair table with the original corrected table after stratum removal. The main source counts and excluded zero-based source IDs are in `public_audit/annotation_sensitivity/provenance.json`.

My exact-manifest assertion applies to the closed assignment only. It is not a claim that receptor-fold assignments, cross-validation AUROCs, source-loader behavior or peptide-distance splits are unchanged. The methods reviewer separately completed the metadata rerun: sorted peptide-plus-HLA AUROC 0.561857, release-order 0.518775, and shuffled/assignment 0.5. The five receptor-fold assignments retained 2,618, 2,625, 2,636, 2,642 and 2,621 rows. These new values support the limited construction conclusion while demonstrating why invariance could not be assumed.

## Biological interpretation checks now satisfied

1. **Label scope.** “Registry contradiction” is a suitable term for an exact beta/peptide/HLA tuple already present in the declared registry. It is not a measured false-negative rate for a full alpha-beta receptor. The draft now explains the omitted alpha chain, V/J genes and full receptor context. Pairs absent from the registry should remain “unobserved pairs of unknown binding status.” The current revision uses that wording.

2. **HLA compatibility and comparator fairness.** The source sampler reassigns peptide and HLA together; assignment and the independent rejection comparator preserve the receiver’s recorded HLA stratum. Consequently, source-versus-assignment differences combine HLA restriction, known-positive exclusion, and donor coordination. They cannot be attributed solely to marginal preservation. The independent within-HLA rejection-versus-assignment comparison isolates independent versus coordinated donor selection under the shared biological constraints. The current draft states this distinction adequately.

3. **Retained population.** Loss of 19 raw annotation strata and exclusion of 15.43% of primary rows are substantive results. Marginal matching is only guaranteed on retained rows; it does not preserve the original full cohort. The coverage discussion and HLA-B*27:05 closed-versus-bipartite example convey this clearly. The two-label sensitivity should be reported alongside, rather than silently replacing, the original counts.

4. **Metadata control interpretation.** The split-first sorted-batch AUROC near 0.556 shows modest label predictability from antigen metadata across held-out exact beta groups within an antigen collection. It does not show unseen-peptide recognition. The exact 0.5 matching control follows from the construction within folds, not learned molecular interaction. The draft now says this. Lower-than-chance globally constructed controls are retained and explained rather than flipped; this is the correct interpretation.

5. **Partition feasibility.** Distance separation is a construction sensitivity, not evidence of improved generalization. No trained binding model or assay-confirmed negative experiment has been run. The draft accurately limits the result.

6. **Implementation scope.** The source audit is one research-code implementation. It should never be called a representative industry or field practice. The current draft makes this explicit and documents identity-tag execution as sequence-selection equivalence, not neural numerical equivalence.

## Prior-art challenge and response

The central danger is rediscovering established negative-data bias with a standard assignment solver. The manuscript should retain explicit acknowledgment of the nearest work:

- Ceder Dens, Kris Laukens, Wout Bittremieux and Pieter Meysman (2023), “The pitfalls of negative data bias for the T-cell epitope specificity challenge,” *Nature Machine Intelligence* 5, 1060–1062. DOI: https://doi.org/10.1038/s42256-023-00727-0. This is the journal version of bioRxiv 2023.04.06.535863. It already establishes how negative-data construction can inflate evaluations and examines shuffled controls.
- Xue Mi et al. (2026), “Mitigating negative data bias to enhance TCR–epitope binding prediction,” *Briefings in Bioinformatics* 27(4), bbag418. DOI: https://doi.org/10.1093/bib/bbag418. Its receptor-junction mutation construction is a different solution to negative-data bias. It must not be omitted because it appeared after the draft’s original citations.
- Castorina et al. (2025), *PLOS ONE* 20(5), e0324011. DOI: https://doi.org/10.1371/journal.pone.0324011. Peptide-distance and unseen-antigen concerns are established prior art, not a new contribution here.

The current manuscript cites these works and presents its novelty cautiously. A defensible response to a skeptical reviewer is: this study traces a particular row-permutation-to-unique-category fallback, measures its dependence on batching, separates known-positive exclusion from marginal control, and gives an executable maximum closed-row-subset contract with actual coverage costs. Partition-placement sensitivity adds a useful finite-data diagnostic. The Hungarian solver and the fact that balancing marginals removes a metadata shortcut are not themselves novel.

Simple shuffling already achieves the main marginal certificate in the raw primary cohort, so the assignment should not be sold as the only effective control or as universally preferable. Its additional guarantee is joint exclusion plus closed marginal preservation, purchased with coverage loss and possible duplicated negative tuples. The current discussion admits these tradeoffs.

## Adversarial reviewer exchange

I challenged the editor’s concern that a simulation-only protocol note would be too slight: the rewritten study now contains a real frozen public cohort, source execution, learned metadata controls, distance-separated construction and independent annotation sensitivities. This supports a substantive implementation audit. I agree with the editor that its general scientific significance remains bounded because one source sampler and one database cannot establish prevalence across unrelated models.

I challenged the methods framing of “same positives” as a sufficient controlled comparison: changing HLA restriction changes the candidate space even on the same positive cohort. The draft now identifies the within-HLA rejection comparator as the coordinated-donor comparison. I also flagged that removing raw HLA labels cannot justify asserting unchanged receptor-fold results; the methods reviewer subsequently reran that sensitivity and obtained the new values reported above.

## Final biological recommendation

Retain the narrow construction-study title and remove “unambiguous” from HLA eligibility. Explicitly report the annotation exclusions, add the official dated source evidence, and incorporate the independently rerun sensitivity results. With those changes, I find no remaining unsupported biological efficacy claim in the current manuscript. The work is scientifically reviewable as an exploratory benchmark-construction audit; stronger claims of broad impact or acceptance at the target journal would exceed the evidence.
