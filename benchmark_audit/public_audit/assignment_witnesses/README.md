# Actual assignment witnesses

These 58 CSV files record the actual donor row chosen by the closed assignment, unlike the older reconstruction manifests, whose antigen column identifies a representative source row with the same peptide/HLA value. Both reconstruct exactly the same decoy tuples.

All row identifiers index the unchanged VDJdb TSV, starting at zero after its header. Within each CSV, receiver and actual donor identifier multisets are identical, every donor has the receiver's source HLA label, and all generated triples lie outside the declared positive registry. `cycle_id` is local to a CSV and `cycle_position` starts at zero; following actual donor IDs traverses each closed cycle. No fixed points are retained.

The coverage is primary, score 3, annotation-filtered primary, five distance-test partitions, and the 25 receptor folds in each of the primary and annotation-filtered split-first protocols. `verification.json` records cardinalities, checksums, and comparisons with the original representative-value manifests. The annotation-filtered cases exclude HLA-A*08:01 and HLA-A*24:01 from both positives and registry.

Regenerate with `python audit_work/export_assignment_witnesses.py` from the repository root. The script passively records solver-local arrays and does not modify the assignment implementation. Alternative software versions can choose different optimal cycles; this exporter requires exact agreement with the archived chosen decoy values. Certificates establish the specified row-assignment contract, not biological nonbinding.
