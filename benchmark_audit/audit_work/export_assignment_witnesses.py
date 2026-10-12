#!/usr/bin/env python3
"""Export actual closed-assignment donor rows without changing the solver.

Existing reconstruction manifests use a representative source row with the
selected peptide/HLA value. Those IDs reproduce decoy values but need not be
the actual donor occurrence chosen by the assignment. This exporter passively
captures the immutable corrected() function's final selected/donors arrays,
checks its row-level certificate, and compares generated values with every
corresponding archived representative-value manifest. No model is fitted.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "public_audit"))
from run_public_audit import corrected, prepare

EXCLUDED = {"HLA-A*08:01", "HLA-A*24:01"}
FOLD_SEEDS = (17, 41, 73, 107, 211)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def captured_assignment(data, known):
    captured = {}
    previous = sys.gettrace()

    def trace(frame, event, arg):
        if frame.f_code is not corrected.__code__:
            return None
        frame.f_trace_lines = False
        if event == "return":
            captured.update({name: frame.f_locals[name].copy()
                             for name in ("selected", "donors")})
        return trace

    sys.settrace(trace)
    try:
        result = corrected(data, known)
    finally:
        sys.settrace(previous)
    if set(captured) != {"selected", "donors"}:
        raise AssertionError("Passive assignment capture failed")
    return result, captured


def cycles(receiver_ids, donor_ids):
    successor = dict(zip(receiver_ids, donor_ids))
    seen = set()
    annotations = {}
    cycle_id = 0
    for start in sorted(receiver_ids):
        if start in seen:
            continue
        current = start
        members = []
        while current not in seen:
            seen.add(current)
            members.append(current)
            current = successor[current]
        assert current == start and len(members) >= 2
        for position, member in enumerate(members):
            annotations[member] = (cycle_id, position)
        cycle_id += 1
    assert seen == set(receiver_ids)
    return annotations, cycle_id


def export_one(name, data, known, manifest_path, raw, out, fold=None):
    (retained, negatives, _), captured = captured_assignment(data, known)
    selected = list(captured["selected"])
    donors = list(captured["donors"])
    assert len(selected) == len(set(selected)) == len(retained) == len(negatives)
    assert Counter(selected) == Counter(donors)
    receiver_ids = [int(data.loc[i, "source_row_0based"]) for i in selected]
    donor_ids = [int(data.loc[j, "source_row_0based"]) for j in donors]
    assert retained.source_row_0based.tolist() == receiver_ids
    assert Counter(receiver_ids) == Counter(donor_ids)
    expected = [(data.loc[i, "beta"], data.loc[j, "peptide"], data.loc[j, "hla"])
                for i, j in zip(selected, donors)]
    assert negatives == expected
    assert all(i != j and data.loc[i, "hla"] == data.loc[j, "hla"]
               for i, j in zip(selected, donors))
    assert all(triple not in known for triple in negatives)
    positives = list(retained[["beta", "peptide", "hla"]].itertuples(index=False, name=None))
    for feature in (lambda t: t[0], lambda t: t[1], lambda t: t[2],
                    lambda t: t[1:], lambda t: (t[0], t[2])):
        assert Counter(map(feature, positives)) == Counter(map(feature, negatives))

    reference = pd.read_csv(manifest_path)
    if fold is not None:
        reference = reference.loc[reference.fold.eq(fold)].copy()
    assert reference.receiver_source_row_0based.tolist() == receiver_ids
    representative_ids = reference.negative_antigen_source_row_0based.tolist()
    assert all(j >= 0 for j in representative_ids)
    archived = [(raw.iloc[i].cdr3, raw.iloc[j]["antigen.epitope"], raw.iloc[j]["mhc.a"])
                for i, j in zip(receiver_ids, representative_ids)]
    assert archived == negatives

    annotations, ncycles = cycles(receiver_ids, donor_ids)
    witness = pd.DataFrame({
        "receiver_source_row_0based": receiver_ids,
        "actual_donor_source_row_0based": donor_ids,
        "cycle_id": [annotations[i][0] for i in receiver_ids],
        "cycle_position": [annotations[i][1] for i in receiver_ids],
    })
    target = out / f"{name}.csv"
    witness.to_csv(target, index=False)
    return {
        "name": name,
        "witness_file": target.name,
        "witness_sha256": digest(target),
        "reference_manifest": str(manifest_path.relative_to(ROOT)),
        "reference_manifest_sha256": digest(manifest_path),
        "fold_filter": fold,
        "original_rows": len(data),
        "retained_rows": len(retained),
        "cycles": ncycles,
        "representative_ids_differing_from_actual_donors":
            sum(i != j for i, j in zip(donor_ids, representative_ids)),
        "actual_donor_closure": True,
        "within_hla": True,
        "zero_reference_contradictions": True,
        "declared_marginals_equal": True,
        "archived_negative_values_equal": True,
    }


def receptor_folds(data, seed):
    # Independently restate the published frozen fold rule.
    groups = np.array(sorted(set(data.beta)), dtype=object)
    np.random.default_rng(seed).shuffle(groups)
    mapping = {b: fold for fold, block in enumerate(np.array_split(groups, 5)) for b in block}
    return data.beta.map(mapping)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=ROOT / "research/data/vdjdb.txt")
    parser.add_argument("--out", type=Path, default=ROOT / "public_audit/assignment_witnesses")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    cohorts, known, _, _ = prepare(args.source)
    raw = pd.read_csv(args.source, sep="\t", keep_default_na=False)
    clean = cohorts["score_ge_1"].loc[~cohorts["score_ge_1"].hla.isin(EXCLUDED)].reset_index(drop=True)
    clean_known = {t for t in known if t[2] not in EXCLUDED}
    records = []

    def add(name, data, registry, manifest, fold=None):
        records.append(export_one(name, data, registry, ROOT / manifest,
                                  raw, args.out, fold))

    add("primary", cohorts["score_ge_1"], known,
        "public_audit/manifests/common_corrected_seed0.csv")
    add("score3", cohorts["score_eq_3"], known,
        "public_audit/manifests/strict_corrected_seed0.csv")
    add("annotation_filtered", clean, clean_known,
        "public_audit/annotation_sensitivity/manifests/clean_corrected_seed0.csv")
    for seed in range(5):
        partition = pd.read_csv(ROOT / f"public_audit/split_sensitivity/manifests/partition_seed{seed}.csv")
        test_ids = set(partition.loc[partition.partition.eq("test"), "source_row_0based"])
        data = cohorts["score_ge_1"].loc[cohorts["score_ge_1"].source_row_0based.isin(test_ids)].reset_index(drop=True)
        add(f"distance_split{seed}", data, known,
            f"public_audit/split_sensitivity/manifests/test_split{seed}_corrected.csv")
    for prefix, folder, data, registry in (
        ("split_first", "split_first_results", cohorts["score_ge_1"], known),
        ("annotation_split_first", "annotation_filtered_results", clean, clean_known),
    ):
        for seed in FOLD_SEEDS:
            folds = receptor_folds(data, seed)
            for fold in range(5):
                subset = data.loc[folds.eq(fold)].reset_index(drop=True)
                add(f"{prefix}_seed{seed}_fold{fold}", subset, registry,
                    f"crossfit_audit/{folder}/manifests/f{seed}_exact_closed_matching_s0.csv", fold)
            print(f"Verified {prefix} fold seed {seed}", flush=True)

    assert len(records) == 58
    result = {
        "purpose": "Actual row-level cycle witnesses complement existing antigen-value reconstruction manifests.",
        "method": "Passively capture final arrays from unmodified corrected(); verify output against archived decoy values.",
        "scope": "Combinatorial certificates only; no biological nonbinding or prediction claim.",
        "source_sha256": digest(args.source),
        "solver_file": "public_audit/run_public_audit.py",
        "solver_file_sha256": digest(ROOT / "public_audit/run_public_audit.py"),
        "export_script_sha256": digest(Path(__file__)),
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "scipy": scipy.__version__},
        "witness_count": len(records),
        "retained_rows_total": sum(r["retained_rows"] for r in records),
        "all_certificates_passed": True,
        "records": records,
    }
    (args.out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    (args.out / "README.md").write_text(
        "# Actual assignment witnesses\n\n"
        "These 58 CSV files record the actual donor row chosen by the closed assignment, "
        "unlike the older reconstruction manifests, whose antigen column identifies a "
        "representative source row with the same peptide/HLA value. Both reconstruct "
        "exactly the same decoy tuples.\n\n"
        "All row identifiers index the unchanged VDJdb TSV, starting at zero after its header. "
        "Within each CSV, receiver and actual donor identifier multisets are identical, "
        "every donor has the receiver's source HLA label, and all generated triples lie "
        "outside the declared positive registry. `cycle_id` is local to a CSV and "
        "`cycle_position` starts at zero; following actual donor IDs traverses each "
        "closed cycle. No fixed points are retained.\n\n"
        "The coverage is primary, score 3, annotation-filtered primary, five distance-test "
        "partitions, and the 25 receptor folds in each of the primary and annotation-filtered "
        "split-first protocols. `verification.json` records cardinalities, checksums, and "
        "comparisons with the original representative-value manifests. The annotation-filtered "
        "cases exclude HLA-A*08:01 and HLA-A*24:01 from both positives and registry.\n\n"
        "Regenerate with `python audit_work/export_assignment_witnesses.py` from the repository root. "
        "The script passively records solver-local arrays and does not modify the assignment "
        "implementation. Alternative software versions can choose different optimal cycles; "
        "this exporter requires exact agreement with the archived chosen decoy values. "
        "Certificates establish the specified row-assignment contract, not biological nonbinding.\n"
    )
    print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
