#!/usr/bin/env python3
"""Fresh-review source portability, numerical replay, and artifact checks.

This intentionally does not import the archived-result verifier. All model
retraining and figure generation occur in temporary directories. It tests a
declared subset (fold seed 17 / sampler seed 0), not the full experiment grid.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import Counter

import pandas as pd


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execute(root, script, *args):
    run = subprocess.run([sys.executable, str(root / script), *map(str, args)],
                         cwd=root, capture_output=True, text=True)
    if run.returncode:
        raise AssertionError(f"{script} failed:\n{run.stdout}\n{run.stderr}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    report = {"scope": "Fresh source portability, selected construction/model replay, and figure/table regeneration; not a full-grid retraining."}
    source_zip = root / "research/data/vdjdb-2023-06-01.zip"
    assert digest(source_zip) == "7c4994ce32eba8e0ff6c23d1ec53b85a982be69651ed1f467ffe64c86c66c1d1"
    with zipfile.ZipFile(source_zip) as archive:
        tsv = archive.read("vdjdb.txt")
        assert hashlib.sha256(tsv).hexdigest() == "2823281ece02f37b15cdd9187be374e99ee7b46e041572fc2cde4b295a9bf52e"
    report["source_zip_and_tsv_hashes"] = "PASS"
    # Reconstruct each new row-level certificate directly from the source TSV,
    # independently of the exporter and its imported assignment implementation.
    witness_dir = root / "public_audit/assignment_witnesses"
    if witness_dir.exists():
        import io
        source_frame = pd.read_csv(io.BytesIO(tsv), sep="\t", keep_default_na=False)
        eligible = source_frame[
            source_frame.species.eq("HomoSapiens") & source_frame.gene.eq("TRB") &
            source_frame["mhc.class"].eq("MHCI") &
            source_frame["mhc.a"].str.fullmatch(r"HLA-[ABC]\*\d{2}:\d{2}") &
            source_frame.cdr3.str.fullmatch(r"C[ACDEFGHIKLMNPQRSTVWY]*[FW]") &
            source_frame["antigen.epitope"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{8,15}")]
        known = set(eligible[["cdr3", "antigen.epitope", "mhc.a"]].itertuples(index=False, name=None))
        witness_metadata = json.loads((witness_dir / "verification.json").read_text())
        total = 0
        for entry in witness_metadata["records"]:
            witness = pd.read_csv(witness_dir / entry["witness_file"])
            receiver = witness.receiver_source_row_0based.tolist()
            donor = witness.actual_donor_source_row_0based.tolist()
            assert len(receiver) == len(set(receiver))
            assert Counter(receiver) == Counter(donor)
            assert all(a != b for a, b in zip(receiver, donor))
            pos = source_frame.iloc[receiver]
            neg_antigen = source_frame.iloc[donor]
            assert pos["mhc.a"].tolist() == neg_antigen["mhc.a"].tolist()
            negative = list(zip(pos.cdr3, neg_antigen["antigen.epitope"], neg_antigen["mhc.a"]))
            assert not any(pair in known for pair in negative)
            representative = pd.read_csv(root / entry["reference_manifest"])
            if entry["fold_filter"] is not None:
                representative = representative[representative.fold.eq(entry["fold_filter"])]
            assert representative.receiver_source_row_0based.tolist() == receiver
            ref_antigen = source_frame.iloc[representative.negative_antigen_source_row_0based]
            assert negative == list(zip(pos.cdr3, ref_antigen["antigen.epitope"], ref_antigen["mhc.a"]))
            total += len(witness)
        assert len(witness_metadata["records"]) == 58 and total == 33839
        report["independent_row_witness_reconstruction"] = {"witnesses": 58, "assignments": total, "status": "PASS"}
    sys.path.insert(0, str(root / "audit_work"))
    from run_source_audit import source_metadata
    provenance = source_metadata(root / "precision_repo")
    assert provenance["commit"] == "0ba2d538065e2d3a8b4642547eb77447ecc89a5b"
    if (root / "precision_repo/.git").exists():
        for name, expected in provenance["files"].items():
            blob = subprocess.check_output(["git", "-C", str(root / "precision_repo"),
                                            "show", f"{provenance['commit']}:{name}"])
            assert hashlib.sha256(blob).hexdigest() == expected
    with tempfile.TemporaryDirectory(prefix="tcr-fresh-review-") as temp:
        work = Path(temp)
        repo = work / "source-without-git"
        repo.mkdir()
        for name in provenance["files"]:
            shutil.copy2(root / "precision_repo" / name, repo / name)
        (repo / "SOURCE_PROVENANCE.json").write_text(json.dumps(provenance))
        assert source_metadata(repo)["files"] == provenance["files"]
        original = (repo / "dataset.py").read_bytes()
        (repo / "dataset.py").write_bytes(original + b"\n# fresh reviewer tamper probe\n")
        try:
            source_metadata(repo)
        except ValueError:
            pass
        else:
            raise AssertionError("Mutated source was accepted")
        (repo / "dataset.py").write_bytes(original)
        report["git_free_source_copy_and_tamper_rejection"] = "PASS"
        execute(root, "audit_work/run_source_audit.py", "--repo", repo,
                "--out", work / "synthetic", "--seeds", 2)
        fresh = pd.read_csv(work / "synthetic/source_audit_seed_results.csv", float_precision="round_trip")
        saved = pd.read_csv(root / "audit_work/source_audit_seed_results.csv", float_precision="round_trip")
        assert fresh.equals(saved[saved.seed.isin([0, 1])].reset_index(drop=True))
        report["synthetic_two_seed_exact_replay"] = {"rows": len(fresh), "status": "PASS"}
        source = work / "vdjdb.txt"
        source.write_bytes(tsv)
        out = work / "split-first"
        execute(root, "crossfit_audit/run_split_first.py", "--repo", repo,
                "--source", source, "--sampling-seeds", 1, "--fold-seeds", 17, "--out", out)
        baseline = root / "crossfit_audit/split_first_results"
        report["fresh_construction_and_retraining"] = {}
        for name in ["construction_metrics.csv", "fold_diagnostics.csv", "oof_metrics.csv",
                     "oof_scores.csv.gz", "retention.csv"]:
            actual = pd.read_csv(out / name, float_precision="round_trip")
            expected = pd.read_csv(baseline / name, float_precision="round_trip")
            if "fold_seed" in expected:
                expected = expected[expected.fold_seed.eq(17)]
            if "sampling_seed" in expected:
                expected = expected[expected.sampling_seed.eq(0)]
            assert actual.equals(expected.reset_index(drop=True)), name
            report["fresh_construction_and_retraining"][name] = {"rows": len(actual), "exact_dataframe_match": True}
        manifests = list((out / "manifests").glob("*.csv"))
        for path in manifests:
            assert path.read_bytes() == (baseline / "manifests" / path.name).read_bytes()
        report["fresh_manifest_files_exact_byte_match"] = len(manifests)
        figure_root = work / "render"
        inputs = ["audit_work/source_audit_results.json", "audit_work/reference_rejection_seed_results.csv",
                  "public_audit/filter_flow.csv", "public_audit/aggregate_results.csv", "public_audit/seed_results.csv",
                  "public_audit/matching_retention_by_hla.csv", "public_audit/split_sensitivity/split_provenance.json",
                  "crossfit_audit/split_first_results/summary.csv", "crossfit_audit/results/summary.csv",
                  "research/results/geometry_control.csv", "paper/analysis/build_tables.py", "paper/analysis/draw_figures.py"]
        for name in inputs:
            target = figure_root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / name, target)
        execute(figure_root, "paper/analysis/build_tables.py")
        execute(figure_root, "paper/analysis/draw_figures.py")
        report["regenerated_artifacts_exact_byte_match"] = {}
        for folder, pattern in [("tables", "*.tex"), ("figures", "*.png")]:
            paths = list((figure_root / "paper" / folder).glob(pattern))
            for path in paths:
                assert path.read_bytes() == (root / "paper" / folder / path.name).read_bytes(), path.name
            report["regenerated_artifacts_exact_byte_match"][folder] = len(paths)
    report["input_sha256"] = {name: digest(root / name) for name in [
        "paper/manuscript.tex", "paper/supplementary.tex", "run_all.py",
        "audit_work/run_source_audit.py", "crossfit_audit/run_crossfit.py",
        "crossfit_audit/run_split_first.py", "public_audit/run_public_audit.py",
        "paper/analysis/draw_figures.py", "paper/analysis/build_tables.py"]}
    report["status"] = "PASS"
    target = args.output or root / "reviews/third_round_reproducibility_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
