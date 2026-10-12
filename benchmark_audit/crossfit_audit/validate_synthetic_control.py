#!/usr/bin/env python3
"""Execute an artificial perfect-peptide signal through the actual CV pipeline.

All strings in this fixture are artificial category identifiers. This tests
whether the implementation learns an available signal and preserves its null
controls; it says nothing about biological recognition or generalization.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import sklearn

from run_crossfit import run_dataset, validate_dataset


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=here / "synthetic_control_results")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for i in range(100):
        for label in [1, 0]:
            rows.append({"tcr_beta": f"SYNTHETIC_RECEPTOR_{i:03d}",
                         "peptide": "POSITIVE_CATEGORY" if label else "NEGATIVE_CATEGORY",
                         "hla": "CONSTANT_SYNTHETIC_HLA", "label": label,
                         "example_id": f"synthetic_{i:03d}_y{label}"})
    frame = pd.DataFrame(rows)
    validate_dataset(frame)
    desc = {"sampler": "artificial_perfect_peptide_signal", "sampling_seed": 0,
            "cohort": "synthetic_100_receptors_two_labels"}
    metrics, folds, scores = run_dataset(frame, desc, args.out, [17, 41])
    metric_frame = pd.DataFrame(metrics)
    fold_frame = pd.DataFrame(folds)
    score_frame = pd.concat(scores, ignore_index=True)
    expected_auc = {"peptide": 1.0, "peptide_hla": 1.0, "hla": 0.5, "receptor": 0.5}
    for feature, expected in expected_auc.items():
        assert metric_frame.loc[metric_frame.feature.eq(feature), "auroc"].eq(expected).all()
        assert fold_frame.loc[fold_frame.feature.eq(feature), "auroc"].eq(expected).all()
    for feature in ["hla", "receptor"]:
        assert np.array_equal(score_frame[f"score_{feature}"], np.full(len(score_frame), 0.5))
    for seed, block in score_frame.groupby("fold_seed"):
        assert block.groupby("receptor_sha256").fold.nunique().eq(1).all()
        assert len(block) == len(frame)
    frame.to_csv(args.out / "synthetic_fixture.csv", index=False)
    metric_frame.to_csv(args.out / "oof_metrics.csv", index=False)
    fold_frame.to_csv(args.out / "fold_diagnostics.csv", index=False)
    score_frame.to_csv(args.out / "oof_scores.csv.gz", index=False,
                       compression={"method": "gzip", "mtime": 0})
    report = {"status": "PASS", "scope": "Entirely artificial implementation check; not biological evidence.",
              "design": {"receptor_categories": 100, "n_pairs": 200,
                         "positive_rows_per_receptor": 1, "negative_rows_per_receptor": 1,
                         "signal": "Peptide category deterministically identifies label; HLA constant.",
                         "fold_seeds": [17, 41], "folds_per_seed": 5,
                         "actual_pipeline": "run_crossfit.run_dataset"},
              "expected_and_observed_auroc_in_every_run_and_fold": expected_auc,
              "null_scores_exactly_0_5": True,
              "receptor_groups_do_not_cross_folds": True,
              "pooled_metric_sets": len(metric_frame), "fold_metric_sets": len(fold_frame),
              "score_rows": len(score_frame),
              "versions": {"python": platform.python_version(), "numpy": np.__version__,
                           "pandas": pd.__version__, "sklearn": sklearn.__version__},
              "source_sha256": {name: hashlib.sha256((here / name).read_bytes()).hexdigest()
                                for name in ["run_crossfit.py", "validate_synthetic_control.py"]}}
    (args.out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
