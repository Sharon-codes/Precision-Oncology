"""Cross-fitted metadata controls for frozen constructed TCR pair datasets.

This is a within-antigen-collection construction diagnostic, NOT an epitope-OOD
binding benchmark. Receptor groups are held out, but peptide/HLA categories may
appear in train and test. No receptor sequence is given to the antigen models.

Input JSON manifest contains a ``datasets`` list; every item has ``path``,
``sampler``, ``sampling_seed`` and optional ``cohort`` (default common_retained).
CSV required columns: tcr_beta, peptide, hla, label. Optional ``example_id``.
All files in the same cohort must have exactly identical positive multisets.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import warnings

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.preprocessing import OneHotEncoder

FEATURES = {
    "peptide_hla": ["peptide", "hla"],
    "peptide": ["peptide"],
    "hla": ["hla"],
    "receptor": ["tcr_beta"],
}


def fold_map(groups, seed, nfolds=5):
    unique = np.array(sorted(set(groups)), dtype=object)
    if len(unique) < nfolds:
        raise ValueError("Fewer receptor groups than folds")
    np.random.default_rng(seed).shuffle(unique)
    return {str(g): f for f, gg in enumerate(np.array_split(unique, nfolds)) for g in gg}


def positive_signature(df):
    return Counter(map(tuple, df.loc[df.label.eq(1), ["tcr_beta", "peptide", "hla"]].values))


def validate_dataset(df):
    required = ["tcr_beta", "peptide", "hla", "label"]
    if set(required) - set(df):
        raise ValueError(f"Missing columns: {set(required)-set(df)}")
    if df[required].isna().any().any():
        raise ValueError("Missing required values")
    if set(df.label) != {0, 1}:
        raise ValueError("Labels must include exactly 0 and 1")
    counts = pd.crosstab(df.tcr_beta, df.label)
    if not counts[0].equals(counts[1]):
        raise ValueError("Each receptor must have balanced positive and decoy multiplicity")
    if not df.example_id.is_unique:
        raise ValueError("Duplicate example IDs")


def run_dataset(df, desc, out_dir, fold_seeds):
    metrics, fold_rows, score_tables = [], [], []
    for fold_seed in fold_seeds:
        if "preassigned_fold" in df:
            if df.groupby("tcr_beta").preassigned_fold.nunique().max() != 1:
                raise AssertionError("One receptor assigned to multiple fixed folds")
            folds = df.preassigned_fold.to_numpy(dtype=int)
            if set(folds) != set(range(5)):
                raise AssertionError("Fixed folds must be 0..4")
        else:
            mapping = fold_map(df.tcr_beta, fold_seed)
            folds = df.tcr_beta.map(mapping).to_numpy()
        # Sequence hashes allow score joins without redistributing source strings.
        receptor_hash = df.tcr_beta.map(lambda x: hashlib.sha256(x.encode()).hexdigest()).to_numpy()
        table = pd.DataFrame({
            "example_id": df.example_id, "label": df.label,
            "receptor_sha256": receptor_hash, "fold": folds,
        })
        for feature_name, columns in FEATURES.items():
            oof = np.full(len(df), np.nan)
            for fold in range(5):
                train = folds != fold
                test = ~train
                if set(df.loc[train, "tcr_beta"]) & set(df.loc[test, "tcr_beta"]):
                    raise AssertionError("Receptor group overlap")
                # Fit every preprocessing transform exclusively within this fold.
                encoder = OneHotEncoder(handle_unknown="ignore", dtype=np.float64)
                x_train = encoder.fit_transform(df.loc[train, columns])
                x_test = encoder.transform(df.loc[test, columns])
                y_train = df.loc[train, "label"].to_numpy()
                y_test = df.loc[test, "label"].to_numpy()
                clf = LogisticRegression(C=1.0, solver="liblinear", tol=1e-8,
                                         max_iter=2000, random_state=1729)
                with warnings.catch_warnings():
                    warnings.simplefilter("error", ConvergenceWarning)
                    clf.fit(x_train, y_train)
                pred = clf.predict_proba(x_test)[:, int(np.where(clf.classes_ == 1)[0][0])]
                oof[test] = pred
                row = {
                    **desc, "fold_seed": fold_seed, "fold": fold,
                    "feature": feature_name, "train_pairs": int(train.sum()),
                    "test_pairs": int(test.sum()),
                    "train_groups": int(df.loc[train, "tcr_beta"].nunique()),
                    "test_groups": int(df.loc[test, "tcr_beta"].nunique()),
                    "auroc": float(roc_auc_score(y_test, pred)),
                    "ap": float(average_precision_score(y_test, pred)),
                }
                for col, cats in zip(columns, encoder.categories_):
                    row[f"unseen_{col}_fraction"] = float((~df.loc[test,col].isin(cats)).mean())
                fold_rows.append(row)
            if np.isnan(oof).any():
                raise AssertionError("OOF score coverage incomplete")
            # Balanced receptor training groups plus all-unknown test receptors
            # must yield the same intercept-only prediction across all rows.
            if feature_name == "receptor" and not np.allclose(oof, 0.5, atol=1e-8, rtol=0):
                raise AssertionError(f"Receptor-null control not 0.5: {oof.min()}, {oof.max()}")
            metrics.append({
                **desc, "fold_seed": fold_seed, "feature": feature_name,
                "n_pairs": len(df), "n_positive": int(df.label.sum()),
                "n_receptor_groups": int(df.tcr_beta.nunique()),
                "auroc": float(roc_auc_score(df.label, oof)),
                "average_precision": float(average_precision_score(df.label, oof)),
                "brier": float(brier_score_loss(df.label, oof)),
                "score_min": float(oof.min()), "score_max": float(oof.max()),
            })
            table[f"score_{feature_name}"] = oof
        table["sampler"] = desc["sampler"]
        table["sampling_seed"] = desc["sampling_seed"]
        table["fold_seed"] = fold_seed
        table["cohort"] = desc["cohort"]
        score_tables.append(table)
    return metrics, fold_rows, score_tables


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fold-seeds", default="17,41,73,107,211")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(args.manifest.read_text())
    fold_seeds = [int(x) for x in args.fold_seeds.split(",")]
    signatures, all_metrics, all_folds, all_scores, inputs = {}, [], [], [], []
    for item in manifest["datasets"]:
        path = Path(item["path"])
        if not path.is_absolute():
            path = args.manifest.parent / path
        df = pd.read_csv(path, dtype={"tcr_beta":str,"beta":str,"peptide":str,"hla":str})
        if "tcr_beta" not in df and "beta" in df:
            df = df.rename(columns={"beta":"tcr_beta"})
        if "example_id" not in df:
            if "source_row_0based" in df:
                df["example_id"] = df.source_row_0based.astype(str) + "_y" + df.label.astype(str)
            else:
                df["example_id"] = [f"e{i:07d}" for i in range(len(df))]
        validate_dataset(df)
        desc = {"sampler":item["sampler"],"sampling_seed":int(item["sampling_seed"]),
                "cohort":item.get("cohort","common_retained")}
        sig = positive_signature(df)
        if desc["cohort"] in signatures and sig != signatures[desc["cohort"]]:
            raise AssertionError(f"Positive cohort mismatch in {path}")
        signatures[desc["cohort"]] = sig
        metrics, folds, scores = run_dataset(df, desc, args.output, fold_seeds)
        all_metrics.extend(metrics)
        all_folds.extend(folds)
        all_scores.extend(scores)
        inputs.append({**desc, "file":path.name,
                       "sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
        print(desc, {r["feature"]:round(r["auroc"],4) for r in metrics if r["fold_seed"]==fold_seeds[0]}, flush=True)
    met = pd.DataFrame(all_metrics)
    met.to_csv(args.output/"oof_metrics.csv", index=False)
    pd.DataFrame(all_folds).to_csv(args.output/"fold_diagnostics.csv", index=False)
    pd.concat(all_scores,ignore_index=True).to_csv(args.output/"oof_scores.csv.gz",index=False,compression="gzip")
    seed_means=met.groupby(["cohort","sampler","feature","sampling_seed"],as_index=False).agg(
        auroc=("auroc","mean"),average_precision=("average_precision","mean"),brier=("brier","mean"))
    seed_means.to_csv(args.output/"sampling_seed_means.csv",index=False)
    summary=seed_means.groupby(["cohort","sampler","feature"],as_index=False).agg(
        sampling_seeds=("sampling_seed","nunique"),auroc_mean=("auroc","mean"),
        auroc_sd_between_sampler_seeds=("auroc","std"),auroc_min=("auroc","min"),
        auroc_max=("auroc","max"),ap_mean=("average_precision","mean"),brier_mean=("brier","mean"))
    summary.to_csv(args.output/"summary.csv",index=False)
    metadata={
        "purpose":"within-antigen-collection label-construction negative control; not epitope-OOD prediction",
        "version":{"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,
                   "scipy":scipy.__version__,"sklearn":sklearn.__version__},
        "fold_seeds":fold_seeds,"n_folds":5,"group":"exact cdr3_beta string",
        "model":{"name":"LogisticRegression","C":1.0,"solver":"liblinear","tol":1e-8,"max_iter":2000},
        "features":FEATURES,"inputs":inputs,
        "uncertainty":"SD and range of fold-averaged metrics across negative-sampling seeds; not biological population confidence intervals",
        "preprocessing":"OneHotEncoder fitted on training fold only; unseen categories ignored",
        "split_rule":"shuffle sorted receptor groups by seed, split into five group blocks; shared across same-cohort samplers",
    }
    (args.output/"run_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(summary.to_string(index=False),flush=True)


if __name__=="__main__":
    main()
