"""Prespecified annotation sensitivity: exclude two flagged HLA source labels.

The original raw-annotation analyses are retained. This wrapper reruns the
identical split-before-construction control after excluding HLA-A*08:01 and
HLA-A*24:01 from all cohorts AND from the known-positive reference. No allele
is recoded. Receptor folds are regenerated, rather than assumed unchanged.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import run_split_first as experiment

EXCLUDED = {"HLA-A*08:01", "HLA-A*24:01"}
original_prepare = experiment.prepare
filter_counts = {}


def annotation_filtered_prepare(source):
    cohorts, known, known_pair, flow = original_prepare(source)
    for name, df in cohorts.items():
        filtered = df.loc[~df.hla.isin(EXCLUDED)].reset_index(drop=True)
        filter_counts[name] = {"before":len(df), "after":len(filtered),
                              "excluded":len(df)-len(filtered),
                              "hla_strata_after":int(filtered.hla.nunique())}
        cohorts[name] = filtered
    old_n = len(known)
    known = {triple for triple in known if triple[2] not in EXCLUDED}
    known_pair = {triple[:2] for triple in known}
    filter_counts["reference"] = {"before":old_n, "after":len(known), "excluded":old_n-len(known)}
    flow.append({"filter":"exclude_flagged_HLA_source_labels", "labels":sorted(EXCLUDED),
                 "primary_n":len(cohorts["score_ge_1"])})
    return cohorts, known, known_pair, flow


if __name__ == "__main__":
    # Resolve this wrapper's own safe default before delegating to the shared
    # experiment parser, whose standalone default belongs to the primary run.
    # Appending the resolved value also preserves --out=PATH/abbreviated forms.
    wrapper_parser = argparse.ArgumentParser(add_help=False)
    wrapper_parser.add_argument("--out", type=Path, default=experiment.HERE/"annotation_filtered_results")
    wrapper_args, _ = wrapper_parser.parse_known_args()
    output = wrapper_args.out
    sys.argv.extend(["--out", str(output)])
    experiment.prepare = annotation_filtered_prepare
    experiment.main()
    meta_path=output/"run_metadata.json"
    meta=json.loads(meta_path.read_text())
    meta["annotation_sensitivity"]={
        "excluded_source_labels":sorted(EXCLUDED),
        "policy":"Remove rows from every cohort and exclusion reference, without recoding; regenerate receptor folds",
        "counts":filter_counts,
        "comparison_scope":"Sensitivity to two flagged source labels, not comprehensive allele nomenclature validation",
    }
    meta["source_code_sha256"]={name:hashlib.sha256((experiment.HERE/name).read_bytes()).hexdigest()
        for name in ["run_crossfit.py","run_split_first.py","run_annotation_filtered.py"]}
    meta_path.write_text(json.dumps(meta,indent=2)+"\n")
    print(json.dumps(meta["annotation_sensitivity"],indent=2),flush=True)
