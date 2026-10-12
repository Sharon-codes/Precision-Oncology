#!/usr/bin/env python3
"""Exclude flagged raw HLA annotations; rerun only closed-assignment construction.

This is an annotation-exclusion sensitivity, not allele recoding or a rerun of
the metadata cross-validation, source collator, or distance partitions.
"""
from pathlib import Path
import hashlib
import json
import sys
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'public_audit'))
from run_public_audit import prepare, corrected, metrics


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'), ensure_ascii=True).encode()).hexdigest()


def hla_tv(original, retained):
    p, q = Counter(original.hla), Counter(retained.hla)
    return .5 * sum(abs(p[h] / len(original) - q[h] / len(retained)) for h in p.keys() | q.keys())


def main():
    source = ROOT / 'research/data/vdjdb.txt'
    cohorts, known, known_pair, _ = prepare(source)
    primary = cohorts['score_ge_1']
    base_retained, base_negative, _ = corrected(primary, known)
    base_rows = list(map(int, base_retained.source_row_0based))
    base_manifest = list(zip(base_rows, base_negative))
    result = {
        'scope': __doc__,
        'source_file_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'construction_code_sha256': hashlib.sha256((ROOT / 'public_audit/run_public_audit.py').read_bytes()).hexdigest(),
        'raw_annotations_not_reassigned': True,
        'official_deleted_alleles_url': 'https://github.com/ANHIG/IMGTHLA/blob/Latest/Deleted_alleles.txt',
        'deleted_allele_observation': 'Registry line HLA00049,A*2401,Sequence shown to contain errors (March 1995). Source HLA-A*24:01 is the colon-form rendering of that legacy annotation.',
        'A08_scope': 'Flagged unresolved raw source annotation: absent from VDJdb bundled HLA nomenclature table and explicitly corrected by later VDJdb curation. No reassignment applied here.',
        'baseline': {'primary_rows': len(primary), 'registry_triples': len(known), 'retained_rows': len(base_retained), 'manifest_sha256': digest(base_manifest)},
        'flagged_label_counts': {},
        'scenarios': {},
        'not_rerun': ['source-collator seed experiments', 'metadata cross-validation', 'peptide-distance partitions'],
    }
    for h in ['HLA-A*08:01', 'HLA-A*24:01']:
        result['flagged_label_counts'][h] = {
            'by_confidence_cohort': {name: int(d.hla.eq(h).sum()) for name, d in cohorts.items()},
            'retained_baseline_rows': int(base_retained.hla.eq(h).sum()),
            'primary_source_rows_0based': list(map(int, primary.loc[primary.hla.eq(h), 'source_row_0based'])),
        }
    for name, excluded in [('exclude_A08', {'HLA-A*08:01'}), ('exclude_A08_and_legacy_A24', {'HLA-A*08:01', 'HLA-A*24:01'})]:
        d = primary.loc[~primary.hla.isin(excluded)].copy()
        k = {t for t in known if t[2] not in excluded}
        kp = {t[:2] for t in k}
        retained, negative, stats = corrected(d, k)
        rows = list(map(int, retained.source_row_0based))
        manifest = list(zip(rows, negative))
        expected = [(r, t) for r, t in base_manifest if t[2] not in excluded]
        m = metrics(retained, negative, k, kp)
        assert manifest == expected, 'Within-stratum exclusion unexpectedly changed remaining assignments'
        result['scenarios'][name] = {
            'excluded_raw_hla_labels': sorted(excluded),
            'primary_rows': len(d), 'registry_triples': len(k),
            'retained_rows': len(retained), 'retention_fraction': len(retained) / len(d),
            'primary_hla_strata': int(d.hla.nunique()), 'retained_hla_strata': int(retained.hla.nunique()),
            'primary_peptides': int(d.peptide.nunique()), 'retained_peptides': int(retained.peptide.nunique()),
            'hla_distribution_tv_original_vs_retained': hla_tv(d, retained),
            'retained_manifest_sha256': digest(manifest),
            'exact_manifest_equals_unexcluded_baseline': manifest == base_manifest,
            'exact_manifest_equals_baseline_with_excluded_strata_removed': manifest == expected,
            'retained_source_row_ids_sha256': digest(rows),
            'negative_tuples_sha256': digest(negative),
            'construction_metrics': m,
        }
    path = ROOT / 'reviews/annotation_sensitivity.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'baseline': result['baseline'], 'scenarios': result['scenarios']}, indent=2))


if __name__ == '__main__':
    main()
