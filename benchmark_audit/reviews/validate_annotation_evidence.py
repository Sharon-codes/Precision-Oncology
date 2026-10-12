#!/usr/bin/env python3
"""Independent nomenclature and source-row reconstruction checks."""
import hashlib
import json
import re
from pathlib import Path
from collections import Counter
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'reviews'
FLAGGED = {'HLA-A*08:01', 'HLA-A*24:01'}
COMMIT = '38398a75e9762ff070d8e9bd714d074332646cd7'


def main():
    nom = OUT / 'IPD_HLA_3_52_0_hla_nom.txt'
    deleted = OUT / 'IPD_HLA_3_52_0_deleted_alleles.txt'
    # Official WMDA fields: locus;allele name;assignment date;deletion date;...
    present = set()
    for line in nom.read_text().splitlines():
        fields = line.split(';')
        if len(fields) < 4 or fields[0] not in {'A*', 'B*', 'C*'} or fields[3]:
            continue
        m = re.match(r'^(\d{2,3}:\d{2,3})(?::|[A-Z]|$)', fields[1])
        if m:
            present.add('HLA-' + fields[0] + m.group(1))
    strata = pd.read_csv(ROOT / 'public_audit/matching_retention_by_hla.csv')
    checks = [{'source_hla_label': row.hla, 'has_active_two_field_prefix': row.hla in present,
               'primary_n': int(row.n_positive), 'retained_n': int(row.n_retained)}
              for row in strata.itertuples()]
    assert {c['source_hla_label'] for c in checks if not c['has_active_two_field_prefix']} == FLAGGED
    source = pd.read_csv(ROOT / 'research/data/vdjdb.txt', sep='\t', keep_default_na=False)
    result = {
        'nomenclature_release': 'IPD-IMGT/HLA 3.52.0, file date 2023-04-17',
        'official_git_commit': COMMIT, 'official_git_tag': 'v3.52.0-alpha',
        'scope': 'Presence of an active official allele with the same two-field prefix; does not validate the biological annotation of a source assay.',
        'nomenclature_sha256': hashlib.sha256(nom.read_bytes()).hexdigest(),
        'deleted_alleles_sha256': hashlib.sha256(deleted.read_bytes()).hexdigest(),
        'official_nomenclature_url': f'https://raw.githubusercontent.com/ANHIG/IMGTHLA/{COMMIT}/wmda/hla_nom.txt',
        'official_deleted_alleles_url': f'https://raw.githubusercontent.com/ANHIG/IMGTHLA/{COMMIT}/Deleted_alleles.txt',
        'legacy_A2401_deletion_record': next(line for line in deleted.read_text().splitlines() if ',A*2401,' in line),
        'checks': checks, 'manifest_reconstructions': [],
    }
    folder = ROOT / 'public_audit/annotation_sensitivity'
    for path in sorted((folder / 'manifests').glob('*.csv')):
        manifest = pd.read_csv(path)
        pairs = pd.read_csv(folder / 'local_pairs' / path.name)
        receivers = source.iloc[manifest.receiver_source_row_0based.to_numpy()]
        donors = source.iloc[manifest.negative_antigen_source_row_0based.to_numpy()]
        positives = pairs[pairs.label.eq(1)]
        negatives = pairs[pairs.label.eq(0)]
        n = len(manifest)
        assert len(positives) == len(negatives) == n
        assert not (manifest.negative_antigen_source_row_0based < 0).any()
        assert positives.beta.tolist() == receivers.cdr3.tolist() == negatives.beta.tolist()
        assert positives.peptide.tolist() == receivers['antigen.epitope'].tolist()
        assert positives.hla.tolist() == receivers['mhc.a'].tolist()
        assert negatives.peptide.tolist() == donors['antigen.epitope'].tolist()
        assert negatives.hla.tolist() == donors['mhc.a'].tolist()
        assert not set(pairs.hla) & FLAGGED
        rec = {'filename': path.name, 'positive_n': n, 'negative_n': n,
               'source_row_reconstruction_passed': True, 'excluded_labels_absent': True,
               'manifest_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        if path.name == 'clean_corrected_seed0.csv':
            for cols in [('beta',), ('peptide',), ('hla',), ('peptide', 'hla')]:
                assert Counter(positives[list(cols)].itertuples(index=False, name=None)) == Counter(negatives[list(cols)].itertuples(index=False, name=None))
            rec['all_declared_marginals_identical'] = True
            raw = pd.read_csv(ROOT / 'public_audit/local_pairs/common_corrected_seed0.csv')
            raw = raw[~raw.hla.isin(FLAGGED)].reset_index(drop=True)
            assert pairs.equals(raw)
            rec['exact_pairs_equal_original_minus_excluded_strata'] = True
        result['manifest_reconstructions'].append(rec)
    target = OUT / 'official_annotation_validation.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'n_hla_labels_checked': len(checks), 'unmatched': [c for c in checks if not c['has_active_two_field_prefix']], 'manifests_passed': len(result['manifest_reconstructions']), 'output': str(target)}, indent=2))


if __name__ == '__main__':
    main()
