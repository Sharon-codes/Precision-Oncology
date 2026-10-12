#!/usr/bin/env python3
"""Reproduce the archived analysis from untouched source data and pinned code.

Default rebuilds every executed main experiment sequentially, writing explicit
per-stage logs. --verify-only recomputes archived numerical summaries from
stored predictions without data extraction, network use, or model fitting.
"""
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import zipfile
from restore_large_artifacts import restore

ROOT=Path(__file__).resolve().parent
TSV_SHA='2823281ece02f37b15cdd9187be374e99ee7b46e041572fc2cde4b295a9bf52e'
ZIP_SHA='7c4994ce32eba8e0ff6c23d1ec53b85a982be69651ed1f467ffe64c86c66c1d1'
def verify_hash(path,expected):
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual!=expected:raise ValueError(f'Unexpected source hash for {path.name}: {actual}')
def run(name,relative,*args):
    logs=ROOT/'reproduction_logs';logs.mkdir(exist_ok=True)
    print(f'Running {name}',flush=True)
    with (logs/f'{name}.log').open('w') as out:
        subprocess.run([sys.executable,str(ROOT/relative),*map(str,args)],cwd=ROOT,
                       check=True,stdout=out,stderr=subprocess.STDOUT)
    print(f'Passed {name}',flush=True)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify-only',action='store_true');ap.add_argument('--compile-pdfs',action='store_true');args=ap.parse_args()
    restore(ROOT)
    if not args.verify_only:
        folder=ROOT/'research/data';archive=folder/'vdjdb-2023-06-01.zip'
        verify_hash(archive,ZIP_SHA)
        with zipfile.ZipFile(archive) as z:
            # Explicit member selection avoids arbitrary path extraction.
            for name in ['vdjdb.txt','vdjdb.meta.txt','LICENSE.txt']:
                (folder/name).write_bytes(z.read(name))
        verify_hash(folder/'vdjdb.txt',TSV_SHA)
        stages=[
          ('source_synthetic','audit_work/run_source_audit.py',[]),
          ('public_grid','public_audit/run_public_audit.py',[]),
          ('reference_sensitivity','public_audit/run_reference_sensitivity.py',[]),
          ('distance_sensitivity','public_audit/run_split_sensitivity.py',[]),
          ('annotation_source','public_audit/run_annotation_sensitivity.py',[]),
          ('manifest_validation','public_audit/validate_audit.py',[]),
          ('closed_matching_exhaustive','audit_work/matching_validation.py',[]),
          ('rejection_reference','audit_work/reference_rejection.py',[]),
          ('metadata_positive_control','crossfit_audit/validate_synthetic_control.py',[]),
          ('crossfit_global','crossfit_audit/run_crossfit.py',['--manifest',ROOT/'crossfit_audit/input_manifest.json','--output',ROOT/'crossfit_audit/results']),
          ('crossfit_split_first','crossfit_audit/run_split_first.py',[]),
          ('crossfit_annotation','crossfit_audit/run_annotation_filtered.py',['--out',ROOT/'crossfit_audit/annotation_filtered_results']),
          ('annotation_assignment','reviews/run_annotation_sensitivity.py',[]),
          ('official_nomenclature','reviews/validate_annotation_evidence.py',[]),
          ('assignment_witnesses','audit_work/export_assignment_witnesses.py',[]),
          ('geometry_counterexample','research/geometry_control.py',[]),
          ('tables','paper/analysis/build_tables.py',[]),
          ('figures','paper/analysis/draw_figures.py',[]),
        ]
        for name,path,params in stages:run(name,path,*params)
    run('archived_result_verification','audit_work/verify_archived_results.py','--root',ROOT)
    if args.compile_pdfs:
        for folder,names in [('paper',['manuscript','supplementary']),('paper/author_only',['review_report','title_page'])]:
            for name in names:
                subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error',name+'.tex'],cwd=ROOT/folder,check=True)
    print('Finished. See reproduction_logs and audit_work/archive_verification_results.json.',flush=True)
if __name__=='__main__':main()
