"""Stage publication artifacts without local sequence materializations or Git history."""
from pathlib import Path
import hashlib,json,shutil,zipfile,tempfile

ROOT=Path(__file__).resolve().parent
STAGE=ROOT/'release/TCR_Benchmark_Audit'
OUT=ROOT/'deliverables'
def copy_file(src,dest):
    dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
def stage():
    global STAGE
    destination=STAGE
    destination.parent.mkdir(parents=True,exist_ok=True)
    STAGE=Path(tempfile.mkdtemp(prefix='TCR_Benchmark_Audit_stage_',dir=destination.parent))
    OUT.mkdir(exist_ok=True)
    for name in ['run_all.py','restore_large_artifacts.py','requirements_reproducibility.txt','package_artifacts.py']:
        copy_file(ROOT/name,STAGE/name)
    copy_file(ROOT/'REPRODUCIBILITY_README.md',STAGE/'README.md')
    copy_file(ROOT/'REPRODUCIBILITY_README.md',STAGE/'REPRODUCIBILITY_README.md')
    # Scientific results, portable code and reviewer records. Private local
    # sequence materializations and transient build files are not distributed.
    for folder in ['audit_work','public_audit','crossfit_audit','reviews']:
        for src in (ROOT/folder).rglob('*'):
            if not src.is_file() or any(x in src.parts for x in ['local_pairs','__pycache__']):continue
            if src.suffix not in {'.py','.json','.csv','.gz','.md','.bib','.txt','.tex','.pdf','.png'}:continue
            if 'tmp_statistics' in src.parts and src.suffix!='.json':continue
            if src.name in ['source.txt','source_clean.txt','journal_special_requirements.png','journal_special_requirements.pdf','source_audit_stdout.json']:continue
            copy_file(src,STAGE/src.relative_to(ROOT))
    for src in (ROOT/'paper').rglob('*'):
        if not src.is_file() or any(x in src.parts for x in ['qa','__pycache__']):continue
        if src.suffix not in {'.py','.tex','.bib','.bbl','.cls','.bst','.md','.pdf','.eps','.png'}:continue
        if src.name=='methods_draft.txt':continue
        copy_file(src,STAGE/src.relative_to(ROOT))
    for name in ['geometry_control.py','results/geometry_control.csv','results/geometry_control.json',
                 'data/vdjdb-2023-06-01.zip','data/VDJdb_LICENSE.txt']:
        copy_file(ROOT/'research'/name,STAGE/'research'/name)
    for src in (ROOT/'precision_repo').rglob('*'):
        if not src.is_file() or '.git' in src.parts or '__pycache__' in src.parts:continue
        if src.suffix not in {'.py','.txt','.md','.json','.sh'} and src.name!='LICENSE':continue
        copy_file(src,STAGE/'precision_repo'/src.relative_to(ROOT/'precision_repo'))
    provenance=json.loads((ROOT/'public_audit/provenance.json').read_text())['source_code']
    (STAGE/'precision_repo/SOURCE_PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
    # The archive carries immutable reference bytes with their license; a
    # validated provenance fallback supplies metadata without Git internals.
    manifest=[]
    for p in sorted(STAGE.rglob('*')):
        if p.is_file() and p.name!='SHA256SUMS.txt':manifest.append(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(STAGE)))
    (STAGE/'SHA256SUMS.txt').write_text('\n'.join(manifest)+'\n')
    # Only this task-owned generated destination is replaced. Rebuilding must
    # never retain locally regenerated sequence tables or stale deleted files.
    if destination.exists():shutil.rmtree(destination)
    STAGE.rename(destination);STAGE=destination
    return manifest
def archives():
    for source,target in [('manuscript.pdf','TCR_Benchmark_Audit.pdf'),('supplementary.pdf','TCR_Supplementary.pdf'),
                          ('author_only/review_report.pdf','Reviewer_Debate_and_Evidence.pdf')]:
        copy_file(ROOT/'paper'/source,OUT/target)
    with zipfile.ZipFile(OUT/'TCR_LaTeX_Source.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted((STAGE/'paper').rglob('*')):
            if p.is_file():z.write(p,Path('TCR_LaTeX_Source')/p.relative_to(STAGE/'paper'))
    with zipfile.ZipFile(OUT/'TCR_Complete_Research_Package.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(STAGE.rglob('*')):
            if p.is_file():z.write(p,Path('TCR_Benchmark_Audit')/p.relative_to(STAGE))
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--stage-only',action='store_true');args=ap.parse_args()
    files=stage()
    if not args.stage_only:archives()
    print(json.dumps({'staged_files':len(files),'staged_root':str(STAGE),'archives_created':not args.stage_only},indent=2))
