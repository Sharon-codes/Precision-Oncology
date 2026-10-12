#!/usr/bin/env python3
"""Sensitivity excluding two exact raw HLA annotations without recoding.

Primary analysis retains source annotations and is a raw-annotation audit.
This secondary analysis removes HLA-A*08:01 (unresolved in the nomenclature
review) and HLA-A*24:01 (officially deleted nomenclature). Annotation absence
or deletion is not an inference about binding biology. No allele is imputed.
"""
import hashlib,json
import pandas as pd
from run_public_audit import HERE,prepare,original,corrected,metrics,save_manifest,summarize

EXCLUDED=('HLA-A*08:01','HLA-A*24:01')

def main():
    out=HERE/'annotation_sensitivity'
    for folder in ['manifests','local_pairs']:(out/folder).mkdir(parents=True,exist_ok=True)
    source=HERE.parent/'research/data/vdjdb.txt';cohorts,known,known_pair,flow=prepare(source)
    # Exclude the two exact labels independently in every cohort and in the
    # broad all-confidence known-positive reference. Never recode to a relative.
    clean={name:d[~d.hla.isin(EXCLUDED)].reset_index(drop=True) for name,d in cohorts.items()}
    known={t for t in known if t[2] not in EXCLUDED};known_pair={t[:2] for t in known}
    full=clean['score_ge_1'];retained,neg,stats=corrected(full,known)
    pd.DataFrame(stats).to_csv(out/'matching_retention_by_hla.csv',index=False)
    save_manifest(out,'clean_corrected_seed0',retained,neg,clean['score_ge_0'])
    rows=[dict(cohort='clean_score_ge_1_retained',sampler='exact_closed_matching',order='fixed_manifest',batch_size=0,seed=0,**metrics(retained,neg,known,known_pair))]
    for population,d in [('clean_score_ge_1',full),('clean_score_ge_1_retained',retained)]:
        for order in ['release_order','peptide_sorted','row_shuffled']:
            for seed in range(30):
                ordered,negative,fallback,hardcoded=original(d,HERE.parent/'precision_repo',order,64,seed,known)
                rows.append(dict(cohort=population,sampler='source_collator',order=order,batch_size=64,seed=seed,**metrics(ordered,negative,known,known_pair,fallback,hardcoded)))
                if seed==0:save_manifest(out,f'{population}_{order}_seed0',ordered,negative,clean['score_ge_0'])
            print('DONE',population,order,flush=True)
    frame,agg=summarize(rows);frame.to_csv(out/'seed_results.csv',index=False);agg.to_csv(out/'aggregate_results.csv',index=False)
    exclusions=[]
    for name,d in cohorts.items():
        for label in EXCLUDED:
            e=d[d.hla.eq(label)]
            exclusions.append({'cohort':name,'excluded_raw_label':label,'n_rows':len(e),'source_rows_0based':[int(x) for x in e.source_row_0based]})
    meta={'scope':'Sensitivity of raw-annotation primary audit to removal of two exact HLA labels; no allele recoding or imputation. Official nomenclature status is not used to infer biological binding.','excluded_labels':list(EXCLUDED),'exclusion_applies_to':'Each analysis cohort AND broad score>=0 known-positive exclusion reference.','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'annotation_evidence':'Authoritative nomenclature evidence is recorded by the biological review; HLA-A*08:01 unresolved, HLA-A*24:01 deleted. No replacement allele inferred.','exclusions':exclusions,'raw_main_n':len(cohorts['score_ge_1']),'clean_main_n':len(full),'clean_main_peptides':full.peptide.nunique(),'clean_main_hlas':full.hla.nunique(),'clean_main_betas':full.beta.nunique(),'clean_known_reference_n':len(known),'corrected_retained_n':len(retained),'corrected_retention_fraction':len(retained)/len(full),'corrected_retained_peptides':retained.peptide.nunique(),'corrected_retained_hlas':retained.hla.nunique(),'corrected_retained_betas':retained.beta.nunique(),'corrected_asserted_invariants':['zero broad-reference-positive collisions','exact beta marginal equality','exact peptide marginal equality','exact HLA marginal equality','exact peptide/HLA joint marginal equality'],'construction_seeds':list(range(30)),'batch_size':64,'uncertainty':'Construction RNG variability, not biological population confidence intervals.','identity_oracle':'Resubstitution marginal diagnostic, not trained held-out performance.'}
    meta['annotation_evidence_sources']={
      'official_deleted_alleles_url':'https://github.com/ANHIG/IMGTHLA/blob/Latest/Deleted_alleles.txt',
      'deleted_allele_record':'HLA00049,A*2401; sequence shown to contain errors; March 1995 (official deletion record supplied by independent biological review).',
      'current_nomenclature_local_sha256':'6620fbd3f7238e9f26947bb0353cc361e68f24dc300b1fe3d17bdf79721e86a5',
      'unresolved_label_policy':'HLA-A*08:01 absent from current two-field nomenclature comparison; excluded without applying a later curatorial replacement.'}
    (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps({k:v for k,v in meta.items() if k not in ['exclusions','construction_seeds']},indent=2),flush=True)
    print(agg[agg.metric.isin(['known_positive_collisions_fraction','peptide_tv','peptide_identity_oracle_auc','fallback_count_fraction'])][['cohort','order','metric','mean','sd','min','max']].to_string(index=False),flush=True)

if __name__=='__main__':main()
