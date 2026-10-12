#!/usr/bin/env python3
"""Verify archived aggregates and every stored out-of-fold prediction.

No external data, networking, model retraining, or repository imports required.
Run: python audit_work/verify_archived_results.py --root <archive-root>
This verifies numeric internal consistency and stored fold/group certificates;
it cannot establish assay labels or reproduce training from scores alone.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss

CHECKS=0

def eq(actual,expected,context,atol=1e-10):
    global CHECKS
    aa=np.asarray(actual,dtype=float);bb=np.asarray(expected,dtype=float)
    if not np.allclose(aa,bb,atol=atol,rtol=0,equal_nan=True):
        raise AssertionError(f'{context}: {actual} != {expected}')
    CHECKS+=1

def truth(value,context):
    global CHECKS
    if not value:raise AssertionError(context)
    CHECKS+=1

def csv(path,**kw):return pd.read_csv(path,float_precision='round_trip',**kw)
def obj(path):return json.loads(path.read_text())

def verify_public(root):
    d=csv(root/'public_audit/seed_results.csv');agg=csv(root/'public_audit/aggregate_results.csv')
    keys=['cohort','sampler','order','batch_size'];measures=[c for c in d if c not in keys+['seed']]
    count=0
    for vals,g in d.groupby(keys,dropna=False):
        select=agg
        for k,v in zip(keys,vals):select=select[select[k].eq(v)]
        for metric in measures:
            row=select[select.metric.eq(metric)];truth(len(row)==1,'public aggregation unique row');r=row.iloc[0];x=g[metric].to_numpy(float)
            actual=[len(x),x.mean(),x.std(ddof=1) if len(x)>1 else 0.,x.min(),x.max(),np.quantile(x,.025),np.quantile(x,.975)]
            eq(actual,[r.n_seeds,r['mean'],r.sd,r['min'],r['max'],r.q025,r.q975],f'public {vals} {metric}');count+=1
    truth(count==len(agg),'public aggregation complete')
    primary=d[(d.cohort=='score_ge_1')&(d.batch_size==64)]
    claim={}
    for order,expected in [('release_order',2.9887),('peptide_sorted',.8820),('row_shuffled',.2284)]:
        x=primary[primary.order==order].known_positive_collisions_fraction.to_numpy()*100
        eq(x.mean(),expected,f'manuscript contradiction {order}',atol=.00005);claim[order]={'mean_percentage':float(x.mean()),'sd_percentage_points':float(x.std(ddof=1))}
    shuffled=d[(d.cohort=='score_ge_1')&(d.order=='row_shuffled')]
    truth(shuffled.peptide_tv.eq(0).all(),'full primary shuffled peptide TV exactly zero')
    ret=csv(root/'public_audit/matching_retention_by_hla.csv')
    eq([ret.n_positive.sum(),ret.n_retained.sum(),len(ret),(ret.n_retained>0).sum()],[3182,2691,41,22],'retention numeric claims')
    tv=.5*sum(abs(ret.n_positive/ret.n_positive.sum()-ret.n_retained/ret.n_retained.sum()))
    eq(tv,.1274,'HLA selection shift rounding',atol=.00005)
    sens=obj(root/'public_audit/additional_sensitivity.json');eq(sens['score3_corrected']['retained_n'],139,'score3 retention');eq(sens['main_confidence_reference_only']['retained_n'],2696,'narrow reference retention');eq(sens['main_confidence_reference_only']['known_collisions_against_all_scores'],8,'narrow reference contradictions')
    c=csv(root/'public_audit/confidence_collision_sensitivity.csv')
    for order,exp in [('release_order',2.6639),('peptide_sorted',.8003),('row_shuffled',.1446)]:
        z=c[c.order==order];eq((100*z.score_ge_1_known_collisions/z.n_negative).mean(),exp,f'narrow reference {order}',atol=.00005)
    splits=obj(root/'public_audit/split_sensitivity/split_provenance.json')['splits']
    truth(len(splits)==5 and all(s['minimum_cross_partition_levenshtein']==3 for s in splits),'split distance certificates recorded')
    return dict(seed_rows=len(d),aggregate_rows=count,primary_contradictions=claim,retained_hla_strata=int((ret.n_retained>0).sum()),hla_selection_tv=float(tv))

def verify_source_and_reference(root):
    ref=obj(root/'audit_work/reference_rejection_results.json');rs=csv(root/'audit_work/reference_rejection_seed_results.csv')
    for col,s in ref['summary'].items():
        x=rs[col].to_numpy(float);eq([x.mean(),x.std(ddof=1),x.min(),x.max()],[s['mean'],s['sd'],s['min'],s['max']],f'rejection {col}')
    truth(rs.known_positive_collisions.eq(0).all() and rs.hla_tv.eq(0).all() and rs.beta_tv.eq(0).all(),'rejection recorded invariant metrics')
    eq(ref['summary']['peptide_tv']['mean'],.1349,'manuscript rejection TV',atol=.00005)
    eq(ref['summary']['peptide_identity_oracle_auc']['mean'],.5990,'manuscript rejection AUROC',atol=.00005)
    eq(ref['closed_matching_same_cohort']['unique_negative_triplets'],2624,'matching unique triplets')
    src=obj(root/'audit_work/source_audit_results.json');ss=csv(root/'audit_work/source_audit_seed_results.csv')
    for order in ['peptide_sorted','row_shuffled']:
        for metric,r in src['synthetic_batch_order_simulation'][order].items():
            if not isinstance(r,dict):continue
            x=ss[ss.order==order][metric].to_numpy(float)
            eq([x.mean(),x.std(ddof=1),x.min(),x.max(),np.quantile(x,.025),np.quantile(x,.975)],[r['mean'],r['sd'],r['min'],r['max'],r['q025'],r['q975']],f'synthetic {order}/{metric}')
    m=obj(root/'audit_work/matching_validation_results.json')
    eq([len(m['fixtures']),sum(f['permutations_examined'] for f in m['fixtures'])],[166,924960],'exhaustive validation numeric claims')
    truth(all(f['retained_actual']==f['retained_exhaustive_optimum'] and f['certificates_passed'] for f in m['fixtures']),'exhaustive optimum certificates')
    return dict(rejection_seed_rows=len(rs),synthetic_seed_rows=len(ss),exhaustive_fixtures=len(m['fixtures']))

def verify_oof(path,split_first,annotation_filtered=False):
    scores=csv(path/'oof_scores.csv.gz');met=csv(path/'oof_metrics.csv');folds=csv(path/'fold_diagnostics.csv')
    seeds=csv(path/'sampling_seed_means.csv');summary=csv(path/'summary.csv')
    keys=['cohort','sampler','sampling_seed','fold_seed'];features=['peptide_hla','peptide','hla','receptor']
    lookup=met.set_index(keys+['feature']);fl=folds.set_index(keys+['feature','fold'])
    truth(not lookup.index.has_duplicates and not fl.index.has_duplicates,'OOF metric keys unique')
    baseline={};n_runs=0;fold_metric_sets=0
    for vals,g in scores.groupby(keys,sort=False):
        truth(g.example_id.is_unique,'OOF example IDs unique per run')
        truth(set(g.label)=={0,1} and set(g.fold)==set(range(5)),'OOF labels/folds complete')
        truth(g.groupby('receptor_sha256').fold.nunique().eq(1).all(),'OOF receptor groups disjoint across folds')
        receptor_counts=g.groupby(['receptor_sha256','label']).size().unstack(fill_value=0)
        truth(receptor_counts[0].equals(receptor_counts[1]),'OOF receptor class multiplicities equal')
        truth(g.example_id.str.endswith('_y1').equals(g.label.eq(1)),'OOF example IDs match labels')
        positive=frozenset(zip(g.loc[g.label.eq(1),'example_id'],g.loc[g.label.eq(1),'receptor_sha256'],g.loc[g.label.eq(1),'fold']))
        cohort_key=(vals[0],vals[3])
        if cohort_key not in baseline:baseline[cohort_key]=positive
        truth(positive==baseline[cohort_key],'same positive cohort and folds across samplers/seeds')
        y=g.label.to_numpy();pcount=int(y.sum());ngroups=g.receptor_sha256.nunique()
        for feature in features:
            pred=g['score_'+feature].to_numpy(float)
            truth(np.isfinite(pred).all() and ((pred>=0)&(pred<=1)).all(),'valid probabilities')
            r=lookup.loc[vals+(feature,)]
            eq([len(g),pcount,ngroups,roc_auc_score(y,pred),average_precision_score(y,pred),brier_score_loss(y,pred),pred.min(),pred.max()],
               [r.n_pairs,r.n_positive,r.n_receptor_groups,r.auroc,r.average_precision,r.brier,r.score_min,r.score_max],f'OOF pooled {vals}/{feature}')
            if feature=='receptor':eq(pred,.5,'receptor-only scores all chance',atol=1e-8)
            for fold in range(5):
                test=g.fold.to_numpy()==fold;z=fl.loc[vals+(feature,fold)]
                eq([(~test).sum(),test.sum(),g.loc[~test,'receptor_sha256'].nunique(),g.loc[test,'receptor_sha256'].nunique(),roc_auc_score(y[test],pred[test]),average_precision_score(y[test],pred[test])],
                   [z.train_pairs,z.test_pairs,z.train_groups,z.test_groups,z.auroc,z.ap],f'OOF fold {vals}/{feature}/{fold}')
                if feature=='receptor':eq(z.unseen_tcr_beta_fraction,1,'all held-out receptors unknown to fold encoder')
                fold_metric_sets+=1
        n_runs+=1
    eq(len(met),n_runs*len(features),'OOF stored metrics fully covered')
    skeys=['sampler','feature','sampling_seed'] if split_first else ['cohort','sampler','feature','sampling_seed']
    recomputed=met.groupby(skeys)[['auroc','average_precision','brier']].mean().sort_index()
    saved=seeds.set_index(skeys).sort_index();truth(recomputed.index.equals(saved.index),'seed mean indices')
    eq(recomputed.to_numpy(),saved[['auroc','average_precision','brier']].to_numpy(),'fold-averaged seed metrics')
    gkeys=['sampler','feature'] if split_first else ['cohort','sampler','feature']
    for vals,g in seeds.groupby(gkeys):
        if not isinstance(vals,tuple):vals=(vals,)
        selected=summary
        for key,value in zip(gkeys,vals):selected=selected[selected[key]==value]
        truth(len(selected)==1,'summary unique');r=selected.iloc[0];x=g.auroc.to_numpy()
        eq([len(g),x.mean(),x.std(ddof=1) if len(x)>1 else np.nan,x.min(),x.max(),g.average_precision.mean(),g.brier.mean()],
           [r.sampling_seeds,r.auroc_mean,r.auroc_sd_between_sampler_seeds,r.auroc_min,r.auroc_max,r.ap_mean,r.brier_mean],f'OOF summary {vals}')
    ph=summary[summary.feature.eq('peptide_hla')].set_index('sampler')
    expected={'source_peptide_sorted':.5619,'source_release_order':.5188,'source_row_shuffled':.5,'exact_closed_matching':.5} if annotation_filtered else {'source_peptide_sorted':.5558,'source_release_order':.5168,'source_row_shuffled':.5,'exact_closed_matching':.5} if split_first else {'source_peptide_sorted':.7146,'source_row_shuffled':.3298,'exact_closed_matching':.3262}
    for sampler,value in expected.items():eq(ph.loc[sampler,'auroc_mean'],value,'main manuscript OOF rounded mean',atol=.00005)
    if split_first:
        totals=csv(path/'retention.csv').groupby('fold_seed').n_retained_positive.sum().tolist()
        eq(totals,[2618,2625,2636,2642,2621] if annotation_filtered else [2667,2659,2652,2670,2658],'split-first retained totals')
        all_sorted=met[(met.feature=='peptide_hla')&(met.sampler=='source_peptide_sorted')].auroc
        eq([all_sorted.min(),all_sorted.max()],[.5532,.5738] if annotation_filtered else [.5435,.5633],'manuscript across-50-run range',atol=.00005)
    return dict(score_rows=len(scores),complete_fivefold_runs=n_runs,pooled_metric_sets=len(met),fold_metric_sets=fold_metric_sets,summary_auroc={str(k):float(v) for k,v in ph.auroc_mean.items()})

def verify_annotation(root):
    folder=root/'public_audit/annotation_sensitivity'
    d=csv(folder/'seed_results.csv');agg=csv(folder/'aggregate_results.csv');meta=obj(folder/'provenance.json')
    keys=['cohort','sampler','order','batch_size'];measures=[c for c in d if c not in keys+['seed']];count=0
    for vals,g in d.groupby(keys,dropna=False):
        select=agg
        for k,v in zip(keys,vals):select=select[select[k].eq(v)]
        for metric in measures:
            row=select[select.metric.eq(metric)];truth(len(row)==1,'annotation unique aggregate');r=row.iloc[0];x=g[metric].to_numpy(float)
            eq([len(x),x.mean(),x.std(ddof=1) if len(x)>1 else 0.,x.min(),x.max(),np.quantile(x,.025),np.quantile(x,.975)],
               [r.n_seeds,r['mean'],r.sd,r['min'],r['max'],r.q025,r.q975],f'annotation aggregate {vals}/{metric}');count+=1
    truth(count==len(agg),'annotation aggregation complete')
    ret=csv(folder/'matching_retention_by_hla.csv')
    eq([ret.n_positive.sum(),ret.n_retained.sum(),len(ret),(ret.n_retained>0).sum()],[3108,2665,39,21],'annotation retained coverage')
    eq([meta['raw_main_n'],meta['clean_main_n'],meta['clean_known_reference_n'],meta['corrected_retained_n']],[3182,3108,29146,2665],'annotation provenance counts')
    primary=d[(d.cohort=='clean_score_ge_1')&(d.batch_size==64)]
    outcomes={}
    for order,contradiction,oracle in [('release_order',2.888,.6423),('peptide_sorted',.933,.7936),('row_shuffled',.245,.5)]:
        z=primary[primary.order==order]
        eq(100*z.known_positive_collisions_fraction.mean(),contradiction,'annotation manuscript contradiction rounding',atol=.0005)
        eq(z.peptide_identity_oracle_auc.mean(),oracle,'annotation manuscript oracle rounding',atol=.00005)
        outcomes[order]={'mean_contradiction_percentage':float(100*z.known_positive_collisions_fraction.mean()),'mean_peptide_oracle_auroc':float(z.peptide_identity_oracle_auc.mean())}
    review=obj(root/'reviews/annotation_sensitivity.json')
    flags=review['flagged_label_counts'];eq([flags['HLA-A*08:01']['by_confidence_cohort']['score_ge_1'],flags['HLA-A*24:01']['by_confidence_cohort']['score_ge_1']],[2,72],'flagged-label primary counts')
    scenario=review['scenarios']['exclude_A08_and_legacy_A24']
    eq([scenario['primary_rows'],scenario['registry_triples'],scenario['retained_rows'],scenario['primary_hla_strata'],scenario['retained_hla_strata']],[3108,29146,2665,39,21],'independent annotation construction review counts')
    truth(review['raw_annotations_not_reassigned'],'annotation review does not recode alleles')
    excluded=set(flags['HLA-A*08:01']['primary_source_rows_0based'])|set(flags['HLA-A*24:01']['primary_source_rows_0based'])
    baseline=csv(root/'public_audit/manifests/common_corrected_seed0.csv')
    clean=csv(folder/'manifests/clean_corrected_seed0.csv')
    expected=baseline[~baseline.receiver_source_row_0based.isin(excluded)].reset_index(drop=True)
    truth(expected.equals(clean),'annotation exact manifest equals baseline with flagged strata removed')
    for metric in ['known_positive_collisions','peptide_tv','hla_tv','beta_tv','peptide_hla_tv']:
        eq(scenario['construction_metrics'][metric],0,'independent annotation certificate '+metric)
    return dict(seed_rows=len(d),aggregate_rows=count,raw_rows=3182,filtered_rows=3108,retained_rows=2665,retained_strata=21,primary_results=outcomes)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent);ap.add_argument('--output',type=Path);a=ap.parse_args();root=a.root.resolve()
    report={'scope':'Archived-result and stored-OOF-score consistency verification; no model retraining or biological label verification.'}
    report['public']=verify_public(root);report['source_and_reference']=verify_source_and_reference(root)
    report['global_first_oof']=verify_oof(root/'crossfit_audit/results',False)
    report['split_first_oof']=verify_oof(root/'crossfit_audit/split_first_results',True)
    report['annotation_public']=verify_annotation(root)
    report['annotation_filtered_oof']=verify_oof(root/'crossfit_audit/annotation_filtered_results',True,annotation_filtered=True)
    report['numeric_and_certificate_assertions']=CHECKS;report['status']='PASS'
    inputs=[root/'crossfit_audit/annotation_filtered_results/oof_scores.csv.gz',root/'public_audit/annotation_sensitivity/seed_results.csv',root/'reviews/annotation_sensitivity.json',root/'paper/manuscript.tex',root/'public_audit/seed_results.csv',root/'crossfit_audit/results/oof_scores.csv.gz',root/'crossfit_audit/split_first_results/oof_scores.csv.gz']
    report['input_sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs if p.exists()}
    out=a.output or root/'audit_work/archive_verification_results.json';out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
