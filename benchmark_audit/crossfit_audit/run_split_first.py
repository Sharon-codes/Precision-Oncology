"""Split receptors BEFORE negative construction to avoid complement imbalance.

Each fold gets a within-HLA cycle-cover-matched retained cohort, independently.
All compared negative samplers use that same per-fold positive population.
The original global-construction-then-crossfit audit is retained unchanged as
an explicit protocol-placement sensitivity, not discarded for its low AUCs.
"""
import argparse
from collections import Counter
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'public_audit'))
from run_public_audit import prepare,corrected,original,metrics
from run_crossfit import fold_map,run_dataset,validate_dataset,positive_signature


def pair_frame(d,neg,fold):
    pos=d[['source_row_0based','beta','peptide','hla']].copy().rename(columns={'beta':'tcr_beta'})
    pos['label']=1
    nd=pd.DataFrame(neg,columns=['tcr_beta','peptide','hla'])
    nd['source_row_0based']=d.source_row_0based.to_numpy()
    nd['label']=0
    z=pd.concat([pos,nd],ignore_index=True)
    z['preassigned_fold']=fold
    z['example_id']=z.source_row_0based.astype(str)+'_y'+z.label.astype(str)
    return z


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',type=Path,default=ROOT/'research/data/vdjdb.txt')
    ap.add_argument('--repo',type=Path,default=ROOT/'precision_repo')
    ap.add_argument('--out',type=Path,default=HERE/'split_first_results')
    ap.add_argument('--sampling-seeds',type=int,default=10)
    ap.add_argument('--fold-seeds',default='17,41,73,107,211')
    args=ap.parse_args();args.out.mkdir(exist_ok=True,parents=True)
    (args.out/'manifests').mkdir(exist_ok=True)
    cohorts,known,known_pair,flow=prepare(args.source)
    d=cohorts['score_ge_1']
    # Identity lookup chooses one source row with the needed peptide-HLA value.
    donor_lookup={}
    for row in cohorts['score_ge_0'].itertuples():
        donor_lookup.setdefault((row.peptide,row.hla),int(row.source_row_0based))
    all_results=[];all_folds=[];all_scores=[];cohort_rows=[];construction=[]
    for fs in map(int,args.fold_seeds.split(',')):
        fm=fold_map(d.beta,fs);fold_ids=d.beta.map(fm)
        fixed_parts=[];parts={};fixed_retained={}
        for fold in range(5):
            raw=d.loc[fold_ids.eq(fold)].reset_index(drop=True)
            ret,neg,_=corrected(raw,known)
            fixed_retained[fold]=ret
            fixed_parts.append(pair_frame(ret,neg,fold))
            cohort_rows.append({'fold_seed':fs,'fold':fold,'n_original_positive':len(raw),
                                'n_retained_positive':len(ret),'retention_fraction':len(ret)/len(raw),
                                'n_receptor_groups_retained':ret.beta.nunique()})
            construction.append({'fold_seed':fs,'fold':fold,'sampler':'exact_closed_matching',
                                 'sampling_seed':0,**metrics(ret,neg,known,known_pair)})
        parts[('exact_closed_matching',0)]=fixed_parts
        for order in ['release_order','peptide_sorted','row_shuffled']:
            for seed in range(args.sampling_seeds):
                blocks=[]
                for fold in range(5):
                    # Independent draws per fold, recorded deterministic seeds.
                    draw_seed=seed+10000*fold
                    ordered,neg,fallback,hardcoded=original(fixed_retained[fold],args.repo,order,64,draw_seed,known)
                    blocks.append(pair_frame(ordered,neg,fold))
                    construction.append({'fold_seed':fs,'fold':fold,'sampler':'source_'+order,
                                         'sampling_seed':seed,**metrics(ordered,neg,known,known_pair,fallback,hardcoded)})
                parts[('source_'+order,seed)]=blocks
        positive_ref=None
        for (sampler,seed),blocks in parts.items():
            df=pd.concat(blocks,ignore_index=True)
            validate_dataset(df)
            signature=positive_signature(df)
            if positive_ref is not None and signature!=positive_ref:raise AssertionError('Common cohort mismatch')
            positive_ref=signature
            desc={'sampler':sampler,'sampling_seed':seed,'cohort':f'split_first_fold_seed_{fs}'}
            met,folds,scores=run_dataset(df,desc,args.out,[fs])
            all_results.extend(met);all_folds.extend(folds);all_scores.extend(scores)
            # Package only source indices and score outputs, not modified sequence text.
            neg_df=df.loc[df.label.eq(0)]
            manifest=pd.DataFrame({'receiver_source_row_0based':neg_df.source_row_0based,
                'negative_antigen_source_row_0based':[donor_lookup.get((r.peptide,r.hla),-1) for r in neg_df.itertuples()],
                'fold':neg_df.preassigned_fold})
            manifest.to_csv(args.out/'manifests'/f'f{fs}_{sampler}_s{seed}.csv',index=False)
        print('Completed fold seed',fs,'retained',sum(len(x) for x in fixed_retained.values()),flush=True)
    frame=pd.DataFrame(all_results)
    frame.to_csv(args.out/'oof_metrics.csv',index=False)
    pd.DataFrame(all_folds).to_csv(args.out/'fold_diagnostics.csv',index=False)
    pd.concat(all_scores,ignore_index=True).to_csv(args.out/'oof_scores.csv.gz',index=False,compression='gzip')
    pd.DataFrame(cohort_rows).to_csv(args.out/'retention.csv',index=False)
    pd.DataFrame(construction).to_csv(args.out/'construction_metrics.csv',index=False)
    # Collapse fold-map seeds before summarizing decoy seed variability.
    seed_means=frame.groupby(['sampler','feature','sampling_seed'],as_index=False).agg(
        auroc=('auroc','mean'),average_precision=('average_precision','mean'),brier=('brier','mean'))
    seed_means.to_csv(args.out/'sampling_seed_means.csv',index=False)
    summary=seed_means.groupby(['sampler','feature'],as_index=False).agg(
        sampling_seeds=('sampling_seed','nunique'),auroc_mean=('auroc','mean'),
        auroc_sd_between_sampler_seeds=('auroc','std'),auroc_min=('auroc','min'),auroc_max=('auroc','max'),
        ap_mean=('average_precision','mean'),brier_mean=('brier','mean'))
    summary.to_csv(args.out/'summary.csv',index=False)
    metadata={'purpose':'Split-before-construction metadata control within an antigen collection; not unseen-peptide binding prediction',
      'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),'cohort':'VDJdb2023-06-01 strict human TRB MHCI score>=1',
      'fold_seeds':list(map(int,args.fold_seeds.split(','))),'n_folds':5,'sampling_seeds':args.sampling_seeds,
      'draw_seed_rule':'sampling_seed + 10000 * fold','retention':'maximum closed within-HLA known-positive-excluding assignment separately in each receptor fold',
      'common_population':'all samplers on same retained positives within each fold seed; cohort differs across fold seeds',
      'training':'onehot fit only training folds; C1 liblinear tol1e-8; no receptor sequence for antigen models',
      'uncertainty':'mean of fold-seed metrics then SD/range across negative-sampling seeds; not population CI',
      'protocol_note':'Original global-construction-then-crossfit retained as separate sensitivity because exact global marginals couple complementary fold label imbalances'}
    (args.out/'run_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(summary.to_string(index=False),flush=True)


if __name__=='__main__':main()
