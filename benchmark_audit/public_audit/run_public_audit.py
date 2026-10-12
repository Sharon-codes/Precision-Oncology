#!/usr/bin/env python3
"""Public-release construction audit. No learned binding model is evaluated.

The unmodified TCRCollate class AST is executed with identity-tag embeddings.
Only sequence-selection logic is exercised. Source-row manifests reconstruct
all released derived data locally from the untouched VDJdb release.
"""
import argparse, hashlib, json, platform, sys
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'audit_work'))
from run_source_audit import make_collate, decode_pairs, item, source_metadata

def prepare(source):
    d=pd.read_csv(source,sep='\t',keep_default_na=False)
    d['source_row_0based']=np.arange(len(d)); flow=[{'filter':'raw_release','n':len(d)}]
    tests=[('human',lambda x:x.species.eq('HomoSapiens')),
      ('TRB',lambda x:x.gene.eq('TRB')),('MHCI',lambda x:x['mhc.class'].eq('MHCI')),
      ('exact_two_field_HLA',lambda x:x['mhc.a'].str.fullmatch(r'HLA-[ABC]\*\d{2}:\d{2}')),
      ('canonical_CDR3',lambda x:x.cdr3.str.fullmatch(r'C[ACDEFGHIKLMNPQRSTVWY]*[FW]')),
      ('peptide_8_to_15',lambda x:x['antigen.epitope'].str.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY]{8,15}'))]
    for name,test in tests:
        d=d[test(d)].copy();flow.append({'filter':name,'n':len(d)})
    d=d.rename(columns={'cdr3':'beta','antigen.epitope':'peptide','mhc.a':'hla','vdjdb.score':'score'})
    d=d[['source_row_0based','beta','peptide','hla','score']]
    cohorts={}
    for name,score in [('score_ge_0',0),('score_ge_1',1),('score_eq_3',3)]:
        e=d[d.score.eq(3) if name=='score_eq_3' else d.score.ge(score)].copy();before=len(e)
        e=e.drop_duplicates(['beta','peptide','hla'],keep='first').reset_index(drop=True)
        cohorts[name]=e
        flow.append({'filter':name,'n_before_dedup':before,'n':len(e),'peptides':e.peptide.nunique(),'hlas':e.hla.nunique(),'betas':e.beta.nunique()})
    known=set(d[['beta','peptide','hla']].itertuples(index=False,name=None))
    known_pair=set(d[['beta','peptide']].itertuples(index=False,name=None))
    return cohorts,known,known_pair,flow

def marginal(a,b):
    p=Counter(a);q=Counter(b);keys=set(p)|set(q);n=len(a)
    tv=.5*sum(abs(p[k]-q[k]) for k in keys)/n
    # Category-frequency resubstitution identity oracle, explicitly not a
    # held-out predictor. Weighted calculation correctly handles tied scores.
    groups={}
    from fractions import Fraction
    for k in keys:
        score=Fraction(p[k],p[k]+q[k]); pp,qq=groups.get(score,(0,0));groups[score]=(pp+p[k],qq+q[k])
    cumneg=0;u=0
    for score,(pp,qq) in sorted(groups.items()):u+=pp*(cumneg+.5*qq);cumneg+=qq
    return float(tv),float(u/n/n)

def metrics(d,neg,known,known_pair,fallback=0,hardcoded=0):
    pos=list(d[['beta','peptide','hla']].itertuples(index=False,name=None));n=len(pos)
    result={'n_positive':n,'n_negative':len(neg),'known_positive_collisions':sum(t in known for t in neg),
      'known_beta_peptide_collisions_any_hla':sum(t[:2] in known_pair for t in neg),
      'fallback_count':fallback,'hardcoded_count':hardcoded}
    for k in ['known_positive_collisions','known_beta_peptide_collisions_any_hla','fallback_count','hardcoded_count']:result[k+'_fraction']=result[k]/n
    for name,fn in [('peptide',lambda t:t[1]),('hla',lambda t:t[2]),('peptide_hla',lambda t:t[1:]),('beta',lambda t:t[0])]:
        tv,auc=marginal([fn(t) for t in pos],[fn(t) for t in neg]);result[name+'_tv']=tv;result[name+'_identity_oracle_auc']=auc
    return result

def original(d,repo,order,batch_size,seed,known):
    ids=np.arange(len(d))
    if order=='peptide_sorted':ids=np.array(sorted(ids,key=lambda i:(d.iloc[i].peptide,d.iloc[i].hla,int(i))))
    elif order=='row_shuffled':ids=np.random.default_rng(100000+seed).permutation(ids)
    ordered=d.iloc[ids].reset_index(drop=True)
    records=[item(b,p,h) for b,p,h in ordered[['beta','peptide','hla']].itertuples(index=False,name=None)]
    pool=list(dict.fromkeys(d[['peptide','hla']].itertuples(index=False,name=None)))
    C,tags,rng=make_collate(repo,seed);collate=C(known,pool)
    calls=[0];load=tags.load
    def counted(x):calls[0]+=1;return load(x)
    collate.__call__.__globals__['load_embedding']=counted
    neg=[];fallback=0
    for k in range(0,len(records),batch_size):
        batch=records[k:k+batch_size];before=calls[0];out=collate(batch)
        fallback+=calls[0]-before-3*len(batch);neg.extend(decode_pairs(batch,out,tags))
    return ordered,neg,fallback,sum(p=='AAAAA' for _,p,_ in neg)

def corrected(d,known):
    """Maximum-cardinality closed within-HLA assignment, including abstention.

    Each valid off-diagonal edge costs zero, a diagonal abstention costs one,
    and a forbidden edge costs N+1. An identity permutation always costs N,
    so no forbidden edge can be selected. Minimizing cost maximizes retained
    rows while ensuring retained source and donor sets are identical.
    """
    selected=[];donors=[];stats=[]
    for h,sub in d.groupby('hla',sort=True):
        ix=sub.index.to_numpy();n=len(ix);p=sub.peptide.to_numpy();b=sub.beta.to_numpy()
        allowed=np.array([[(bb,pp,h) not in known for pp in p] for bb in b],dtype=bool)
        np.fill_diagonal(allowed,False)
        maxmatch=maximum_bipartite_matching(csr_matrix(allowed),perm_type='column')
        cost=np.full((n,n),n+1,dtype=np.int32);cost[allowed]=0;np.fill_diagonal(cost,1)
        r,c=linear_sum_assignment(cost);retain=r[c!=r];dst=c[c!=r]
        assert np.all(allowed[retain,dst]);assert set(retain)==set(dst)
        selected.extend(ix[retain]);donors.extend(ix[dst])
        stats.append({'hla':h,'n_positive':n,'n_retained':len(retain),'dominant_peptide_n':int(sub.peptide.value_counts().max()),'maximum_bipartite_matching_size':int((maxmatch>=0).sum()),'n_distinct_peptides':sub.peptide.nunique()})
    perm=np.argsort(selected);selected=np.asarray(selected)[perm];donors=np.asarray(donors)[perm]
    retained=d.loc[selected].reset_index(drop=True)
    neg=[(d.loc[i,'beta'],d.loc[j,'peptide'],d.loc[j,'hla']) for i,j in zip(selected,donors)]
    assert not any(t in known for t in neg)
    for col,idx in [('beta',0),('peptide',1),('hla',2)]:assert Counter(retained[col])==Counter(t[idx] for t in neg)
    assert Counter(retained[['peptide','hla']].itertuples(index=False,name=None))==Counter(t[1:] for t in neg)
    return retained,neg,stats

def save_manifest(out,name,d,neg,all_source):
    lookup={}
    for r,p,h in all_source[['source_row_0based','peptide','hla']].itertuples(index=False,name=None):lookup.setdefault((p,h),int(r))
    manifest=pd.DataFrame({'receiver_source_row_0based':d.source_row_0based.to_numpy(),'negative_antigen_source_row_0based':[lookup.get(t[1:],-1) for t in neg]})
    # -1 means the source collator's hardcoded AAAAA / A*02:01 sentinel.
    manifest.to_csv(out/'manifests'/f'{name}.csv',index=False)
    # Local-only sequence materialization for independent analysis, not included
    # in redistributable source package; users regenerate it from original data.
    pos=d[['source_row_0based','beta','peptide','hla']].copy();pos['label']=1
    nd=pd.DataFrame(neg,columns=['beta','peptide','hla']);nd['source_row_0based']=d.source_row_0based.to_numpy();nd['label']=0
    pd.concat([pos,nd],ignore_index=True).to_csv(out/'local_pairs'/f'{name}.csv',index=False)

def summarize(rows):
    frame=pd.DataFrame(rows);group=['cohort','sampler','order','batch_size'];metrics_cols=[c for c in frame if c not in group+['seed']]
    agg=[]
    for keys,g in frame.groupby(group,dropna=False):
        for metric in metrics_cols:
            a=g[metric].to_numpy(dtype=float)
            agg.append(dict(zip(group,keys),metric=metric,n_seeds=len(a),mean=float(a.mean()),sd=float(a.std(ddof=1)) if len(a)>1 else 0.,min=float(a.min()),max=float(a.max()),q025=float(np.quantile(a,.025)),q975=float(np.quantile(a,.975))))
    return frame,pd.DataFrame(agg)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,default=HERE.parent/'research/data/vdjdb.txt');ap.add_argument('--repo',type=Path,default=HERE.parent/'precision_repo');ap.add_argument('--out',type=Path,default=HERE);ap.add_argument('--seeds',type=int,default=30);args=ap.parse_args();out=args.out
    for folder in ['manifests','local_pairs']:(out/folder).mkdir(parents=True,exist_ok=True)
    cohorts,known,known_pair,flow=prepare(args.source)
    meta={'scope':'Public-release data-construction audit; no PLM inference or original-cohort replication.','source_release':'VDJdb 2023-06-01','source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),'source_rows':92771,'filter_flow':flow,'known_positive_reference':'All eligible human TRB MHCI canonical sequence exact two-field-HLA records, score >=0; exact beta/peptide/HLA triples and beta/peptide pairs (any HLA) separately.','identity_oracle':'Resubstitution diagnostic of category class marginals on the same generated dataset, not predictive generalization.','uncertainty':'Variation across negative-sampler RNG seeds; not biological or population confidence intervals.','order_rng':'numpy.default_rng(100000+seed), separate from source collator Python random.Random(seed); global fallback pool preserves fixed source order regardless of loader order.','versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},'source_code':source_metadata(args.repo)}
    (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n');pd.DataFrame(flow).to_csv(out/'filter_flow.csv',index=False)
    rows=[]
    # Correction and common-cohort comparisons first: independent reviewer can
    # immediately begin cross-fitting on these fixed generated pair datasets.
    d=cohorts['score_ge_1'];retained,neg,stats=corrected(d,known)
    pd.DataFrame(stats).to_csv(out/'matching_retention_by_hla.csv',index=False)
    save_manifest(out,'common_corrected_seed0',retained,neg,cohorts['score_ge_0'])
    rows.append(dict(cohort='score_ge_1_retained',sampler='exact_closed_matching',order='fixed_manifest',batch_size=0,seed=0,**metrics(retained,neg,known,known_pair)))
    print('CORRECTION',len(d),len(retained),flush=True)
    for order in ['release_order','peptide_sorted','row_shuffled']:
        for seed in range(args.seeds):
            ordered,neg,fallback,hardcoded=original(retained,args.repo,order,64,seed,known)
            rows.append(dict(cohort='score_ge_1_retained',sampler='source_collator',order=order,batch_size=64,seed=seed,**metrics(ordered,neg,known,known_pair,fallback,hardcoded)))
            if seed<10:save_manifest(out,f'common_original_{order}_seed{seed}',ordered,neg,cohorts['score_ge_0'])
        print('DONE common',order,flush=True)
    frame,agg=summarize(rows);frame.to_csv(out/'seed_results.csv',index=False);agg.to_csv(out/'aggregate_results.csv',index=False)
    for cohort,d in cohorts.items():
        batches=[16,64,256] if cohort=='score_ge_1' else [64]
        for bs in batches:
            for order in ['release_order','peptide_sorted','row_shuffled']:
                for seed in range(args.seeds):
                    ordered,neg,fallback,hardcoded=original(d,args.repo,order,bs,seed,known)
                    rows.append(dict(cohort=cohort,sampler='source_collator',order=order,batch_size=bs,seed=seed,**metrics(ordered,neg,known,known_pair,fallback,hardcoded)))
                    if seed==0:save_manifest(out,f'{cohort}_{order}_B{bs}_seed0',ordered,neg,cohorts['score_ge_0'])
                frame,agg=summarize(rows);frame.to_csv(out/'seed_results.csv',index=False);agg.to_csv(out/'aggregate_results.csv',index=False)
                print('DONE',cohort,bs,order,flush=True)
    meta['retention']={'original_n':len(cohorts['score_ge_1']),'retained_n':len(retained),'retention_fraction':len(retained)/len(cohorts['score_ge_1']),'algorithm':'Minimum-cost full assignment allowing diagonal abstention, independently within each HLA; maximum cardinality closed subset of rows admitting a known-positive-excluding permutation.','known_positive_collisions_corrected':0,'exact_equal_marginals':['beta','peptide','hla','peptide_hla']}
    (out/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(agg[agg.metric.isin(['peptide_tv','peptide_identity_oracle_auc','known_positive_collisions_fraction','fallback_count_fraction'])].to_string(index=False),flush=True)

if __name__=='__main__':main()
