#!/usr/bin/env python3
"""Verify decoy diagnostics after strict peptide-component partitioning.

Each distance<=2 connected component belongs to exactly one partition.
This is a construction audit, not a trained unseen-peptide prediction test.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from run_public_audit import prepare,original,corrected,metrics,save_manifest,summarize,HERE
from run_source_audit import lev

def main():
    cohorts,known,kp,flow=prepare(HERE.parent/'research/data/vdjdb.txt');d=cohorts['score_ge_1']
    peps=sorted(d.peptide.unique());n=len(peps);dist=np.zeros((n,n),dtype=np.int16)
    for i in range(n):
        for j in range(i):dist[i,j]=dist[j,i]=lev(peps[i],peps[j])
    ncomp,labels=connected_components(csr_matrix(dist<=2),directed=False)
    pcomp=dict(zip(peps,labels));rows=[];splits=[];out=HERE/'split_sensitivity'
    for f in ['manifests','local_pairs']:(out/f).mkdir(parents=True,exist_ok=True)
    for splitseed in range(5):
        order=np.random.default_rng(200000+splitseed).permutation(ncomp)
        # Fixed fractions of components, with true row fractions reported.
        groups={'train':order[:int(.6*ncomp)],'validation':order[int(.6*ncomp):int(.8*ncomp)],'test':order[int(.8*ncomp):]}
        partmap={int(c):name for name,cs in groups.items() for c in cs}
        partition=d.peptide.map(lambda p:partmap[int(pcomp[p])]);test=d[partition.eq('test')].reset_index(drop=True)
        ppart=np.array([partmap[int(c)] for c in labels]);cross=ppart[:,None]!=ppart[None,:]
        mincross=int(dist[cross].min());assert mincross>2
        splitrows=d[['source_row_0based']].copy();splitrows['component_id']=d.peptide.map(pcomp);splitrows['partition']=partition;splitrows.to_csv(out/'manifests'/f'partition_seed{splitseed}.csv',index=False)
        retained,neg,stats=corrected(test,known)
        pd.DataFrame(stats).to_csv(out/f'retention_split{splitseed}.csv',index=False)
        save_manifest(out,f'test_split{splitseed}_corrected',retained,neg,cohorts['score_ge_0'])
        splits.append({'split_seed':splitseed,'n_components':ncomp,'partition_component_counts':{k:len(v) for k,v in groups.items()},'partition_row_counts':partition.value_counts().to_dict(),'minimum_cross_partition_levenshtein':mincross,'test_retained_n':len(retained),'test_retention_fraction':len(retained)/len(test),'test_unique_peptides':test.peptide.nunique(),'retained_unique_peptides':retained.peptide.nunique()})
        rows.append(dict(cohort=f'test_split{splitseed}_retained',sampler='exact_closed_matching',order='fixed_manifest',batch_size=0,seed=0,**metrics(retained,neg,known,kp)))
        for subset_name,sub in [('full',test),('retained',retained)]:
            for loader_order in ['release_order','peptide_sorted','row_shuffled']:
                for seed in range(10):
                    ordered,neg,fall,hard=original(sub,HERE.parent/'precision_repo',loader_order,64,seed,known)
                    rows.append(dict(cohort=f'test_split{splitseed}_{subset_name}',sampler='source_collator',order=loader_order,batch_size=64,seed=seed,**metrics(ordered,neg,known,kp,fall,hard)))
                    if seed==0:save_manifest(out,f'test_split{splitseed}_{subset_name}_{loader_order}_seed0',ordered,neg,cohorts['score_ge_0'])
        print('DONE split',splitseed,'test',len(test),'retained',len(retained),'mindist',mincross,flush=True)
    frame,agg=summarize(rows);frame.to_csv(out/'seed_results.csv',index=False);agg.to_csv(out/'aggregate_results.csv',index=False)
    (out/'split_provenance.json').write_text(json.dumps({'design':'Five deterministic peptide-component partitions; ten collator RNG seeds per test partition/order; partition weights60/20/20 by number of connected components, actual row counts reported. No model fitting or OOD performance claims.','split_rng':'numpy.default_rng(200000+split_seed)','distance':'unit-cost Levenshtein edit distance; connected components at distance<=2','splits':splits},indent=2)+'\n')
    print(json.dumps(splits,indent=2))
if __name__=='__main__':main()
