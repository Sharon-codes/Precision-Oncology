#!/usr/bin/env python3
"""Meaningful reproducibility checks for source-index manifests and assignment."""
import itertools,json
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from run_public_audit import HERE,prepare,corrected

def main():
    source=pd.read_csv(HERE.parent/'research/data/vdjdb.txt',sep='\t',keep_default_na=False)
    checks=[]
    for path in sorted(HERE.rglob('manifests/*.csv')):
        if path.name.startswith('partition_'):continue
        material=path.parent.parent/'local_pairs'/path.name
        d=pd.read_csv(path);p=pd.read_csv(material,keep_default_na=False);neg=p[p.label.eq(0)].reset_index(drop=True)
        receiver=d.receiver_source_row_0based.to_numpy();donor=d.negative_antigen_source_row_0based.to_numpy()
        assert len(d)==len(neg)
        assert np.all(neg.beta.to_numpy()==source.iloc[receiver].cdr3.to_numpy())
        regular=donor>=0
        assert np.all(neg.peptide.to_numpy()[regular]==source.iloc[donor[regular]]['antigen.epitope'].to_numpy())
        assert np.all(neg.hla.to_numpy()[regular]==source.iloc[donor[regular]]['mhc.a'].to_numpy())
        assert np.all(neg.peptide.to_numpy()[~regular]=='AAAAA')
        assert np.all(neg.hla.to_numpy()[~regular]=='A*02:01')
        checks.append({'manifest':str(path.relative_to(HERE)),'n_negative':len(neg),'reconstruction_exact':True})
    rng=np.random.default_rng(420);count=0
    for n in range(2,7):
        for k in range(20):
            allowed=rng.random((n,n))<.4;np.fill_diagonal(allowed,False)
            cost=np.full((n,n),n+1,dtype=np.int32);cost[allowed]=0;np.fill_diagonal(cost,1)
            r,c=linear_sum_assignment(cost);retained=r[r!=c];donors=c[r!=c]
            assert set(retained)==set(donors)
            assert np.all(allowed[retained,donors])
            best=0
            for perm in itertools.permutations(range(n)):
                if all(i==j or allowed[i,j] for i,j in enumerate(perm)):
                    best=max(best,sum(i!=j for i,j in enumerate(perm)))
            assert len(retained)==best;count+=1
    cohorts,known,_,_=prepare(HERE.parent/'research/data/vdjdb.txt')
    a,b,_=corrected(cohorts['score_ge_1'],known);saved=pd.read_csv(HERE/'local_pairs/common_corrected_seed0.csv')
    assert b==list(saved[saved.label.eq(0)][['beta','peptide','hla']].itertuples(index=False,name=None))
    result={'all_passed':True,'source_row_indexing':'0-based data rows excluding TSV header','checks':checks,'exhaustive_small_graph_maximality_checks':count,'deterministic_main_assignment_replay':True}
    (HERE/'manifest_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'manifest_files':len(checks),'reconstructed_rows':sum(c['n_negative'] for c in checks),'exhaustive_graph_checks':count,'deterministic_replay':True}))
if __name__=='__main__':main()
