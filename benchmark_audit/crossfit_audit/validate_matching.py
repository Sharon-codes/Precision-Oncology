"""Independent exhaustive validation of the public audit's cycle-cover result."""
import itertools,json,sys
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'public_audit'))
from run_public_audit import corrected,prepare

rng=np.random.default_rng(771)
checks=0
for n in range(2,8):
    for rep in range(20):
        d=pd.DataFrame({'beta':[f'R{i}' for i in range(n)],'peptide':[f'P{i}' for i in range(n)],
                        'hla':['H']*n,'source_row_0based':range(n)})
        forbidden=rng.random((n,n))<.5
        np.fill_diagonal(forbidden,True)
        known={(f'R{i}',f'P{j}','H') for i in range(n) for j in range(n) if forbidden[i,j]}
        kept,neg,_=corrected(d,known)
        best=0
        for perm in itertools.permutations(range(n)):
            if all(j==i or not forbidden[i,j] for i,j in enumerate(perm)):
                best=max(best,sum(i!=j for i,j in enumerate(perm)))
        assert len(kept)==best,(n,rep,len(kept),best)
        checks+=1
cohorts,known,known_pair,flow=prepare(ROOT/'research/data/vdjdb.txt')
pairs=pd.read_csv(ROOT/'public_audit/local_pairs/common_corrected_seed0.csv')
pos=pairs[pairs.label.eq(1)];neg=pairs[pairs.label.eq(0)]
for cols in [['beta'],['peptide'],['hla'],['peptide','hla']]:
    assert Counter(map(tuple,pos[cols].values))==Counter(map(tuple,neg[cols].values))
collisions=sum(tuple(x) in known for x in neg[['beta','peptide','hla']].values)
assert collisions==0
result={'independent_bruteforce_graph_checks':checks,'graph_size_range':[2,7],
    'public_retained_positive_rows':len(pos),'public_original_positive_rows':len(cohorts['score_ge_1']),
    'exact_marginals_independently_confirmed':['beta','peptide','hla','peptide_hla'],
    'known_positive_triple_collisions':collisions,
    'verdict':'cycle-cover correction passed independent finite exhaustive tests and public CSV invariants'}
(ROOT/'crossfit_audit/independent_matching_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
