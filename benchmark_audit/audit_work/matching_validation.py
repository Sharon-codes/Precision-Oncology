#!/usr/bin/env python3
"""Independent exhaustive validation of corrected closed within-HLA assignment.

Enumerates every permutation for tiny synthetic fixtures with repeated receptor
and peptide identities and cross-reactivity. Executes the imported original
corrected() function and passively records its return-frame variables; does
not alter public_audit source. This tests a combinatorial contract, not binding.
"""
import argparse, hashlib, importlib.util, itertools, json, math, sys
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

ROOT=Path(__file__).resolve().parent

def load_public(path):
    spec=importlib.util.spec_from_file_location('validated_public_audit',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def allowed_matrix(d,known):
    rows=list(d[['beta','peptide','hla']].itertuples(index=False,name=None))
    return np.array([[i!=j and h==hh and (b,pp,h) not in known
                      for j,(_,pp,hh) in enumerate(rows)]
                     for i,(b,_,h) in enumerate(rows)],dtype=bool)

def brute_force(allowed):
    n=len(allowed);best=-1;witness=None;feasible=0
    for perm in itertools.permutations(range(n)):
        if all(i==j or allowed[i,j] for i,j in enumerate(perm)):
            feasible+=1
            retained=sum(i!=j for i,j in enumerate(perm))
            if retained>best:best=retained;witness=list(perm)
    return best,witness,feasible,math.factorial(n)

def recorded_call(function,*args):
    captured={};old=sys.gettrace()
    def tracer(frame,event,arg):
        if frame.f_code is function.__code__ and event=='return':
            captured.update({k:frame.f_locals[k].copy() for k in ['selected','donors']})
        return tracer
    sys.settrace(tracer)
    try:result=function(*args)
    finally:sys.settrace(old)
    return result,captured

def validate_fixture(module,name,d,known):
    a=allowed_matrix(d,known);opt,witness,n_feasible,n_examined=brute_force(a)
    (retained,neg,stats),captured=recorded_call(module.corrected,d,known)
    selected=list(captured['selected']);donors=list(captured['donors'])
    assert len(retained)==len(neg)==opt,(name,len(retained),opt)
    assert len(set(selected))==len(selected) and set(selected)==set(donors),name
    assert all(i!=j and d.loc[i,'hla']==d.loc[j,'hla'] and
               (d.loc[i,'beta'],d.loc[j,'peptide'],d.loc[j,'hla']) not in known
               for i,j in zip(selected,donors)),name
    assert all(t not in known for t in neg),name
    assert Counter(retained.beta)==Counter(t[0] for t in neg),name
    assert Counter(retained.peptide)==Counter(t[1] for t in neg),name
    assert Counter(retained.hla)==Counter(t[2] for t in neg),name
    assert Counter(zip(retained.peptide,retained.hla))==Counter(t[1:] for t in neg),name
    assert all(h==t[2] for h,t in zip(retained.hla,neg)),name
    assert sum(s['n_retained'] for s in stats)==opt,name
    unrestricted=int((maximum_bipartite_matching(csr_matrix(a),perm_type='column')>=0).sum())
    assert unrestricted>=opt,name
    return dict(name=name,n=len(d),unique_betas=d.beta.nunique(),unique_peptides=d.peptide.nunique(),n_hlas=d.hla.nunique(),
                n_known=len(known),retained_actual=len(retained),retained_exhaustive_optimum=opt,
                unconstrained_bipartite_upper_bound=unrestricted,permutations_examined=n_examined,
                feasible_permutations=n_feasible,witness_permutation=witness,certificates_passed=True)

def fixture_from_edges(n,edges):
    d=pd.DataFrame([dict(beta=f'B{i}',peptide=f'P{i}',hla='H0',source_row_0based=i) for i in range(n)])
    known={(f'B{i}',f'P{j}','H0') for i in range(n) for j in range(n) if (i,j) not in edges or i==j}
    return d,known

def randomized_fixture(n,seed):
    rng=np.random.default_rng(seed)
    # Several fixtures deliberately contain repeated beta/peptide/HLA records;
    # the row-level contract remains well defined, including identical rows.
    n_beta=max(1,int(rng.integers(1,n+1)));n_pep=max(1,int(rng.integers(1,n+1)))
    n_hla=1+(seed%3==0)
    d=pd.DataFrame([dict(beta=f'B{rng.integers(n_beta)}',peptide=f'P{rng.integers(n_pep)}',hla=f'H{rng.integers(n_hla)}',source_row_0based=i) for i in range(n)])
    known=set(d[['beta','peptide','hla']].itertuples(index=False,name=None))
    for b in d.beta.unique():
        for p in d.peptide.unique():
            for h in d.hla.unique():
                if rng.random()<.23:known.add((b,p,h))
    # Exercise arbitrary unique DataFrame index labels and out-of-order rows.
    d.index=np.arange(n)*7+11
    d=d.iloc[rng.permutation(n)]
    return d,known

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--public-source',type=Path,default=ROOT.parent/'public_audit/run_public_audit.py');ap.add_argument('--out',type=Path,default=ROOT/'matching_validation_results.json');args=ap.parse_args()
    module=load_public(args.public_source);records=[]
    fixtures=[('directed_path_3',3,{(0,1),(1,2)}),('directed_cycle_3',3,{(0,1),(1,2),(2,0)}),
      ('cycle_plus_unreturnable_path',5,{(0,1),(1,0),(2,3),(3,4)}),
      ('disjoint_cycles_5',5,{(0,1),(1,0),(2,3),(3,4),(4,2)}),
      ('no_edges_4',4,set()),('all_offdiagonal_4',4,{(i,j) for i in range(4) for j in range(4) if i!=j})]
    for name,n,edges in fixtures:
        d,known=fixture_from_edges(n,edges);records.append(validate_fixture(module,name,d,known))
    for n in range(1,9):
        for k in range(20):
            seed=10000+n*100+k;d,known=randomized_fixture(n,seed)
            records.append(validate_fixture(module,f'repeated_identity_n{n}_seed{seed}',d,known))
    assert records[0]['unconstrained_bipartite_upper_bound']==2 and records[0]['retained_actual']==0
    result=dict(scope='Independent exhaustive combinatorial verification on synthetic tiny fixtures; no molecular or predictive validation.',
      public_source_sha256=hashlib.sha256(args.public_source.read_bytes()).hexdigest(),
      fixtures_tested=len(records),permutations_examined=sum(r['permutations_examined'] for r in records),
      all_certificates_passed=all(r['certificates_passed'] for r in records),
      repeated_beta_fixtures=sum(r['unique_betas']<r['n'] for r in records),
      repeated_peptide_fixtures=sum(r['unique_peptides']<r['n'] for r in records),
      multi_hla_fixtures=sum(r['n_hlas']>1 for r in records),
      fixtures_with_strict_upper_bound_gap=sum(r['unconstrained_bipartite_upper_bound']>r['retained_actual'] for r in records),
      proof='An identity permutation costs n. Any forbidden edge alone costs n+1, hence cannot occur in an optimum. Every allowed off-diagonal costs zero and each fixed point costs one, so minimum cost equals n minus maximum retained rows. Nonfixed cycles have identical source and donor row sets. The result optimizes closed row subsets within each HLA under the supplied known-positive reference only.',
      limitations=['This is maximum cardinality under the stated closed within-HLA row-assignment contract, not maximum cardinality among every possible distribution-balanced resampling design.',
      'Known-positive exclusion does not establish that unobserved pairs are biologically negative.',
      'Matching fixes the retained-cohort distribution; it does not preserve the original population if rows are excluded.',
      'Matching globally before train/test or CV splitting may create cross-fold dependence; construct independently inside frozen partitions.',
      'At most eight rows per fixture; exhaustive tests complement the mathematical proof and do not benchmark large-cohort memory/time.'],fixtures=records)
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='fixtures'},indent=2))
if __name__=='__main__':main()
