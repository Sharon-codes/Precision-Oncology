#!/usr/bin/env python3
"""Executable audit of source logic, not neural-model inference.

Uses the original, unmodified AST class body of TCRCollate. A minimal NumPy
shim implements only stack/cat/tensor; synthetic identity-tag embeddings let
us recover which sequence the original sampler selected. Nothing here
executes ESM, Mamba, learned checkpoints, or original study cohorts.
"""
import argparse, ast, csv, hashlib, json, logging, random, subprocess, sys
from collections import defaultdict, Counter
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

class TorchShim:
    float = np.float64
    @staticmethod
    def stack(x): return np.stack(x)
    @staticmethod
    def cat(x, dim=0): return np.concatenate(x, axis=dim)
    @staticmethod
    def tensor(x, dtype=None): return np.array(x, dtype=dtype)

class SequenceTags:
    def __init__(self): self.keys = {}; self.reverse = {}
    def load(self, value):
        if value not in self.keys:
            tag = len(self.keys) + 1
            self.keys[value] = tag; self.reverse[tag] = value
        return np.array([self.keys[value]], dtype=np.int64)
    def decode(self, tensor): return [self.reverse[int(row[0])] for row in tensor]

class AuditRandom:
    def __init__(self, seed): self.rng = random.Random(seed); self.choice_calls = 0
    def shuffle(self, x): self.rng.shuffle(x)
    def choice(self, x): self.choice_calls += 1; return self.rng.choice(x)

def load_ast(path, class_name=None, functions=None, methods=None, namespace=None):
    source = Path(path).read_text()
    tree = ast.parse(source)
    if class_name:
        nodes = [n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == class_name]
        if methods: nodes[0].body = [n for n in nodes[0].body if isinstance(n, ast.FunctionDef) and n.name in methods]
    else:
        nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in functions]
    ns = {} if namespace is None else namespace
    exec(compile(ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[])), str(path), 'exec'), ns)
    return ns

def lev(a,b):
    r=list(range(len(b)+1))
    for i,x in enumerate(a):
        n=[i+1]
        for j,y in enumerate(b): n.append(min(n[-1]+1,r[j+1]+1,r[j]+(x!=y)))
        r=n
    return r[-1]

def item(beta, peptide, hla='A*02:01'):
    return dict(raw_cdr3_beta=beta, raw_cdr3_alpha='', raw_peptide=peptide, raw_hla_allele=hla, disease='synthetic')

def source_metadata(repo):
    names=['dataset.py','preprocess.py','train.py','run_official_sota_replications.py']
    hashes={n:hashlib.sha256((repo/n).read_bytes()).hexdigest() for n in names}
    if (repo/'.git').exists():
        commit=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
    else:
        # Archive portability only: selection logic still executes the original
        # class AST. Validate copied source bytes against recorded Git provenance.
        manifest=json.loads((repo/'SOURCE_PROVENANCE.json').read_text())
        if hashes!=manifest['files']:
            raise ValueError('Archived reference source does not match pinned fingerprints')
        commit=manifest['commit']
    return dict(commit=commit,files=hashes,
                python=sys.version, numpy=np.__version__, pandas=pd.__version__)

def make_collate(repo, seed):
    tags=SequenceTags(); rng=AuditRandom(seed)
    ns=load_ast(repo/'dataset.py',class_name='TCRCollate',namespace={
      'torch':TorchShim,'load_embedding':tags.load,'get_hla_pseudo':lambda h:'|HLA='+h,
      'random':rng,'print':lambda *a,**kw:None})
    return ns['TCRCollate'], tags, rng

def decode_pairs(batch, out, tags):
    n=len(batch); keys=tags.decode(out['peptide_plus_hla'][n:])
    return [(batch[i]['raw_cdr3_beta'],*key.split('|HLA=')) for i,key in enumerate(keys)]

def counterexamples(repo):
    p0,p1,p2='AAAAAAAAA','CCAAAAAAA','CCCCAAAAA'
    preliminary_train=[p0]; preliminary_test=[p1,p2]
    retained=[p for p in preliminary_test if all(lev(p,q)>2 for q in preliminary_train)]
    final_train=[p for p in preliminary_train+preliminary_test if p not in retained]
    violations=[dict(test=p,train=q,distance=lev(p,q)) for p in retained for q in final_train if lev(p,q)<=2]
    assert retained == [p2] and violations
    logger=logging.getLogger('audit'); logger.addHandler(logging.NullHandler())
    ns=load_ast(repo/'preprocess.py',class_name='TCRDataPipeline',methods={'levenshtein_distance','cluster_peptides','split_data'},namespace={
      'np':np,'defaultdict':defaultdict,'tqdm':lambda x,**kw:x,'logger':logger})
    pipeline=ns['TCRDataPipeline'](); pipeline.shapes_tracker={}
    components=pipeline.cluster_peptides([p0,p1,p2])
    assert len(components)==1
    # Synthetic alphabetically distinct extra components allow three actual partitions.
    peptides=[p0,p1,p2]+[aa*9 for aa in 'DEFGHIKLMNPQRSTVWY']
    df=pd.DataFrame([dict(cdr3_beta=f'CASS{i}F',cdr3_alpha='',peptide=p,hla_allele='A*02:01',disease='synthetic') for i,p in enumerate(peptides)])
    partitions=pipeline.split_data(df)
    cross_distances=[lev(a,b) for i in range(3) for j in range(i+1,3) for a in partitions[i].peptide for b in partitions[j].peptide]
    assert all(d>2 for d in cross_distances)
    # Run the exact leave_one_disease branch and exact helper functions from train.py.
    train_ns=load_ast(repo/'train.py',functions={'levenshtein_distance','cluster_peptides','split_by_peptide_clusters'},namespace={
      'np':np,'pd':pd,'defaultdict':defaultdict})
    tree=ast.parse((repo/'train.py').read_text())
    branch=next(n for n in ast.walk(tree) if isinstance(n,ast.If) and n.lineno < 300 and isinstance(n.test,ast.Compare)
                and any(isinstance(c,ast.Constant) and c.value=='leave_one_disease' for c in n.test.comparators))
    df2=pd.DataFrame([
      dict(cdr3_beta='CASSAF',peptide=p0,hla_allele='A*02:01',disease='background'),
      dict(cdr3_beta='CASSCF',peptide='CAAAAAAAA',hla_allele='A*02:01',disease='HeldOut')])
    train_ns.update(df_all=df2,args=SimpleNamespace(held_out_disease='HeldOut'),logger=logger)
    exec(compile(ast.Module(body=branch.body,type_ignores=[]),str(repo/'train.py'),'exec'),train_ns)
    disease_min=min(lev(a,b) for a in train_ns['train_df'].peptide for b in train_ns['test_df'].peptide)
    assert disease_min==1
    # A two-item derangement necessarily maps every beta to the other peptide.
    Collate,tags,rng=make_collate(repo,42)
    batch=[item('CASSAF',p0),item('CASSCF',p2)]
    known={(b,p,'A*02:01') for b in ['CASSAF','CASSCF'] for p in [p0,p2]}
    pool=[(p0,'A*02:01'),(p2,'A*02:01')]
    out=Collate(known,pool)(batch); pairs=decode_pairs(batch,out,tags)
    assert all(p in known for p in pairs)
    # Same RNG and batch yield exactly equal labels/tensors regardless of supplied positive set.
    C2,tags2,rng2=make_collate(repo,42); out2=C2(set(),pool)(batch)
    assert decode_pairs(batch,out2,tags2)==pairs
    # Singleton/one-peptide pool: hardcoded peptide fallback (not a biological negative).
    C3,tags3,rng3=make_collate(repo,17); singleton=[item('CASSAF',p0)]
    singleout=C3({('CASSAF',p0,'A*02:01')},[(p0,'A*02:01')])(singleton)
    singlepairs=decode_pairs(singleton,singleout,tags3)
    assert singlepairs[0][1]=='AAAAA' and rng3.choice_calls==20
    tree=ast.parse((repo/'dataset.py').read_text())
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='TCRCollate')
    attrs=[n for n in ast.walk(cls) if isinstance(n,ast.Attribute) and n.attr=='pos_set']
    reads=[n.lineno for n in attrs if isinstance(n.ctx,ast.Load)]
    writes=[n.lineno for n in attrs if isinstance(n.ctx,ast.Store)]
    assert reads==[] and len(writes)==1
    return {
      'manuscript_algorithm_counterexample':dict(preliminary_train=preliminary_train,preliminary_test=preliminary_test,retained_test=retained,final_train=final_train,violations=violations),
      'repo_component_split_control':dict(chain_components=components,partition_sizes=[len(d) for d in partitions],minimum_cross_partition_distance=min(cross_distances),cross_partition_pairs_checked=len(cross_distances)),
      'repo_disease_holdout_counterexample':dict(train_peptides=train_ns['train_df'].peptide.tolist(),test_peptides=train_ns['test_df'].peptide.tolist(),minimum_distance=disease_min,source_branch_lines=[branch.lineno,branch.end_lineno]),
      'known_positive_collision_counterexample':dict(generated_negatives=pairs,known_positive_collisions=len(pairs),negatives=len(pairs),same_output_with_empty_positive_set=True,pos_set_load_lines=reads,pos_set_store_lines=writes),
      'single_peptide_pool_fallback':dict(generated_negative=singlepairs[0],global_choices=rng3.choice_calls)
    }

def run_batch_order_simulation(repo,seeds):
    # Predeclared design: 4096 synthetic positive records; 9 peptide categories;
    # one at 50%, eight at 6.25% each. All beta IDs unique; one HLA.
    # Prespecified nuisance score is 1 iff peptide is the majority category.
    peptides=[aa*9 for aa in 'ACDEFGHIK']; counts=[2048]+[256]*8
    records=[]
    for p,n in zip(peptides,counts):
        for _ in range(n): records.append(item(f'CASS{len(records):05d}F',p))
    known={(r['raw_cdr3_beta'],r['raw_peptide'],r['raw_hla_allele']) for r in records}
    pool=[(p,'A*02:01') for p in peptides]; n=len(records); rows=[]
    for seed in range(seeds):
        for order in ['peptide_sorted','row_shuffled']:
            ordered=records.copy()
            if order=='row_shuffled': random.Random(seed+100000).shuffle(ordered)
            C,tags,rng=make_collate(repo,seed); collate=C(known,pool)
            positives=[]; negatives=[]; collisions=0; hardcoded=0; global_fallbacks=0
            # The original class calls load_embedding an extra time per fallback;
            # recover their number with call counting, without changing its AST.
            original_load=tags.load
            calls=[0]
            def counted_load(x): calls[0]+=1; return original_load(x)
            collate.__call__.__globals__['load_embedding']=counted_load
            for k in range(0,n,64):
                batch=ordered[k:k+64]; start=calls[0]
                out=collate(batch)
                global_fallbacks+=calls[0]-start-3*len(batch)
                pairs=decode_pairs(batch,out,tags)
                positives.extend(r['raw_peptide'] for r in batch)
                negatives.extend(p[1] for p in pairs)
                collisions+=sum(p in known for p in pairs)
                hardcoded+=sum(p[1]=='AAAAA' for p in pairs)
            y=np.r_[np.ones(n),np.zeros(n)]
            score=np.array([p==peptides[0] for p in positives+negatives],dtype=float)
            posfreq=Counter(positives); negfreq=Counter(negatives)
            tv=.5*sum(abs(posfreq[p]/n-negfreq[p]/n) for p in set(posfreq)|set(negfreq))
            rows.append(dict(seed=seed,order=order,n_positive=n,n_negative=n,majority_score_auroc=float(roc_auc_score(y,score)),
              majority_negative_fraction=negfreq[peptides[0]]/n,global_fallback_fraction=global_fallbacks/n,
              fallback_choice_calls=rng.choice_calls,hardcoded_fallbacks=hardcoded,known_positive_collisions=collisions,total_variation_peptide_marginal=tv))
    # Deterministic within-split cyclic matching preserves the complete peptide
    # marginal exactly for this feasible synthetic design. Not a general solver.
    shift=max(counts); cyclic=records[shift:]+records[:shift]
    control_neg=[r['raw_peptide'] for r in cyclic]
    assert all(r['raw_peptide']!=q['raw_peptide'] for r,q in zip(records,cyclic))
    assert Counter(control_neg)==Counter(r['raw_peptide'] for r in records)
    control_auc=roc_auc_score(np.r_[np.ones(n),np.zeros(n)],np.array([r['raw_peptide']==peptides[0] for r in records]+[p==peptides[0] for p in control_neg]))
    summary={}
    for order in ['peptide_sorted','row_shuffled']:
        sub=[r for r in rows if r['order']==order]
        summary[order]={}
        for metric in ['majority_score_auroc','global_fallback_fraction','majority_negative_fraction','total_variation_peptide_marginal']:
            a=np.array([r[metric] for r in sub]); summary[order][metric]=dict(mean=float(a.mean()),sd=float(a.std(ddof=1)),min=float(a.min()),max=float(a.max()),q025=float(np.quantile(a,.025)),q975=float(np.quantile(a,.975)))
        summary[order]['hardcoded_fallbacks_total']=sum(r['hardcoded_fallbacks'] for r in sub)
        summary[order]['known_positive_collisions_total']=sum(r['known_positive_collisions'] for r in sub)
    a=np.array([r['majority_score_auroc'] for r in rows if r['order']=='peptide_sorted'])
    b=np.array([r['majority_score_auroc'] for r in rows if r['order']=='row_shuffled'])
    summary['paired_order_difference']=dict(mean=float((a-b).mean()),sd=float((a-b).std(ddof=1)),q025=float(np.quantile(a-b,.025)),q975=float(np.quantile(a-b,.975)))
    summary['marginal_preserving_cyclic_control']=dict(majority_score_auroc=float(control_auc),known_positive_collisions=0,peptide_marginal_total_variation=0.0)
    summary['design']=dict(synthetic=True,n_records=n,n_peptides=len(peptides),counts=counts,batch_size=64,seeds=list(range(seeds)),score='indicator(peptide == majority peptide)',uncertainty='Monte Carlo distribution across sampler seeds; not a population confidence interval')
    return rows,summary

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo',type=Path,default=Path(__file__).resolve().parent.parent/'precision_repo'); ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parent); ap.add_argument('--seeds',type=int,default=100); args=ap.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    results=dict(scope='Source-logic execution and explicitly synthetic simulations only; no original-cohort reconstruction and no neural-model evaluation.',provenance=source_metadata(args.repo),counterexamples=counterexamples(args.repo))
    rows,summary=run_batch_order_simulation(args.repo,args.seeds); results['synthetic_batch_order_simulation']=summary
    (args.out/'source_audit_results.json').write_text(json.dumps(results,indent=2)+'\n')
    with (args.out/'source_audit_seed_results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(results,indent=2))
if __name__=='__main__': main()
