#!/usr/bin/env python3
"""Rebuilds vector audit figures from measured source-audit outputs."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
df=pd.read_csv(ROOT/'source_audit_seed_results.csv')
colors=['#276E8E','#D77D37','#32816C']; labels=['Peptide sorted','Rows shuffled','Fixed matching']
fig,axs=plt.subplots(1,3,figsize=(10.2,3.35),layout='constrained')
for i,(metric,title,ylabel,scale) in enumerate([
 ('majority_score_auroc','A  Peptide-only nuisance score','AUROC',1),
 ('global_fallback_fraction','B  Resorting changes the sampler','Global fallback (%)',100),
 ('majority_negative_fraction','C  Negative peptide distribution','Majority peptide among negatives (%)',100)]):
 ax=axs[i]
 for j,order in enumerate(['peptide_sorted','row_shuffled']):
  a=df[df.order==order][metric].to_numpy()*scale
  jitter=np.random.default_rng(246+j).uniform(-.1,.1,len(a))
  ax.scatter(j+jitter,a,s=9,color=colors[j],alpha=.35,linewidths=0)
  ax.plot([j-.14,j+.14],[a.mean(),a.mean()],color=colors[j],lw=2)
  ax.annotate(f'{a.mean():.3f}' if scale==1 else f'{a.mean():.2f}%',(j,a.mean()),xytext=(0,13),textcoords='offset points',ha='center',weight='bold',fontsize=9)
 control={'majority_score_auroc':.5,'global_fallback_fraction':0,'majority_negative_fraction':.5}[metric]*scale
 ax.scatter([2],[control],marker='D',s=32,color=colors[2],zorder=3)
 ax.annotate(f'{control:.3f}' if scale==1 else f'{control:.2f}%',(2,control),xytext=(0,13),textcoords='offset points',ha='center',weight='bold',fontsize=9)
 ax.set_xticks(range(3),labels,rotation=18,ha='right');ax.set_ylabel(ylabel);ax.set_title(title,loc='left',fontsize=10,weight='bold');ax.grid(axis='y',alpha=.15);ax.set_xlim(-.4,2.5)
 if i==0: ax.set_ylim(.46,.8);ax.axhline(.5,color='#777',lw=.8,ls='--',zorder=0)
 elif i==1:ax.set_ylim(-5,119)
 else:ax.set_ylim(0,60);ax.axhline(50,color='#777',lw=.8,ls='--',zorder=0)
fig.suptitle('Synthetic control: 4,096 pairs · 9 peptide categories · batch size 64 · 100 seeds',fontsize=10,y=1.035)
fig.savefig(ROOT/'synthetic_sampler_audit.pdf',bbox_inches='tight');fig.savefig(ROOT/'synthetic_sampler_audit.png',dpi=220,bbox_inches='tight');plt.close(fig)

fig,ax=plt.subplots(figsize=(9.5,3.6));ax.set_xlim(0,10);ax.set_ylim(0,4);ax.axis('off')
xs=[1.8,5,8.2]; seqs=['AAAAAAAAA','CCAAAAAAA','CCCCAAAAA']; roles=['Initial training','Candidate test → rejected','Candidate test → retained']
for i,(x,s,r) in enumerate(zip(xs,seqs,roles)):
 face='#DAE9EE' if i<2 else '#FBE5D7'
 rect=FancyBboxPatch((x-1.35,1.6),2.7,.8,boxstyle='round,pad=0.05,rounding_size=.08',facecolor=face,edgecolor=colors[0] if i<2 else colors[1])
 ax.add_patch(rect);ax.text(x,2,s,ha='center',va='center',family='DejaVu Sans Mono',fontsize=11,weight='bold');ax.text(x,1.3,r,ha='center',va='center',fontsize=9)
ax.annotate('',xy=(3.6,2),xytext=(3.15,2),arrowprops=dict(arrowstyle='<->',color='#444'));ax.text(3.4,2.45,'d = 2',ha='center')
ax.annotate('',xy=(6.8,2),xytext=(6.35,2),arrowprops=dict(arrowstyle='<->',color=colors[1],lw=2));ax.text(6.6,2.45,'d = 2',ha='center',color=colors[1],weight='bold')
ax.annotate('',xy=(8.2,3.1),xytext=(1.8,3.1),arrowprops=dict(arrowstyle='<->',color='#777'));ax.text(5,3.35,'d = 4: passes the preliminary-training check',ha='center',fontsize=9)
ax.text(5,.6,'Final split: rejected candidate returns to training, creating a train–test distance of 2.',ha='center',fontsize=10)
ax.text(5,.13,'Connected-component assignment keeps all three peptides in one partition.',ha='center',fontsize=10,color=colors[2])
fig.savefig(ROOT/'split_counterexample.pdf',bbox_inches='tight');fig.savefig(ROOT/'split_counterexample.png',dpi=220,bbox_inches='tight');plt.close(fig)
