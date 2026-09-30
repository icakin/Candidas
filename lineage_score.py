import os, warnings, collections
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
AA=set('ACDEFGHIKLMNPQRSTVWY')
STAB={('K','R'),('L','I'),('V','I'),('T','I'),('M','I'),('G','A'),('S','A'),
      ('N','D'),('Q','E')}
DESTAB={(b,a) for a,b in STAB}
LIM={'auris':44,'haemulonii':38,'duobushaemulonii':38,'parapsilosis':40}
def fasta(p):
    h,s=None,[]
    for l in open(p):
        if l.startswith('>'):
            if h: yield h.split()[0],''.join(s)
            h,s=l[1:].strip(),[]
        else: s.append(l.strip())
    if h: yield h.split()[0],''.join(s)
AU={i:s for i,s in fasta('phylo/proteomes/auris_cladeI.faa')}
SEQD={'haemulonii':{i:s for i,s in fasta('phylo/proteomes/haemulonii.faa')},
      'duobushaemulonii':{i:s for i,s in fasta('phylo/proteomes/duobushaemulonii.faa')},
      'lusitaniae':{i:s for i,s in fasta('phylo/proteomes/lusitaniae.faa')}}
pu=pd.read_csv('gem/inputs/parap_uniprot.tsv',sep='\t',header=None,usecols=[1,2],names=['l','s'])
SEQD['parapsilosis']=dict(zip(pu.l,pu.s))
rd=lambda s: dict(pd.read_csv(f'gem/tables/rbh/pid_{s}.tsv',sep='\t',header=None,
                              names=['a','b','p'])[['a','b']].values)
O={s:rd(s) for s in ['haemulonii','duobushaemulonii','parapsilosis','lusitaniae']}
RIB=set(h.split('|')[1] for h,_ in fasta('gem/tables/machinery_ribosome.faa')
        if h.startswith('auris|'))
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score,al.extend_gap_score=-11,-1; al.mode='global'
def proj(x,y):
    a=al.align(x,y)[0]; o=['-']*len(x)
    for (s0,e0),(s1,e1) in zip(a.aligned[0],a.aligned[1]):
        for k in range(e0-s0): o[s0+k]=y[s1+k]
    return o
genes=[g for g in AU if all(g in O[s] and O[s][g] in SEQD[s] for s in O)]
print(f"{len(genes)} genes in all five ({len(set(genes)&RIB)} ribosomal)",flush=True)
SP=['auris','haemulonii','duobushaemulonii','parapsilosis']
C={ (s,c):collections.Counter() for s in SP for c in ('rib','all') }
for n,g in enumerate(genes):
    s=AU[g]
    if not (60<=len(s)<=2500): continue
    try:
        seq={'auris':list(s)}
        for k in SEQD: seq[k]=proj(s,SEQD[k][O[k][g]])
    except Exception: continue
    cls=('rib' if g in RIB else 'all')
    for i in range(len(s)):
        v={k:seq[k][i] for k in seq}
        if any(x not in AA for x in v.values()): continue
        for X in SP:
            others=[v[k] for k in list(SP)+['lusitaniae'] if k!=X]
            if len(set(others))==1 and others[0]!=v[X]:
                anc,der=others[0],v[X]
                for c in ({cls,'all'} if cls=='rib' else {'all'}):
                    C[(X,c)]['n']+=1
                    if (anc,der) in STAB: C[(X,c)]['s']+=1
                    if (anc,der) in DESTAB: C[(X,c)]['d']+=1
    if n%800==0: print(f"  {n}/{len(genes)}",flush=True)
print(f"\n{'species':>18}{'limit':>7}{'subs':>8}{'stab':>7}{'destab':>8}{'net score':>11}")
rows=[]
for X in SP:
    c=C[(X,'all')]; sc=(c['s']-c['d'])/max(c['n'],1)
    rows.append((X,LIM[X],sc))
    print(f"{X:>18}{LIM[X]:>7}{c['n']:>8}{c['s']:>7}{c['d']:>8}{sc:>11.4f}")
R=pd.DataFrame(rows,columns=['sp','lim','score'])
print(f"\n  score vs limit: r = {np.corrcoef(R.score,R.lim)[0,1]:+.3f}")
h,d=R[R.sp=='haemulonii'].score.iloc[0],R[R.sp=='duobushaemulonii'].score.iloc[0]
print(f"  NULL (hae - duo, same limit) = {h-d:+.4f}   [should be ~0]")
print(f"\n{'species':>18}{'ribosomal net score':>22}")
for X in SP:
    c=C[(X,'rib')]
    if c['n']>50: print(f"{X:>18}{(c['s']-c['d'])/c['n']:>22.4f}   (n={c['n']})")
