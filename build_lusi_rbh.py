import os, collections, time
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
K,TOP=5,6
def fasta(p):
    h,s=None,[]
    for l in open(p):
        if l.startswith('>'):
            if h: yield h.split()[0],''.join(s)
            h,s=l[1:].strip(),[]
        else: s.append(l.strip())
    if h: yield h.split()[0],''.join(s)
A={i:s for i,s in fasta('phylo/proteomes/auris_cladeI.faa') if 60<=len(s)<=2000}
B={i:s for i,s in fasta('phylo/proteomes/lusitaniae.faa')   if 60<=len(s)<=2000}
print(len(A),len(B),flush=True)
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score,al.extend_gap_score=-11,-1; al.mode='global'
def index(P):
    ix=collections.defaultdict(set)
    for i,s in P.items():
        for j in range(0,len(s)-K,2): ix[s[j:j+K]].add(i)
    return ix
def short(seq,ix):
    c=collections.Counter()
    for j in range(0,len(seq)-K,2): c.update(ix.get(seq[j:j+K],()))
    return [i for i,_ in c.most_common(TOP)]
def ident(x,y):
    try: a=al.align(x,y)[0]
    except Exception: return 0.0
    m=sum(sum(1 for k in range(e0-s0) if x[s0+k]==y[s1+k])
          for (s0,e0),(s1,e1) in zip(a.aligned[0],a.aligned[1]))
    return 100.0*m/min(len(x),len(y))
def bh(P,Q,ixQ,tag):
    o={}; t0=time.time()
    for n,(i,s) in enumerate(P.items()):
        best=(None,0.0)
        for j in short(s,ixQ):
            v=ident(s,Q[j])
            if v>best[1]: best=(j,v)
        if best[0] and best[1]>=30: o[i]=best
        if n%1000==0: print(f'{tag} {n}/{len(P)} {time.time()-t0:.0f}s',flush=True)
    return o
AB=bh(A,B,index(B),'A->B'); BA=bh(B,A,index(A),'B->A')
r=[(i,j,v) for i,(j,v) in AB.items() if BA.get(j,(None,))[0]==i]
pd.DataFrame(r).to_csv('gem/tables/rbh/pid_lusitaniae.tsv',sep='\t',header=False,index=False)
print(f'{len(r)} RBH, median identity {pd.DataFrame(r)[2].median():.1f}%')
