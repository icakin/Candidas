import argparse, itertools, os, subprocess, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
AA=set('ACDEFGHIKLMNPQRSTVWY'); SLOPE=937.0; PROT='phylo/proteomes'
LIMIT={'auris_cladeI':44,'auris_cladeII':44,'auris_cladeIII':44,'auris_cladeIV':44,
       'haemulonii':38,'duobushaemulonii':38,'parapsilosis':40,'pseudohaemulonii':None}
def fasta(p):
    h,s=None,[]
    for l in open(p):
        if l.startswith('>'):
            if h: yield h.split()[0],''.join(s)
            h,s=l[1:].strip(),[]
        else: s.append(l.strip())
    if h: yield h.split()[0],''.join(s)
def ivywrel(seq):
    n=k=0
    for c in seq:
        if c in AA:
            n+=1
            if c in 'IVYWREL': k+=1
    return k/n if n>=50 else None
def blast(q,s,cache):
    out=os.path.join(cache,f'{q}__{s}.tsv')
    if not os.path.exists(out):
        db=os.path.join(cache,s)
        if not os.path.exists(db+'.dmnd'):
            subprocess.run(['diamond','makedb','--in',f'{PROT}/{s}.faa','-d',db,'--quiet'],check=True)
        subprocess.run(['diamond','blastp','-q',f'{PROT}/{q}.faa','-d',db,'-o',out,'--quiet',
            '--sensitive','--max-target-seqs','1','--evalue','1e-10','--outfmt','6',
            'qseqid','sseqid','pident'],check=True)
    d=pd.read_csv(out,sep='\t',header=None,names=['q','s','pid']).drop_duplicates('q')
    return dict(zip(d.q,zip(d.s,d.pid)))
def rbh(a,b,c):
    ab,ba=blast(a,b,c),blast(b,a,c)
    return [(x,y,p) for x,(y,p) in ab.items() if ba.get(y,(None,))[0]==x]
ap=argparse.ArgumentParser(); ap.add_argument('--cache',default='rbh0'); a=ap.parse_args()
os.makedirs(a.cache,exist_ok=True)
for g in LIMIT:
    if not os.path.exists(f'{PROT}/{g}.faa'): sys.exit(f'missing {PROT}/{g}.faa')
F={g:{i:ivywrel(s) for i,s in fasta(f'{PROT}/{g}.faa')} for g in LIMIT}
print(f"{len(F)} proteomes indexed\n"); rng=np.random.default_rng(0); rows=[]
for x,y in itertools.combinations(list(LIMIT),2):
    pr=rbh(x,y,a.cache)
    d=np.array([F[x][p]-F[y][q] for p,q,_ in pr if F[x].get(p) is not None and F[y].get(q) is not None])
    pid=np.median([p for _,_,p in pr])
    bs=np.array([d[rng.integers(0,len(d),len(d))].mean() for _ in range(2000)])
    lo,hi=np.percentile(bs,[2.5,97.5]); lx,ly=LIMIT[x],LIMIT[y]
    dlim=(lx-ly) if (lx is not None and ly is not None) else np.nan
    kind='NULL' if dlim==0 else ('test' if dlim==dlim else 'distance')
    rows.append(dict(a=x,b=y,n=len(d),pid=pid,dT=SLOPE*d.mean(),lo=SLOPE*lo,hi=SLOPE*hi,
                     dlim=dlim,kind=kind))
    print(f"  {x:16s} vs {y:16s} n={len(d):5d} id={pid:5.1f}%  dT={SLOPE*d.mean():+6.2f} "
          f"[{SLOPE*lo:+6.2f},{SLOPE*hi:+6.2f}]  {kind}")
R=pd.DataFrame(rows); R.to_csv('stage0_null_control.csv',index=False)
nul=R[R.kind=='NULL']; tst=R[R.kind=='test']
print("\n"+"="*72)
print(f"NULL (true answer 0): n={len(nul)}  |dT| mean {nul.dT.abs().mean():.2f}, max {nul.dT.abs().max():.2f}")
print(f"TEST (auris vs relatives): |dT| mean {tst.dT.abs().mean():.2f}, max {tst.dT.abs().max():.2f}")
both=R[R.dlim.notna()]
print(f"\n  |dT| vs identity  {np.corrcoef(both.pid,both.dT.abs())[0,1]:+.3f}  (negative = time)")
print(f"  |dT| vs limit gap {np.corrcoef(both.dlim.abs(),both.dT.abs())[0,1]:+.3f}  (positive = temperature)")
print("\n  identity-matched comparison:")
for _,t in tst.iterrows():
    near=nul.iloc[(nul.pid-t.pid).abs().argsort()[:2]]
    print(f"    {t.a} vs {t.b} (id {t.pid:.1f}%) dT {t.dT:+.2f}  |  nulls: "+", ".join(
        f"{n.a[:12]}/{n.b[:12]} id {n.pid:.1f}% dT {n.dT:+.2f}" for _,n in near.iterrows()))
