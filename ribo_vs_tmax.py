import os, re, glob, numpy as np, pandas as pd
KEEP=re.compile(r'(40S|60S) ribosomal protein|ribosomal protein [SL]\d', re.I)
DROP=re.compile(r'mitochondrial|37S|54S|chloroplast|ribosome biogenesis|assembly', re.I)
def fasta(p):
    h,s=None,[]
    for l in open(p):
        if l.startswith('>'):
            if h: yield h,''.join(s)
            h,s=l[1:].strip(),[]
        else: s.append(l.strip())
    if h: yield h,''.join(s)
rows=[]
for f in glob.glob('phylo/ribo/*.faa'):
    sp=os.path.basename(f)[:-4].replace('_',' ')
    R=K=n=aR=aK=aN=0
    for h,s in fasta(f):
        aR+=s.count('R'); aK+=s.count('K'); aN+=1
        if KEEP.search(h) and not DROP.search(h) and 50<=len(s)<=1000:
            R+=s.count('R'); K+=s.count('K'); n+=1
    if n>=30: rows.append(dict(sp=sp,n_cyto=n,n_all=aN,rib_RRK=R/(R+K),
                               all_RRK=aR/(aR+aK)))
D=pd.DataFrame(rows)
M=pd.read_csv('phylo/inputs/thermal_manifest.csv')
M['sp']=M.species.astype(str).str.split().str[:2].str.join(' ')
M['tmax']=pd.to_numeric(M.tmax,errors='coerce')
X=D.merge(M[['sp','tmax','division','order']].drop_duplicates('sp'),on='sp')
X.to_csv('phylo/ribo_vs_tmax.csv',index=False)
print(f"{len(D)} species passed the filter, {len(X)} matched to a Tmax")
print(f"Tmax: min {X.tmax.min():.0f} max {X.tmax.max():.0f} median {X.tmax.median():.0f}")
print(f"cytosolic ribosomal proteins per species: median {X.n_cyto.median():.0f}\n")
def rep(lab,d):
    if len(d)<12: return
    r=np.corrcoef(d.rib_RRK,d.tmax)[0,1]; r2=np.corrcoef(d.all_RRK,d.tmax)[0,1]
    b=np.polyfit(d.rib_RRK,d.tmax,1)[0]
    print(f"{lab:<26} n={len(d):<4} r(cyto)={r:+.3f}  r(unfiltered)={r2:+.3f}  "
          f"slope={b:6.0f} C per unit")
rep('ALL',X)
for k,d in X.groupby(X.division.astype(str).str.strip()):
    rep(f'  division {k}',d)
for k,d in X.groupby(X.order.astype(str).str.strip()):
    rep(f'  order {k}',d)
