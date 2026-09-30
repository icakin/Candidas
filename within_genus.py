import numpy as np, pandas as pd
from scipy import stats
X=pd.read_csv('phylo/ribo_vs_tmax.csv')
X['genus']=X.sp.str.split().str[0]
rows=[]
for g,d in X.groupby('genus'):
    d=d.dropna(subset=['rib_RRK','tmax'])
    if len(d)<3 or d.tmax.nunique()<2 or (d.tmax.max()-d.tmax.min())<3: continue
    r,p=stats.pearsonr(d.rib_RRK,d.tmax)
    b=np.polyfit(d.rib_RRK,d.tmax,1)[0]
    rows.append(dict(genus=g,n=len(d),span=d.tmax.max()-d.tmax.min(),r=r,p=p,slope=b))
R=pd.DataFrame(rows).sort_values('n',ascending=False)
print(R.round(3).to_string(index=False))
pos=(R.r>0).sum(); tot=len(R)
print(f"\n{tot} genera tested;  positive slopes {pos}, negative {tot-pos}")
print(f"  sign test p = {stats.binomtest(pos,tot,0.5).pvalue:.4f}")
print(f"  mean r = {R.r.mean():+.3f}   median r = {R.r.median():+.3f}")
w=R[R.n>=5]
if len(w): print(f"  restricted to genera with n>=5 ({len(w)}): "
                 f"mean r = {w.r.mean():+.3f}, positive {(w.r>0).sum()}/{len(w)}")
