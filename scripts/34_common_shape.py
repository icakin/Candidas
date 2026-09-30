"""Is the temperature dependence of respiration-per-unit-growth ONE shape, shifted between
taxa, or do taxa differ in shape as well?

y = log(K/r) per growing well.  Biophysical form (same as the pipeline):
    y = a_taxon + dE*CINV*(1/TREF - 1/TK) + log1p(exp(Eh*CINV*(1/Th - 1/TK)))
  dE = E_R - E  (low-temperature slope: shape)
  Eh              (steepness of the high-temperature collapse: width)
  Th              (position on the temperature axis: the shift)

  M1 translation : dE, Eh SHARED ; Th, a per taxon        -> pure horizontal shift
  M2 shift+width : dE SHARED     ; Eh, Th, a per taxon
  M3 free        : dE, Eh, Th, a per taxon

Model comparison by LEAVE-ONE-ISOLATE-OUT cross-validation: every well of the held-out
isolate, at every temperature, is withheld together (no pseudoreplication). Scored by mean
log predictive density and RMSE. Retain 'translation' only if M1 predicts held-out isolates
essentially as well as M3.
"""
import os, numpy as np, pandas as pd
from scipy.optimize import least_squares
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15
cd=pd.read_csv(f"{C}/results/tables/fit_coefficients_wide.csv")
on=pd.read_csv(f"{C}/results/tables/otu_names.csv")
d=cd.merge(on,on="OTU"); d=d[(d.keep)&(d.fit_valid)&(d.r>0)&(d.K>0)].copy()
d["y"]=np.log(d.K/d.r); d["TK"]=d["T"]+273.15
TAX=sorted(d.group.unique()); ISO=sorted(d.otu_name.unique())
print(f"{len(d)} growing wells | {len(ISO)} isolates | {len(TAX)} taxa")
ti={t:i for i,t in enumerate(TAX)}
def predict(p,df,model):
    k=df.group.map(ti).values; TK=df.TK.values; n=len(TAX)
    a=p[:n]
    if model==1:   dE=np.full(n,p[n]); Eh=np.full(n,p[n+1]); Th=p[n+2:n+2+n]
    elif model==2: dE=np.full(n,p[n]); Eh=p[n+1:n+1+n];      Th=p[n+1+n:n+1+2*n]
    else:          dE=p[n:2*n];        Eh=p[2*n:3*n];        Th=p[3*n:4*n]
    u=np.clip(Eh[k]*CINV*(1/Th[k]-1/TK),-50,50)
    return a[k]+dE[k]*CINV*(1/TREF-1/TK)+np.log1p(np.exp(u))
def start(model):
    n=len(TAX); a=[np.log(d[d.group==t].K.median()/d[d.group==t].r.median()) for t in TAX]
    if model==1:   return np.array(a+[-0.4,3.0]+[309.0]*n)
    elif model==2: return np.array(a+[-0.4]+[3.0]*n+[309.0]*n)
    else:          return np.array(a+[-0.4]*n+[3.0]*n+[309.0]*n)
def bounds(model):
    n=len(TAX)
    lo=[-20]*n; hi=[20]*n
    if model==1:   lo+=[-4,0.3]+[295]*n;      hi+=[4,40]+[330]*n
    elif model==2: lo+=[-4]+[0.3]*n+[295]*n;  hi+=[4]+[40]*n+[330]*n
    else:          lo+=[-4]*n+[0.3]*n+[295]*n; hi+=[4]*n+[40]*n+[330]*n
    return (lo,hi)
def fit(df,model):
    return least_squares(lambda p: predict(p,df,model)-df.y.values,
                         start(model),bounds=bounds(model),max_nfev=20000)
res={}
for m in (1,2,3):
    f=fit(d,m); n=len(d); rss=np.sum(f.fun**2); k=len(f.x)
    res[m]=dict(rmse=np.sqrt(rss/n),k=k,aic=n*np.log(rss/n)+2*k)
    print(f"  in-sample M{m}: RMSE {res[m]['rmse']:.4f}  params {k}  AIC {res[m]['aic']:.1f}")
print("\nleave-one-ISOLATE-out (all temperatures and replicates of the isolate held out together):")
rows=[]
for iso in ISO:
    tr=d[d.otu_name!=iso]; te=d[d.otu_name==iso]
    if te.empty or tr[tr.group==te.group.iloc[0]].otu_name.nunique()<1: continue
    r={}
    for m in (1,2,3):
        f=fit(tr,m); resid=predict(f.x,tr,m)-tr.y.values; sd=np.sqrt(np.mean(resid**2))
        e=predict(f.x,te,m)-te.y.values
        r[f"rmse{m}"]=np.sqrt(np.mean(e**2))
        r[f"lpd{m}"]=np.mean(-0.5*np.log(2*np.pi*sd**2)-e**2/(2*sd**2))
    rows.append(dict(isolate=iso,taxon=te.group.iloc[0],n=len(te),**r))
cv=pd.DataFrame(rows); cv.to_csv(f"{C}/results/tables/common_shape_loio.csv",index=False)
print(cv.round(3).to_string(index=False)); print()
for m,nm in ((1,"M1 translation"),(2,"M2 shift+width"),(3,"M3 free shape")):
    print(f"  {nm:16s}  mean LPD {cv[f'lpd{m}'].mean():+.4f}   mean RMSE {cv[f'rmse{m}'].mean():.4f}")
dl31=cv.lpd3-cv.lpd1; dl21=cv.lpd2-cv.lpd1
se=lambda v: v.std(ddof=1)/np.sqrt(len(v))
print(f"\n  M3 - M1 : dLPD {dl31.mean():+.4f} +/- {se(dl31):.4f} (SE over isolates)  isolates favouring M3: {(dl31>0).sum()}/{len(dl31)}")
print(f"  M2 - M1 : dLPD {dl21.mean():+.4f} +/- {se(dl21):.4f}                       isolates favouring M2: {(dl21>0).sum()}/{len(dl21)}")
