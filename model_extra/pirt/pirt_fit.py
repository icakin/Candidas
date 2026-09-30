"""Pirt decomposition across temperature: R(T) = Y*G(T) + m0*exp(Em/k*(1/Tref-1/T)).
Per-cell G, R in fg C cell^-1 h^-1 from derived_N0_R_results_with_carbon.csv.
Fit in log space per isolate and per taxon (pooled wells). Compare with Y=0 (pure Arrhenius)."""
import pandas as pd, numpy as np, sys
from scipy.optimize import least_squares
k=8.617e-5; Tref=273.15+30
d=pd.read_csv('results/tables/derived_N0_R_results_with_carbon.csv')
o=pd.read_csv('results/tables/otu_names.csv'); d=d.merge(o[['OTU','group']],on='OTU')
d=d[d.group.isin(['Clade1','Clade2','Clade3','Clade4','para'])].copy()
d['Tk']=d['T']+273.15
def model(p,G,Tk):
    Y,lm0,Em=p
    return Y*G+np.exp(lm0)*np.exp(Em/k*(1/Tref-1/Tk))
def fit(s,fixY=None):
    G=s.growth_fgC_h.values; R=s.respiration_fgC_h.values; Tk=s.Tk.values
    if fixY is None:
        f=lambda p: np.log(model(p,G,Tk))-np.log(R)
        r=least_squares(f,[0.3,np.log(300),0.5],bounds=([0,-5,-1],[5,15,3]))
        p=r.x
    else:
        f=lambda p: np.log(model([fixY,*p],G,Tk))-np.log(R)
        r=least_squares(f,[np.log(300),0.5],bounds=([-5,-1],[15,3]))
        p=[fixY,*r.x]
    n=len(R); kk=3 if fixY is None else 2
    sse=np.sum(r.fun**2); aicc=n*np.log(sse/n)+2*kk+2*kk*(kk+1)/(n-kk-1)
    return p,sse,aicc,n
rows=[]
for lvl,key in [('taxon','group'),('isolate','otu_name')]:
    for g,s in d.groupby(key):
        p,sse,a1,n=fit(s); p0,sse0,a0,_=fit(s,fixY=0)
        grp=s.group.iloc[0]
        # maintenance fraction at 37 and 40 from taxon-median G at those T (interpolate medians)
        med=s.groupby('T').agg(G=('growth_fgC_h','median')).reset_index()
        out=dict(level=lvl,name=g,group=grp,n=n,Y=p[0],m0=np.exp(p[1]),Em=p[2],dAICc_vs_noY=a1-a0,Em_noY=p0[2])
        for T in (30,37,40):
            Gt=np.interp(T,med['T'],med.G)
            m=np.exp(p[1])*np.exp(p[2]/k*(1/Tref-1/(T+273.15)))
            out[f'mfrac{T}']=m/(p[0]*Gt+m)
        rows.append(out)
res=pd.DataFrame(rows)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(res.round(3).to_string())
res.to_csv('model_extra/pirt/pirt_fits.csv',index=False)
