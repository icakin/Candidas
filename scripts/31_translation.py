"""Thermotolerance as a rigid translation of the growth-respiration relationship.
Per posterior draw, per isolate: T_opt (growth) and T_carbon (carbon-economy optimum,
= Sharpe-Schoolfield optimum with E replaced by E - E_R). Regress T_carbon on T_opt
across the 20 isolates within each draw; slope ~ 1 with a conserved gap means warm
adaptation slides the whole curve along the temperature axis rather than reshaping it.
Needs results/tables/{growth,resp}_draws_isolate.csv from 30_extract_draws.R.
Scale-free: both optima come from r and K only (no cell number, quota or RQ)."""
import pandas as pd, numpy as np
CINV=11604.51812
G=pd.read_csv("results/tables/growth_draws_isolate.csv")
R=pd.read_csv("results/tables/resp_draws_isolate.csv")
ISO=['Clade1_2068','Clade1_2069','Clade1_2070','Clade2_2071','Clade2_2072','Clade2_2073',
'Clade3_2074','Clade3_2075','Clade3_2076','Clade4_2077','Clade4_2078','Clade4_2079',
'Duo_1770','Duo_1771','Hae_1724','Hae_1768','Hae_1769','para_2051','para_2052','para_2053']
grp=lambda i:i.split('_')[0]
def topt(Es,Eh,ThK):
    s=Es/Eh
    with np.errstate(all='ignore'): t=1/(1/ThK - np.log(s/(1-s))/(Eh*CINV))-273.15
    return np.where((s>0)&(s<1)&np.isfinite(t),t,np.nan)
n=len(G); Topt=np.zeros((n,20)); Tcar=np.zeros((n,20))
for j,iso in enumerate(ISO):
    g=grp(iso)
    E=G[f"b_E_Group{g}"].values+G[f"r_Isolate__E[{iso},Intercept]"].values
    Th=G[f"b_Th_Group{g}"].values+G[f"r_Isolate__Th[{iso},Intercept]"].values
    Eh=G[f"b_Eh_Group{g}"].values
    ER=R[f"b_E_Group{g}"].values+R[f"r_Isolate__E[{iso},Intercept]"].values
    Topt[:,j]=topt(E,Eh,Th); Tcar[:,j]=topt(E-ER,Eh,Th)
slp=np.full(n,np.nan); dgap=np.full(n,np.nan)
for d in range(n):
    x,y=Topt[d],Tcar[d]; m=np.isfinite(x)&np.isfinite(y)
    if m.sum()<10: continue
    b,a=np.polyfit(x[m],y[m],1); slp[d]=b
    xlo,xhi=np.percentile(x[m],[10,90]); dgap[d]=(xlo-(a+b*xlo))-(xhi-(a+b*xhi))
q=lambda v:np.percentile(v[np.isfinite(v)],[2.5,50,97.5])
lo,md,hi=q(slp); print(f"slope T_carbon~T_opt: {md:.2f} [{lo:.2f},{hi:.2f}]  P(<1)={np.mean(slp[np.isfinite(slp)]<1):.2f}")
lo,md,hi=q(dgap); print(f"gap widening low-minus-high Topt: {md:.1f} [{lo:.1f},{hi:.1f}]  P(>0)={np.mean(dgap[np.isfinite(dgap)]>0):.2f}")
out=pd.DataFrame({"isolate":ISO,"group":[grp(i) for i in ISO],
    "T_opt":np.nanmedian(Topt,0),"T_carbon":np.nanmedian(Tcar,0)})
out["gap"]=out.T_opt-out.T_carbon
out.to_csv("results/tables/isolate_translation.csv",index=False)
print(f"T_opt spans {out.T_opt.max()-out.T_opt.min():.1f} C; gap median {out.gap.median():.1f} C, sd {out.gap.std():.2f}")
