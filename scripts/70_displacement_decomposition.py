"""Idea 1: which curve parameter carries the C. auris displacement?
Group-level posterior draws (68_isolate_params_dump.R). For every draw and every
C. auris clade vs each relative, recompute the 3x crossing (and the growth half-peak T50)
after swapping ONE parameter from the relative into the C. auris clade, holding the others
at C. auris values. The change in the crossing is that parameter's contribution to the
displacement (Shapley-style average over both directions to remove order dependence).
Parameters: E (growth activation energy), Eh (deactivation energy), Th (deactivation
midpoint), ER (respiration activation energy)."""
import os, itertools, math, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); CINV=11604.51812; TREF=293.15
pp=pd.read_csv(f"{C}/results/tables/isolate_params_draws.csv")
G=pp[pp.Isolate.str.startswith("GROUP_")].copy(); G["Group"]=G.Isolate.str.replace("GROUP_","")
def ln_cost(TC,E,ER,Eh,Th):
    TK=TC+273.15; u=Eh*CINV*(1/Th-1/TK); return (ER-E)*CINV*(1/TREF-1/TK)+np.where(u>30,u,np.log1p(np.exp(np.minimum(u,30))))
def topt(Es,Eh,Th):
    s=Es/Eh
    with np.errstate(invalid="ignore",divide="ignore"):
        return np.where((s>0)&(s<1),1/(1/Th-np.log(s/(1-s))/(Eh*CINV))-273.15,np.nan)
GRID=np.linspace(20,70,1001)
def cross3(E,ER,Eh,Th,M=3.0):
    """vectorised over draws: first T > Tmin where ln_cost - c0 >= ln M"""
    Tmin=topt(E-ER,Eh,Th); out=np.full(len(E),np.nan)
    for i in range(len(E)):
        if not np.isfinite(Tmin[i]): continue
        c0=ln_cost(Tmin[i],E[i],ER[i],Eh[i],Th[i]); lc=ln_cost(GRID,E[i],ER[i],Eh[i],Th[i])-c0-np.log(M)
        ok=(GRID>Tmin[i])&(lc>=0); out[i]=GRID[ok][0] if ok.any() else np.nan
    return out
PARS=["E","Eh","Th","ER"]
def params(g,n): x=G[G.Group==g].iloc[:n]; return {p:x[p].values for p in PARS}
N=600  # draws per comparison
rows=[]
for a in ["Clade1","Clade2","Clade3","Clade4"]:
    A=params(a,N)
    for r in ["para","Hae","Duo"]:
        R=params(r,N)
        full=cross3(A["E"],A["ER"],A["Eh"],A["Th"]); base=cross3(R["E"],R["ER"],R["Eh"],R["Th"])
        # Shapley values over the 4 parameters
        def val(S):  # crossing with parameters in S taken from C. auris, others from relative
            P={p:(A[p] if p in S else R[p]) for p in PARS}; return cross3(P["E"],P["ER"],P["Eh"],P["Th"])
        cache={}
        def v(S):
            k=tuple(sorted(S))
            if k not in cache: cache[k]=val(set(k))
            return cache[k]
        shap={p:np.zeros(N) for p in PARS}
        for p in PARS:
            others=[q for q in PARS if q!=p]
            for k in range(len(others)+1):
                for S in itertools.combinations(others,k):
                    w=math.factorial(k)*math.factorial(len(PARS)-k-1)/math.factorial(len(PARS))
                    shap[p]+=w*(v(set(S)|{p})-v(set(S)))
        tot=full-base
        row=dict(auris=a,relative=r,displacement_med=np.nanmedian(tot))
        for p in PARS: row[f"{p}_med"]=np.nanmedian(shap[p]); row[f"{p}_share"]=np.nanmedian(shap[p]/tot)
        rows.append(row)
        print(f"{a} vs {r}: displacement {np.nanmedian(tot):+.2f} C = "+" + ".join(f"{p} {np.nanmedian(shap[p]):+.2f}" for p in PARS),flush=True)
D=pd.DataFrame(rows); D.to_csv(f"{C}/results/tables/displacement_decomposition.csv",index=False)
print("\nmean share of the displacement (median over draws, averaged over the 12 pairs):")
print(D[[f"{p}_share" for p in PARS]].mean().round(2).to_string())
print("\nmean contribution in C:"); print(D[[f"{p}_med" for p in PARS]].mean().round(2).to_string())
