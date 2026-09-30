#!/usr/bin/env python3
"""tmr_uncertainty.py -- confidence intervals on the fitted TmR of each species
and on the between-species differences (task 2).

Two independent routes, same background as model_v5 (z_auris.npy from
model_compare.py; recomputed if missing):

  A. Two-stage bootstrap.  Within each species resample isolates with
     replacement, then wells with replacement within each (isolate, T) cell,
     take the median per T (the v5 species-level definition) and refit
     TmR + Ea.  B replicates -> percentile CIs on TmR and on dTmR.
     For the single-isolate species (Hae_1768, Duo, para) the first stage is
     trivial and this is a plain well bootstrap.

  B. Profile likelihood.  TmR fixed on a grid, Ea re-optimised at each point.
     The 95% interval is where the scaled RSS stays below
     RSS_min * (1 + F_{0.95}(1, N-2) / (N-2)).  This asks a different
     question from A: how sharply does *this* dataset pin TmR given the model.

Usage:  python3 tmr_uncertainty.py [B=300] [n_jobs=all cores]
Output: tmr_uncertainty.csv, tmr_boot_draws.csv, tmr_uncertainty.png
"""
import sys, os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import f as fdist
from multiprocessing import get_context, cpu_count

HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model_v5_core import (SPS,COLS,SP_LABEL,load_species,measured,unfolded,predict,
                           fit_auris_background,sp_map,_raw)

B=int(sys.argv[1]) if len(sys.argv)>1 else 300
NJ=int(sys.argv[2]) if len(sys.argv)>2 else max(1,cpu_count()-1)
GROUPS=SPS+['hae_1724']

# ── data ──────────────────────────────────────────────────────────────────────
print("Loading ...",flush=True)
prot=load_species(); M=measured()
Ts=np.array(sorted(M['T'].unique()),float)
U={s:unfolded(*prot[s],Ts) for s in GROUPS}
W=_raw()                                   # well-level rows, column sp already mapped
W=W[W.sp.isin(GROUPS)].copy()

zfile=os.path.join(HERE,'z_auris.npy')
if os.path.exists(zfile): z_auris=np.load(zfile); print("  background: z_auris.npy")
else:
    a=M[M.sp=='auris'].sort_values('T')
    z_auris,_=fit_auris_background(U['auris'],Ts,a.mu.values,a.resp.values); np.save(zfile,z_auris)

def to_vectors(sub):
    """well rows of one species -> (ym, yr) on the common T grid, v5 convention"""
    g=sub.groupby('T')[['growth_C_per_C_h','respiration_C_per_C_h']].median()
    ym=np.zeros(len(Ts)); yr=np.full(len(Ts),np.nan)
    for i,T in enumerate(Ts):
        if T in g.index: ym[i]=g.loc[T,'growth_C_per_C_h']; yr[i]=g.loc[T,'respiration_C_per_C_h']
    return ym,yr

def loss_sp(x,u,ym,yr):
    z=z_auris.copy(); z[7]=x[0]; z[2]=np.log(max(x[1],1000.))
    mu,re=predict(z,u,Ts)
    if not(np.all(np.isfinite(mu)) and np.all(np.isfinite(re))): return 1e6
    ok=~np.isnan(yr)
    gl=np.mean((mu-ym)**2)/max(np.var(ym),1e-12)
    rl=np.mean((re[ok]-yr[ok])**2)/max(np.var(yr[ok]),1e-12) if ok.any() else 0
    return gl+rl

def fit_sp(u,ym,yr,starts):
    best=None
    for s in starts:
        r=minimize(loss_sp,s,args=(u,ym,yr),method='Nelder-Mead',
                   options=dict(maxiter=1500,fatol=1e-10,xatol=1e-6))
        if best is None or r.fun<best.fun: best=r
    return best.x[0],best.x[1],best.fun

# ── point estimates (full grid, as v5) ────────────────────────────────────────
print("Point estimates ...",flush=True)
FULL=[[t,e] for t in np.arange(30,50,2.) for e in (30e3,50e3,80e3,120e3,160e3)]
POINT={}
for sp in GROUPS:
    ym,yr=to_vectors(W[W.sp==sp]); POINT[sp]=fit_sp(U[sp],ym,yr,FULL)
    print(f"  {sp:18} TmR={POINT[sp][0]:6.2f}  Ea={POINT[sp][1]/1e3:5.0f} kJ/mol")

# ── A. two-stage bootstrap ────────────────────────────────────────────────────
def resample(sub,rng):
    isos=sub.otu_name.unique()
    pick=rng.choice(isos,size=len(isos),replace=True)
    parts=[]
    for k,iso in enumerate(pick):
        si=sub[sub.otu_name==iso]
        for T,cell in si.groupby('T'):
            idx=rng.integers(0,len(cell),len(cell))
            c=cell.iloc[idx].copy(); c['otu_name']=f'{iso}#{k}'; parts.append(c)
    return pd.concat(parts)

# shared Ea (model_compare M2: TmR per species, one Ea) for the fixed-Ea variants
EA_SHARED=83e3
mc=os.path.join(HERE,'model_compare.csv')
if os.path.exists(mc):
    t=pd.read_csv(mc); m2=t[t.model.str.startswith('M2')]
    if len(m2): EA_SHARED=float(m2.Ea_kJ_auris.iloc[0])*1e3
print(f"  shared Ea for fixed-Ea variants: {EA_SHARED/1e3:.0f} kJ/mol")

def fit_tmr_only(u,ym,yr,t0):
    best=None
    for s in (t0,t0-4,t0+4,t0-8,t0+8):
        r=minimize(lambda x: loss_sp([x[0],EA_SHARED],u,ym,yr),[s],method='Nelder-Mead',
                   options=dict(maxiter=600,fatol=1e-11,xatol=1e-5))
        if best is None or r.fun<best.fun: best=r
    return best.x[0]

POINT_FIX={sp:fit_tmr_only(U[sp],*to_vectors(W[W.sp==sp]),POINT[sp][0]) for sp in GROUPS}

def one_rep(seed):
    rng=np.random.default_rng(seed); out={}
    for sp in GROUPS:
        ym,yr=to_vectors(resample(W[W.sp==sp],rng))
        t0,e0,_=POINT[sp]
        starts=[[t0,e0],[t0-3,e0],[t0+3,e0],[t0,e0*0.6],[t0,e0*1.6],[t0-6,e0*1.3],[t0+6,e0*0.7]]
        tm,ea,_=fit_sp(U[sp],ym,yr,starts); out[f'TmR_{sp}']=tm; out[f'Ea_{sp}']=ea
        out[f'TmRfix_{sp}']=fit_tmr_only(U[sp],ym,yr,POINT_FIX[sp])
    return out

print(f"Bootstrap: B={B} on {NJ} cores ...",flush=True)
with get_context("fork").Pool(NJ) as pool: draws=pool.map(one_rep,range(B),chunksize=max(1,B//(NJ*8)))
D=pd.DataFrame(draws)
for a,b in [('auris','haemulonii'),('auris','duobushaemulonii'),('auris','parapsilosis'),('hae_1724','haemulonii')]:
    D[f'd_{a}-{b}']=D[f'TmR_{a}']-D[f'TmR_{b}']
    D[f'dfix_{a}-{b}']=D[f'TmRfix_{a}']-D[f'TmRfix_{b}']
D.to_csv(os.path.join(HERE,'tmr_boot_draws.csv'),index=False)

def ci(v): return np.percentile(v,[2.5,97.5])

# ── B. profile likelihood ────────────────────────────────────────────────────
print("Profile likelihood ...",flush=True)
grid=np.arange(28,50.01,0.25); PROF={}; PCI={}; PROFF={}; PCIF={}
for sp in GROUPS:
    ym,yr=to_vectors(W[W.sp==sp]); u=U[sp]
    N=len(ym)+int((~np.isnan(yr)).sum())
    prof=[]
    for t in grid:
        best=None
        for e0 in (30e3,60e3,100e3,150e3):
            r=minimize(lambda e: loss_sp([t,e[0]],u,ym,yr),[e0],method='Nelder-Mead',
                       options=dict(maxiter=800,fatol=1e-11,xatol=1e-3))
            if best is None or r.fun<best.fun: best=r
        prof.append(best.fun)
    prof=np.array(prof); PROF[sp]=prof
    thr=prof.min()*(1+fdist.ppf(0.95,1,N-2)/(N-2))
    inside=grid[prof<=thr]
    PCI[sp]=(inside.min(),inside.max(),grid[np.argmin(prof)])
    # same profile with Ea fixed at the shared value (1 free parameter)
    pf=np.array([loss_sp([t,EA_SHARED],u,ym,yr) for t in grid]); PROFF[sp]=pf
    thrf=pf.min()*(1+fdist.ppf(0.95,1,N-1)/(N-1)); ins=grid[pf<=thrf]
    PCIF[sp]=(ins.min(),ins.max(),grid[np.argmin(pf)])

# ── table ─────────────────────────────────────────────────────────────────────
rows=[]
for sp in GROUPS:
    v=D[f'TmR_{sp}']; lo,hi=ci(v); plo,phi,pmin=PCI[sp]
    rows.append(dict(quantity=f'TmR_{sp}',variant='Ea free',point=POINT[sp][0],boot_median=v.median(),boot_sd=v.std(),
                     boot_lo95=lo,boot_hi95=hi,profile_lo95=plo,profile_hi95=phi,profile_argmin=pmin))
    v=D[f'TmRfix_{sp}']; lo,hi=ci(v); plo,phi,pmin=PCIF[sp]
    rows.append(dict(quantity=f'TmR_{sp}',variant=f'Ea fixed {EA_SHARED/1e3:.0f}',point=POINT_FIX[sp],boot_median=v.median(),boot_sd=v.std(),
                     boot_lo95=lo,boot_hi95=hi,profile_lo95=plo,profile_hi95=phi,profile_argmin=pmin))
for a,b in [('auris','haemulonii'),('auris','duobushaemulonii'),('auris','parapsilosis'),('hae_1724','haemulonii')]:
    v=D[f'd_{a}-{b}']; lo,hi=ci(v)
    rows.append(dict(quantity=f'dTmR_{a}-{b}',variant='Ea free',point=POINT[a][0]-POINT[b][0],boot_median=v.median(),
                     boot_sd=v.std(),boot_lo95=lo,boot_hi95=hi,P_le_0=float((v<=0).mean())))
    v=D[f'dfix_{a}-{b}']; lo,hi=ci(v)
    rows.append(dict(quantity=f'dTmR_{a}-{b}',variant=f'Ea fixed {EA_SHARED/1e3:.0f}',point=POINT_FIX[a]-POINT_FIX[b],boot_median=v.median(),
                     boot_sd=v.std(),boot_lo95=lo,boot_hi95=hi,P_le_0=float((v<=0).mean())))
T=pd.DataFrame(rows); T.to_csv(os.path.join(HERE,'tmr_uncertainty.csv'),index=False)

print("\n══ TmR UNCERTAINTY (background fixed, v5) ══")
print(f"{'quantity':30} {'variant':>12} {'point':>7} {'boot 95% CI':>18} {'profile 95% CI':>18}")
for _,r in T.iterrows():
    prof=f"[{r.profile_lo95:5.1f}, {r.profile_hi95:5.1f}]" if 'profile_lo95' in r and pd.notna(r.get('profile_lo95')) else ''
    extra=f"   P(Δ≤0)={r.P_le_0:.3f}" if 'P_le_0' in r and pd.notna(r.get('P_le_0')) else ''
    print(f"  {r.quantity:28} {r.variant:>12} {r.point:7.2f}   [{r.boot_lo95:5.1f}, {r.boot_hi95:5.1f}]   {prof:>18}{extra}")

# ── figure ────────────────────────────────────────────────────────────────────
fig,AX=plt.subplots(2,3,figsize=(16,9))
for row,(PR,PC,tk,dk,lab) in enumerate([(PROF,PCI,'TmR_','d_','Ea free (v5)'),
                                        (PROFF,PCIF,'TmRfix_','dfix_',f'Ea fixed at {EA_SHARED/1e3:.0f} kJ/mol (model-comparison M2)')]):
    ax=AX[row]
    for sp in GROUPS:
        p=PR[sp]; ax[0].plot(grid,p/p.min(),color=COLS[sp],lw=2,label=SP_LABEL[sp].split('\n')[0])
        ax[0].axvspan(PC[sp][0],PC[sp][1],color=COLS[sp],alpha=.08)
    ax[0].axhline(1,color='k',lw=.5); ax[0].set_ylim(0.95,2.5)
    ax[0].set_xlabel('TmR (°C)'); ax[0].set_ylabel('loss / min loss'); ax[0].set_title(f'Profile likelihood on TmR, {lab}',fontsize=9.5)
    if row==0: ax[0].legend(fontsize=7.5)
    for sp in GROUPS:
        ax[1].hist(D[f'{tk}{sp}'],bins=40,color=COLS[sp],alpha=.55)
        ax[1].axvline((POINT[sp][0] if row==0 else POINT_FIX[sp]),color=COLS[sp],lw=1.2,ls='--')
    ax[1].set_xlabel('bootstrap TmR (°C)'); ax[1].set_title(f'Two-stage bootstrap, B={B}, {lab}',fontsize=9.5)
    for key,c in [('auris-haemulonii','#377eb8'),('auris-duobushaemulonii','#4daf4a'),('auris-parapsilosis','#984ea3')]:
        v=D[f'{dk}{key}']; lo,hi=ci(v)
        ax[2].hist(v,bins=40,color=c,alpha=.55,label=f"{key.replace('-',' − ')}: {v.median():+.1f} [{lo:+.1f}, {hi:+.1f}]")
    ax[2].axvline(0,color='k',lw=1)
    ax[2].set_xlabel('ΔTmR (°C)'); ax[2].set_title(f'Between-species differences, {lab}',fontsize=9.5); ax[2].legend(fontsize=7.5)
fig.suptitle('How well is TmR determined by the growth + respiration data?  (v5 background fixed)',fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(HERE,'tmr_uncertainty.png'),dpi=150,bbox_inches='tight')
print("\nSaved: tmr_uncertainty.csv, tmr_boot_draws.csv, tmr_uncertainty.png")
