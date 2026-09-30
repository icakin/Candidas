#!/usr/bin/env python3
"""isolate_tmr.py -- TmR and Ea fitted for every isolate separately (task 3).

Same background as model_v5 (z_auris.npy), same per-species proteome u(T)
(all auris clades share the auris etcGEM; all Hae isolates share the
haemulonii one, etc.).  Only TmR and Ea are free per isolate.  Nothing is
excluded here: Hae_1724 and Hae_1769 are fitted like everyone else so the
reader can see where they fall.

Questions this answers
  * within-species spread of TmR versus the between-species difference
  * whether the 12 auris isolates (4 clades x 3) cluster above the relatives
  * where Hae_1724 and Hae_1769 sit, with the fit quality (loss) alongside,
    so 'excluded because it is a poor TPC' is visible rather than asserted
  * one-way ANOVA / Kruskal-Wallis on isolate TmR by species, and a
    permutation test of auris vs relatives on isolate TmR (the 40 degC
    Fig. 3 contrast, but on the model parameter)

Output: isolate_tmr.csv, isolate_tmr.png
"""
import sys, os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import kruskal, f_oneway, mannwhitneyu

HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model_v5_core import (SPS,COLS,load_species,measured,measured_isolates,unfolded,
                           predict,fit_auris_background)

# ── data ──────────────────────────────────────────────────────────────────────
print("Loading ...",flush=True)
prot=load_species(); M=measured(); I=measured_isolates()
Ts=np.array(sorted(M['T'].unique()),float)
U={s:unfolded(*prot[s],Ts) for s in SPS}
zfile=os.path.join(HERE,'z_auris.npy')
if os.path.exists(zfile): z_auris=np.load(zfile); print("  background: z_auris.npy")
else:
    a=M[M.sp=='auris'].sort_values('T')
    z_auris,_=fit_auris_background(U['auris'],Ts,a.mu.values,a.resp.values); np.save(zfile,z_auris)

def loss_sp(x,u,ym,yr):
    z=z_auris.copy(); z[7]=x[0]; z[2]=np.log(max(x[1],1000.))
    mu,re=predict(z,u,Ts)
    if not(np.all(np.isfinite(mu)) and np.all(np.isfinite(re))): return 1e6
    ok=~np.isnan(yr)
    gl=np.mean((mu-ym)**2)/max(np.var(ym),1e-12)
    rl=np.mean((re[ok]-yr[ok])**2)/max(np.var(yr[ok]),1e-12) if ok.any() else 0
    return gl+rl

def fit(u,ym,yr):
    best=None
    for t0 in np.arange(30,50,2.):
        for e0 in (30e3,50e3,80e3,120e3,160e3):
            r=minimize(loss_sp,[t0,e0],args=(u,ym,yr),method='Nelder-Mead',
                       options=dict(maxiter=1500,fatol=1e-10,xatol=1e-6))
            if best is None or r.fun<best.fun: best=r
    return best.x[0],best.x[1],best.fun

# ── shared Ea from the model comparison (M2 was preferred: TmR per species, one Ea) ──
EA_SHARED=83e3
mc=os.path.join(HERE,'model_compare.csv')
if os.path.exists(mc):
    t=pd.read_csv(mc); m2=t[t.model.str.startswith('M2')]
    if len(m2): EA_SHARED=float(m2.Ea_kJ_auris.iloc[0])*1e3
print(f"  shared Ea for the TmR-only fits: {EA_SHARED/1e3:.0f} kJ/mol (from model_compare M2)")

def fit_tmr_only(u,ym,yr):
    best=None
    for t0 in np.arange(28,50,1.):
        r=minimize(lambda x: loss_sp([x[0],EA_SHARED],u,ym,yr),[t0],method='Nelder-Mead',
                   options=dict(maxiter=600,fatol=1e-11,xatol=1e-5))
        if best is None or r.fun<best.fun: best=r
    return best.x[0],best.fun

# ── fit each isolate ──────────────────────────────────────────────────────────
print("Fitting isolates ...",flush=True)
rows=[]
for iso,d in I.groupby('iso'):
    d=d.sort_values('T'); sp=d.sp.iloc[0]
    ym=d.mu.values; yr=d.resp.values
    Tmax=d[d.alive]['T'].max() if d.alive.any() else np.nan
    tm1,lo1=fit_tmr_only(U[sp],ym,yr)           # primary: TmR only, Ea shared
    tm,ea,lo=fit(U[sp],ym,yr)                   # secondary: TmR + Ea free
    clade=iso.split('_')[0] if sp=='auris' else sp
    rows.append(dict(iso=iso,sp=sp,clade=clade,TmR=tm1,loss=lo1,
                     TmR_freeEa=tm,Ea_kJ_freeEa=ea/1e3,loss_freeEa=lo,
                     Tmax_obs=Tmax,n_T=int(d.alive.sum()),
                     excluded_in_model=iso in ('Hae_1724','Hae_1769')))
    print(f"  {iso:14} {sp:17} TmR={tm1:6.2f} (loss {lo1:.3f})   | free-Ea: TmR={tm:6.2f} Ea={ea/1e3:5.0f} (loss {lo:.3f})   Tmax={Tmax}")
R=pd.DataFrame(rows).sort_values(['sp','TmR']); R.to_csv(os.path.join(HERE,'isolate_tmr.csv'),index=False)

# ── statistics ────────────────────────────────────────────────────────────────
print("\n══ ISOLATE-LEVEL TmR ══")
for sp in SPS:
    v=R[R.sp==sp].TmR
    print(f"  {sp:17} n={len(v):2d}  mean={v.mean():6.2f}  sd={v.std(ddof=1) if len(v)>1 else float('nan'):5.2f}"
          f"  range=[{v.min():.1f}, {v.max():.1f}]")
au=R[R.sp=='auris'].TmR.values
rel_all=R[R.sp!='auris'].TmR.values
rel_clean=R[(R.sp!='auris')&(~R.excluded_in_model)].TmR.values
def perm_p(a,b,n=20000,seed=1):
    rng=np.random.default_rng(seed); obs=a.mean()-b.mean(); pool=np.concatenate([a,b]); k=len(a); c=0
    for _ in range(n):
        rng.shuffle(pool); c+= (pool[:k].mean()-pool[k:].mean())>=obs
    return obs,(c+1)/(n+1)
for lab,rel in [('all 8 relatives (incl. Hae_1724, Hae_1769)',rel_all),('6 relatives used in the model',rel_clean)]:
    obs,p=perm_p(au,rel); _,pmw=mannwhitneyu(au,rel,alternative='greater')
    print(f"\n  auris (n={len(au)}) vs {lab} (n={len(rel)}):")
    print(f"    mean ΔTmR = {obs:+.2f}°C   permutation p = {p:.4f}   Mann-Whitney p = {pmw:.4f}")
    print(f"    auris min {au.min():.1f}  vs  relatives max {rel.max():.1f}  ->  "
          +("no overlap" if au.min()>rel.max() else "OVERLAP"))
groups=[R[R.sp==s].TmR.values for s in SPS if (R.sp==s).sum()>1]
print(f"\n  Kruskal-Wallis across species: H={kruskal(*groups).statistic:.2f}, p={kruskal(*groups).pvalue:.2e}")
print(f"  one-way ANOVA:                 F={f_oneway(*groups).statistic:.2f}, p={f_oneway(*groups).pvalue:.2e}")
sw=R[R.sp=='auris'].groupby('clade').TmR.agg(['mean','std','count'])
print("\n  auris by clade:\n"+sw.to_string())
print(f"\n  within-auris sd = {au.std(ddof=1):.2f}°C   between-species (auris − hae_1768) = "
      f"{au.mean()-R[R.iso=='Hae_1768'].TmR.iloc[0]:+.2f}°C")

# ── figure ────────────────────────────────────────────────────────────────────
fig,ax=plt.subplots(1,2,figsize=(13,5),gridspec_kw=dict(width_ratios=[1.3,1]))
order=SPS; xs={s:i for i,s in enumerate(order)}
rng=np.random.default_rng(0)
for _,r in R.iterrows():
    x=xs[r.sp]+rng.uniform(-.18,.18)
    mk='x' if r.excluded_in_model else 'o'
    ax[0].scatter(x,r.TmR,color=COLS[r.sp],s=70 if not r.excluded_in_model else 90,marker=mk,
                  zorder=5,edgecolor='k' if not r.excluded_in_model else COLS[r.sp],lw=.6)
    if r.sp!='auris' or r.excluded_in_model:
        ax[0].annotate(r.iso,(x,r.TmR),xytext=(6,0),textcoords='offset points',fontsize=7,va='center')
for s in order:
    v=R[(R.sp==s)&(~R.excluded_in_model)].TmR
    ax[0].hlines(v.mean(),xs[s]-.3,xs[s]+.3,color=COLS[s],lw=2.5)
ax[0].set_xticks(range(len(order))); ax[0].set_xticklabels([f"C. {s}" for s in order],fontsize=8.5,style='italic')
ax[0].set_ylabel('fitted TmR (°C), Ea shared'); ax[0].set_title(f'TmR per isolate, Ea fixed at {EA_SHARED/1e3:.0f} kJ/mol (bar = species mean of isolates used in the model;\n× = excluded from the species-level model)',fontsize=9)
sc=ax[1].scatter(R.Tmax_obs+rng.uniform(-.3,.3,len(R)),R.TmR,c=[COLS[s] for s in R.sp],s=60,edgecolor='k',lw=.5,zorder=5)
for _,r in R[R.excluded_in_model].iterrows(): ax[1].annotate(r.iso,(r.Tmax_obs,r.TmR),xytext=(6,-3),textcoords='offset points',fontsize=7)
ax[1].set_xlabel('highest temperature with growth (°C, observed)'); ax[1].set_ylabel('fitted TmR (°C)')
ax[1].set_title('Fitted TmR against the raw phenotype it summarises',fontsize=9)
fig.suptitle('Isolate-level fits: is the auris shift larger than the within-species spread?',fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(HERE,'isolate_tmr.png'),dpi=150,bbox_inches='tight')
print("\nSaved: isolate_tmr.csv, isolate_tmr.png")
