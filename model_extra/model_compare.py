#!/usr/bin/env python3
"""model_compare.py -- does the data require a species-specific TmR?

Nested comparison on the same four species and same background as model_v5
(7 background parameters fixed from the auris fit; only TmR and Ea vary):

  M0  null       : 1 TmR + 1 Ea shared by all species           (2 free)
  M1  shared-TmR : 1 TmR shared, Ea free per species             (5 free)
  M2  shared-Ea  : TmR free per species, 1 Ea shared             (5 free)
  M3  full (v5)  : TmR + Ea free per species                     (8 free)

Fit criterion is the v5 loss (normalised MSE of growth + respiration).  For
model selection the residuals are put on one scale (divided by the per-species
SD of each channel) and treated as Gaussian with a common variance, which
gives the least-squares form of AIC:  AIC = N ln(RSS/N) + 2k, plus the
small-sample correction (AICc) and nested F-tests.

Output: model_compare.csv, model_compare.png, z_auris.npy (cached background).
"""
import sys, os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import f as fdist

HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,HERE)
from model_v5_core import (SPS,COLS,SP_LABEL,load_species,measured,unfolded,
                           predict,fit_auris_background,compute_fivywrel,ZELDOVICH_SLOPE)

# ── data ──────────────────────────────────────────────────────────────────────
print("Loading ...",flush=True)
prot=load_species(); M=measured()
Ts=np.array(sorted(M['T'].unique()),float)
U={s:unfolded(*prot[s],Ts) for s in SPS}
T_fine=np.linspace(Ts.min()-2,Ts.max()+2,200)
U_fine={s:unfolded(*prot[s],T_fine) for s in SPS}

# per-species observed vectors and normalisation constants
OBS={}
for sp in SPS:
    d=M[M.sp==sp].sort_values('T')
    ym=d.mu.values; yr=d.resp.values; ok=~np.isnan(yr)
    OBS[sp]=dict(ym=ym,yr=yr,ok=ok,sg=np.sqrt(max(np.var(ym),1e-12)),
                 sr=np.sqrt(max(np.var(yr[ok]),1e-12)) if ok.any() else 1.0)
N_total=sum(len(o['ym'])+int(o['ok'].sum()) for o in OBS.values())
print(f"  data points (growth + respiration, 4 species): N = {N_total}")

# ── background (cached) ───────────────────────────────────────────────────────
zfile=os.path.join(HERE,'z_auris.npy')
if os.path.exists(zfile):
    z_auris=np.load(zfile); print("  background: loaded z_auris.npy")
else:
    print("Fitting auris background (v5 procedure) ...",flush=True)
    z_auris,_=fit_auris_background(U['auris'],Ts,OBS['auris']['ym'],OBS['auris']['yr'])
    np.save(zfile,z_auris); print("  background: saved z_auris.npy")

def species_resid(sp,TmR,Ea):
    """scaled residual vector (growth then respiration) for one species"""
    z=z_auris.copy(); z[7]=TmR; z[2]=np.log(max(Ea,1000.))
    mu,re=predict(z,U[sp],Ts)
    if not(np.all(np.isfinite(mu)) and np.all(np.isfinite(re))): return None
    o=OBS[sp]
    return np.concatenate([(mu-o['ym'])/o['sg'],(re[o['ok']]-o['yr'][o['ok']])/o['sr']])

def rss(TmR_d,Ea_d):
    tot=0.
    for sp in SPS:
        r=species_resid(sp,TmR_d[sp],Ea_d[sp])
        if r is None: return 1e9
        tot+=np.sum(r**2)
    return tot

# ── IVYWREL fractions (sequence input for the v6-type models) ─────────────────
FIV=compute_fivywrel(); F=np.array([FIV[s] for s in SPS])
Fc=F-F[SPS.index('auris')]              # centred on auris: c0 is then auris' TmR
print("  F(IVYWREL): "+"  ".join(f"{s[:4]}={FIV[s]:.5f}" for s in SPS))

# ── model definitions: parameter vector -> (TmR dict, Ea dict) ────────────────
def m_null(x):   return {s:x[0] for s in SPS},{s:x[1] for s in SPS}
def m_shTmR(x):  return {s:x[0] for s in SPS},{s:x[1+i] for i,s in enumerate(SPS)}
def m_shEa(x):   return {s:x[1+i] for i,s in enumerate(SPS)},{s:x[0] for s in SPS}
def m_full(x):   return {s:x[i] for i,s in enumerate(SPS)},{s:x[4+i] for i,s in enumerate(SPS)}
# v6 family: TmR(sp) = c0 + slope * (F_sp - F_auris)   (c0 = auris TmR; slope fixed 937 or free)
def m_v6(x):     return {s:x[0]+ZELDOVICH_SLOPE*Fc[i] for i,s in enumerate(SPS)},{s:x[1+i] for i,s in enumerate(SPS)}
def m_v6sh(x):   return {s:x[0]+ZELDOVICH_SLOPE*Fc[i] for i,s in enumerate(SPS)},{s:x[1] for s in SPS}
def m_v6free(x): return {s:x[0]+x[2]*Fc[i] for i,s in enumerate(SPS)},{s:x[1] for s in SPS}
MODELS={'M0 null (1 TmR, 1 Ea)':                 (m_null,  2),
        'M1 shared TmR, Ea per species':         (m_shTmR, 5),
        'M2 TmR per species, shared Ea':         (m_shEa,  5),
        'M3 full v5 (TmR+Ea per species)':       (m_full,  8),
        'M4 v6: TmR=937F+c0, Ea per species':    (m_v6,    5),
        'M5 v6: TmR=937F+c0, shared Ea':         (m_v6sh,  2),
        'M6 TmR=b*F+c0 (slope free), shared Ea': (m_v6free,3)}

def fit(name):
    unpackf,k=MODELS[name]
    def obj(x): return rss(*unpackf(x))
    starts=[]
    TmR0s=np.arange(30,50,2.); Ea0s=[30e3,50e3,80e3,120e3,160e3]
    rng=np.random.default_rng(0)
    for t0 in TmR0s:
        for e0 in Ea0s:
            if name.startswith(('M0','M5')): starts.append([t0,e0])
            elif name.startswith('M1'): starts.append([t0]+[e0*rng.uniform(.7,1.4) for _ in SPS])
            elif name.startswith('M2'): starts.append([e0]+[t0+rng.uniform(-3,3) for _ in SPS])
            elif name.startswith('M4'): starts.append([t0]+[e0*rng.uniform(.7,1.4) for _ in SPS])
            elif name.startswith('M6'): starts.append([t0,e0,rng.choice([0.,300.,937.,2000.])])
            else:      starts.append([t0+rng.uniform(-3,3) for _ in SPS]+[e0*rng.uniform(.7,1.4) for _ in SPS])
    best=None
    for s in starts:
        r=minimize(obj,s,method='Nelder-Mead',options=dict(maxiter=6000,maxfev=6000,fatol=1e-10,xatol=1e-8))
        if best is None or r.fun<best.fun: best=r
    # polish
    r=minimize(obj,best.x,method='Nelder-Mead',options=dict(maxiter=20000,maxfev=20000,fatol=1e-12,xatol=1e-10))
    if r.fun<best.fun: best=r
    return best.x,best.fun,k,unpackf

# ── fit all ───────────────────────────────────────────────────────────────────
res={}
for name in MODELS:
    print(f"Fitting {name} ...",flush=True)
    x,RSS,k,unpackf=fit(name)
    TmR_d,Ea_d=unpackf(x)
    N=N_total
    AIC=N*np.log(RSS/N)+2*k
    AICc=AIC+2*k*(k+1)/(N-k-1)
    BIC=N*np.log(RSS/N)+k*np.log(N)
    res[name]=dict(x=x,RSS=RSS,k=k,AIC=AIC,AICc=AICc,BIC=BIC,TmR=TmR_d,Ea=Ea_d)
    print(f"  RSS={RSS:.3f}  k={k}  AICc={AICc:.2f}")
    print("  TmR: "+"  ".join(f"{s[:4]}={TmR_d[s]:.2f}" for s in SPS))
    print("  Ea : "+"  ".join(f"{s[:4]}={Ea_d[s]/1e3:.0f}" for s in SPS))
    if name.startswith('M6'): print(f"  fitted IVYWREL slope b = {x[2]:.0f} °C per unit   (Zeldovich 937)")

# ── table ─────────────────────────────────────────────────────────────────────
best_aicc=min(r['AICc'] for r in res.values())
rows=[]
for name,r in res.items():
    rows.append(dict(model=name,k=r['k'],RSS=r['RSS'],AIC=r['AIC'],AICc=r['AICc'],
                     dAICc=r['AICc']-best_aicc,BIC=r['BIC'],
                     **{f'TmR_{s}':r['TmR'][s] for s in SPS},
                     **{f'Ea_kJ_{s}':r['Ea'][s]/1e3 for s in SPS}))
tab=pd.DataFrame(rows)
w=np.exp(-0.5*tab.dAICc); tab['Akaike_weight']=w/w.sum()
tab.to_csv(os.path.join(HERE,'model_compare.csv'),index=False)

def ftest(red,full):
    a,b=res[red],res[full]; N=N_total
    df1=b['k']-a['k']; df2=N-b['k']
    F=((a['RSS']-b['RSS'])/df1)/(b['RSS']/df2)
    return F,df1,df2,fdist.sf(F,df1,df2)

names=list(MODELS)
print("\n══ MODEL COMPARISON (same background, same 4 species) ══")
print(f"{'model':36} {'k':>2} {'RSS':>8} {'AICc':>8} {'ΔAICc':>7} {'w':>6}")
for _,t in tab.iterrows():
    print(f"  {t.model:34} {t.k:>2} {t.RSS:>8.3f} {t.AICc:>8.2f} {t.dAICc:>7.2f} {t.Akaike_weight:>6.3f}")
print("\nNested F-tests (reduced vs fuller):")
for red,full in [(names[0],names[1]),(names[0],names[2]),(names[1],names[3]),(names[2],names[3]),
                 (names[5],names[6]),(names[6],names[2]),(names[5],names[2]),(names[4],names[3])]:
    F_,d1,d2,p=ftest(red,full)
    print(f"  {red[:2]} -> {full[:2]}: F({d1},{d2}) = {F_:6.2f}   p = {p:.2e}")

r3=res[names[3]]; r2=res[names[2]]
print(f"\n  ΔTmR(auris-hae)  full: {r3['TmR']['auris']-r3['TmR']['haemulonii']:+.2f}°C"
      f"   shared-Ea: {r2['TmR']['auris']-r2['TmR']['haemulonii']:+.2f}°C")
print("\nReading: if M1 (one TmR for everyone) is far worse than M3 while M2 (one Ea)\n"
      "is close to M3, the species difference is carried by TmR, not by Ea.")

# ── joint refit: background free under every model (fairness check) ──────────
# The fixed background was co-optimised with auris at TmR~39, which could
# handicap models that force auris lower.  Here the 7 background parameters
# are re-optimised together with each model's own parameters.
BG_IDX=[0,1,3,4,5,6,8]          # kappa,g0,m0,cg,cd,PHI,Nr  (z[2]=Ea, z[7]=TmR are model params)
def rss_joint(zbg,TmR_d,Ea_d):
    tot=0.
    for sp in SPS:
        z=z_auris.copy(); z[BG_IDX]=zbg; z[7]=TmR_d[sp]; z[2]=np.log(max(Ea_d[sp],1000.))
        mu,re=predict(z,U[sp],Ts)
        if not(np.all(np.isfinite(mu)) and np.all(np.isfinite(re))): return 1e9
        o=OBS[sp]
        tot+=np.sum(((mu-o['ym'])/o['sg'])**2)+np.sum(((re[o['ok']]-o['yr'][o['ok']])/o['sr'])**2)
    return tot

print("\nJoint refit (background free under each model) ...",flush=True)
resJ={}; rng=np.random.default_rng(7)
for name,(unpackf,k) in MODELS.items():
    x0=res[name]['x']; zb0=z_auris[BG_IDX]
    def obj(v): return rss_joint(v[:7],*unpackf(v[7:]))
    best=None
    for i in range(6):
        s=np.concatenate([zb0+(rng.normal(0,.3,7) if i else 0), x0+(rng.normal(0,1,len(x0))*np.where(np.abs(x0)>1e3,0.15*np.abs(x0),1.5) if i else 0)])
        r=minimize(obj,s,method='Nelder-Mead',options=dict(maxiter=8000,maxfev=8000,fatol=1e-10,xatol=1e-8,adaptive=True))
        if best is None or r.fun<best.fun: best=r
    TmR_d,Ea_d=unpackf(best.x[7:]); kk=k+7; N=N_total
    AICc=N*np.log(best.fun/N)+2*kk+2*kk*(kk+1)/(N-kk-1)
    resJ[name]=dict(RSS=best.fun,k=kk,AICc=AICc,TmR=TmR_d,Ea=Ea_d,zbg=best.x[:7])
    print(f"  {name:34} RSS={best.fun:.3f} k={kk} AICc={AICc:.2f}  "
          +" ".join(f"TmR_{s[:4]}={TmR_d[s]:.1f}" for s in SPS)+"  "
          +" ".join(f"Ea_{s[:4]}={Ea_d[s]/1e3:.0f}" for s in SPS))
bJ=min(r['AICc'] for r in resJ.values())
print("\n══ JOINT (background free) ══")
print(f"{'model':36} {'k':>2} {'RSS':>8} {'AICc':>8} {'ΔAICc':>7}")
for name,r in resJ.items():
    print(f"  {name:34} {r['k']:>2} {r['RSS']:>8.3f} {r['AICc']:>8.2f} {r['AICc']-bJ:>7.2f}")
def ftestJ(red,full):
    a,b=resJ[red],resJ[full]; d1=b['k']-a['k']; d2=N_total-b['k']
    F=((a['RSS']-b['RSS'])/d1)/(b['RSS']/d2); return F,d1,d2,fdist.sf(F,d1,d2)
for red,full in [(names[1],names[3]),(names[2],names[3]),(names[0],names[1]),(names[0],names[2]),(names[5],names[2]),(names[4],names[3])]:
    F_,d1,d2,p=ftestJ(red,full); print(f"  {red[:2]} -> {full[:2]}: F({d1},{d2}) = {F_:6.2f}   p = {p:.2e}")
pd.DataFrame([dict(model=n,fit='joint',k=r['k'],RSS=r['RSS'],AICc=r['AICc'],dAICc=r['AICc']-bJ,
                   **{f'TmR_{s}':r['TmR'][s] for s in SPS},**{f'Ea_kJ_{s}':r['Ea'][s]/1e3 for s in SPS})
              for n,r in resJ.items()]).to_csv(os.path.join(HERE,'model_compare_joint.csv'),index=False)

# ── plot ─────────────────────────────────────────────────────────────────────
STY={names[0]:dict(ls=':',lw=1.2,alpha=.7),names[1]:dict(ls='--',lw=1.5,alpha=.9),
     names[2]:dict(ls='-.',lw=1.5,alpha=.9),names[3]:dict(ls='-',lw=2.2,alpha=1),
     names[5]:dict(ls=(0,(1,1)),lw=2.0,alpha=1,color='k')}
PLOT_NAMES=[names[i] for i in (0,1,2,3,5)]
fig,axes=plt.subplots(2,4,figsize=(18,7.5),sharey='row')
for c,sp in enumerate(SPS):
    d=M[M.sp==sp]
    for row,(key,yl) in enumerate([('mu','Growth (C/C/h)'),('resp','Respiration (C/C/h)')]):
        ax=axes[row,c]
        pts=d[d.alive] if row==0 else d.dropna(subset=['resp'])
        ax.scatter(pts['T'],pts[key],color=COLS[sp],s=38,zorder=6,label='data')
        for name in PLOT_NAMES:
            r=res[name]; z=z_auris.copy(); z[7]=r['TmR'][sp]; z[2]=np.log(max(r['Ea'][sp],1000.))
            mu,re=predict(z,U_fine[sp],T_fine); y=mu if row==0 else re
            st=dict(STY[name]); c=st.pop('color',COLS[sp])
            ax.plot(T_fine,y,color=c,label=name.split(' (')[0].split(':')[0],**st)
        ax.axhline(0,color='k',lw=.4,ls=':')
        if row==0:
            ax.set_ylim(-0.01,None)
            ax.set_title(f"{SP_LABEL[sp].split(chr(10))[0]}\n"
                         f"TmR full={r3['TmR'][sp]:.1f}  sharedEa={r2['TmR'][sp]:.1f}  "
                         f"sharedTmR={res[names[1]]['TmR'][sp]:.1f}",fontsize=8.5)
        if c==0: ax.set_ylabel(yl)
        if row==1: ax.set_xlabel('Temperature (°C)')
axes[0,0].legend(fontsize=7,loc='upper left')
sub="  |  ".join(f"{n.split(' ')[0]}: k={res[n]['k']}, ΔAICc={res[n]['AICc']-best_aicc:.1f}" for n in names)
sub=sub.replace('M4','\nM4')
fig.suptitle("Does the data require a species-specific TmR?  Nested models on the v5 background\n"+sub,fontsize=9.5)
plt.tight_layout(); plt.savefig(os.path.join(HERE,'model_compare.png'),dpi=150,bbox_inches='tight')
print("\nSaved: model_compare.csv, model_compare.png")
