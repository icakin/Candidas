#!/usr/bin/env python3
"""model_v5_core.py -- the model_v5 machinery as an importable module.

Everything here is copied from model_v5.py unchanged (two-state unfolding,
u(T) from the etcGEM enzyme list, the allocation model, the auris background
fit).  model_compare.py, tmr_uncertainty.py, isolate_tmr.py and
ribo_seq2tm.py import from here so all analyses share one definition.

Species-level data: median of growth / respiration across isolates at each T.
  haemulonii  = Hae_1768 only; Hae_1724 kept as its own group 'hae_1724';
  Hae_1769 excluded (as in model_v5).
"""
import sys, os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from scipy.optimize import minimize

def _find_repo():
    """CANDIDAS_ROOT env var, else walk up from this file until gem/gempaths.py is found"""
    if os.environ.get('CANDIDAS_ROOT'): return os.environ['CANDIDAS_ROOT']
    p=os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        if os.path.exists(os.path.join(p,'gem','gempaths.py')): return p
        p=os.path.dirname(p)
    return '/home/claude/cand'
REPO=_find_repo(); os.chdir(REPO); sys.path.insert(0,os.path.join(REPO,'gem'))
import gempaths as G
print(f"  repo: {REPO}")
R_gas=8.314462618; TH,TS=373.15,385.0
SPS=['auris','haemulonii','duobushaemulonii','parapsilosis']
COLS={'auris':'#e41a1c','haemulonii':'#377eb8','duobushaemulonii':'#4daf4a',
      'parapsilosis':'#984ea3','hae_1724':'#ff7f00'}
SP_LABEL={'auris':'C. auris','haemulonii':'C. haemulonii\n(Hae_1768 only)',
          'duobushaemulonii':'C. duobushaemulonii','parapsilosis':'C. parapsilosis',
          'hae_1724':'Hae_1724\n(thermotolerant outlier)'}

def two_state(N,TmC,T_C):
    N=np.atleast_1d(np.asarray(N,float)); TmK=np.atleast_1d(np.asarray(TmC,float))+273.15
    dH=(4.*N+143.)*1000.; dS=13.27*N+448.
    den=(TmK-TH)-TmK*np.log(TmK/TS)
    dCp=np.where(np.abs(den)>1e-9,(TmK*dS-dH)/den,0.)
    out=[]
    for T in np.atleast_1d(np.asarray(T_C,float))+273.15:
        dG=dH+dCp*(T-TH)-T*dS-T*dCp*np.log(T/TS)
        out.append(1./(1.+np.exp(-np.clip(dG/(R_gas*T),-50,50))))
    return np.array(out)

def unfolded(L,tm,Ts): return 1.-two_state(L,tm,Ts).mean(axis=1)

def load_species(tm_sd=5.9):
    tmd=pd.read_csv(G.TABLES/'thermal_tm.csv')
    tmd['sp']=tmd.id.str.split('|').str[0]; tmd['g']=tmd.id.str.split('|').str[1]
    tmd['use']=(51.9+(tmd.pred_tm-tmd.pred_tm.mean())*(tm_sd/tmd.pred_tm.std(ddof=1))
                if tm_sd>0 else tmd.pred_tm)
    out={}
    for sp in SPS:
        mwt=pd.read_csv(G.TABLES/f'enzyme_mw_{sp}.csv'); L=dict(zip(mwt.gene,mwt.length_aa))
        rows=[(L[r.g],r.use) for r in tmd[tmd.sp==sp].itertuples() if r.g in L]
        out[sp]=(np.array([a for a,_ in rows]),np.array([b for _,b in rows]))
    out['hae_1724']=out['haemulonii']
    return out

def sp_map(o):
    o=str(o)
    if o.startswith('Clade'): return 'auris'
    if o=='Hae_1768': return 'haemulonii'
    if o=='Hae_1724': return 'hae_1724'
    if o=='Hae_1769': return None
    if o.startswith('Duo'): return 'duobushaemulonii'
    if o.startswith('para'): return 'parapsilosis'
    return None

def _raw():
    d=pd.read_csv('results/tables/derived_N0_R_results_with_carbon.csv')
    d=d[d.fit_valid.astype(str).str.lower().isin(['true','1'])]
    d=d.copy(); d['sp']=d.otu_name.map(sp_map)
    return d[d.sp.notna()].dropna(subset=['growth_C_per_C_h','respiration_C_per_C_h'])

def _tidy(g,groups,key,allT):
    rows=[]
    for s in groups:
        have=set(g[g[key]==s]['T'])
        for T in allT:
            if T in have:
                rr=g[(g[key]==s)&(g['T']==T)].iloc[0]
                rows.append((s,T,rr.growth_C_per_C_h,rr.respiration_C_per_C_h,True))
            else: rows.append((s,T,0.,np.nan,False))
    return pd.DataFrame(rows,columns=[key,'T','mu','resp','alive'])

def measured():
    """species-level medians (v5 definition)"""
    d=_raw()
    g=d.groupby(['sp','T'])[['growth_C_per_C_h','respiration_C_per_C_h']].median().reset_index()
    return _tidy(g,SPS+['hae_1724'],'sp',sorted(d['T'].unique()))

def measured_isolates():
    """isolate-level medians across wells; returns (frame with 'iso','sp', ...)"""
    d=pd.read_csv('results/tables/derived_N0_R_results_with_carbon.csv')
    d=d[d.fit_valid.astype(str).str.lower().isin(['true','1'])].copy()
    d=d.dropna(subset=['growth_C_per_C_h','respiration_C_per_C_h'])
    d['iso']=d.otu_name.astype(str)
    def iso_sp(o):
        if o.startswith('Clade'): return 'auris'
        if o.startswith('Hae'): return 'haemulonii'
        if o.startswith('Duo'): return 'duobushaemulonii'
        if o.startswith('para'): return 'parapsilosis'
        return None
    d['sp']=d.iso.map(iso_sp); d=d[d.sp.notna()]
    allT=sorted(d['T'].unique())
    g=d.groupby(['iso','T'])[['growth_C_per_C_h','respiration_C_per_C_h']].median().reset_index()
    out=_tidy(g,sorted(d.iso.unique()),'iso',allT)
    out['sp']=out.iso.map(iso_sp)
    return out

def raw_wells():
    """well-level rows for bootstrap: columns otu_name, sp, T, Replicate?, growth, resp"""
    d=_raw()
    return d

def unpack(z):
    return (np.exp(z[0]),np.exp(z[1]),np.exp(z[2]),np.exp(z[3]),np.exp(z[4]),np.exp(z[5]),
            0.15+0.55/(1+np.exp(-z[6])),z[7],np.exp(z[8]))

def predict(z,u,Ts):
    kappa,g0,Ea,m0,cg,cd,PHI,TmR,Nr=unpack(z)
    TK=np.asarray(Ts,float)+273.15
    phiD=kappa*u/np.clip(1-u,1e-9,None); phiG=np.clip(PHI-phiD,0.,None)
    fR=two_state(Nr,TmR,Ts)[:,0]
    mu=g0*np.exp(-Ea/R_gas*(1/TK-1/303.15))*fR*phiG
    resp=m0+cg*mu+cd*np.minimum(phiD,PHI)
    return mu,resp

def loss_species(z,u,Ts,ym,yr):
    ok=~np.isnan(yr)
    mu,re=predict(z,u,Ts)
    if not(np.all(np.isfinite(mu)) and np.all(np.isfinite(re))): return 1e6
    gl=np.mean((mu-ym)**2)/max(np.var(ym),1e-12)
    rl=np.mean((re[ok]-yr[ok])**2)/max(np.var(yr[ok]),1e-12) if ok.any() else 0
    return gl+rl

def fit_auris_background(u,Ts,ym,yr,n_restart=60,seed=3):
    """the v5 auris fit: all 9 parameters free, multi-start Nelder-Mead"""
    def loss_au(z): return loss_species(z,u,Ts,ym,yr)
    z0=np.array([np.log(7.5),np.log(1.1),np.log(60000.),np.log(0.20),
                 np.log(0.15),np.log(2.0),0.,41.,np.log(400.)])
    rng=np.random.default_rng(seed); best=None
    for i in range(n_restart):
        s=z0.copy()
        if i:
            s+=rng.normal(0,.5,9)*np.array([1,1,.3,1,1,1,1,0,1]); s[7]=rng.uniform(36,50)
        r=minimize(loss_au,s,method='Nelder-Mead',options=dict(maxiter=4000,maxfev=4000,fatol=1e-12,xatol=1e-10))
        if best is None or r.fun<best.fun: best=r
    return best.x,best.fun

def fit_tmr_ea(z_bg,u,Ts,ym,yr,grid_tmr=np.arange(30,50,2),grid_ea=(30e3,50e3,80e3,120e3,160e3)):
    """the v5 per-species fit: TmR and Ea free, everything else from z_bg"""
    def loss2(x):
        z=z_bg.copy(); z[7]=x[0]; z[2]=np.log(max(x[1],1000.))
        return loss_species(z,u,Ts,ym,yr)
    best=None
    for t0 in grid_tmr:
        for e0 in grid_ea:
            r=minimize(loss2,[t0,e0],method='Nelder-Mead',options=dict(maxiter=1000,fatol=1e-10))
            if best is None or r.fun<best.fun: best=r
    return best.x[0],best.x[1],best.fun

# ── IVYWREL composition of ribosomal proteins (from model_v6) ─────────────────
ZELDOVICH_SLOPE=937.0
def _fasta(p):
    h,s=None,[]
    for l in open(p):
        if l.startswith('>'):
            if h: yield h.split()[0],''.join(s)
            h,s=l[1:].strip(),[]
        else: s.append(l.strip())
    if h: yield h.split()[0],''.join(s)
def compute_fivywrel(path='gem/tables/machinery_ribosome.faa'):
    acc={}
    for h,seq in _fasta(path):
        sp=h.split('|')[0]; acc.setdefault(sp,[0,0])
        acc[sp][0]+=sum(1 for aa in seq if aa in set('IVYWREL')); acc[sp][1]+=len(seq)
    return {sp:v[0]/v[1] for sp,v in acc.items()}
