#!/usr/bin/env python3
"""Counterfactual sufficiency test for the etcGEM.
Question (two-sided, per ChatGPT correction): what MINIMUM interspecies difference in a
model parameter is required to reproduce BOTH sides of the phenotype —
  (i)  C. auris still grows (detectably) at 44 C, AND
  (ii) the median close-relative fails (below detection) by 40 C,
while (iii) permissive-temperature growth stays approximately correct?
Then compare the model-required separation with the separation predicted from sequence.
It is a model-implied minimum, not a biological Tm estimate, and 44 C is assay-censored.
"""
import cobra, pandas as pd, numpy as np, json, statistics
from pathlib import Path
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
G = GEM
SP={'auris':('auris_iRV973_rekeyed.xml','medium_iRV973_auris.csv'),
    'haemulonii':('haemulonii_draft.xml','medium_iRV973_auris.csv'),
    'duobushaemulonii':('duobushaemulonii_draft.xml','medium_iRV973_auris.csv'),
    'parapsilosis':('parapsilosis_iDC1003.xml','medium_iDC1003_parapsilosis.csv')}
RELS=['haemulonii','duobushaemulonii','parapsilosis']
calib=json.load(open(TABLES / 'etcgem_calib.json')); SIG,W,P0,SCALE=calib['sig'],calib['w'],calib['P'],calib['SCALE']
topt=pd.read_csv(TABLES / 'thermal_topt.csv'); tm=pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (topt,tm):
    df['g']=df['id'].str.split('|').str[1]; df['sp']=df['id'].str.split('|').str[0]

def therm_maps(sp):
    tov={r.g:r.pred_topt for r in topt[topt.sp==sp].itertuples()}
    tmv={r.g:r.pred_tm for r in tm[tm.sp==sp].itertuples()}
    return tov,tmv,(statistics.median(tov.values()) if tov else 35.0),(statistics.median(tmv.values()) if tmv else 54.0)

def setup_pool(sp):
    xml,medf=SP[sp]; m=cobra.io.read_sbml_model(str(MODELS / xml)); med=pd.read_csv(INPUTS / medf)
    EX={r.id for r in m.reactions if r.id.startswith(("EX_","Drain")) or r.boundary}
    for r in m.reactions:
        if r.id in EX: r.lower_bound=0.0
    for _,row in med.iterrows():
        rid=str(row['exchange_id'])
        if rid in EX:
            if row['setting']=='OPEN': m.reactions.get_by_id(rid).lower_bound=-1000
            elif row['setting']=='FIT': m.reactions.get_by_id(rid).lower_bound=-10
    for _,row in med[med.setting=='AA_POOL'].iterrows():
        if str(row['exchange_id']) in EX: m.reactions.get_by_id(str(row['exchange_id'])).lower_bound=-3
    # ATP maintenance: irreversible and forced, in every species.
    # The three KEGG-based models ship it pinned at 3.89, but iDC1003 (C. parapsilosis)
    # ships it REVERSIBLE at lb=-3.9, and the solver takes that direction at every
    # temperature -- synthesising ATP from ADP + Pi. C. parapsilosis was therefore paying
    # no maintenance AND collecting 3.9 mmol/gDW/h of free ATP, a 7.8 swing against the
    # other three, which inflated its predicted mu by 11-15% and made it look unreachable
    # by any uniform Tm shift. Corrected, it is reachable at 32.6 C -- next to
    # C. duobushaemulonii's 32.4 C.
    for _rid in ('ATP_maintenance_NGAM__cyto','ATP_Maintenance__cyto'):
        if _rid in m.reactions:
            _r=m.reactions.get_by_id(_rid); _v=abs(_r.upper_bound); _r.bounds=(_v,_v)
    mw={r['gene']:float(r['MW_kDa']) for _,r in pd.read_csv(TABLES / f'enzyme_mw_{sp}.csv').iterrows()}
    kc=pd.read_csv(TABLES / f'kcat_reaction_{sp}.csv')
    kmap={r['reaction']:(float(r['kcat_per_s']),r['best_gene']) for _,r in kc.iterrows()}
    return m,EX,mw,kmap

def precompute(m,EX,mw,kmap,tov,tmv,med_to,med_tm):
    terms=[]
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws=[mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        kcat_s,bg=kmap.get(r.id,(13.7,None))
        base_coef=statistics.mean(mws)/(kcat_s*3600.0)
        terms.append((r.forward_variable,r.reverse_variable,base_coef,tov.get(bg,med_to),tmv.get(bg,med_tm)))
    con=m.problem.Constraint(0.0,ub=1.0,name='pp'); m.add_cons_vars(con); m.solver.update()
    return terms,con

def act(T,Topt,Tm):
    return np.exp(-((T-Topt)**2)/(2.0*SIG*SIG))*(1.0/(1.0+np.exp((T-Tm)/W)))

def mu(S,T,dTopt=0.0,dTm=0.0,Pfac=1.0):
    m,terms,con=S['m'],S['terms'],S['con']; coeffs={}
    for fv,rv,bc,Topt,Tm in terms:
        a=max(act(T,Topt+dTopt,Tm+dTm),1e-6); c=bc/a; coeffs[fv]=c; coeffs[rv]=c
    con.ub=P0*Pfac; con.set_linear_coefficients(coeffs)
    try: g=m.slim_optimize()
    except Exception: g=float('nan')
    return (g if g==g else 0.0)*SCALE

# ---- setup all species ----
STATE={}
for sp in SP:
    tov,tmv,mto,mtm=therm_maps(sp); m,EX,mw,kmap=setup_pool(sp)
    terms,con=precompute(m,EX,mw,kmap,tov,tmv,mto,mtm)
    STATE[sp]=dict(m=m,terms=terms,con=con,med_to=mto,med_tm=mtm)
    print(f"setup {sp}: {len(terms)} enzyme-costed reactions",flush=True)

# ---- sequence-predicted separations (auris vs each relative) ----
def med(df,sp,col):
    v=df[df.sp==sp][col]; return float(np.median(v)) if len(v) else np.nan
def pct1(df,sp,col):
    v=df[df.sp==sp][col]; return float(np.percentile(v,1)) if len(v) else np.nan
seq={}
for rel in RELS:
    seq[rel]=dict(dTopt_med=med(topt,'auris','pred_topt')-med(topt,rel,'pred_topt'),
                  dTm_med=med(tm,'auris','pred_tm')-med(tm,rel,'pred_tm'),
                  dTm_p1=pct1(tm,'auris','pred_tm')-pct1(tm,rel,'pred_tm'))
print("sequence-predicted auris-minus-relative separations:")
for rel in RELS: print(f"  {rel}: dTopt_med={seq[rel]['dTopt_med']:.2f}  dTm_med={seq[rel]['dTm_med']:.2f}  dTm_1pct={seq[rel]['dTm_p1']:.2f}")

# ---- baseline sanity ----
base={sp:{T:mu(STATE[sp],T) for T in (34,40,42,44)} for sp in SP}
print("baseline predicted mu:", {sp:{T:round(v,3) for T,v in d.items()} for sp,d in base.items()})

def bisect_offset(sp,T,param,thr,lo=-35.0,hi=0.0,tol=0.1):
    """smallest downward offset (most negative to 0) on param making mu(sp,T,offset) < thr.
    returns offset (<=0) at the crossing, or lo if even lo doesn't reach thr."""
    f=lambda d: mu(STATE[sp],T, dTopt=(d if param=='Topt' else 0.0), dTm=(d if param=='Tm' else 0.0))
    if f(0.0) < thr: return 0.0                     # already fails at baseline
    if f(lo) >= thr: return None                    # cannot fail even at extreme
    a,b=lo,0.0
    for _ in range(40):
        mid=(a+b)/2
        if f(mid) < thr: a=mid
        else: b=mid
        if b-a<tol: break
    return (a+b)/2

RESULTS=dict(seq=seq, baseline={sp:{str(T):round(v,4) for T,v in d.items()} for sp,d in base.items()},
            calib=calib, mechanisms={})
THRS=[0.03,0.05,0.10]  # detection thresholds (measured mu units); 0.05 primary
PERM_FRAC=0.7          # permissive (34C) growth must stay >= 0.7x baseline

for thr in THRS:
    key=f"thr={thr}"
    RESULTS['mechanisms'][key]={}
    # auris must still grow at 44 unshifted:
    auris_ok = base['auris'][44] >= thr
    for mech,param in (('Tm','Tm'),('Topt','Topt')):
        per_rel={}
        for rel in RELS:
            off=bisect_offset(rel,40,param,thr)     # offset to make rel fail at 40C
            if off is None: per_rel[rel]=dict(required=None,note="cannot fail at 40C within +/-35C"); continue
            # permissive preserved? check mu(rel,34) at this offset
            perm=mu(STATE[rel],34, dTopt=(off if param=='Topt' else 0.0), dTm=(off if param=='Tm' else 0.0))
            perm_ok = perm >= PERM_FRAC*base[rel][34]
            # model-required interspecies separation = auris(0) - rel(off) = -off  (auris unshifted, ok at 44)
            per_rel[rel]=dict(required=round(-off,2), perm_growth=round(perm,3),
                              perm_baseline=round(base[rel][34],3), permissive_preserved=bool(perm_ok))
        reqs=[v['required'] for v in per_rel.values() if v.get('required') is not None]
        RESULTS['mechanisms'][key][mech]=dict(auris_grows_at_44=bool(auris_ok), per_rel=per_rel,
            median_required=(round(float(np.median(reqs)),2) if reqs else None))
    # flat pool mechanism: can lowering P make rel fail at 40 while keeping 34? (test selectivity)
    poolres={}
    for rel in RELS:
        # find Pfac making mu(rel,40)<thr
        lo,hi=0.01,1.0
        if mu(STATE[rel],40,Pfac=1.0)<thr: pf=1.0
        else:
            a,b=lo,hi
            for _ in range(30):
                mid=(a+b)/2
                if mu(STATE[rel],40,Pfac=mid)<thr: a=mid
                else: b=mid
            pf=(a+b)/2
        perm=mu(STATE[rel],34,Pfac=pf)
        poolres[rel]=dict(Pfac=round(pf,3),growth40=round(mu(STATE[rel],40,Pfac=pf),3),
                          growth34=round(perm,3),growth34_baseline=round(base[rel][34],3),
                          still_grows_permissive=bool(perm>=PERM_FRAC*base[rel][34]))
    RESULTS['mechanisms'][key]['pool_flat']=poolres

# ---- sweep for figure panel B: mu vs imposed Tm shift (auris + representative relative) ----
sweep=[]
for sp in ['auris','duobushaemulonii']:
    for d in np.arange(-35,6,1.0):
        for T in (34,40,44):
            sweep.append(dict(species=sp,dTm=float(d),T=T,mu=round(mu(STATE[sp],T,dTm=d),4)))
pd.DataFrame(sweep).to_csv(TABLES / 'counterfactual_sweep.csv',index=False)

json.dump(RESULTS,open(TABLES / 'counterfactual_results.json','w'),indent=2)
print("\n=== PRIMARY (thr=0.05) ===")
pr=RESULTS['mechanisms']['thr=0.05']
for mech in ('Tm','Topt'):
    print(f"\n{mech}: auris grows@44={pr[mech]['auris_grows_at_44']}, median required separation={pr[mech]['median_required']} C")
    for rel,v in pr[mech]['per_rel'].items():
        print(f"   {rel}: required {v.get('required')} C | permissive preserved={v.get('permissive_preserved')} (mu34 {v.get('perm_growth')} vs base {v.get('perm_baseline')})")
    sseq=np.median([seq[r]['dTm_med' if mech=='Tm' else 'dTopt_med'] for r in RELS])
    print(f"   sequence-predicted median separation: {sseq:.2f} C")
print("\npool_flat (selectivity test):")
for rel,v in pr['pool_flat'].items():
    print(f"   {rel}: Pfac={v['Pfac']} -> mu40={v['growth40']}, mu34={v['growth34']} (base {v['growth34_baseline']}), permissive kept={v['still_grows_permissive']}")
print("\nwrote counterfactual_results.json")
