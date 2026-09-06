import cobra,pandas as pd,numpy as np,json,statistics,sys
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
MODE=sys.argv[1]; N=int(sys.argv[2])
S_TOPT, S_TM = 12.26, 7.57          # Seq2Topt / Seq2Tm published RMSE (predictor README)
cal=json.load(open(TABLES / 'etcgem_calib.json')); SIG,W,P0,SC=cal['sig'],cal['w'],cal['P'],cal['SCALE']
to=pd.read_csv(TABLES / 'thermal_topt.csv'); tmv=pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (to,tmv): df['g']=df['id'].str.split('|').str[1]; df['sp']=df['id'].str.split('|').str[0]
SPX={'auris':('auris_iRV973_rekeyed.xml','medium_iRV973_auris.csv'),
     'haemulonii':('haemulonii_draft.xml','medium_iRV973_auris.csv'),
     'duobushaemulonii':('duobushaemulonii_draft.xml','medium_iRV973_auris.csv'),
     'parapsilosis':('parapsilosis_iDC1003.xml','medium_iDC1003_parapsilosis.csv')}
REL=['haemulonii','duobushaemulonii','parapsilosis']; DET=0.05
def build(sp):
    xml,mf=SPX[sp]; m=cobra.io.read_sbml_model(str(MODELS / xml)); med=pd.read_csv(INPUTS / mf)
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
    mw={r['gene']:float(r['MW_kDa']) for _,r in pd.read_csv(TABLES / f'enzyme_mw_{sp}.csv').iterrows()}
    kc=pd.read_csv(TABLES / f'kcat_reaction_{sp}.csv'); km={r['reaction']:(float(r['kcat_per_s']),r['best_gene']) for _,r in kc.iterrows()}
    tov={r.g:r.pred_topt for r in to[to.sp==sp].itertuples()}; tmm={r.g:r.pred_tm for r in tmv[tmv.sp==sp].itertuples()}
    mt=statistics.median(tov.values()); mm=statistics.median(tmm.values())
    fv=[];rv=[];bc=[];To=[];Tm=[]
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        ws=[mw[g.id] for g in r.genes if g.id in mw]
        if not ws: continue
        ks,bg=km.get(r.id,(13.7,None))
        fv.append(r.forward_variable); rv.append(r.reverse_variable)
        bc.append(statistics.mean(ws)/(ks*3600.0)); To.append(tov.get(bg,mt)); Tm.append(tmm.get(bg,mm))
    con=m.problem.Constraint(0.0,ub=P0,name='pp'); m.add_cons_vars(con); m.solver.update()
    return dict(m=m,con=con,fv=fv,rv=rv,bc=np.array(bc),To=np.array(To),Tm=np.array(Tm))
def mu(st,T,dTo,dTm):
    a=np.exp(-((T-(st['To']+dTo))**2)/(2*SIG*SIG))/(1+np.exp((T-(st['Tm']+dTm))/W))
    c=st['bc']/np.maximum(a,1e-6)
    st['con'].set_linear_coefficients(dict(zip(st['fv'],c))|dict(zip(st['rv'],c)))
    g=st['m'].slim_optimize(); return (g if g==g else 0.0)*SC
S={sp:build(sp) for sp in SPX}
# sanity: zero offset must reproduce the pipeline
base={sp:{T:mu(S[sp],T,0.0,0.0) for T in (34,40,44)} for sp in SPX}
P=pd.read_csv(TABLES / 'etcgem_tpc_pred.csv').pivot_table(index='T',columns='species',values='pred_mu')
ok=all(abs(base[sp][T]-P[sp].loc[T])<1e-3 for sp in SPX for T in (34,40,44))
print(f"sanity: zero-offset reproduces pipeline = {ok}")
rng=np.random.default_rng(2024)
rows=[]
for d in range(N):
    off={}; v={}
    for sp in SPX:
        n=len(S[sp]['To'])
        if MODE=='systematic':
            a,b=rng.normal(0,S_TOPT), rng.normal(0,S_TM)
            jo=np.full(n,a); jm=np.full(n,b)
        else:
            a=b=np.nan; jo=rng.normal(0,S_TOPT,n); jm=rng.normal(0,S_TM,n)
        off[sp]=(a,b); v[sp]={T:mu(S[sp],T,jo,jm) for T in (34,40,44)}
    rows.append(dict(draw=d,
        **{f'mu34_{s}':v[s][34] for s in SPX}, **{f'mu40_{s}':v[s][40] for s in SPX},
        **{f'mu44_{s}':v[s][44] for s in SPX},
        **{f'offTopt_{s}':off[s][0] for s in SPX}, **{f'offTm_{s}':off[s][1] for s in SPX}))
D=pd.DataFrame(rows); D.to_csv(TABLES / f'ens_final_{MODE}.csv',index=False)
def wilson(k,n):
    if n==0: return (np.nan,np.nan)
    z=1.96; p=k/n; d=1+z*z/n
    c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return max(0,c-h),min(1,c+h)
print(f"\n=== MODE={MODE} | N={N} | sigma Topt={S_TOPT} Tm={S_TM} ===")
pheno=(D.mu44_auris>=DET)&(D.mu40_haemulonii<DET)&(D.mu40_duobushaemulonii<DET)&(D.mu40_parapsilosis<DET)
rank=D.mu44_auris>D[[f'mu44_{s}' for s in REL]].max(axis=1)
print(f"  P(auris ranks top at 44C) = {rank.mean():.3f}  95%CI {tuple(round(x,3) for x in wilson(rank.sum(),N))}  (chance 0.25)")
for floor in (0.20,0.30,0.40):
    viable=(D[[f'mu34_{s}' for s in SPX]]>=floor).all(axis=1)
    nv=int(viable.sum()); k=int((pheno&viable).sum())
    ci=wilson(k,nv)
    print(f"  permissive floor mu34>={floor:.2f}: viable draws {nv}/{N}"
          f" | phenotype among viable {k}/{nv}"
          + (f" = {k/nv:.3f}  95%CI {tuple(round(x,3) for x in ci)}" if nv else ""))
sub=D[pheno]
if len(sub):
    print(f"\n  in the {len(sub)} phenotype draws, realised offsets (auris minus mean relative):")
    dTo=sub.offTopt_auris-sub[[f'offTopt_{s}' for s in REL]].mean(axis=1)
    dTm=sub.offTm_auris-sub[[f'offTm_{s}' for s in REL]].mean(axis=1)
    if np.isfinite(dTo).any():
        print(f"    dTopt median {np.nanmedian(dTo):+.2f} C | dTm median {np.nanmedian(dTm):+.2f} C")
