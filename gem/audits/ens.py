import cobra, pandas as pd, numpy as np, json, statistics, sys
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
SIGMA=float(sys.argv[1]); KSIG=float(sys.argv[2]); N=int(sys.argv[3]); TAG=sys.argv[4]
cal=json.load(open(TABLES / 'etcgem_calib.json')); SIG,W,P0,SC=cal['sig'],cal['w'],cal['P'],cal['SCALE']
to=pd.read_csv(TABLES / 'thermal_topt.csv'); tmv=pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (to,tmv): df['g']=df['id'].str.split('|').str[1]; df['sp']=df['id'].str.split('|').str[0]
SPX={'auris':('auris_iRV973_rekeyed.xml','medium_iRV973_auris.csv'),
     'haemulonii':('haemulonii_draft.xml','medium_iRV973_auris.csv'),
     'duobushaemulonii':('duobushaemulonii_draft.xml','medium_iRV973_auris.csv'),
     'parapsilosis':('parapsilosis_iDC1003.xml','medium_iDC1003_parapsilosis.csv')}
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
    con=m.problem.Constraint(0.0,ub=1.0,name='pp'); m.add_cons_vars(con); m.solver.update()
    return dict(m=m,con=con,fv=fv,rv=rv,bc=np.array(bc),To=np.array(To),Tm=np.array(Tm))
def mu(st,T,dTo=0.0,dTm=0.0,klog=None):
    a=np.exp(-((T-(st['To']+dTo))**2)/(2*SIG*SIG))/(1+np.exp((T-(st['Tm']+dTm))/W))
    bc=st['bc'] if klog is None else st['bc']/(10.0**klog)
    c=bc/np.maximum(a,1e-6)
    st['con'].set_linear_coefficients(dict(zip(st['fv'],c))|dict(zip(st['rv'],c)))
    g=st['m'].slim_optimize(); return (g if g==g else 0.0)*SC
S={sp:build(sp) for sp in SPX}
for sp in SPX: S[sp]['con'].ub=P0
REL=['haemulonii','duobushaemulonii','parapsilosis']; THR=0.05
rng=np.random.default_rng(7); rec={sp:{34:[],40:[],44:[]} for sp in SPX}; joint=0; joint_ok=0; rank=0
for d in range(N):
    vals={}
    for sp in SPX:
        n=len(S[sp]['To'])
        jo=rng.normal(0,SIGMA,n); jm=rng.normal(0,SIGMA,n)
        kl=rng.normal(0,KSIG,n) if KSIG>0 else None
        vals[sp]={T:mu(S[sp],T,jo,jm,kl) for T in (34,40,44)}
        for T in (34,40,44): rec[sp][T].append(vals[sp][T])
    ok34=all(vals[sp][34]>=0.30 for sp in SPX)          # permissive growth must be retained
    if vals['auris'][44]>=THR and all(vals[r][40]<THR for r in REL):
        joint+=1
        if ok34: joint_ok+=1
    if vals['auris'][44]>max(vals[r][44] for r in REL): rank+=1
print(f"=== {TAG} | N={N} | INDEPENDENT per-species errors (most favourable to the model) ===")
res={}
for sp in SPX:
    a={T:np.array(rec[sp][T]) for T in (34,40,44)}
    print(f"  {sp:18s} mu44 med {np.median(a[44]):.4f} [95% {np.percentile(a[44],2.5):.4f}-{np.percentile(a[44],97.5):.4f}]"
          f" | P(mu44>=det) {np.mean(a[44]>=THR):.2f} | P(mu40<det) {np.mean(a[40]<THR):.2f} | P(mu34>=det) {np.mean(a[34]>=THR):.2f}")
    res[sp]={str(T):[float(np.median(a[T])),float(np.percentile(a[T],2.5)),float(np.percentile(a[T],97.5))] for T in (34,40,44)}
print(f"  --> P(auris grows@44 AND all relatives fail@40)                    = {joint/N:.3f}")
print(f"  --> ...AND all species retain permissive growth at 34C (mu>=0.30)  = {joint_ok/N:.3f}")
print(f"  --> P(auris ranks above every relative at 44C)                     = {rank/N:.3f}   (chance 0.25)")
json.dump({'tag':TAG,'sigma':SIGMA,'ksig':KSIG,'N':N,'joint':joint/N,'joint_permissive_ok':joint_ok/N,'rank':rank/N,'quantiles':res},
          open(TABLES / f'ensemble_{TAG}.json','w'),indent=2)
