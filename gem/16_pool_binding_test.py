import cobra, csv, re, statistics, pandas as pd
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
GEMU = GEM
KCAT_CLASS={"1":13.7,"2":13.7,"3":79.0,"4":15.0,"5":6.9,"6":10.0}
def kcat_h(ecs):
    if not ecs: return 13.7*3600
    return statistics.median([KCAT_CLASS.get(e.split(".")[0],13.7) for e in ecs])*3600.0
mw={}
for row in csv.DictReader(open(f"{TABLES}/enzyme_mw_auris.csv")): mw[row["gene"]]=float(row["MW_kDa"])
def rx_ecs(r):
    a=r.annotation.get("ec-code",[]); a=[a] if isinstance(a,str) else a
    return {x for e in a for x in re.split(r"[;, ]+",e) if re.match(r"\d+\.\d+\.\d+\.\d+",x)}
m=cobra.io.read_sbml_model(f"{MODELS}/auris_iRV973_rekeyed.xml")
med=pd.read_csv(INPUTS / "medium_iRV973_auris.csv")
EX={r.id for r in m.reactions if r.id.startswith(("EX_","Drain")) or r.boundary}
for r in m.reactions:
    if r.id in EX: r.lower_bound=0.0
for _,row in med.iterrows():
    rid=str(row["exchange_id"])
    if rid in EX:
        if row["setting"]=="OPEN": m.reactions.get_by_id(rid).lower_bound=-1000
        elif row["setting"]=="FIT": m.reactions.get_by_id(rid).lower_bound=-10
for _,row in med[med.setting=="AA_POOL"].iterrows():
    if str(row["exchange_id"]) in EX: m.reactions.get_by_id(str(row["exchange_id"])).lower_bound=-3
g0=m.slim_optimize()
coeff={}
for r in m.reactions:
    if r.id in EX or "iomass" in r.id or not r.genes: continue
    mws=[mw[g.id] for g in r.genes if g.id in mw]
    if not mws: continue
    coeff[r]=statistics.mean(mws)/kcat_h(rx_ecs(r))
expr=sum((c*r.forward_variable + c*r.reverse_variable) for r,c in coeff.items())
for P in (0.25,0.10,0.05):
    con=m.problem.Constraint(expr, ub=P, name="pp")
    m.add_cons_vars(con); m.solver.update()
    g=m.slim_optimize()
    print(f"P={P:5} g/gDW  ecGEM growth={g:.3f}/h  (plain={g0:.3f})", flush=True)
    m.remove_cons_vars([con]); m.solver.update()
