"""Model-free cross-kingdom table. Per taxon: T_opt = temperature of maximum median r,
T_carbon = temperature of minimum median K/r (respiration per unit growth, scale-free),
both on the measured grid, with a within-temperature well bootstrap. Only wells with a
detected growth fit enter (as in the manuscript)."""
import os, numpy as np, pandas as pd
rng=np.random.default_rng(11); NB=2000
P=os.path.expanduser("~/mnt/Projects"); C=os.path.expanduser("~/mnt/Candidas"); OUT=f"{C}/model_extra/cross_kingdom"
sets=[]   # (taxon, kingdom, dataset, df[T,r,K])
# Candida (manuscript wells; r,K per minute)
cd=pd.read_csv(f"{C}/results/tables/fit_coefficients_wide.csv"); on=pd.read_csv(f"{C}/results/tables/otu_names.csv")
cd=cd.merge(on,on="OTU"); cd=cd[(cd.keep)&(cd.fit_valid)]
names={"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV","para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
for g,s in cd.groupby("group"): sets.append((names[g],"Fungi","Candidas (this study)",s[["T","r","K"]]))
# R2A bacteria + E. coli (this session; r,K per hour)
r2=pd.read_csv(f"{OUT}/r2a_wells.csv"); r2=r2[r2.state=="growing"]
for st,s in r2.groupby("strain"): sets.append((st if st.startswith("E.") else f"Soil {st}","Bacteria","R2A series" if not st.startswith("E.") else "E. coli R2A series",s[["T","r","K"]]))
# LB soil OTUs (Parsa; per minute)
pb=pd.read_csv(f"{P}/Parsa_bacteria/tables/fit_coefficients_wide.csv")
for o,s in pb.groupby("OTU"): sets.append((f"LB OTU{o}","Bacteria","LB series",s[["T","r","K"]]))
# MicroAdapt Pseudomonas, E. coli M9
q=pd.read_csv(f"{P}/MicroAdapt/Experiments/Temperature_23_05/oxygen_model_results_good_only.csv")
sets.append(("Pseudomonas sp.","Bacteria","MicroAdapt",q.rename(columns={"Temperature":"T","r_per_hour":"r","resp_rate":"K"})[["T","r","K"]]))
e=pd.read_csv(f"{P}/MicroAdapt/Experiments/Ecoli_extended_temp/oxygen_model_results_good_only.csv")
sets.append(("E. coli (M9)","Bacteria","MicroAdapt",e.rename(columns={"Temperature":"T","r_per_hour":"r","resp_rate":"K"})[["T","r","K"]]))
# Zymoseptoria, two experiments batch-scaled on shared temperatures
def zymo(path):
    d=pd.read_csv(path); w=d.pivot_table(index=["T","Condition","Replicate"],columns="parameter",values="Estimate").reset_index()
    return w[w.Condition=="Control"].rename(columns={"resp_rate":"K"})[["T","r","K"]]
za=zymo(f"{P}/Zymo/Zymo_3days_redone/o2_model_from_filtered/o2_model_coefficients_from_filtered.csv").assign(b="A")
zb=zymo(f"{P}/Zymo/Zymo_3_days_TPC/o2_model_outputs/o2_model_coefficients.csv").assign(b="B")
z=pd.concat([za,zb]); z=z[(z.r>0)&(z.K>0)]; ov=sorted(set(za["T"])&set(zb["T"]))
for col in ("r","K"):
    f=np.exp(np.mean(np.log(z[(z.b=="A")&z["T"].isin(ov)].groupby("T")[col].median()/z[(z.b=="B")&z["T"].isin(ov)].groupby("T")[col].median())))
    z.loc[z.b=="B",col]*=f
sets.append(("Zymoseptoria tritici","Fungi","Zymo series",z[["T","r","K"]]))

def one(df):
    df=df[(df.r>0)&(df.K>0)].copy(); df["q"]=df.K/df.r
    n=df.groupby("T").size(); Ts=n[n>=3].index.values; df=df[df["T"].isin(Ts)]
    med=df.groupby("T").agg(r=("r","median"),q=("q","median"))
    Topt=med.r.idxmax(); Tc=med.q.idxmin()
    groups={T:s for T,s in df.groupby("T")}
    bo=[];bc=[]
    for _ in range(NB):
        m={T:(s.r.values[i:=rng.integers(0,len(s),len(s))],s.q.values[i]) for T,s in groups.items()}
        rr=pd.Series({T:np.median(v[0]) for T,v in m.items()}); qq=pd.Series({T:np.median(v[1]) for T,v in m.items()})
        bo.append(rr.idxmax()); bc.append(qq.idxmin())
    bo=np.array(bo); bc=np.array(bc)
    return dict(nT=len(Ts),Tlo=Ts.min(),Thi=Ts.max(),step=int(np.min(np.diff(np.sort(Ts)))),n_wells=len(df),
        T_opt=Topt,T_opt_lo=np.quantile(bo,.025),T_opt_hi=np.quantile(bo,.975),
        T_carbon=Tc,T_carbon_lo=np.quantile(bc,.025),T_carbon_hi=np.quantile(bc,.975),
        P_below=np.mean(bc<bo),P_below_or_equal=np.mean(bc<=bo),gap=Topt-Tc,
        q_rise_past_opt=float(med.q.loc[med.index>=Topt].max()/med.q.loc[Topt]) if (med.index>Topt).any() else np.nan,
        curve=med)
rows=[];curves=[]
for tax,kg,ds,df in sets:
    o=one(df); cv=o.pop("curve"); rows.append(dict(taxon=tax,kingdom=kg,dataset=ds,**o))
    curves.append(cv.reset_index().assign(taxon=tax,kingdom=kg,q_rel=lambda d:d.q/d.q.min(),r_rel=lambda d:d.r/d.r.max(),dT=lambda d:d["T"]-o["T_opt"]))
t=pd.DataFrame(rows); t.to_csv(f"{OUT}/cross_kingdom_modelfree.csv",index=False)
pd.concat(curves).to_csv(f"{OUT}/cross_kingdom_curves.csv",index=False)
pd.set_option("display.width",250); print(t.drop(columns=["dataset"]).round(2).to_string())
print("\nT_carbon < T_opt strictly:",(t.gap>0).sum(),"of",len(t),"| ties:",(t.gap==0).sum(),"| above:",(t.gap<0).sum())
print("median P(T_carbon<T_opt):",t.P_below.median().round(3))
