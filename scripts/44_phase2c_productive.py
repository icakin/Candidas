"""PHASE 2C - productive-respiration fraction.
For each isolate x temperature: the fraction of observed oxygen consumption (fixed
1.5-6.5 h window, Phase 2B) that occurs in CONFIRMED-GROWTH wells.
  LOWER bound: ambiguous consumption counted as NON-productive
  UPPER bound: ambiguous consumption counted as productive
NOTE ON CENSORING: growing wells hit the oxygen floor inside the fixed window far more
often (45%) than respiration-only wells (7%), so their consumption is right-censored and
the productive fraction reported here is itself a LOWER bound wherever that occurs.
"""
import os, numpy as np, pandas as pd
C=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
c=pd.read_csv(f"{C}/results/tables/phase2b_consumption.csv")
c=c[np.isfinite(c.cons)].copy()
rows=[]
for (g,iso,T),s in c.groupby(["group","otu_name","T"]):
    tot=s.cons.sum()
    if tot<=0: continue
    prod_lo=s[s.S=="growth"].cons.sum()
    prod_hi=s[s.S.isin(["growth","ambig"])].cons.sum()
    rows.append(dict(group=g,isolate=iso,T=T,tot=tot,
                     frac_lo=prod_lo/tot, frac_hi=prod_hi/tot,
                     n_wells=len(s), cens_frac=s.censored.mean()))
p=pd.DataFrame(rows)
p["clade"]=np.where(p.group.str.startswith("Clade"),"C. auris","relatives")
p.to_csv(f"{C}/results/tables/phase2c_productive.csv",index=False)
pd.set_option("display.width",200)
print("=== productive-respiration fraction, LOWER bound (ambiguous = non-productive) ===")
print(p.pivot_table(index=["clade","isolate"],columns="T",values="frac_lo").round(2)[[34,36,38,40,42,44]].to_string())
print("\n=== UPPER bound (ambiguous = productive) ===")
print(p.pivot_table(index=["clade","isolate"],columns="T",values="frac_hi").round(2)[[34,36,38,40,42,44]].to_string())
def collapse(df,col):
    w=df.pivot_table(index="isolate",columns="T",values=col)
    return (w[38]-w[40])
print("\n=== COLLAPSE between 38 and 40 °C (positive = productive fraction falls) ===")
for col,lab in (("frac_lo","LOWER bound"),("frac_hi","UPPER bound")):
    print(f"\n  {lab}:")
    for grp,sub in (("C. auris",p[p.clade=="C. auris"]),("relatives",p[p.clade=="relatives"]),
                    ("relatives, no Hae_1724",p[(p.clade=="relatives")&(p.isolate!="Hae_1724")])):
        d=collapse(sub,col).dropna()
        print(f"    {grp:24s} median drop {d.median():+.2f}   isolates dropping >0.5: {(d>0.5).sum()}/{len(d)}")
print("\n=== separate contrasts requested ===")
for lab,sel in (("C. auris vs C. parapsilosis",["para"]),
                ("C. auris vs C. haemulonii + C. duobushaemulonii",["Hae","Duo"]),
                ("...the same, excluding Hae_1724",["Hae","Duo"])):
    print(f"\n  {lab}")
    for col in ("frac_lo","frac_hi"):
        a=p[(p.clade=="C. auris")&(p["T"]==40)][col]
        r=p[(p.group.isin(sel))&(p["T"]==40)]
        if "excluding" in lab: r=r[r.isolate!="Hae_1724"]
        r=r[col]
        print(f"    {col} at 40 °C:  C. auris median {a.median():.2f} (n={len(a)})   "
              f"comparator median {r.median():.2f} (n={len(r)})")
