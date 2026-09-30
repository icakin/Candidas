"""Carbon optimum vs niche temperature, one footing for every taxon with paired G and R.

For each taxon: Sharpe-Schoolfield on the per-well specific growth rate r, Arrhenius on the
per-well volumetric consumption rate K. Both scale-free: no N0, carbon quota, O2-to-C factor
or RQ enters, so nothing here depends on the assumptions flagged in the Candidas manuscript.

    T_opt(growth)  = SS optimum of the fitted curve
    T_opt(carbon)  = same formula with E replaced by E - E_R   (argmin of K/r)

Uncertainty by non-parametric bootstrap over wells, stratified by temperature.
Reports the median and 95% percentile interval, and P(T_carbon < niche temperature).
"""
import os, numpy as np, pandas as pd
from scipy.optimize import least_squares
k = 8.617e-5; TREF = 293.15; NBOOT = 600
rng = np.random.default_rng(7)
P = os.path.expanduser("~/mnt/Projects")
C = os.path.expanduser("~/mnt/Candidas")

def ss_topt(E, Eh, Th):
    if not (0 < E < Eh): return np.nan
    return 1/(1/Th - np.log(E/(Eh-E))/(Eh/k)) - 273.15

def fit_once(T, r, K):
    x = 1/(k*TREF) - 1/(k*(T+273.15)); Tk = T + 273.15
    def resid(p):
        lnB0, E, Eh, Th = p
        return lnB0 + E*x - np.log1p(np.exp(np.clip(Eh*(1/(k*Th)-1/(k*Tk)), -50, 50))) - np.log(r)
    g = least_squares(resid, [np.log(np.median(r)), 0.8, 3.0, np.median(Tk)+5],
                      bounds=([-40, 0.05, 0.3, 285], [30, 4, 25, 345]))
    lnB0, E, Eh, Th = g.x
    a = least_squares(lambda p: (p[0] + p[1]*x) - np.log(K), [np.log(np.median(K)), 0.5],
                      bounds=([-80, 0.0], [30, 3]))
    ER = a.x[1]
    return E, ER, Eh, Th, ss_topt(E, Eh, Th), ss_topt(E-ER, Eh, Th)

def analyse(name, df, niche, group, note=""):
    """df needs columns T, r, K (per well). r and K in any consistent units."""
    d = df.dropna(subset=["T","r","K"]); d = d[(d.r > 0) & (d.K > 0)]
    if d["T"].nunique() < 5 or len(d) < 12: 
        return dict(taxon=name, n=len(d), note="too few points")
    base = fit_once(d["T"].values.astype(float), d.r.values, d.K.values)
    boots = []
    for _ in range(NBOOT):
        idx = d.groupby("T").apply(lambda s: s.sample(len(s), replace=True, random_state=rng.integers(1e9)).index)
        b = d.loc[np.concatenate([np.asarray(i) for i in idx])]
        try: boots.append(fit_once(b["T"].values.astype(float), b.r.values, b.K.values))
        except Exception: pass
    B = np.array(boots)
    q = lambda col, p: np.nanpercentile(B[:, col], p)
    tc = B[:, 5]
    return dict(taxon=name, group=group, niche_C=niche, n=len(d),
                Trange=f"{d['T'].min():.0f}-{d['T'].max():.0f}", nT=d["T"].nunique(),
                E=base[0], ER=base[1], Topt=base[4], Topt_lo=q(4,2.5), Topt_hi=q(4,97.5),
                Tcarbon=base[5], Tc_lo=q(5,2.5), Tc_hi=q(5,97.5),
                margin=niche-base[5] if niche==niche else np.nan,
                P_below_niche=float(np.mean(tc[np.isfinite(tc)] < niche)) if niche==niche else np.nan,
                boot_ok=int(np.isfinite(tc).sum()), boot_n=len(B),
                note=note)

rows = []

# --- Candida (this study): per-well r and K from the pipeline -----------------
cd = pd.read_csv(f"{C}/results/tables/fit_coefficients_wide.csv")
on = pd.read_csv(f"{C}/results/tables/otu_names.csv")
cd = cd.merge(on[["OTU","group"]], on="OTU")
LAB = {"Clade1":"C. auris I","Clade2":"C. auris II","Clade3":"C. auris III","Clade4":"C. auris IV",
       "para":"C. parapsilosis","Hae":"C. haemulonii","Duo":"C. duobushaemulonii"}
for g, s in cd.groupby("group"):
    rows.append(analyse(LAB[g], s.rename(columns={"K":"K"})[["T","r","K"]], 37.0,
                        "mammal, opportunist"))

# --- E. coli ------------------------------------------------------------------
e = pd.read_csv(f"{P}/MicroAdapt/Experiments/Ecoli_extended_temp/oxygen_model_results_good_only.csv")
e = e.assign(T=e.Temperature, r=e.r_per_minute, K=e.resp_rate*e.N0)
rows.append(analyse("E. coli", e[["T","r","K"]], 37.0, "mammal, adapted", "6 temperatures, 6 C steps"))

# --- Pseudomonas (species unknown) -------------------------------------------
q = pd.read_csv(f"{P}/MicroAdapt/Experiments/Temperature_23_05/oxygen_model_results_good_only.csv")
q = q.assign(T=q.Temperature, r=q.r_per_minute, K=q.resp_rate*q.N0)
rows.append(analyse("Pseudomonas sp.", q[["T","r","K"]], np.nan, "unknown", "species not recorded"))

# --- Parsa OTUs (identities unknown) -----------------------------------------
pb = pd.read_csv(f"{P}/Parsa_bacteria/tables/fit_coefficients_wide.csv")
for o, s in pb.groupby("OTU"):
    rows.append(analyse(f"OTU{o}", s[["T","r","K"]], np.nan, "environmental?", "identity not recorded"))

# --- Zymoseptoria: two experiments, batch-corrected onto one scale ------------
def zymo(path, lab):
    d = pd.read_csv(path)
    w = d.pivot_table(index=["T","Condition","Replicate"], columns="parameter", values="Estimate").reset_index()
    w = w[w.Condition == "Control"]
    return w.assign(batch=lab).rename(columns={"resp_rate":"K"})[["T","r","K","batch"]]
z = pd.concat([zymo(f"{P}/Zymo/Zymo_3days_redone/o2_model_from_filtered/o2_model_coefficients_from_filtered.csv","A"),
               zymo(f"{P}/Zymo/Zymo_3_days_TPC/o2_model_outputs/o2_model_coefficients.csv","B")])
z = z[(z.r > 0) & (z.K > 0)]
ov = sorted(set(z[z.batch=="A"]["T"]) & set(z[z.batch=="B"]["T"]))
for col in ("r","K"):                                   # scale batch B onto batch A at shared temperatures
    fa = z[(z.batch=="A") & z["T"].isin(ov)].groupby("T")[col].median()
    fb = z[(z.batch=="B") & z["T"].isin(ov)].groupby("T")[col].median()
    f = np.exp(np.mean(np.log(fa/fb)))
    z.loc[z.batch=="B", col] *= f
rows.append(analyse("Zymoseptoria tritici", z[["T","r","K"]], 18.0, "plant, adapted",
                    "two experiments, batch-scaled on shared temperatures; niche = wheat canopy"))

out = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(out.round(2).to_string(index=False))
out.to_csv(f"{C}/model_extra/cross_taxon/cross_taxon_optima.csv", index=False)
