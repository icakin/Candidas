#!/usr/bin/env python3
"""carbon_budget_check.py — does the etcGEM reproduce the MEASURED carbon budget?

The gate check that has to pass before any maintenance-cost extension is worth fitting.

WHY THIS RATIO. Carbon use efficiency, CUE = C_biomass / (C_biomass + C_respired), is
dimensionless. It is therefore independent of SCALE, the constant in etcgem_calib.json
that maps model mu onto measured mu. So this compares the model's ENERGY BUDGET against
the data without borrowing anything from the growth-rate calibration — the model cannot
pass by having been tuned to the growth rates.

Measured side: results/tables/derived_N0_R_results_with_carbon.csv already carries CUE
and resp_over_growth per well, from the oxygen traces. Live wells only (drawdown >= 2.0
mg/L), the same dead-well filter as Fig 3 and gem/21_counts_at_44.py.

Model side: the machinery of gem/19_etcgem_counterfactual.py verbatim — same medium, same
enzyme-pool constraint, same activation function, same calibration — so what comes out is
a property of their etcGEM, not of a model rebuilt here. Reproducing their committed
baseline mu to 3 dp (auris 34/40/42/44 = 0.799/0.664/0.560/0.447) is the check on that.

TWO ACCOUNTING TRAPS, both of which silently give a wrong answer rather than an error:

  1. The biomass reaction consumes formula-less lumped pseudo-metabolites (e_Protein,
     e_Lipid, e_DNA, ...), so its carbon content CANNOT be read off its own stoichiometry
     — that route returns 0 and CUE comes out NaN. Carbon into biomass is obtained here
     by balance across the exchanges, C_bio = C_uptake - C_CO2 - C_organic_secreted,
     which is exact for a carbon-balanced model. Organic secretion is tracked separately
     (overflow_frac) rather than being folded into either term.

  2. CO2 must be identified by FORMULA, not by an id substring. These models use KEGG
     ids, where CO2 is EX_C00011__extr: a `'co2' in r.id.lower()` test matches nothing,
     scores respiration as exactly zero, and reports a flawless CUE of 1.000.

Usage:  python3 gem/carbon_budget_check.py [repo_root]
Writes: gem/model_carbon_budget.csv
"""
import sys, json, statistics
from pathlib import Path
import cobra, numpy as np, pandas as pd
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

G = GEM
SP = {'auris': ('auris_iRV973_rekeyed.xml', 'medium_iRV973_auris.csv'),
      'haemulonii': ('haemulonii_draft.xml', 'medium_iRV973_auris.csv'),
      'duobushaemulonii': ('duobushaemulonii_draft.xml', 'medium_iRV973_auris.csv'),
      'parapsilosis': ('parapsilosis_iDC1003.xml', 'medium_iDC1003_parapsilosis.csv')}
TEMPS = [22, 26, 30, 34, 36, 38, 40, 42, 44]

calib = json.load(open(TABLES / 'etcgem_calib.json'))
SIG, W, P0, SCALE = calib['sig'], calib['w'], calib['P'], calib['SCALE']
topt = pd.read_csv(TABLES / 'thermal_topt.csv'); tm = pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (topt, tm):
    df['g'] = df['id'].str.split('|').str[1]; df['sp'] = df['id'].str.split('|').str[0]


def therm_maps(sp):
    tov = {r.g: r.pred_topt for r in topt[topt.sp == sp].itertuples()}
    tmv = {r.g: r.pred_tm for r in tm[tm.sp == sp].itertuples()}
    return (tov, tmv,
            statistics.median(tov.values()) if tov else 35.0,
            statistics.median(tmv.values()) if tmv else 54.0)


def setup_pool(sp):                                   # identical to 19_etcgem_counterfactual.py
    xml, medf = SP[sp]
    m = cobra.io.read_sbml_model(str(MODELS /  xml)); med = pd.read_csv(INPUTS / medf)
    EX = {r.id for r in m.reactions if r.id.startswith(("EX_", "Drain")) or r.boundary}
    for r in m.reactions:
        if r.id in EX: r.lower_bound = 0.0
    for _, row in med.iterrows():
        rid = str(row['exchange_id'])
        if rid in EX:
            if row['setting'] == 'OPEN':  m.reactions.get_by_id(rid).lower_bound = -1000
            elif row['setting'] == 'FIT': m.reactions.get_by_id(rid).lower_bound = -10
    for _, row in med[med.setting == 'AA_POOL'].iterrows():
        if str(row['exchange_id']) in EX:
            m.reactions.get_by_id(str(row['exchange_id'])).lower_bound = -3
    mw = {r['gene']: float(r['MW_kDa']) for _, r in pd.read_csv(TABLES / f'enzyme_mw_{sp}.csv').iterrows()}
    kc = pd.read_csv(TABLES / f'kcat_reaction_{sp}.csv')
    kmap = {r['reaction']: (float(r['kcat_per_s']), r['best_gene']) for _, r in kc.iterrows()}
    return m, EX, mw, kmap


def precompute(m, EX, mw, kmap, tov, tmv, med_to, med_tm):
    terms = []
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws = [mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        kcat_s, bg = kmap.get(r.id, (13.7, None))
        terms.append((r.forward_variable, r.reverse_variable,
                      statistics.mean(mws) / (kcat_s * 3600.0),
                      tov.get(bg, med_to), tmv.get(bg, med_tm)))
    con = m.problem.Constraint(0.0, ub=1.0, name='pp'); m.add_cons_vars(con); m.solver.update()
    return terms, con


def act(T, Topt, Tm):
    return np.exp(-((T - Topt) ** 2) / (2.0 * SIG * SIG)) / (1.0 + np.exp((T - Tm) / W))


def nC(met):
    try:    return float((met.elements or {}).get('C', 0) or 0)
    except Exception: return 0.0


def budget(S, T):
    coeffs = {}
    for fv, rv, bc, Topt, Tm in S['terms']:
        c = bc / max(act(T, Topt, Tm), 1e-6); coeffs[fv] = c; coeffs[rv] = c
    S['con'].ub = P0; S['con'].set_linear_coefficients(coeffs)
    sol = S['m'].optimize()
    if sol.status != 'optimal' or sol.objective_value <= 1e-9:
        return None
    f = sol.fluxes
    C_co2 = C_in = C_org = O2_in = 0.0
    for r in S['exch']:
        v = f[r.id]
        if abs(v) < 1e-12: continue
        forms = {(mt.formula or '').replace(' ', '') for mt in r.metabolites}
        if 'O2' in forms and any(mt.id.startswith('C00007') for mt in r.metabolites):
            if v < 0: O2_in += -v
            continue
        n = sum(nC(mt) * abs(c) for mt, c in r.metabolites.items())
        if n == 0: continue
        if v < 0:                       C_in  += -v * n     # uptake
        elif 'CO2' in forms:            C_co2 +=  v * n     # respired
        else:                           C_org +=  v * n     # overflow secretion
    C_bio = C_in - C_co2 - C_org
    return dict(T=T, mu=sol.objective_value * SCALE,
                CUE=C_bio / (C_bio + C_co2) if C_bio + C_co2 > 0 else np.nan,
                resp_over_growth=C_co2 / C_bio if C_bio > 0 else np.nan,
                C_bio=C_bio, C_co2=C_co2, C_in=C_in, C_org=C_org, O2_in=O2_in,
                overflow_frac=C_org / C_in if C_in > 0 else np.nan)


rows = []
for sp in SP:
    tov, tmv, mto, mtm = therm_maps(sp)
    m, EX, mw, kmap = setup_pool(sp)
    terms, con = precompute(m, EX, mw, kmap, tov, tmv, mto, mtm)
    S = dict(m=m, terms=terms, con=con,
             exch=[r for r in m.reactions if r.id.startswith(('EX_', 'Drain')) or r.boundary])
    print(f'setup {sp}: {len(terms)} enzyme-costed reactions', flush=True)
    for T in TEMPS:
        b = budget(S, T)
        if b is None:
            print(f'  {sp} {T}C infeasible', flush=True); continue
        b['sp'] = sp; rows.append(b)
        print(f"  {sp:<18} {T}C  mu={b['mu']:.3f}  CUE={b['CUE']:.3f}"
              f"  resp/growth={b['resp_over_growth']:.3f}"
              f"  overflow={b['overflow_frac']:.3f}", flush=True)

D = pd.DataFrame(rows)
D.to_csv(G / 'model_carbon_budget.csv', index=False)
print('\nwrote', G / 'model_carbon_budget.csv')
print('\nmodel CUE across 34-44 C (span = how much thermal signal the model carries):')
print(D[D['T'] >= 34].groupby('sp').CUE.agg(
    lambda v: f'{v.min():.3f}-{v.max():.3f}   span {v.max() - v.min():.3f}').to_string())
print('\nCompare with the measured spans over the same window (live wells, median):')
print('  auris 0.323   haemulonii 0.364   duobushaemulonii 0.406   parapsilosis 0.311')
