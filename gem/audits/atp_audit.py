#!/usr/bin/env python3
"""atp_audit.py — energy-conservation audit of the four etcGEMs.

Motivated by a bound discrepancy: three models pin ATP maintenance irreversibly at
3.89 mmol/gDW/h, while the C. parapsilosis model has it REVERSIBLE at lb = -3.9, so the
solver may run it backwards and synthesise ATP from ADP + Pi.

Three tests:
  1. Is the reverse direction actually used at the optimum? (flux < 0 = free ATP taken)
  2. Energy-generating cycle test: with every uptake closed, can the model still make ATP?
     A correct model returns exactly 0. Anything above 0 is a thermodynamically
     infeasible cycle and inflates every growth prediction that model has ever made.
  3. What happens to baseline mu when C. parapsilosis's maintenance is made irreversible
     and forced, matching the other three?

Run:  python3 atp_audit.py [repo_root]
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
MAINT = {'auris': 'ATP_maintenance_NGAM__cyto', 'haemulonii': 'ATP_maintenance_NGAM__cyto',
         'duobushaemulonii': 'ATP_maintenance_NGAM__cyto', 'parapsilosis': 'ATP_Maintenance__cyto'}

calib = json.load(open(TABLES / 'etcgem_calib.json'))
SIG, W, P0, SCALE = calib['sig'], calib['w'], calib['P'], calib['SCALE']
topt = pd.read_csv(TABLES / 'thermal_topt.csv'); tm = pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (topt, tm):
    df['g'] = df['id'].str.split('|').str[1]; df['sp'] = df['id'].str.split('|').str[0]


def build(sp):
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
    tov = {r.g: r.pred_topt for r in topt[topt.sp == sp].itertuples()}
    tmv = {r.g: r.pred_tm for r in tm[tm.sp == sp].itertuples()}
    mto = statistics.median(tov.values()) if tov else 35.0
    mtm = statistics.median(tmv.values()) if tmv else 54.0
    terms = []
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws = [mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        kcat_s, bg = kmap.get(r.id, (13.7, None))
        terms.append((r.forward_variable, r.reverse_variable,
                      statistics.mean(mws) / (kcat_s * 3600.0),
                      tov.get(bg, mto), tmv.get(bg, mtm)))
    con = m.problem.Constraint(0.0, ub=1.0, name='pp'); m.add_cons_vars(con); m.solver.update()
    return dict(m=m, terms=terms, con=con, EX=EX)


def act(T, To, Tm_): return np.exp(-((T - To) ** 2) / (2 * SIG * SIG)) / (1 + np.exp((T - Tm_) / W))


def set_T(S, T):
    c = {}
    for fv, rv, bc, To, Tm_ in S['terms']:
        v = bc / max(act(T, To, Tm_), 1e-6); c[fv] = v; c[rv] = v
    S['con'].ub = P0; S['con'].set_linear_coefficients(c)


print('=' * 78)
print('TEST 1 & 3 — maintenance flux at the optimum, and the effect of fixing it')
print('=' * 78)
TS = [30, 34, 40, 42, 44]
rows = []
for sp in SP:
    S = build(sp); m = S['m']; mid = MAINT[sp]
    rx = m.reactions.get_by_id(mid)
    lb0, ub0 = rx.lower_bound, rx.upper_bound
    print(f'\n{sp}   {mid}  bounds as shipped: [{lb0:g}, {ub0:g}]'
          f'{"   <-- REVERSIBLE" if lb0 < 0 else ""}')
    for T in TS:
        set_T(S, T)
        sol = m.optimize()
        mu0 = sol.objective_value * SCALE
        vmaint = sol.fluxes[mid]
        # counterfactual: force it forward at the same magnitude the others use
        rx.lower_bound, rx.upper_bound = abs(ub0), abs(ub0)
        set_T(S, T); s2 = m.optimize()
        mu1 = (s2.objective_value * SCALE) if s2.status == 'optimal' else float('nan')
        rx.lower_bound, rx.upper_bound = lb0, ub0
        rows.append(dict(sp=sp, T=T, maint_flux=vmaint, mu_shipped=mu0, mu_forced=mu1,
                         d_pct=100 * (mu1 - mu0) / mu0 if mu0 else np.nan))
        print(f'   T={T}  maintenance flux={vmaint:+8.3f}   mu shipped={mu0:.4f}'
              f'   mu with maintenance forced forward={mu1:.4f}'
              f'   ({100*(mu1-mu0)/mu0:+.1f}%)')

print()
print('=' * 78)
print('TEST 2 — energy-generating cycle test (all uptake closed, maximise ATP hydrolysis)')
print('        a correct model returns 0.000; anything else is free energy')
print('=' * 78)
for sp in SP:
    S = build(sp); m = S['m']; mid = MAINT[sp]
    with m:
        for rid in S['EX']:
            m.reactions.get_by_id(rid).lower_bound = 0.0     # nothing may enter
        rx = m.reactions.get_by_id(mid)
        rx.lower_bound, rx.upper_bound = 0.0, 1000.0
        m.objective = rx
        set_T(S, 34)
        v = m.slim_optimize()
    flag = 'OK' if (v is None or not (v > 1e-6)) else '<-- FREE ATP'
    print(f'  {sp:<18} max ATP hydrolysis with no uptake = {0.0 if v!=v else v:8.4f}   {flag}')

pd.DataFrame(rows).to_csv(TABLES / 'atp_audit.csv', index=False)
print('\nwrote', TABLES / 'atp_audit.csv')
