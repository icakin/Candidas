#!/usr/bin/env python3
"""thermal_sensitivity.py — the model's thermal sensitivity against nature's.

Every earlier test asked how far a parameter must move to kill growth at ONE temperature.
That answer (32.5 C of Tm) is hard to compare with anything measured, because no
experiment reports "the Tm shift needed to abolish growth at 40 C".

This asks the question in a currency that IS measured. Define the model's thermal limit as
the temperature at which predicted mu crosses detection, and sweep an imposed uniform Tm
shift. The slope d(limit)/d(Tm) is how many degrees of growth limit the model buys per
degree of enzyme thermostability.

Nature's value for the same quantity is now known. Walunjkar et al. 2025 (Mol Biol Evol
42:msaf137) measured a mean ortholog Tm difference of 1.6 C between S. cerevisiae and
S. uvarum, whose growth thermal limits (IT50) differ by 8 C. That is 8/1.6 = 5.0 degrees of
thermal limit per degree of Tm.

If the model's slope is far below 5, the model is quantitatively too insensitive to enzyme
thermostability, by a factor that can be stated -- and stated against a measurement rather
than against a prediction. That is the same conclusion the counterfactual reaches, in a
unit a reader can check against the literature.
"""
import sys, json, statistics
from pathlib import Path
import cobra, numpy as np, pandas as pd
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

G = GEM
SP = {'auris': ('auris_iRV973_rekeyed.xml', 'medium_iRV973_auris.csv'),
      'haemulonii': ('haemulonii_draft.xml', 'medium_iRV973_auris.csv'),
      'duobushaemulonii': ('duobushaemulonii_draft.xml', 'medium_iRV973_auris.csv'),
      'parapsilosis': ('parapsilosis_iDC1003.xml', 'medium_iDC1003_parapsilosis.csv')}
MAINT = ('ATP_maintenance_NGAM__cyto', 'ATP_Maintenance__cyto')
THR = 0.05
NATURE_DLIMIT, NATURE_DTM = 8.0, 1.6      # Walunjkar et al. 2025

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
    for rid in MAINT:                                   # parapsilosis fix, as everywhere else
        if rid in m.reactions:
            r = m.reactions.get_by_id(rid); v = abs(r.upper_bound); r.bounds = (v, v)
    mw = {r['gene']: float(r['MW_kDa']) for _, r in pd.read_csv(TABLES / f'enzyme_mw_{sp}.csv').iterrows()}
    kc = pd.read_csv(TABLES / f'kcat_reaction_{sp}.csv')
    kmap = {r['reaction']: (float(r['kcat_per_s']), r['best_gene']) for _, r in kc.iterrows()}
    tov = {r.g: r.pred_topt for r in topt[topt.sp == sp].itertuples()}
    tmv = {r.g: r.pred_tm for r in tm[tm.sp == sp].itertuples()}
    mto = statistics.median(tov.values()); mtm = statistics.median(tmv.values())
    terms = []
    for r in m.reactions:
        if r.id in EX or 'iomass' in r.id or not r.genes: continue
        mws = [mw[g.id] for g in r.genes if g.id in mw]
        if not mws: continue
        ks, bg = kmap.get(r.id, (13.7, None))
        terms.append((r.forward_variable, r.reverse_variable,
                      statistics.mean(mws) / (ks * 3600.0), tov.get(bg, mto), tmv.get(bg, mtm)))
    con = m.problem.Constraint(0.0, ub=1.0, name='pp'); m.add_cons_vars(con); m.solver.update()
    return dict(m=m, terms=terms, con=con)


def mu(S, T, dTm=0.0, w=None):
    w = w or W
    c = {}
    for fv, rv, bc, To, Tm_ in S['terms']:
        a = np.exp(-((T - To) ** 2) / (2 * SIG * SIG)) / (1 + np.exp((T - (Tm_ + dTm)) / w))
        v = bc / max(a, 1e-6); c[fv] = v; c[rv] = v
    S['con'].ub = P0; S['con'].set_linear_coefficients(c)
    g = S['m'].slim_optimize()
    return (g if g == g else 0.0) * SCALE


def limit(S, dTm=0.0, w=None, lo=20.0, hi=90.0, tol=0.02):
    """Thermal limit: highest T at which predicted mu still reaches detection."""
    if mu(S, lo, dTm, w) < THR: return float('nan')
    if mu(S, hi, dTm, w) >= THR: return float('inf')
    a, b = lo, hi
    while b - a > tol:
        mid = (a + b) / 2
        if mu(S, mid, dTm, w) >= THR: a = mid
        else: b = mid
    return (a + b) / 2


print('Thermal-limit sensitivity to enzyme thermostability')
print(f'  limit = highest T with predicted mu >= {THR} h-1')
print(f'  nature: {NATURE_DLIMIT:.0f} C of growth limit per {NATURE_DTM:.1f} C of ortholog Tm')
print(f'          = {NATURE_DLIMIT / NATURE_DTM:.1f} C limit per C Tm '
      f'(S. cerevisiae vs S. uvarum, Walunjkar et al. 2025)\n')

SHIFTS = [0.0, -2.0, -5.0, -10.0, -20.0]
for wname, wval in (('as calibrated (w=8.79)', None), ('physical w=2.2', 2.2)):
    print(f'-- {wname}')
    print(f'{"species":<18}' + ''.join(f'{f"dTm{d:+.0f}":>10}' for d in SHIFTS) + f'{"slope":>12}')
    for sp in SP:
        S = build(sp)
        L = [limit(S, d, wval) for d in SHIFTS]
        good = [(d, x) for d, x in zip(SHIFTS, L) if np.isfinite(x)]
        if len(good) >= 2:
            xs = np.array([g[0] for g in good]); ys = np.array([g[1] for g in good])
            slope = float(np.polyfit(xs, ys, 1)[0])
        else:
            slope = float('nan')
        print(f'{sp:<18}' + ''.join(f'{x:>10.2f}' if np.isfinite(x) else f'{"inf":>10}'
                                    for x in L) + f'{slope:>11.3f} ')
    print()

S = build('auris')
for wname, wval in (('as calibrated', None), ('physical w=2.2', 2.2)):
    L0, L1 = limit(S, 0.0, wval), limit(S, -1.0, wval)
    sl = L0 - L1
    if np.isfinite(sl) and sl > 0:
        print(f'{wname:<16} model gives {sl:.3f} C of thermal limit per C of Tm; '
              f'nature gives {NATURE_DLIMIT / NATURE_DTM:.1f} '
              f'-> model is {(NATURE_DLIMIT / NATURE_DTM) / sl:.0f}x too insensitive')
