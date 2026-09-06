#!/usr/bin/env python3
"""ngam_falsification.py — the two-sided sufficiency test on the MAINTENANCE axis.

Exactly the structure of gem/etcgem_counterfactual.py, applied to the one mechanism that
is not already in the model and already dead. The etcGEM pins non-growth-associated
maintenance at a single constant (3.89 mmol ATP/gDW/h), identical for every species and
identical at every temperature: the model assumes that being at 44 C costs a cell no more
to stay alive than being at 22 C, and that C. auris and its relatives pay the same.

Question, two-sided as before: what interspecies difference in maintenance demand makes
the relative fail at 40 C (below detection 0.05) while its permissive 34 C growth stays
at or above 0.7x baseline?

The multiplier is applied at BOTH temperatures, not only at 40 C. Raising the cost only
where you want failure would trivially satisfy the permissive constraint and prove
nothing; a species that pays more maintenance at 40 C pays more at 34 C too. This is the
same discipline as the uniform Tm shift, which also applies at all temperatures.

The comparator is measured, not assumed: from the oxygen traces, median respiration per
unit growth at 40 C is 0.80 in C. auris against 1.07 in the relatives, a ratio of 1.34
(1.67 at 38 C). That is the interspecies separation that actually exists on this axis.

C. parapsilosis's maintenance reaction ships REVERSIBLE (lb = -3.9), letting the solver
synthesise ATP from ADP + Pi; it is corrected to irreversible here, matching the other
three models.
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
RELS = ['haemulonii', 'duobushaemulonii', 'parapsilosis']
THR, PERM_FRAC = 0.05, 0.7

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
    rx = m.reactions.get_by_id(MAINT[sp])
    base_ngam = abs(rx.upper_bound)
    rx.bounds = (base_ngam, base_ngam)           # irreversible, forced (fixes parapsilosis)
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
        ks, bg = kmap.get(r.id, (13.7, None))
        terms.append((r.forward_variable, r.reverse_variable,
                      statistics.mean(mws) / (ks * 3600.0), tov.get(bg, mto), tmv.get(bg, mtm)))
    con = m.problem.Constraint(0.0, ub=1.0, name='pp'); m.add_cons_vars(con); m.solver.update()
    return dict(m=m, terms=terms, con=con, rx=rx, ngam=base_ngam)


def act(T, To, Tm_): return np.exp(-((T - To) ** 2) / (2 * SIG * SIG)) / (1 + np.exp((T - Tm_) / W))


T0 = 32.0   # anchor: the temperature at which the measured maintenance burden bottoms out


def mu(S, T, beta=0.0):
    """Growth at T with a HEAT-ASSOCIATED maintenance demand.

        NGAM(T) = NGAM0 * (1 + beta * max(0, T - T0))

    A uniform multiplier is the wrong shape for a heat cost and gives a false negative:
    it charges the cell as much at 22 C as at 44 C, so any value large enough to kill
    40 C growth also destroys permissive growth, and the mechanism is rejected for a
    reason that has nothing to do with it. The measured burden is not uniform either --
    the relatives are CHEAPER than C. auris at 30-32 C (0.375 vs 0.447) and dearer at
    40 C (1.074 vs 0.801). What separates them is the slope, so the slope is what is
    tested. Measured beta: 0.080 per C in C. auris, 0.241 per C in the relatives.
    """
    k = 1.0 + beta * max(0.0, T - T0)
    S['rx'].bounds = (S['ngam'] * k, S['ngam'] * k)
    c = {}
    for fv, rv, bc, To, Tm_ in S['terms']:
        v = bc / max(act(T, To, Tm_), 1e-6); c[fv] = v; c[rv] = v
    S['con'].ub = P0; S['con'].set_linear_coefficients(c)
    g = S['m'].slim_optimize()
    return (g if g == g else 0.0) * SCALE


def required_beta(S, hi=50.0, tol=1e-3):
    """Smallest heat-cost slope that drives 40 C growth below detection."""
    if mu(S, 40, 0.0) < THR: return 0.0
    if mu(S, 40, hi) >= THR: return None
    lo, up = 0.0, hi
    while up - lo > tol:
        mid = (lo + up) / 2
        if mu(S, 40, mid) < THR: up = mid
        else: lo = mid
    return (lo + up) / 2


MEAS = {'auris': 0.080, 'relative': 0.241}     # measured slope, anchored at T0 = 32 C
print('Two-sided sufficiency test on the HEAT-ASSOCIATED MAINTENANCE axis')
print(f'  NGAM(T) = NGAM0 * (1 + beta * max(0, T - {T0:.0f}))')
print(f'  detection {THR} h-1 at 40 C; permissive 34 C growth must stay >= {PERM_FRAC}x baseline')
print(f'  measured beta: {MEAS["auris"]:.3f}/C in C. auris, {MEAS["relative"]:.3f}/C in the relatives'
      f'  (separation {MEAS["relative"]-MEAS["auris"]:.3f}/C, ratio {MEAS["relative"]/MEAS["auris"]:.1f}x)\n')
print(f'{"species":<18}{"mu34":>8}{"mu40":>8}  {"required beta":>14}  {"NGAM x at 40C":>14}'
      f'  {"mu34 there":>11}  {"floor":>7}  verdict')
out = {}
for sp in SP:
    S = build(sp)
    b34, b40 = mu(S, 34, 0.0), mu(S, 40, 0.0)
    beta = required_beta(S)
    if beta is None:
        print(f'{sp:<18}{b34:>8.3f}{b40:>8.3f}  {"none < 50/C":>14}'
              f'  {"-":>14}  {"-":>11}  {"-":>7}  unreachable')
        out[sp] = dict(baseline34=b34, baseline40=b40, required_beta=None); continue
    m34 = mu(S, 34, beta); floor = PERM_FRAC * b34
    ok = m34 >= floor
    print(f'{sp:<18}{b34:>8.3f}{b40:>8.3f}  {beta:>13.3f}/C  {1+beta*(40-T0):>13.1f}x'
          f'  {m34:>11.3f}  {floor:>7.3f}  {"BOTH SATISFIED" if ok else "permissive growth destroyed"}')
    out[sp] = dict(baseline34=b34, baseline40=b40, required_beta=beta, mu34_at_beta=m34,
                   perm_floor=floor, both_satisfied=bool(ok))

fin = [v['required_beta'] for s, v in out.items() if s in RELS and v.get('required_beta') is not None]
print()
if fin:
    md = float(np.median(fin))
    print(f'  median required beta across relatives : {md:.3f} per C')
    print(f'  measured beta in the relatives        : {MEAS["relative"]:.3f} per C')
    print(f'  fold gap                              : {md/MEAS["relative"]:.0f}x')
json.dump(out, open(TABLES / 'ngam_falsification.json', 'w'), indent=2)
print('\nwrote', TABLES / 'ngam_falsification.json')
