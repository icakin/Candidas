#!/usr/bin/env python3
"""aa_pool_test.py — is the model's thermal blindness caused by its MEDIUM?

The diagnosis this follows from: at the growth optimum the etcGEM takes up ~1.5 mmol
O2/gDW/h against ~35 mmol C/gDW/h, giving a respiratory quotient of 6.7-7.3. A physical
RQ is 0.7-1.3. The model is not respiring. It is assembling biomass from the 17 free
amino acids the medium supplies at 3 mmol/gDW/h each -- 291 mmol C/gDW/h of pre-made
precursors against 60 from glucose. Protein synthesis from free amino acids needs
polymerisation ATP and little else: no carbon-skeleton biosynthesis, no TCA flux, no
oxidative phosphorylation.

That would explain every negative result at once. If growth does not run through
respiration, then nothing that happens to respiratory enzymes with temperature can reach
growth, and the enzyme-pool constraint is the only channel temperature has -- which is
exactly the near-scalar throttle already diagnosed. It would also explain why capping
oxygen changes nothing, and why the predicted carbon-use efficiency is flat: the carbon
budget is fixed by assimilation stoichiometry, and assimilation stoichiometry has no
temperature dependence.

The test: tighten the amino-acid pool towards zero, forcing the cell to build its own
precursors from glucose, and watch whether respiration, oxygen dependence and thermal
sensitivity come back. The quantity that matters is the last column -- the interspecies
Tm separation the model requires to reproduce the phenotype. If the medium is the cause,
that requirement should fall sharply as the pool closes.
"""
import sys, json, statistics
from pathlib import Path
import cobra, numpy as np, pandas as pd
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

G = GEM
SP = {'auris': ('auris_iRV973_rekeyed.xml', 'medium_iRV973_auris.csv'),
      'duobushaemulonii': ('duobushaemulonii_draft.xml', 'medium_iRV973_auris.csv')}
MAINT = 'ATP_maintenance_NGAM__cyto'
THR, PERM = 0.05, 0.7
calib = json.load(open(TABLES / 'etcgem_calib.json'))
SIG, W, P0, SCALE = calib['sig'], calib['w'], calib['P'], calib['SCALE']
topt = pd.read_csv(TABLES / 'thermal_topt.csv'); tm = pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (topt, tm):
    df['g'] = df['id'].str.split('|').str[1]; df['sp'] = df['id'].str.split('|').str[0]


def build(sp, aa):
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
            m.reactions.get_by_id(str(row['exchange_id'])).lower_bound = -aa      # <- swept
    rx = m.reactions.get_by_id(MAINT); ng = abs(rx.upper_bound); rx.bounds = (ng, ng)
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
    return dict(m=m, terms=terms, con=con, EX=[m.reactions.get_by_id(x) for x in EX])


def setT(S, T, dTm=0.0):
    c = {}
    for fv, rv, bc, To, Tm_ in S['terms']:
        a = np.exp(-((T - To) ** 2) / (2 * SIG * SIG)) / (1 + np.exp((T - (Tm_ + dTm)) / W))
        v = bc / max(a, 1e-6); c[fv] = v; c[rv] = v
    S['con'].ub = P0; S['con'].set_linear_coefficients(c)


def mu(S, T, dTm=0.0):
    setT(S, T, dTm); g = S['m'].slim_optimize()
    return (g if g == g else 0.0) * SCALE


def budget(S, T):
    setT(S, T); sol = S['m'].optimize()
    if sol.status != 'optimal': return None
    f = sol.fluxes; co2 = o2 = cin = 0.0
    for r in S['EX']:
        v = f[r.id]
        if abs(v) < 1e-12: continue
        forms = {(mt.formula or '').replace(' ', '') for mt in r.metabolites}
        if 'O2' in forms and any(mt.id.startswith('C00007') for mt in r.metabolites):
            if v < 0: o2 += -v
            continue
        n = sum(float((mt.elements or {}).get('C', 0) or 0) * abs(c) for mt, c in r.metabolites.items())
        if n == 0: continue
        if v < 0: cin += -v * n
        elif 'CO2' in forms: co2 += v * n
    return o2, co2, cin


def req_tm(S, lo=-60.0, tol=0.1):
    if mu(S, 40) < THR: return 0.0
    if mu(S, 40, lo) >= THR: return None
    a, b = lo, 0.0
    while b - a > tol:
        mid = (a + b) / 2
        if mu(S, 40, mid) < THR: a = mid
        else: b = mid
    return -(a + b) / 2


print('Sweeping the amino-acid pool bound (mmol/gDW/h per amino acid, 17 of them).')
print('3.0 is as shipped; 0 forces full de novo biosynthesis from glucose.\n')
print(f'{"AA pool":>8} | {"mu34":>6}{"mu40":>7}{"mu44":>7} | {"O2 34C":>7}{"RQ":>6} | '
      f'{"auris mu44":>10} | {"required Tm sep":>16}')
for aa in (3.0, 1.0, 0.3, 0.1, 0.03, 0.0):
    try:
        SA = build('auris', aa); SD = build('duobushaemulonii', aa)
    except Exception as e:
        print(f'{aa:>8} | build failed: {e}'); continue
    a34 = mu(SA, 34)
    if a34 < 1e-6:
        print(f'{aa:>8} | auris cannot grow at all -- model lacks the biosynthesis to close the gap')
        continue
    d34, d40, d44 = mu(SD, 34), mu(SD, 40), mu(SD, 44)
    bg = budget(SD, 34)
    o2, co2, cin = bg if bg else (float('nan'),) * 3
    rq = co2 / o2 if o2 > 1e-9 else float('inf')
    r = req_tm(SD)
    print(f'{aa:>8} | {d34:>6.3f}{d40:>7.3f}{d44:>7.3f} | {o2:>7.2f}{rq:>6.2f} | '
          f'{mu(SA, 44):>10.3f} | '
          f'{("unreachable" if r is None else f"{r:.1f} C"):>16}')
print('\n  (available separation from sequence: 0.53 C.  A physical RQ is 0.7-1.3.)')
