#!/usr/bin/env python3
"""allocation_and_trehalose.py — the two things the catalysis tests do not cover.

(1) PROTEOME ALLOCATION, TEMPERATURE-DEPENDENT.

The existing mechanism table has a "flat protein-pool capacity" row, and that row is what
made me tell the user that enzyme abundance was already covered. It is not. Flat means
temperature-independent, and a temperature-independent capacity cut fails for a trivial
reason: it removes as much capacity at 34 C as at 40 C, so permissive growth dies with the
restrictive. That is the same error made with the maintenance axis on the first attempt.

The mechanism actually proposed -- proteome diverted to heat-shock and repair functions as
temperature rises, leaving less for metabolism -- is SLOPED:

    P(T) = P0 * (1 - gamma * max(0, T - T0))

which is near-full capacity at permissive temperature and tightening only where the
phenotype appears. That is a genuinely different test and it has not been run.

The comparator is what is known about heat-shock proteome reallocation: chaperones and
stress proteins rise to roughly 10-30% of total protein under strong heat shock in yeast,
so a mechanism demanding that most of the metabolic proteome disappear is not physical.

(2) TREHALOSE. A tps2 deletion in C. auris fails to grow at 42 C (Nat Commun 2025), so
trehalose is a live candidate. But trehalose protects by acting as a chemical chaperone --
stabilising proteins and membranes physically -- and a genome-scale model represents a
metabolite only as a flux. There is no way to encode "this metabolite stabilises that
enzyme" in this formulation. What CAN be checked is whether the pathway is even present and
carries flux in each model, which bounds whether the difference could be metabolic capacity
rather than a physical protective role.
"""
import sys, json, statistics, re
from pathlib import Path
import cobra, numpy as np, pandas as pd
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

G = GEM
SP = {'auris': ('auris_iRV973_rekeyed.xml', 'medium_iRV973_auris.csv'),
      'haemulonii': ('haemulonii_draft.xml', 'medium_iRV973_auris.csv'),
      'duobushaemulonii': ('duobushaemulonii_draft.xml', 'medium_iRV973_auris.csv'),
      'parapsilosis': ('parapsilosis_iDC1003.xml', 'medium_iDC1003_parapsilosis.csv')}
MAINT = ('ATP_maintenance_NGAM__cyto', 'ATP_Maintenance__cyto')
RELS = ['haemulonii', 'duobushaemulonii', 'parapsilosis']
THR, PERM, T0 = 0.05, 0.7, 32.0

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
    for rid in MAINT:
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


def act(T, To, Tm_): return np.exp(-((T - To) ** 2) / (2 * SIG * SIG)) / (1 + np.exp((T - Tm_) / W))


def mu(S, T, gamma=0.0):
    c = {}
    for fv, rv, bc, To, Tm_ in S['terms']:
        v = bc / max(act(T, To, Tm_), 1e-6); c[fv] = v; c[rv] = v
    frac = max(1.0 - gamma * max(0.0, T - T0), 1e-6)
    S['con'].ub = P0 * frac
    S['con'].set_linear_coefficients(c)
    g = S['m'].slim_optimize()
    return (g if g == g else 0.0) * SCALE


def required_gamma(S, hi=0.125, tol=1e-4):
    """Smallest slope of proteome loss that drives 40 C growth below detection.
    hi is capped so that the pool is not driven negative at 40 C (1 - 8*gamma > 0)."""
    if mu(S, 40, 0.0) < THR: return 0.0
    if mu(S, 40, hi) >= THR: return None
    lo, up = 0.0, hi
    while up - lo > tol:
        mid = (lo + up) / 2
        if mu(S, 40, mid) < THR: up = mid
        else: lo = mid
    return (lo + up) / 2


print('(1) TEMPERATURE-SLOPED PROTEOME ALLOCATION')
print(f'    P(T) = P0 * (1 - gamma * max(0, T - {T0:.0f}))')
print(f'    40 C must fall below {THR}; 34 C must stay >= {PERM}x baseline\n')
print(f'{"species":<18}{"mu34":>8}{"mu40":>8}  {"required gamma":>15}  {"pool left at 40C":>17}'
      f'  {"mu34 there":>11}  {"floor":>7}  verdict')
res = {}
for sp in SP:
    S = build(sp)
    b34, b40 = mu(S, 34), mu(S, 40)
    g = required_gamma(S)
    if g is None:
        print(f'{sp:<18}{b34:>8.3f}{b40:>8.3f}  {"unreachable":>15}'
              f'  {"-":>17}  {"-":>11}  {"-":>7}  unreachable'); continue
    m34 = mu(S, 34, g); floor = PERM * b34
    left = 1.0 - g * (40 - T0)
    ok = m34 >= floor
    res[sp] = (g, left, ok)
    print(f'{sp:<18}{b34:>8.3f}{b40:>8.3f}  {g:>14.4f}/C  {left*100:>16.1f}%'
          f'  {m34:>11.3f}  {floor:>7.3f}  '
          f'{"BOTH SATISFIED" if ok else "permissive growth destroyed"}')
fin = [v for k, v in res.items() if k in RELS]
if fin:
    med_left = float(np.median([v[1] for v in fin]))
    print(f'\n  median metabolic proteome remaining at 40 C: {med_left*100:.1f}%')
    print(f'  i.e. {100-med_left*100:.0f}% of the metabolic proteome must be diverted by 40 C')
    print(f'  heat-shock proteins reach ~10-30% of total protein under strong heat shock')

print('\n\n(2) TREHALOSE PATHWAY IN THE MODELS')
PAT = re.compile(r'trehalose', re.I)
for sp in SP:
    xml, _ = SP[sp]
    m = cobra.io.read_sbml_model(str(MODELS /  xml))
    rx = [r for r in m.reactions
          if PAT.search(r.name or '') or any(PAT.search(mt.name or '') for mt in r.metabolites)]
    mets = {mt.id for r in rx for mt in r.metabolites if PAT.search(mt.name or '')}
    print(f'  {sp:<18} {len(rx)} trehalose-associated reactions, {len(mets)} trehalose metabolites')
    for r in rx[:6]:
        print(f'      {r.id:<22} {(r.name or "")[:58]}')
print('\n  Presence of the pathway is not the question the phenotype poses. Trehalose acts as')
print('  a chemical chaperone, and a flux model has no term for a metabolite stabilising a')
print('  protein, so this mechanism is outside what the etcGEM can express either way.')
