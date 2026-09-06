#!/usr/bin/env python3
"""more_mechanisms.py — four further mechanisms, beyond the five already falsified.

Already ruled out: uniform Tm, uniform Topt, flat protein pool, targeted sparse enzyme
sets, heat-associated maintenance. Each of those moves ONE axis. What follows tests the
gaps that leaves.

M1  JOINT Tm + Topt.  The Fig 4 caption currently concedes that the search "does not
    exhaustively enumerate ... arbitrarily large mixed Tm-plus-Topt changes". This closes
    that: a 2-D grid over both offsets at once, reporting the frontier of combinations
    that kill 40 C, and whether any of them spares permissive growth.

M2  WEAKEST LINK.  Every shift so far moves the whole proteome, and the targeted search
    picks enzymes by flux. But the enzyme-cost term is 1/a(T), so it is the LEAST stable
    enzymes that blow up first, whatever their flux. This shifts only the bottom decile
    of the Tm distribution -- a tail difference rather than a mean difference, which is
    also the more realistic way two proteomes actually differ.

M3  OXYGEN.  The model takes oxygen from an unbounded exchange (lb = -1e6). The assay it
    is being compared against measures oxygen drawdown in sealed wells, and O2 solubility
    falls ~14% between 34 and 44 C. So the model is given free access to the one resource
    the experiment shows is actually being exhausted. Capped here by Henry's law.

M4  OXYGEN + MEASURED MAINTENANCE, together.  The most promising remaining candidate, and
    the only one whose terms are all measured or physical rather than fitted. Heat raises
    the ATP demand (measured: beta = 0.080/C in C. auris, 0.241/C in the relatives) while
    lowering O2 availability (physical: Henry's law). Respiration is the ATP source, so
    the two squeeze from opposite sides and compound. A shared constraint that tightens
    with temperature is exactly what amplifies a small interspecies difference into a
    sharp boundary -- which no single-axis test can produce. One free parameter (the O2
    cap at reference); everything else is measured or physical.
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
O2EX = {'auris': 'EX_C00007__extr', 'haemulonii': 'EX_C00007__extr',
        'duobushaemulonii': 'EX_C00007__extr', 'parapsilosis': None}   # resolved by formula
RELS = ['haemulonii', 'duobushaemulonii', 'parapsilosis']
THR, PERM = 0.05, 0.7
BETA = {'auris': 0.080, 'haemulonii': 0.241, 'duobushaemulonii': 0.241, 'parapsilosis': 0.241}
T0 = 32.0

calib = json.load(open(TABLES / 'etcgem_calib.json'))
SIG, W, P0, SCALE = calib['sig'], calib['w'], calib['P'], calib['SCALE']
topt = pd.read_csv(TABLES / 'thermal_topt.csv'); tm = pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (topt, tm):
    df['g'] = df['id'].str.split('|').str[1]; df['sp'] = df['id'].str.split('|').str[0]

# dissolved O2 at air saturation, mg/L (standard table), normalised to 34 C
_T = np.array([20, 25, 30, 35, 40, 45]); _S = np.array([9.08, 8.26, 7.56, 6.95, 6.41, 5.93])
def o2_rel(T): return float(np.interp(T, _T, _S) / np.interp(34.0, _T, _S))


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
    rx = m.reactions.get_by_id(MAINT[sp]); ng = abs(rx.upper_bound); rx.bounds = (ng, ng)
    o2 = None
    for r in m.reactions:
        if r.id in EX and any((mt.formula or '') == 'O2' and mt.id.startswith('C00007')
                              for mt in r.metabolites):
            o2 = r; break
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
        terms.append([r.forward_variable, r.reverse_variable,
                      statistics.mean(mws) / (ks * 3600.0), tov.get(bg, mto), tmv.get(bg, mtm)])
    tms = np.array([t[4] for t in terms])
    con = m.problem.Constraint(0.0, ub=1.0, name='pp'); m.add_cons_vars(con); m.solver.update()
    return dict(m=m, terms=terms, con=con, rx=rx, ngam=ng, o2=o2,
                weak=np.percentile(tms, 10))


def mu(S, T, dTm=0.0, dTopt=0.0, weak_only=False, beta=None, o2cap=None):
    c = {}
    for fv, rv, bc, To, Tm_ in S['terms']:
        shift = dTm if (not weak_only or Tm_ <= S['weak']) else 0.0
        a = np.exp(-((T - (To + dTopt)) ** 2) / (2 * SIG * SIG)) / (1 + np.exp((T - (Tm_ + shift)) / W))
        v = bc / max(a, 1e-6); c[fv] = v; c[rv] = v
    S['con'].ub = P0; S['con'].set_linear_coefficients(c)
    k = 1.0 + beta * max(0.0, T - T0) if beta else 1.0
    S['rx'].bounds = (S['ngam'] * k, S['ngam'] * k)
    if S['o2'] is not None:
        S['o2'].lower_bound = -(o2cap * o2_rel(T)) if o2cap else -1000.0
    g = S['m'].slim_optimize()
    return (g if g == g else 0.0) * SCALE


def bisect(f, lo, tol=0.05):
    """smallest |offset| in [lo,0] driving f below THR; None if unreachable"""
    if f(0.0) < THR: return 0.0
    if f(lo) >= THR: return None
    a, b = lo, 0.0
    while b - a > tol:
        mid = (a + b) / 2
        if f(mid) < THR: a = mid
        else: b = mid
    return -(a + b) / 2


print('=' * 86)
print('M1  JOINT Tm + Topt  (closes the caption\'s stated limitation)')
print('=' * 86)
S = build('duobushaemulonii'); b34 = mu(S, 34); floor = PERM * b34
print(f'C. duobushaemulonii: baseline mu34={b34:.3f}, permissive floor={floor:.3f}, detection={THR}')
print(f'{"dTopt":>7} | minimum extra dTm needed to kill 40 C | mu34 there | permissive OK?')
best = None
for dto in [0, -2, -4, -6, -8, -10, -12, -14, -16]:
    r = bisect(lambda d: mu(S, 40, dTm=d, dTopt=dto), -40.0)
    if r is None:
        print(f'{dto:>7} | {"unreachable within 40 C":>37} |          - | -'); continue
    m34 = mu(S, 34, dTm=-r, dTopt=dto); ok = m34 >= floor
    if ok and (best is None or r + abs(dto) < best[0]): best = (r + abs(dto), r, dto)
    print(f'{dto:>7} | {r:>36.1f} | {m34:>10.3f} | {"YES" if ok else "no"}')
print('  -> no combination satisfies both sides.' if best is None else
      f'  -> feasible at dTopt={best[2]}, dTm={best[1]:.1f} (total {best[0]:.1f} C)')

print()
print('=' * 86)
print('M2  WEAKEST LINK  (shift only the least-stable decile of each proteome)')
print('=' * 86)
print(f'{"species":<18}{"Tm 10th pct":>12}   {"required shift, bottom decile only":>36}')
for sp in RELS:
    S = build(sp)
    r = bisect(lambda d: mu(S, 40, dTm=d, weak_only=True), -60.0)
    print(f'{sp:<18}{S["weak"]:>12.2f}   '
          f'{("unreachable within 60 C" if r is None else f"{r:.1f} C"):>36}')

print()
print('=' * 86)
print('M3 / M4  OXYGEN CAP, alone and with the MEASURED maintenance slope')
print('=' * 86)
print(f'  O2 availability by Henry\'s law, relative to 34 C: '
      f'40 C {o2_rel(40):.3f}, 44 C {o2_rel(44):.3f}')
print(f'  measured maintenance slope beta: auris {BETA["auris"]}/C, relatives {BETA["haemulonii"]}/C\n')
ST = {sp: build(sp) for sp in SP}
unc = {sp: mu(ST[sp], 34) for sp in SP}
print(f'{"O2 cap":>7} | ' + ' | '.join(f'{s[:12]:>12}' for s in SP) + ' |  verdict')
print(f'{"(mmol)":>7} | ' + ' | '.join(f'{"mu40 / mu44":>12}' for s in SP) + ' |')
for cap in [None, 20, 10, 6, 4, 3, 2.5, 2, 1.5, 1.0]:
    cells, ok34 = [], True
    a44 = r40 = None
    for sp in SP:
        S = ST[sp]; b = BETA[sp]
        m34 = mu(S, 34, beta=b, o2cap=cap)
        m40 = mu(S, 40, beta=b, o2cap=cap); m44 = mu(S, 44, beta=b, o2cap=cap)
        cells.append(f'{m40:>5.2f} /{m44:>5.2f}')
        if m34 < PERM * unc[sp]: ok34 = False
        if sp == 'auris': a44 = m44
        else: r40 = max(r40 or 0, m40)
    good = (a44 is not None and a44 >= THR) and (r40 is not None and r40 < THR) and ok34
    v = 'PHENOTYPE REPRODUCED' if good else ('auris dies too' if (a44 or 0) < THR
         else ('permissive lost' if not ok34 else 'relatives still grow'))
    print(f'{str(cap):>7} | ' + ' | '.join(cells) + f' |  {v}')
print('\n  target: auris mu44 >= 0.05, every relative mu40 < 0.05, all mu34 >= 0.7x unconstrained')
