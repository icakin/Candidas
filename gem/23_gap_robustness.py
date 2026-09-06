#!/usr/bin/env python3
"""23_gap_robustness.py — does the model's thermal-limit gap survive correcting its defects?

THE CLAIM THE GRANT MAKES
The etcGEM's own predicted failure temperature (the T at which predicted growth falls
below the 0.05 h^-1 detection floor) sits well above where the cells actually die. That
gap is the motivation for the whole proposal. A reviewer will ask the obvious question:
the model has three documented defects, all of which INFLATE capacity, so is the gap just
an artefact of a too-generous model? This settles it by correcting all three and asking
whether the gap closes.

THE THREE DEFECTS, EACH INFLATES CAPACITY
  1. Respiratory quotient 6.7-7.3. The medium supplies 17 free amino acids at 3 mmol/gDW/h
     each, ~291 mmol pre-made precursor carbon against glucose's 60, so biomass is
     assembled rather than respired. CORRECTION: tighten the amino-acid pool bound.
  2. Energy-generating cycle. Guanylate cyclase + cGMP phosphodiesterase both reversible
     give net ADP+Pi -> ATP with no uptake. CORRECTION: block the reverse of the cycle
     (and C. parapsilosis maintenance is already forced irreversible in build()).
  3. Unfolding width w = 8.79 C implies a van 't Hoff enthalpy of ~101 kJ/mol against
     250-600 for cooperative globular proteins. CORRECTION: w = 2.2 C, the physical value,
     which lowers the thermal limit the most of the three.

Every correction can only LOWER the model's predicted limit. If the limit still sits well
above the observed death temperature with ALL THREE applied together, the gap is a
property of how temperature enters the model (a near-scalar capacity throttle), not of the
defects, and the grant's central premise stands.

OBSERVED LIMITS (this study): C. auris ~44-45 C (grows at 44, pilot no-growth at 45);
relatives ~40 C. The model limit is compared against these.

Run:  python3 23_gap_robustness.py [repo_root]     (needs cobra + the models; run on device)
"""
import sys, json, statistics
from pathlib import Path
import cobra, numpy as np, pandas as pd
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

G = GEM
THR = 0.05
SP = {'auris': ('auris_iRV973_rekeyed.xml', 'medium_iRV973_auris.csv'),
      'haemulonii': ('haemulonii_draft.xml', 'medium_iRV973_auris.csv'),
      'duobushaemulonii': ('duobushaemulonii_draft.xml', 'medium_iRV973_auris.csv'),
      'parapsilosis': ('parapsilosis_iDC1003.xml', 'medium_iDC1003_parapsilosis.csv')}
MAINT = {'auris': 'ATP_maintenance_NGAM__cyto', 'haemulonii': 'ATP_maintenance_NGAM__cyto',
         'duobushaemulonii': 'ATP_maintenance_NGAM__cyto', 'parapsilosis': 'ATP_Maintenance__cyto'}
OBS = {'auris': 45.0, 'haemulonii': 40.0, 'duobushaemulonii': 40.0, 'parapsilosis': 40.0}

calib = json.load(open(TABLES / 'etcgem_calib.json'))
SIG, W0, P0, SCALE = calib['sig'], calib['w'], calib['P'], calib['SCALE']
topt = pd.read_csv(TABLES / 'thermal_topt.csv'); tm = pd.read_csv(TABLES / 'thermal_tm.csv')
for df in (topt, tm):
    df['g'] = df['id'].str.split('|').str[1]; df['sp'] = df['id'].str.split('|').str[0]


def build(sp, aa_bound=3.0, block_egc=False):
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
    for _, row in med[med.setting == 'AA_POOL'].iterrows():         # DEFECT 1 knob
        if str(row['exchange_id']) in EX:
            m.reactions.get_by_id(str(row['exchange_id'])).lower_bound = -aa_bound
    for rid in MAINT.values():                                     # C.para maintenance fix
        if rid in m.reactions:
            r = m.reactions.get_by_id(rid); v = abs(r.upper_bound); r.bounds = (v, v)
    if block_egc:                                                  # DEFECT 2
        for rid in ('R00434', 'R01234'):
            for r in [x for x in m.reactions if rid in x.id]:
                if r.lower_bound < 0: r.lower_bound = 0.0
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


def mu(S, T, w, scale):
    c = {}
    for fv, rv, bc, To, Tm_ in S['terms']:
        a = np.exp(-((T - To) ** 2) / (2 * SIG * SIG)) / (1 + np.exp((T - Tm_) / w))
        v = bc / max(a, 1e-6); c[fv] = v; c[rv] = v
    S['con'].ub = P0; S['con'].set_linear_coefficients(c)
    g = S['m'].slim_optimize()
    return (g if g == g else 0.0) * scale


def calib_scale(S, w):
    """Refit the growth scale so C. auris peak matches ~0.80 h^-1 under this w.
    Without this, narrowing w would lower growth everywhere and the limit shift would
    conflate 'enzymes now unaffordable' with 'curve simply scaled down'."""
    peak = max(mu(S, T, w, 1.0) for T in np.arange(28, 40, 1.0))
    return (0.80 / peak) if peak > 0 else SCALE


def limit(S, w, scale, lo=20.0, hi=90.0, tol=0.02):
    if mu(S, lo, w, scale) < THR: return float('nan')
    if mu(S, hi, w, scale) >= THR: return float('inf')
    a, b = lo, hi
    while b - a > tol:
        mid = (a + b) / 2
        if mu(S, mid, w, scale) >= THR: a = mid
        else: b = mid
    return (a + b) / 2


SCENARIOS = [
    ('baseline (as published)',              dict(aa=3.0, w=W0,  egc=False, recal=False)),
    ('+ RQ fix (AA pool -3 -> -0.3)',        dict(aa=0.3, w=W0,  egc=False, recal=False)),
    ('+ EGC blocked',                        dict(aa=3.0, w=W0,  egc=True,  recal=False)),
    ('+ physical width w=2.2 (rescaled)',    dict(aa=3.0, w=2.2, egc=False, recal=True)),
    ('ALL THREE corrected',                  dict(aa=0.3, w=2.2, egc=True,  recal=True)),
]

print('MODEL THERMAL LIMIT (highest T with predicted mu >= 0.05 h^-1), per correction')
print('observed death temperatures: C. auris ~45, relatives ~40\n')
hdr = f'{"scenario":<34}' + ''.join(f'{sp[:11]:>13}' for sp in SP)
print(hdr); print('-' * len(hdr))
rows = []
for name, cfg in SCENARIOS:
    line = f'{name:<34}'
    rec = {'scenario': name}
    for sp in SP:
        S = build(sp, aa_bound=cfg['aa'], block_egc=cfg['egc'])
        scale = calib_scale(S, cfg['w']) if cfg['recal'] else SCALE
        L = limit(S, cfg['w'], scale)
        rec[sp] = L
        line += f'{("dead" if L!=L else f"{L:.1f}"):>13}'
    print(line); rows.append(rec)

print('\nGAP = model limit minus observed limit (positive = model too generous)')
base = rows[-1]
gl = f'{"ALL THREE corrected":<34}'
for sp in SP:
    L = base[sp]; gap = (L - OBS[sp]) if L == L else float('nan')
    gl += f'{("n/a" if gap!=gap else f"+{gap:.1f}"):>13}'
print(gl)
pd.DataFrame(rows).to_csv(TABLES / 'gap_robustness_results.csv', index=False)
print('\nwrote gem/gap_robustness_results.csv')
print('\nREADING: if the ALL-THREE row still sits well above the observed death')
print('temperatures, the thermal-limit gap is not an artefact of the model defects.')
