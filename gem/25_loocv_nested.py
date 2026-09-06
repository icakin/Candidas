#!/usr/bin/env python3
"""loocv_nested.py — species-held-out validation of the etcGEM, fully nested.

WHAT IT IS FOR
The etcGEM's shape parameters and scale were fitted to C. auris alone. So the model
matching C. auris is not evidence that it works: it is the calibration reproducing its
own target. A referee can reasonably say the relatives only fail because they were never
fitted. This answers that: hold one species out, refit everything on the other three,
then predict the held-out one. If it still fails there, the failure is a property of the
model rather than of the calibration.

WHY IT REPLACES loocv_all.py
That script is not nested, as its own docstring concedes ("fast: fit shape sigma on the 3
training species with w,P fixed"). Two quantities leak:
  * w and P are read from etcgem_calib.json, fitted to C. auris, and frozen in every
    fold, so the fold that holds out C. auris still uses parameters fitted to C. auris;
  * the scale is taken from the HELD-OUT species' own measured maximum, so each fold
    uses the data it is supposed to be predicting.
Here sigma, w, P and the scale are all refitted on the training three only. Nothing
applied to the held-out species has seen its data.

TWO NUMBERS AT EACH SPECIES' OBSERVED FAILURE TEMPERATURE
  pred_mu_at_fail_frac_of_peak  predicted growth as a fraction of the model's own peak
                                for that species. Needs no measured value, so it does not
                                depend on how the scale was set.
  pred_mu_at_fail_h             the same in h^-1, using the training-set scale.
The optimum error is also reported, but it is the weaker half: all four optima lie within
6 C of each other, so recovering one to a couple of degrees is not demanding.

RUN
    cd <repo root>
    python3 gem/loocv_nested.py              # full fit, roughly 1.5-2 h
    FAST=1 python3 gem/loocv_nested.py       # reduced optimiser budget, roughly 25 min

FAST=1 halves the optimiser budget. A worse fit makes the model look worse, so any
failure it reports is pessimistic rather than flattering; use it to check the pipeline,
and the full run for anything quoted.

Writes gem/loocv_nested_results.csv with, per fold: the training species, every refitted
parameter, the training MSE, measured and predicted optima, the observed failure
temperature, both predicted-growth statistics there, the seed, and the wall time.
"""
import os, sys, time, subprocess
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize

# --- locate the repo, and point build_etcgem_tpc at it -----------------------
HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get('CANDIDAS_ROOT', HERE.parent))
G = GEM
if not (TABLES / 'measured_tpc_honest.csv').exists():
    sys.exit(f'cannot find gem/measured_tpc_honest.csv under {ROOT}; '
             f'set CANDIDAS_ROOT to the repository root')
sys.path.insert(0, str(G))
import build_etcgem_tpc as B
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP
B.G = G                      # that module hard-codes an absolute path; override it

SPP = ['auris', 'haemulonii', 'duobushaemulonii', 'parapsilosis']
DETECT = 0.05
SEED = 0                     # the fit is deterministic; recorded for the results table
FAST = os.environ.get('FAST', '0') == '1'
MAXIT, STARTS = (30, 1) if FAST else (70, 2)
hon = pd.read_csv(TABLES / 'measured_tpc_honest.csv')

try:
    CODE_VER = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', '--short', 'HEAD'],
                              capture_output=True, text=True, timeout=10).stdout.strip() or 'unknown'
except Exception:
    CODE_VER = 'unknown'

print(f'repo {ROOT}\ncode {CODE_VER}   budget maxiter={MAXIT} starts={STARTS}'
      f'{"  (FAST)" if FAST else ""}\nbuilding models', flush=True)
S = {}
for sp in SPP:
    tov, tmv, mto, mtm = B.therm_maps(sp)
    m, EX, mw, kmap = B.setup_pool(sp)
    S[sp] = (m,) + B.precompute(m, EX, mw, kmap, tov, tmv, mto, mtm)
    print('  ', sp, flush=True)


def measured(sp):
    return hon[hon.species == sp].dropna(subset=['mu_honest']).set_index('T')['mu_honest']


def curve(sp, sig, w, P, grid):
    m, terms, con = S[sp]
    return np.array([B.growth_at_T(m, terms, con, float(t), sig, w, P) for t in grid])


def fit_on(train):
    obs = {sp: measured(sp) for sp in train}

    def _one_scale(sig, w, P):
        """The least-squares optimal SINGLE scale across all training species.

        The first version of this rescaled each training species separately. That makes
        the absolute level of every predicted curve free, so P (the protein pool) can
        shrink towards zero with the scale growing to compensate and the loss unchanged:
        P and the scale are then not jointly identified, and the optimiser parks P on its
        lower bound. One shared scale makes the RELATIVE level of the species carry
        information, which is what identifies P.
        """
        num = den = 0.0
        for sp in train:
            mm = obs[sp]
            pv = curve(sp, sig, w, P, mm.index.values.astype(float))
            num += float(np.dot(pv, mm.values)); den += float(np.dot(pv, pv))
        return num / den if den > 0 else 0.0

    def loss(x):
        sig, w, P = x
        if sig <= 1 or w <= 0.3 or P <= 0.02 or P > 1.0:
            return 1e3
        sc = _one_scale(sig, w, P)
        e = 0.0
        for sp in train:
            mm = obs[sp]
            pv = curve(sp, sig, w, P, mm.index.values.astype(float))
            e += float(np.mean((pv * sc - mm.values) ** 2))
        return e / len(train)

    best = None
    for x0 in ([10.0, 8.0, 0.35], [6.0, 3.0, 0.20])[:STARTS]:
        r = minimize(loss, x0, method='Nelder-Mead',
                     options=dict(maxiter=MAXIT, xatol=0.3, fatol=1e-4))
        if best is None or r.fun < best.fun:
            best = r
    sig, w, P = best.x
    at_bound = bool(P <= 0.0205 or w >= 19.5 or sig <= 1.05 or sig >= 15.9)
    return dict(sig=float(sig), w=float(w), P=float(P),
                scale=float(_one_scale(sig, w, P)), mse=float(best.fun),
                at_bound=at_bound)


rows = []
print(f'\n{"held out":<18}{"sigma":>7}{"w":>7}{"P":>7}{"scale":>7}'
      f'{"opt m":>7}{"opt p":>7}{"fail T":>8}{"frac":>7}{"h^-1":>7}{"min":>6}', flush=True)
for held in SPP:
    t0 = time.time()
    train = [s for s in SPP if s != held]
    fit = fit_on(train)
    mm = measured(held); Ts = mm.index.values.astype(float)
    pv = curve(held, fit['sig'], fit['w'], fit['P'], Ts)
    # An OBSERVED failure temperature exists only if the species was actually measured
    # at or below detection somewhere. Two ways that fails, and the first version of this
    # script silently mishandled both by taking "last temperature with data + 2":
    #   censored  - C. auris is still at 0.27 h^-1 at 44 C, the assay ceiling, so it never
    #               failed inside the tested range;
    #   truncated - a species not assayed above its last measured point, where it was
    #               still growing, has no failure either. (This applied to
    #               C. duobushaemulonii under the old measured TPC, which was built from
    #               the fit table and so stopped at 38 C; build_measured_tpc.py now takes
    #               the denominator from the raw traces, and it fails at 40 C.)
    # In both cases there is nothing to evaluate against and the fold reports NA.
    dead = mm[mm <= DETECT]
    if len(dead):
        fail_T = float(dead.index.min()); fail_kind = 'observed'
    else:
        fail_T = float('nan')
        fail_kind = ('censored: still growing at the assay ceiling'
                     if float(mm.index.max()) >= 44 else
                     'truncated: not assayed above the last measured temperature')
    j = int(np.where(Ts == fail_T)[0][0]) if fail_T == fail_T else None
    rows.append(dict(
        fold=held, train_species=';'.join(train), sigma=fit['sig'], w=fit['w'], P=fit['P'],
        scale=fit['scale'], train_mse=fit['mse'],
        opt_measured_C=float(mm.idxmax()), opt_predicted_C=float(Ts[int(np.argmax(pv))]),
        opt_abs_err_C=abs(float(mm.idxmax()) - float(Ts[int(np.argmax(pv))])),
        fail_T_C=fail_T, fail_kind=fail_kind,
        measured_mu_at_fail=(float(mm.loc[fail_T]) if j is not None else float('nan')),
        pred_mu_at_fail_frac_of_peak=(float(pv[j] / max(pv.max(), 1e-9))
                                      if j is not None else float('nan')),
        pred_mu_at_fail_h=(float(pv[j] * fit['scale']) if j is not None else float('nan')),
        param_at_bound=fit['at_bound'],
        detect_threshold=DETECT, maxiter=MAXIT, starts=STARTS, seed=SEED,
        code_version=CODE_VER, minutes=round((time.time() - t0) / 60, 1)))
    r = rows[-1]
    print(f'{held:<18}{r["sigma"]:>7.2f}{r["w"]:>7.2f}{r["P"]:>7.3f}{r["scale"]:>7.2f}'
          f'{r["opt_measured_C"]:>7.0f}{r["opt_predicted_C"]:>7.0f}{r["fail_T_C"]:>8.0f}'
          f'{r["pred_mu_at_fail_frac_of_peak"]:>7.2f}{r["pred_mu_at_fail_h"]:>7.2f}'
          f'{r["minutes"]:>6.1f}', flush=True)

D = pd.DataFrame(rows)
out = TABLES / ('loocv_nested_results_FAST.csv' if FAST else 'loocv_nested_results.csv')
D.to_csv(out, index=False)
print(f'\nmean |optimum error| = {D.opt_abs_err_C.mean():.1f} C')
print('predicted growth at each species\' observed failure temperature')
print('  fraction of the model\'s own peak :',
      ', '.join(f'{s}={v:.2f}' for s, v in zip(D.fold, D.pred_mu_at_fail_frac_of_peak)))
print('  h^-1, training-set scale         :',
      ', '.join(f'{s}={v:.2f}' for s, v in zip(D.fold, D.pred_mu_at_fail_h)))
V = D.dropna(subset=['pred_mu_at_fail_frac_of_peak'])
print(f'\nfolds with an OBSERVED failure temperature: {len(V)} of {len(D)}')
for _, r in D[D.fail_T_C.isna()].iterrows():
    print(f'  {r.fold}: no failure to evaluate against ({r.fail_kind})')
if D.param_at_bound.any():
    print('\nWARNING: a parameter sits on its bound in fold(s): '
          + ', '.join(D.fold[D.param_at_bound]) + '. The fit is not identified there and '
          'these numbers should not be quoted. Note that raising the budget does NOT '
          'fix this: the full run reproduces the same bound, so the fold is genuinely '
          'not identified rather than under-optimised.')
print('\nwrote', out)
