#!/usr/bin/env python3
"""proteostasis_burden.py -- does the damaged-protein burden carry a species signal?

THE HYPOTHESIS BEING TESTED
Route temperature through a proteostasis cost rather than through catalytic turnover:

    f_i(T) = 1 / (1 + exp((T - Tm_i)/s))          folded fraction of protein i
    D_s(T) = mean_i a_i * (1 - f_i(T))            damaged-protein burden

and charge D_s(T) as ATP maintenance (m0 + c*D) and/or lost proteome capacity (P0 - q*D),
with c, q and s SHARED across species so that every interspecies difference comes from
the sequence-predicted Tm distributions.

This differs from the sloped-maintenance and sloped-allocation rows already in the
mechanism table, which had no species specificity at all -- their required slopes came
out at 0.108-0.114 /C across the three relatives. Here the slope is derived per species.

WHY THIS SCRIPT COMES BEFORE THE MODEL
The whole idea rests on one quantity differing between species: D(40) against D(34). If
the burden at the restrictive temperature is the same in all four proteomes, then no
value of c or q can kill a relative at 40 C without killing C. auris too, and no amount
of model building changes that. So compute the burden first. It costs a minute; the model
costs months.

MEAN, NOT SUM
The proteomes differ in size (5,173 to 5,830 proteins; 662 to 997 in the metabolic
subset), so a raw sum would report annotation depth and reconstruction size as if it were
biology. Everything here is a mean over proteins, i.e. the expected unfolded fraction of
a randomly chosen protein.

WHAT WOULD SUPPORT THE HYPOTHESIS
  * a low tail that differs between species -- C. auris with materially fewer proteins
    predicted to unfold near 40 C than the relatives;
  * D(40) for the relatives well above C. auris, in the right direction for ALL THREE.
A ratio near 1, or a relative that falls BELOW C. auris, kills it at any c and q.

    python3 proteostasis_burden.py <dir of *_tm.csv> [--metabolic gem/thermal_tm.csv]
"""
import argparse, glob, os, sys
import numpy as np, pandas as pd
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

SPP = ['auris', 'haemulonii', 'duobushaemulonii', 'parapsilosis']
WIDTHS = [1.0, 2.0, 4.0, 8.79]        # 8.79 is the calibrated w in the paper's etcGEM


def burden(tm, T, s):
    return float(np.mean(1.0 / (1.0 + np.exp((tm - T) / s))))


def report(tag, byspp):
    print(f'\n{"="*78}\n{tag}\n{"="*78}')
    print(f'{"n proteins":<24}' + ''.join(f'{len(byspp[s]):>13}' for s in SPP))
    print(f'{"median Tm":<24}' + ''.join(f'{np.median(byspp[s]):>13.2f}' for s in SPP))
    print(f'{"5th percentile Tm":<24}' + ''.join(f'{np.quantile(byspp[s],.05):>13.2f}' for s in SPP))
    print(f'{"1st percentile Tm":<24}' + ''.join(f'{np.quantile(byspp[s],.01):>13.2f}' for s in SPP))
    print(f'{"min Tm":<24}' + ''.join(f'{np.min(byspp[s]):>13.2f}' for s in SPP))
    print('\nthe vulnerable tail -- fraction of the proteome predicted to melt below:')
    for lim in (38, 40, 42, 45, 48, 50):
        row = ''.join(f'{(byspp[s] < lim).mean():>13.4f}' for s in SPP)
        print(f'{f"  Tm < {lim} C":<24}{row}')

    for s in WIDTHS:
        d34 = {sp: burden(byspp[sp], 34.0, s) for sp in SPP}
        d40 = {sp: burden(byspp[sp], 40.0, s) for sp in SPP}
        base = d40['auris']
        print(f'\n--- unfolding width s = {s} C ---')
        print(f'{"":<24}' + ''.join(f'{sp[:12]:>13}' for sp in SPP))
        print(f'{"D(34)":<24}' + ''.join(f'{d34[sp]:>13.3e}' for sp in SPP))
        print(f'{"D(40)":<24}' + ''.join(f'{d40[sp]:>13.3e}' for sp in SPP))
        print(f'{"D(40) - D(34)":<24}' + ''.join(f'{d40[sp]-d34[sp]:>13.3e}' for sp in SPP))
        print(f'{"D(40) / auris":<24}' + ''.join(f'{d40[sp]/base:>13.3f}' for sp in SPP))
        rel = {sp: d40[sp] / base for sp in SPP[1:]}
        worst, best = min(rel.values()), max(rel.values())
        print(f'  species signal at 40 C: {best:.3f}x at most, {worst:.3f}x at least')
        if worst < 1.0:
            wrong = [sp for sp, v in rel.items() if v < 1.0]
            print(f'  WRONG DIRECTION for: {", ".join(wrong)} '
                  f'(burden BELOW C. auris, so no c or q can make it fail first)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir', help='directory of <species>_tm.csv files with id,pred_tm')
    ap.add_argument('--metabolic', help='gem/thermal_tm.csv, for the side-by-side')
    a = ap.parse_args()

    if a.metabolic and os.path.exists(a.metabolic):
        m = pd.read_csv(a.metabolic)
        m['sp'] = m.id.str.split('|').str[0]
        report('METABOLIC ENZYMES ONLY (the set the etcGEM actually uses)',
               {sp: m[m.sp == sp].pred_tm.values for sp in SPP})

    byspp = {}
    for sp in SPP:
        hits = glob.glob(os.path.join(a.dir, f'*{sp}*_tm.csv'))
        hits = [h for h in hits if 'pseudo' not in os.path.basename(h)]
        if not hits:
            sys.exit(f'no *{sp}*_tm.csv in {a.dir}')
        byspp[sp] = pd.concat([pd.read_csv(h) for h in sorted(hits)]).pred_tm.values
    report('COMPLETE PROTEOMES (the set the proteostasis hypothesis needs)', byspp)

    print(f'\n{"="*78}\nREADING\n{"="*78}')
    print('The burden mechanism needs D(40) to be substantially larger in all three')
    print('relatives than in C. auris. A ratio near 1.0, or any relative below 1.0,')
    print('means the mechanism cannot separate the species at any shared c or q, and')
    print('the model does not need to be built to know that.')


if __name__ == '__main__':
    main()
