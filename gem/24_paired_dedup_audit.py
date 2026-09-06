#!/usr/bin/env python3
"""paired_dedup_audit.py -- is the paired thermal-parameter difference pseudo-replicated?

Figure 4B/C pair enzymes by REACTION: 19_fig4_etcgem.py merges the two species' kcat
tables on 'reaction' and reads each side's best_gene. A gene catalysing several reactions
therefore contributes several rows, so n counts reactions, not protein pairs. With only
704 C. auris enzymes carrying a prediction, the reported n = 1041 cannot be one row per
pair.

Two consequences, and the second is the one that matters:
  * the standard error is too small, because independent observations are over-counted;
  * the MEAN is weighted by promiscuity. A gene appearing under 23 reactions counts 23
    times, and highly connected enzymes are central-metabolism enzymes, not a random
    sample of the proteome. That is a bias, not just a precision problem.

Reported per relative and per parameter: the published reaction-level statistic, the same
after collapsing to unique (auris gene, relative gene) pairs, and a cluster bootstrap that
resamples PAIRS rather than rows. The unique-pair version is the estimator to quote; the
bootstrap is there to show the interval is not an artefact of how ties were collapsed.
"""
import numpy as np, pandas as pd
from scipy import stats as st
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

SP4 = ['auris', 'haemulonii', 'duobushaemulonii', 'parapsilosis']
RELS = SP4[1:]
PAR = {'Tm': ('thermal_tm.csv', 'pred_tm'), 'Topt': ('thermal_topt.csv', 'pred_topt')}
REQ = {'Tm': 32.54, 'Topt': 15.62}          # model-required separation, from the sweep
RNG = np.random.default_rng(0)

kc = {s: pd.read_csv(TABLES / f'kcat_reaction_{s}.csv')[['reaction', 'best_gene']].dropna()
      for s in SP4}


def load(param):
    f, col = PAR[param]
    d = pd.read_csv(TABLES / f)
    d['sp'] = d.id.str.split('|').str[0]; d['g'] = d.id.str.split('|').str[1]
    return {s: dict(zip(d[d.sp == s].g, d[d.sp == s][col])) for s in SP4}


def boot(pairs, vals, n=20000):
    u = sorted(set(pairs)); idx = {p: [] for p in u}
    for p, v in zip(pairs, vals): idx[p].append(v)
    cell = {p: float(np.mean(idx[p])) for p in u}
    m = np.array([np.mean([cell[u[j]] for j in RNG.integers(0, len(u), len(u))])
                  for _ in range(n)])
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def ci(v):
    se = v.std(ddof=1) / np.sqrt(len(v))
    return v.mean() - 1.96 * se, v.mean() + 1.96 * se


summary = {}
for param in PAR:
    V = load(param)
    print(f'\n{"="*76}\n{param}   (model requires {REQ[param]:.2f} C of separation)\n{"="*76}')
    print(f'{"":<22}{"n":>7}{"mean":>8}{"95% CI":>18}{"p":>11}{"fold gap":>10}')
    best = {}
    for r in RELS:
        m = kc['auris'].merge(kc[r], on='reaction', suffixes=('_a', '_b'))
        rows = [(x, y, V['auris'][x] - V[r][y])
                for x, y in zip(m.best_gene_a, m.best_gene_b)
                if x in V['auris'] and y in V[r]]
        v = np.array([z for _, _, z in rows])
        keys = [(x, y) for x, y, _ in rows]
        uq = {}
        for k, z in zip(keys, v): uq.setdefault(k, []).append(z)
        u = np.array([np.mean(z) for z in uq.values()])

        print(f'\n{r}')
        for tag, arr in (('  reaction-level', v), ('  unique pairs', u)):
            lo, hi = ci(arr)
            print(f'{tag:<22}{len(arr):>7}{arr.mean():>8.3f}{f"[{lo:.2f}, {hi:.2f}]":>18}'
                  f'{st.ttest_1samp(arr,0)[1]:>11.1e}'
                  f'{(REQ[param]/arr.mean() if arr.mean()>0 else float("nan")):>10.0f}')
        blo, bhi = boot(keys, v)
        print(f'{"  cluster bootstrap":<22}{len(uq):>7}{v.mean():>8.3f}'
              f'{f"[{blo:.2f}, {bhi:.2f}]":>18}')
        ga = pd.Series([x for x, _, _ in rows]).value_counts()
        print(f'    {len(v)} rows from {ga.size} distinct C. auris genes, '
              f'{len(uq)} unique pairs, max {ga.max()} reactions per gene')
        best[r] = (float(u.mean()), len(u), ci(u), float(st.ttest_1samp(u, 0)[1]))

    top = max(best, key=lambda r: best[r][0])
    mean, n, (lo, hi), p = best[top]
    summary[param] = (top, mean, n, lo, hi, p, REQ[param] / mean)
    print(f'\n  largest deduplicated difference: {top}, {mean:.3f} C '
          f'[{lo:.2f}, {hi:.2f}], n = {n}, p = {p:.1e}  ->  {REQ[param]/mean:.0f}x')

print(f'\n{"="*76}\nFIGURE VALUES TO USE\n{"="*76}')
for param, (top, mean, n, lo, hi, p, fold) in summary.items():
    print(f'{param:<6} {mean:.2f} C [{lo:.2f}, {hi:.2f}]  n = {n}  p = {p:.1e}  '
          f'fold gap {fold:.0f}x   ({top})')
