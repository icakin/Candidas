#!/usr/bin/env python3
"""26_fig4_tables.py -- every number Figure 4 draws, as flat tables for scripts/19_fig4.R.

The figure is drawn in R with the rest of the manuscript figures. R has no JSON reader in
base, and the paired ortholog statistic needs the kcat tables and a t-test, so this step
turns the etcGEM outputs into four small CSVs that the R script reads verbatim. Nothing is
computed in the R script that is not in these tables; nothing is computed here that is not
reported in gem/FIG4_LOCKED.md.

Reads   gem/tables/counterfactual_results.json      19_etcgem_counterfactual.py
        gem/tables/dyn_sparse_{Tm,Topt}.json        20_dyn_sparse.py
        gem/tables/thermal_{tm,topt}.csv            14_run_seq2topt_seq2tm.sh
        gem/tables/kcat_reaction_<species>.csv      11_aggregate_kcat.py
Writes  gem/tables/fig4/requirements.csv    param, threshold, median_required (per detection
                                            threshold; None where no shift reaches failure)
        gem/tables/fig4/paired.csv          param, relative, n, mean, se, lo, hi, p
                                            (mean paired difference auris - relative over
                                            UNIQUE ortholog pairs; see 24_paired_dedup_audit.py)
        gem/tables/fig4/tm_predicted.csv    species, gene, pred_tm       (panel B densities)
        gem/tables/fig4/dyn_sparse.csv      param, relative, n, rid, mu40, mu34  (panel D)
Panel A reads gem/tables/measured_tpc_honest.csv and etcgem_tpc_pred_fine.csv directly.

Run:  python3 gem/26_fig4_tables.py   then   Rscript scripts/19_fig4.R
"""
import json
import numpy as np, pandas as pd
from scipy import stats
from gempaths import *  # TABLES, SP

SP4 = list(SP); RELS = SP4[1:]
OUT = TABLES / 'fig4'; OUT.mkdir(exist_ok=True)

def paired(param):
    f, col = (('thermal_tm.csv', 'pred_tm') if param == 'Tm' else ('thermal_topt.csv', 'pred_topt'))
    d = pd.read_csv(TABLES / f)
    d['sp'] = d.id.str.split('|').str[0]; d['g'] = d.id.str.split('|').str[1]
    v = {s: dict(zip(d[d.sp == s].g, d[d.sp == s][col])) for s in SP4}
    kc = {s: pd.read_csv(TABLES / f'kcat_reaction_{s}.csv')[['reaction', 'best_gene']].dropna() for s in SP4}
    rows = []
    for r in RELS:
        m = kc['auris'].merge(kc[r], on='reaction', suffixes=('_a', '_b'))
        uq = {(x, y): v['auris'][x] - v[r][y] for x, y in zip(m.best_gene_a, m.best_gene_b)
              if x in v['auris'] and y in v[r]}
        dd = np.array(list(uq.values()))
        t, p = stats.ttest_1samp(dd, 0)
        se = dd.std(ddof=1) / np.sqrt(len(dd))
        rows.append(dict(param=param, relative=r, n=len(dd), mean=dd.mean(), se=se,
                         lo=dd.mean() - 1.96 * se, hi=dd.mean() + 1.96 * se, p=p))
    return rows

def main():
    R = json.load(open(TABLES / 'counterfactual_results.json'))
    req = [dict(param=k, threshold=float(t.split('=')[1]), median_required=R['mechanisms'][t][k]['median_required'])
           for t in R['mechanisms'] for k in ('Tm', 'Topt')]
    pd.DataFrame(req).to_csv(OUT / 'requirements.csv', index=False)

    pd.DataFrame(paired('Tm') + paired('Topt')).to_csv(OUT / 'paired.csv', index=False, float_format='%.6g')

    tm = pd.read_csv(TABLES / 'thermal_tm.csv')
    pd.DataFrame({'species': tm.id.str.split('|').str[0], 'gene': tm.id.str.split('|').str[1],
                  'pred_tm': tm.pred_tm}).to_csv(OUT / 'tm_predicted.csv', index=False)

    rows = []
    for k in ('Tm', 'Topt'):
        D = json.load(open(TABLES / f'dyn_sparse_{k}.json'))
        for rel in RELS:
            for t in D['relatives'][rel]['greedy']['trajectory']:
                rows.append(dict(param=k, relative=rel, n=t['n'], rid=t['rid'], mu40=t['mu40'], mu34=t['mu34']))
    pd.DataFrame(rows).to_csv(OUT / 'dyn_sparse.csv', index=False)

    p = pd.read_csv(OUT / 'paired.csv'); r5 = {x['param']: x['median_required'] for x in req if x['threshold'] == 0.05}
    for k in ('Tm', 'Topt'):
        b = p[p.param == k].sort_values('mean').iloc[-1]
        print(f'{k:4s} required {r5[k]:6.2f} C   largest paired {b["mean"]:.3f} C ({b.relative}, n={b.n}, p={b.p:.1e})'
              f'   fold {r5[k] / b["mean"]:.0f}x')
    print('wrote', OUT)

if __name__ == '__main__':
    main()
