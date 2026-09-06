#!/usr/bin/env python3
"""17_build_measured_tpc.py — the measured thermal performance curve, with dead wells kept.

REPLACES the hand-built gem/measured_tpc_honest.csv, which was honest in the middle of
the range and not at the top.

The bug it fixes. A well whose culture is dead produces a flat oxygen trace, which cannot
be fitted, so it never reaches results/tables/derived_N0_R_results_with_carbon.csv. The
previous file was built from that fit table, so at temperatures where EVERY well was dead
it had no rows at all and the species curve simply stopped. C. duobushaemulonii ended at
38 C and C. parapsilosis at 40 C, which reads as "not tested" when in fact both were
tested to 44 C and everything died. The curve truncated exactly where the collapse it is
meant to show actually happened.

The denominator therefore comes from the RAW traces (Oxygen_Data_Filtered.csv), which is
the record of what was assayed: every species has wells at all twelve temperatures,
22-44 C. A well is scored

    mu = r * 60          if it has a valid fit AND consumed >= 2 mg/L of oxygen
    mu = 0               otherwise

so a dead or unfittable well is an observed zero, not a missing observation. r is stored
per minute (it enters the oxygen model as exp(r * T_end_min)), hence the factor 60; the
result is a specific growth rate in h^-1, involving no carbon-quota assumption, and is
directly comparable with the etcGEM's mu.

Writes gem/measured_tpc_honest.csv with species, T, mu_honest, n (wells assayed), and
n_alive (wells contributing a non-zero value).
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

TAB, G = RESULTS_TABLES, GEM
MIN_O2_DRAWDOWN = 2.0                      # mg/L, as in scripts/15_fig3.R

O = pd.read_csv(TAB / 'Oxygen_Data_Filtered.csv')
D = pd.read_csv(TAB / 'derived_N0_R_results_with_carbon.csv')
DD = pd.read_csv(TAB / 'well_drawdown_peak.csv')
for x in (D, DD):
    x['Replicate'] = x.Replicate.str.upper()
O['Replicate'] = O.Replicate.str.upper()

lut = D[['OTU', 'otu_name']].drop_duplicates().set_index('OTU').otu_name.to_dict()


def species_of(n):
    n = str(n)
    return ('auris' if n.startswith('Clade') else
            'haemulonii' if n.startswith('Hae') else
            'duobushaemulonii' if n.startswith('Duo') else
            'parapsilosis' if n.startswith('para') else 'other')


# every well that was actually run: the assay record, not the fit record
W = O[['T', 'OTU', 'Replicate']].drop_duplicates()
W['otu_name'] = W.OTU.map(lut)
W['sp'] = W.otu_name.map(species_of)
W = W[W.sp != 'other'].copy()

F = (D.merge(DD, on=['T', 'OTU', 'Replicate'], how='left')
       [['T', 'OTU', 'Replicate', 'r', 'fit_valid', 'drawdown']])
W = W.merge(F, on=['T', 'OTU', 'Replicate'], how='left')

alive = (W.fit_valid == True) & (W.drawdown >= MIN_O2_DRAWDOWN) & W.r.notna()
W['mu'] = np.where(alive, W.r * 60.0, 0.0)
W['alive'] = alive

out = (W.groupby(['sp', 'T'])
         .agg(mu_honest=('mu', 'median'), n=('mu', 'size'), n_alive=('alive', 'sum'))
         .reset_index().rename(columns={'sp': 'species'}))
out['mu_honest'] = out.mu_honest.round(4)
out = out.sort_values(['species', 'T'])
out.to_csv(TABLES / 'measured_tpc_honest.csv', index=False)

print(out.pivot(index='species', columns='T', values='mu_honest').round(2).to_string())
print('\nwells assayed / of which alive, at the top of the range:')
for sp in sorted(out.species.unique()):
    t = out[(out.species == sp) & (out['T'] >= 38)]
    print(f'  {sp:<18}' + '  '.join(f'{int(r.T)}: {int(r.n_alive)}/{int(r.n)}'
                                    for r in t.itertuples()))
print('\nwrote', TABLES / 'measured_tpc_honest.csv')
