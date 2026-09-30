#!/usr/bin/env python3
"""counts_at_44.json — observed isolate-level detection at 44 C for Fig 4 panel A.

Uses the SAME detection criterion as Fig 3 (scripts/15_fig3.R), so the two figures
cannot drift apart: a well is positive when it consumed >= MIN_O2_DRAWDOWN mg/L of O2
over the trace AND its curvature rt = r * fit-window length reaches CURV_MIN_RT. An
isolate is detected at a temperature if at least one of its wells is positive; a well
with no metabolic signal is a negative, not a missing observation.
"""
import json, os, numpy as np, pandas as pd
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

TAB = str(RESULTS_TABLES)
OUT = str(TABLES / 'counts_at_44.json')
MIN_O2_DRAWDOWN, CURV_MIN_RT, R2_MIN, N_MIN, TFOCUS = 2.0, 0.5, 0.90, 50, 44

L   = pd.read_csv(f'{TAB}/fit_coefficients_long.csv')
MET = pd.read_csv(f'{TAB}/fit_metrics.csv')
DER = pd.read_csv(f'{TAB}/derived_N0_R_results_with_carbon.csv')
DD  = pd.read_csv(f'{TAB}/well_drawdown_peak.csv')      # O2 peak minus trace minimum
lut = DER[['OTU', 'otu_name']].drop_duplicates().set_index('OTU').otu_name.to_dict()

r = L[L.parameter == 'r'].copy()
r['otu_name'] = r.OTU.map(lut); r['Replicate'] = r.Replicate.str.upper()
MET['Replicate'] = MET.Replicate.str.upper(); DD['Replicate'] = DD.Replicate.str.upper()
q = MET[['T', 'OTU', 'Replicate', 'r2', 'n', 'r_at_bound', 'K_at_bound']].copy()
q['ok_fit'] = MET.fit_valid
r = (r.merge(q, on=['T', 'OTU', 'Replicate'], how='left')
       .merge(DD, on=['T', 'OTU', 'Replicate'], how='left'))
r = r[(r.ok_fit == True) & (~r.r_at_bound.astype(bool)) & (~r.K_at_bound.astype(bool))
      & (r.r2 >= R2_MIN) & (r.n >= N_MIN)].copy()
r['rt'] = r.Estimate * r.T_end_min

def spof(nm):
    nm = str(nm)
    return ('auris' if nm.startswith('Clade') else 'haemulonii' if nm.startswith('Hae')
            else 'duobushaemulonii' if nm.startswith('Duo')
            else 'parapsilosis' if nm.startswith('para') else 'other')
r['sp'] = r.otu_name.map(spof); r = r[r.sp != 'other']
r['live'] = r.drawdown.fillna(0) >= MIN_O2_DRAWDOWN
r['rt_live'] = np.where(r.live, r.rt, 0.0)

iso_all = r.groupby(['sp', 'otu_name']).size().reset_index()          # every sampled isolate
f = r[r['T'] == TFOCUS].groupby(['sp', 'otu_name']).rt_live.max().reset_index()
counts = {}
for sp in ['auris', 'haemulonii', 'duobushaemulonii', 'parapsilosis']:
    n = int((iso_all.sp == sp).sum())                                  # denominator = all isolates
    b = f[f.sp == sp]
    k = int((b.rt_live >= CURV_MIN_RT).sum())
    counts[sp] = [k, n]
    print(f'{sp:<18} {k}/{n}')

json.dump(counts, open(OUT, 'w'), indent=2)
print('wrote', OUT)
