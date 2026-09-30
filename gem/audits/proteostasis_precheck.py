#!/usr/bin/env python3
"""proteostasis_precheck.py -- is there ANY species signal in the damaged-protein burden?

The proposal is to route temperature through a proteostasis burden

    D_s(T) = sum_i a_i * [1 - f_i(T)],   f_i(T) = 1/(1+exp((T - Tm_i)/s))

and charge it as ATP maintenance and/or lost proteome capacity, with c, q and s SHARED
across species so that all interspecies difference comes from the Tm distributions.

That is a real proposal and it differs from the sloped-maintenance and sloped-allocation
rows already in the mechanism table, because those had no species specificity at all
(the required slopes were 0.108-0.114 /C across the three relatives). Here the slope is
derived per species from sequence.

But the whole idea rests on one quantity being different between species: D_s(40) versus
D_s(34). If the burden at the restrictive temperature is essentially the same in all four
proteomes, then no value of c or q can kill a relative at 40 C without killing C. auris
too, exactly as before -- and no amount of model building changes that. So compute it
first. It costs a minute; the model costs months.

Reported per species, at s = 1, 2, 4 and 8.79 (the calibrated width):
  D(34), D(40)          unweighted mean unfolded fraction, so model size cannot confound it
  D(40) - D(34)         the slope the mechanism would have to exploit
  ratio to C. auris     the species signal, which is what has to be large
  frac Tm < 40, < 45    the low tail, since that is where the burden comes from
"""
import numpy as np, pandas as pd
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from gempaths import *  # GEM, ROOT, INPUTS, MODELS, TABLES, EXTERNAL, RESULTS_TABLES, SP

d = pd.read_csv(TABLES / 'thermal_tm.csv')
d['sp'] = d.id.str.split('|').str[0]
SPP = ['auris', 'haemulonii', 'duobushaemulonii', 'parapsilosis']

def D(tm, T, s):
    # MEAN not sum: the four models carry different numbers of enzymes (704/679/662/997),
    # so a raw sum would report reconstruction size as if it were biology.
    return float(np.mean(1.0 / (1.0 + np.exp((tm - T) / s))))

print(f'{"n enzymes":<22}' + ''.join(f'{len(d[d.sp==s]):>12}' for s in SPP))
print(f'{"median Tm":<22}' + ''.join(f'{d[d.sp==s].pred_tm.median():>12.2f}' for s in SPP))
print(f'{"5th pct Tm":<22}' + ''.join(f'{d[d.sp==s].pred_tm.quantile(.05):>12.2f}' for s in SPP))
print(f'{"min Tm":<22}' + ''.join(f'{d[d.sp==s].pred_tm.min():>12.2f}' for s in SPP))
for lim in (40, 45, 48):
    print(f'{f"frac Tm < {lim} C":<22}'
          + ''.join(f'{(d[d.sp==s].pred_tm < lim).mean():>12.4f}' for s in SPP))

for s in (1.0, 2.0, 4.0, 8.79):
    print(f'\n--- unfolding width s = {s} C ---')
    print(f'{"":<22}' + ''.join(f'{sp[:11]:>12}' for sp in SPP))
    d34 = {sp: D(d[d.sp == sp].pred_tm.values, 34.0, s) for sp in SPP}
    d40 = {sp: D(d[d.sp == sp].pred_tm.values, 40.0, s) for sp in SPP}
    print(f'{"D(34)":<22}' + ''.join(f'{d34[sp]:>12.3e}' for sp in SPP))
    print(f'{"D(40)":<22}' + ''.join(f'{d40[sp]:>12.3e}' for sp in SPP))
    print(f'{"D(40) - D(34)":<22}' + ''.join(f'{d40[sp]-d34[sp]:>12.3e}' for sp in SPP))
    base = d40['auris']
    print(f'{"D(40) / auris":<22}' + ''.join(f'{d40[sp]/base:>12.3f}' for sp in SPP))
    rel = [d40[sp] / base for sp in SPP[1:]]
    print(f'  largest relative excess over C. auris at 40 C: {max(rel):.3f}x')
    print(f'  D(40)/D(34) within C. auris (the temperature signal): {d40["auris"]/max(d34["auris"],1e-12):.1f}x')
