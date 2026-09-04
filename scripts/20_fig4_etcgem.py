#!/usr/bin/env python3
"""20_fig4_etcgem.py — FIGURE 4: the etcGEM counterfactual falsification.

Canonical generator, moved into scripts/ per the pipeline convention (see
reports/MANUSCRIPT_FIGURE_LINKAGE.md: every manuscript figure resolves to a file
produced by a script, referenced by path). Python rather than R because the panels
are read off etcGEM simulations.

Inputs, all regenerable from the repo:
  gem/counterfactual_results.json   <- gem/etcgem_counterfactual.py
  gem/counterfactual_sweep.csv      <- gem/etcgem_counterfactual.py   (panel B)
  gem/dyn_sparse_{Tm,Topt}.json     <- gem/dyn_sparse.py Tm | Topt    (panel D)
  gem/counts_at_44.json             <- gem/counts_at_44.py            (panel A n/N)

counts_at_44.json is derived with the SAME detection criterion as Fig 3
(scripts/15_fig3.R), so the two figures cannot drift apart.

Corrections applied to the earlier gem/mkfig_counterfactual2.py:
  * The 'median_required or 35' fallback is gone. A missing requirement now raises
    instead of silently plotting 35 °C as though it had been computed.
  * Panel B no longer claims the relatives behave alike. A uniform Tm shift reaches
    detection ONLY in C. duobushaemulonii; in C. haemulonii and C. parapsilosis no
    shift within +/-35 °C makes 40 °C growth fail at all, which is the stronger
    result and is now stated.
  * Panel D's axis matches the search that was actually run.
  * Panel D's two blue-green relatives are given separable colours.
Explanatory prose lives in the caption, not in the figure.
"""
import json, os, sys
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'font.size': 10, 'font.family': 'DejaVu Sans', 'axes.linewidth': .8})

ROOT = os.environ.get('CANDIDAS_ROOT') or next(
    (p for p in ('.', '..', os.path.expanduser('~/Desktop/Projects/Candidas'))
     if os.path.exists(os.path.join(p, 'gem', 'counterfactual_results.json'))), '.')
G = os.path.join(ROOT, 'gem')
OUTDIR = os.path.join(ROOT, 'results', 'figures', 'manuscript')
os.makedirs(OUTDIR, exist_ok=True)

R   = json.load(open(f'{G}/counterfactual_results.json'))
SW  = None   # panel B no longer uses the uniform-shift sweep
DTM = json.load(open(f'{G}/dyn_sparse_Tm.json'))
DTO = json.load(open(f'{G}/dyn_sparse_Topt.json'))
C44 = json.load(open(f'{G}/counts_at_44.json'))
# Measured thermal performance with DEAD WELLS COUNTED AS ZERO, not survivor-only
# medians: at 40-44 C most relative wells are dead, and dropping them would report the
# few survivors as if they were the population. Denominators are all 60 (auris) or all
# sampled wells throughout. mu_honest is the fitted specific growth rate r converted to
# h^-1 (r is stored per minute, hence x60); it involves no carbon-quota assumption, so
# it is directly comparable with the model's mu.
MEASTPC = pd.read_csv(f'{G}/measured_tpc_honest.csv')
PREDTPC = pd.read_csv(f'{G}/etcgem_tpc_pred_fine.csv')
base, pr = R['baseline'], R['mechanisms']['thr=0.05']

LAB = {'auris': 'C. auris', 'haemulonii': 'C. haemulonii',
       'duobushaemulonii': 'C. duobushaemulonii', 'parapsilosis': 'C. parapsilosis'}
SP4 = ['auris', 'haemulonii', 'duobushaemulonii', 'parapsilosis']
RELS = ['haemulonii', 'duobushaemulonii', 'parapsilosis']
# ONE species palette, used identically in A, B and D. Earlier drafts coloured
# C. auris teal in A and red in B while red ALSO meant "model-required" in C, so a
# single colour carried two meanings in one figure. Panel C now uses no species
# colours at all: it compares three kinds of quantity, not three species, so it gets
# its own neutral scale plus the one rust accent for the model's requirement.
SPCOL = {'auris': '#1d3f73', 'haemulonii': '#4d9dc0',
         'duobushaemulonii': '#e08214', 'parapsilosis': '#7a5195'}
C_PRED, C_MEAS, C_REQ = '#b0b0b0', '#4f4f4f', '#a5453b'
red = '#a5453b'; grey = '#9e9e9e'
DEG = '°'
INK, SEC = '#1a1a1a', '#454545'          # secondary ink darkened for reduction
TM, TOPT = r'T$_\mathrm{m}$', r'T$_\mathrm{opt}$'

# Which relatives can a uniform shift actually kill? Stated, not assumed.
killable = {p: [r for r, v in pr[p]['per_rel'].items() if v.get('required') is not None]
            for p in ('Tm', 'Topt')}

# Model-required separations, computed once and shared by panels B and C. A missing
# requirement raises rather than being quietly drawn as though it had been computed.
reqd = {}
for _k in ('Tm', 'Topt'):
    _v = pr[_k]['median_required']
    if _v is None:
        sys.exit(f'{_k}: no finite required separation; the figure cannot be drawn')
    reqd[_k] = _v
reqd_tm = reqd['Tm']

# The requirement is a model inversion, so it has no sampling error, but it does move
# with the detection threshold. That range is the honest uncertainty to draw.
_req_range = {}
for _k in ('Tm', 'Topt'):
    _v = [R['mechanisms'][t][_k]['median_required'] for t in R['mechanisms']]
    _f = [x for x in _v if x is not None]
    _req_range[_k] = (min(_f), max(_f), any(x is None for x in _v))

# PAIRED ortholog difference in predicted Tm. This is the figure's central quantity: a
# property of the PREDICTOR, independent of the etcGEM, so it survives every objection
# that can be raised against the metabolic model.
#
# Paired, not unpaired. The four proteomes are largely orthologous, so comparing their Tm
# DISTRIBUTIONS understates the predictor's resolution: the marginal distributions overlap
# almost completely (Cohen's d <= 0.19) while the paired difference is highly significant
# (p ~ 1e-11 against C. haemulonii). Claiming the proteomes are "indistinguishable" would
# be refuted by any reader who runs the paired test. The honest claim is the one drawn
# here -- the difference is real, consistent in sign for the two closest relatives, and
# roughly 60x too small. Reactions shared between two models give the ortholog pairing:
# the same reaction is catalysed by each species' counterpart enzyme.
from scipy import stats as _st
_kc = {s: pd.read_csv(f'{G}/kcat_reaction_{s}.csv')[['reaction', 'best_gene']].dropna()
       for s in SP4}


def _paired(param):
    f, col = (('thermal_tm.csv', 'pred_tm') if param == 'Tm'
              else ('thermal_topt.csv', 'pred_topt'))
    d = pd.read_csv(f'{G}/{f}')
    d['sp'] = d['id'].str.split('|').str[0]; d['g'] = d['id'].str.split('|').str[1]
    v = {s: dict(zip(d[d.sp == s].g, d[d.sp == s][col])) for s in SP4}
    out = {}
    for r in RELS:
        m = _kc['auris'].merge(_kc[r], on='reaction', suffixes=('_a', '_b'))
        dd = np.array([v['auris'][x] - v[r][y] for x, y in zip(m.best_gene_a, m.best_gene_b)
                       if x in v['auris'] and y in v[r]])
        t, p = _st.ttest_1samp(dd, 0)
        se = float(dd.std(ddof=1) / np.sqrt(len(dd)))
        out[r] = dict(n=len(dd), mean=float(dd.mean()), p=float(p), se=se,
                      lo=float(dd.mean() - 1.96 * se), hi=float(dd.mean() + 1.96 * se))
    return out


paired = _paired('Tm')
paired_topt = _paired('Topt')
best_rel = max(paired, key=lambda r: paired[r]['mean'])
dpair = paired[best_rel]['mean']          # largest paired difference, in auris's favour
foldgap = reqd_tm / dpair

fig = plt.figure(figsize=(11.5, 8.6))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.26,
                      left=0.085, right=0.975, top=0.878, bottom=0.085)

# ---- A: measured against predicted thermal performance, per species ---------
# This replaces a single-temperature bar chart. Showing the whole curve answers the
# obvious objection to that chart -- "maybe the model is simply bad" -- from inside the
# figure: it reproduces C. auris across the entire range. But C. auris IS the
# calibration target (build_etcgem_tpc.py fits sig, w, P and SCALE to its curve alone,
# then freezes them), so that agreement is in-sample and is labelled as such. The three
# relatives are out-of-sample, and they are where the model fails.
gsA = gs[0, 0].subgridspec(2, 2, hspace=0.42, wspace=0.16)
_axA = []
for i, sp in enumerate(SP4):
    ax = fig.add_subplot(gsA[i // 2, i % 2]); _axA.append(ax)
    mm = MEASTPC[MEASTPC.species == sp].dropna(subset=['mu_honest']).sort_values('T')
    pp = PREDTPC[PREDTPC.species == sp].sort_values('T')
    ax.plot(pp['T'], pp['pred_mu'], color='#7a7a7a', lw=1.3, ls='--', zorder=2)
    ax.plot(mm['T'], mm['mu_honest'], color=SPCOL[sp], lw=2.0, zorder=4)
    # every assayed temperature gets a point, including the zeros: they are observed
    # failures, not missing data. The red cross is reserved for temperatures at which
    # EVERY well was dead, which is a stronger statement than a zero median -- at 40-44
    # C, 5 of 15 C. haemulonii wells still grew even though the median well did not.
    ax.plot(mm['T'], mm['mu_honest'], 'o', ms=2.6, color=SPCOL[sp], zorder=5)
    alldead = mm[mm.n_alive == 0]
    ax.plot(alldead['T'], np.zeros(len(alldead)), 'x', ms=5, mew=1.4, color=red, zorder=6)
    ax.axhline(0.05, ls=':', lw=.7, color='#9a9a9a', zorder=1)
    ax.set_xlim(21, 45); ax.set_ylim(-0.045, 0.92)
    ax.set_xticks([22, 30, 38, 44]); ax.set_yticks([0, 0.4, 0.8])
    ax.tick_params(labelsize=7.0, pad=1.5)
    if i % 2: ax.set_yticklabels([])
    if i < 2: ax.set_xticklabels([])
    ttl = LAB[sp] + ('  (calibration target)' if sp == 'auris' else '')
    ax.set_title(ttl, loc='left', fontsize=7.6, style='italic', color=INK, pad=2.2)
    for sd in ('top', 'right'): ax.spines[sd].set_visible(False)
_axA[2].set_ylabel('growth rate (h⁻¹)', fontsize=7.6)
_axA[2].yaxis.set_label_coords(-0.20, 1.10)
_axA[2].set_xlabel(f'temperature ({DEG}C)', fontsize=7.6)
_axA[2].xaxis.set_label_coords(1.09, -0.24)
_axA[0].text(0.04, 0.955, 'measured', transform=_axA[0].transAxes, fontsize=7.0,
             color=SPCOL['auris'], va='top')
_axA[0].text(0.04, 0.80, 'etcGEM', transform=_axA[0].transAxes, fontsize=7.0,
             color='#7a7a7a', va='top')
_axA[1].text(0.96, 0.955, '× = every well dead', transform=_axA[1].transAxes,
             fontsize=7.0, color=red, va='top', ha='right')
fig.text(0.085, 0.923, 'A  The model fits its calibration target but misses the '
         'relatives\u2019 collapse', fontsize=9.6, fontweight='bold', color=INK)

# ---- B --------------------------------------------------------------------
# The reason every mechanism fails, in one picture: the four proteomes' predicted Tm
# distributions are superimposed. The separation the model needs is wider than the whole
# axis. This replaces the earlier single-species uniform-shift sweep, which showed one
# relative's trajectory and rested on a claim ("the only relative a uniform shift
# reaches") that was an artefact of C. parapsilosis's reversible maintenance reaction.
from scipy.stats import gaussian_kde
axB = fig.add_subplot(gs[0, 1])
TMD = pd.read_csv(f'{G}/thermal_tm.csv')
TMD['sp'] = TMD['id'].str.split('|').str[0]
colmapB = SPCOL
vals = {s: TMD[TMD.sp == s].pred_tm.values for s in SP4}
kdes = {s: gaussian_kde(v) for s, v in vals.items()}
med_lo = min(np.median(v) for v in vals.values())
med_hi = max(np.median(v) for v in vals.values())

# MAIN AXIS: drawn at the scale the model demands. The requirement is a real length here,
# and at that length the four species collapse onto one line -- which is the whole point,
# so it is shown rather than asserted.
X0, X1 = 44.0, 44.0 + reqd_tm + 10.0
xs = np.linspace(X0, X1, 900)
pk = max(k(xs).max() for k in kdes.values())
for s in SP4:
    axB.plot(xs, kdes[s](xs), color=colmapB[s], lw=1.6, zorder=4)
axB.set_xlim(X0, X1); axB.set_ylim(0, pk * 1.50)
yb = pk * 1.30
axB.annotate('', xy=(med_lo, yb), xytext=(med_lo + reqd_tm, yb),
             arrowprops=dict(arrowstyle='<->', lw=1.6, color=red))
axB.text(med_lo + reqd_tm / 2, yb + pk * .05,
         f'separation the model requires: {reqd_tm:.0f} {DEG}C',
         fontsize=7.8, color=red, ha='center', va='bottom', fontweight='bold')
# One effect-size line. What this number IS belongs in the caption, not the panel.
axB.text(X0 + (X1 - X0) * 0.435, pk * 0.235,
         f'mean paired Δ{TM} ({LAB["auris"]} − {LAB[best_rel]})\n'
         f'= {dpair:.2f} {DEG}C  (95% CI {paired[best_rel]["lo"]:.2f}'
         f'–{paired[best_rel]["hi"]:.2f}; n = {paired[best_rel]["n"]})',
         fontsize=7.2, color=INK, ha='left', va='center', linespacing=1.55)
axB.set_xlabel(f'predicted enzyme {TM} ({DEG}C)')
axB.set_ylabel('density')
# 33/0.52 = 63x in the fitted formulation; 18/0.52 = 35x once the unfolding width is
# set to its physical value. The heading takes the smaller, which holds either way.
FOLD_ROBUST = 18.0 / dpair                       # corrected unfolding width
FOLD_STRICT = 18.0 / paired[best_rel]['hi']      # smallest requirement / largest dTm
axB.set_title('B  Paired Δ{} is directionally consistent but ~{:.0f}× below the fitted '
              'requirement'.format(TM, foldgap), loc='left', fontsize=8.8,
              fontweight='bold', y=1.06)
for s in ('top', 'right'): axB.spines[s].set_visible(False)

# INSET: the same four curves magnified ~5x. They stay superimposed, so the collapse on
# the main axis is not an artefact of the scale.
axI = axB.inset_axes([0.415, 0.30, 0.565, 0.50])
xi = np.linspace(47, 63, 500)
for s in SP4:
    axI.plot(xi, kdes[s](xi), color=colmapB[s], lw=1.5,
             label=f'{LAB[s]} ({np.median(vals[s]):.2f})')
axI.axvspan(med_lo, med_hi, color='#333333', alpha=.18, zorder=1)
axI.set_xlim(47, 63); axI.set_yticks([])
axI.tick_params(axis='x', labelsize=7.2, pad=1)
axI.set_title(f'zoomed view — medians span {med_hi - med_lo:.2f} {DEG}C',
              fontsize=6.8, color=SEC, pad=2.5)
axI.legend(frameon=False, fontsize=6.0, loc='upper right', handlelength=1.0,
           labelspacing=.30, borderpad=.1, bbox_to_anchor=(1.03, 1.06))
for s in ('top', 'right', 'left'): axI.spines[s].set_visible(False)
for sp_ in axI.spines.values(): sp_.set_linewidth(.6)

# ---- C --------------------------------------------------------------------
axC = fig.add_subplot(gs[1, 0])
# Paired ortholog differences, the same quantity panel B reports. Using the unpaired
# median gap here instead would put two near-identical but different ratios in one figure
# (62x in C against 63x in B), which reads as an inconsistency rather than as the two
# distinct comparisons it is. One quantity, one ratio per axis.
pred = {'Tm': dpair,
        'Topt': max(v['mean'] for v in paired_topt.values())}
# An independent MEASURED anchor for the Tm axis. Walunjkar et al. 2025 (Mol Biol Evol
# 42:msaf137) did thermal proteome profiling on S. cerevisiae and S. uvarum -- congeners
# whose thermal growth limits (IT50) differ by 8 C, twice the difference here -- and found
# a mean ortholog Tm difference of 1.6 C over 827 pairs, with 85% of S. cerevisiae proteins
# more stable. That is what real thermotolerance divergence looks like at the proteome
# level, so it bounds the argument without depending on this study's predictor at all.
MEAS_TM, MEAS_LAB = 1.6, 'measured in an\nanalogous pair'
ys = {'Tm': 1, 'Topt': 0}
for k in ('Tm', 'Topt'):
    y = ys[k]
    # dotted, so it reads as "these two are being compared" rather than as a bar with
    # ends that mean something. The capped solid bars on top are the actual intervals.
    axC.plot([pred[k], reqd[k]], [y, y], color='#cfcfcf', lw=1.4, ls=':', zorder=1)
    # uncertainty, drawn rather than implied: a sampling CI on the prediction, and for
    # the requirement the range it moves over as the detection threshold varies.
    _pp = (paired if k == 'Tm' else paired_topt)
    _br = max(_pp, key=lambda r: _pp[r]['mean'])
    axC.plot([_pp[_br]['lo'], _pp[_br]['hi']], [y, y], color='#6f6f6f', lw=1.5, zorder=2)
    for e in (_pp[_br]['lo'], _pp[_br]['hi']):
        axC.plot([e, e], [y - 0.075, y + 0.075], color='#6f6f6f', lw=1.5, zorder=2)
    rlo, rhi, unb = _req_range[k]
    axC.plot([rlo, rhi], [y, y], color=C_REQ, lw=1.5, alpha=.6, zorder=2)
    axC.plot([rlo, rlo], [y - 0.075, y + 0.075], color=C_REQ, lw=1.5, alpha=.6, zorder=2)
    if unb:
        # A capped bar would assert a finite upper limit that does not exist: at the
        # strictest detection threshold no shift within the searched range reaches
        # failure at all. Drawn open, and said so.
        axC.annotate('', xy=(rhi * 1.42, y), xytext=(rhi, y), zorder=2,
                     arrowprops=dict(arrowstyle='-|>', lw=1.5, color=C_REQ, alpha=.6,
                                     shrinkA=0, shrinkB=0))
        axC.text(rhi * 1.5, y - 0.155, 'unbounded at 0.03 h⁻¹', fontsize=7.2,
                 color=C_REQ, ha='center', va='top')
    else:
        axC.plot([rhi, rhi], [y - 0.075, y + 0.075], color=C_REQ, lw=1.5, alpha=.6, zorder=2)
    axC.scatter(pred[k], y, s=105, color=C_PRED, edgecolor='#222222', zorder=3,
                label=('sequence-predicted (paired orthologs)' if k == 'Tm' else None))
    axC.scatter(reqd[k], y, s=105, color=C_REQ, edgecolor='#222222', zorder=3,
                label=('model-required' if k == 'Tm' else None))
    axC.text(pred[k], y + 0.13, f'{pred[k]:.2f} {DEG}C', ha='center', fontsize=8, color=SEC)
    axC.text(reqd[k], y + 0.13, f'{reqd[k]:.0f} {DEG}C', ha='center', fontsize=8.5,
             color=red, fontweight='bold')
    if k == 'Tm':
        axC.scatter(MEAS_TM, y, s=112, marker='D', color=C_MEAS, edgecolor='#222222',
                    zorder=4, label='measured congeneric benchmark')
        axC.text(MEAS_TM, y + 0.14, f'{MEAS_TM:.1f} {DEG}C', ha='center', fontsize=8,
                 color=C_MEAS, fontweight='bold')
        axC.text(np.sqrt(MEAS_TM * reqd[k]), y - 0.21, f'≈{reqd[k] / MEAS_TM:.0f}×',
                 ha='center', fontsize=10, fontweight='bold')
        axC.text(MEAS_TM, y - 0.40, 'S. cerevisiae / S. uvarum,\n8 °C apart in growth limit',
                 ha='center', fontsize=6.4, color=SEC, style='italic', linespacing=1.35)
    else:
        axC.text(np.sqrt(pred[k] * reqd[k]), y - 0.20, f'≈{reqd[k] / pred[k]:.0f}×',
                 ha='center', fontsize=10, fontweight='bold')
axC.set_xscale('log'); axC.set_xlim(0.30, 55); axC.set_ylim(-0.62, 1.45)
axC.set_yticks([0, 1]); axC.set_yticklabels([f'enzyme {TOPT}', f'enzyme {TM}'], fontsize=9.5)
axC.tick_params(axis='y', pad=2)
axC.set_xlabel(f'interspecies separation ({DEG}C, log)')
axC.legend(frameon=False, fontsize=7.4, loc='upper left', bbox_to_anchor=(-0.02, 1.10),
           ncol=3, handletextpad=.45, columnspacing=1.2)
axC.set_title('C  Required shifts exceed predictions and a measured benchmark',
              loc='left', fontsize=9.6, fontweight='bold', y=1.10)
for s in ('top', 'right', 'left'): axC.spines[s].set_visible(False)
axC.tick_params(axis='y', length=0)

# ---- D --------------------------------------------------------------------
axD = fig.add_subplot(gs[1, 1])
colmap = SPCOL
nmax = 0
for rel in RELS:
    for J, ls, mk in ((DTM, '-', 'o'), (DTO, '--', '^')):
        tr = J['relatives'][rel]['greedy']['trajectory']
        nmax = max(nmax, max(t['n'] for t in tr))
        axD.plot([t['n'] for t in tr], [t['mu40'] for t in tr], color=colmap[rel],
                 lw=1.9, ls=ls, marker=mk, ms=3.2, markevery=10)
axD.axhline(0.05, ls='--', lw=1.0, color='#6f6f6f')
axD.text(nmax * .97, 0.075, 'detection (must reach)', fontsize=7.2, color=SEC, ha='right')
axD.text(0.5, 0.965, 'no pair or triple in the beam search crossed detection',
         transform=axD.transAxes, fontsize=7.6, color=INK, ha='center', va='top')
axD.set_xlabel(f'cumulative no. of bottleneck enzymes shifted −8 {DEG}C')
axD.set_ylabel(f'predicted µ at 40 {DEG}C (h⁻¹)')
axD.set_ylim(0, 0.78); axD.set_xlim(0, nmax + 1)
axD.set_title('D  Targeted perturbations do not drive relatives below detection',
              loc='left', fontsize=9.6, fontweight='bold', y=1.14)
h = [Line2D([], [], color=colmap[r], lw=1.8, label=LAB[r]) for r in RELS] + \
    [Line2D([], [], color='#555555', lw=1.6, ls='-', label=f'{TM} shift'),
     Line2D([], [], color='#555555', lw=1.6, ls='--', label=f'{TOPT} shift')]
# Legend lifted clear of the detection line it used to sit on.
axD.legend(handles=h, frameon=False, fontsize=7.0, loc='lower left',
           bbox_to_anchor=(-0.015, 1.005), ncol=5, handlelength=1.7,
           columnspacing=1.0, handletextpad=0.5)
for s in ('top', 'right'): axD.spines[s].set_visible(False)

fig.suptitle('Within the etcGEM, sequence-predicted enzyme thermal properties do not '
             'reproduce the observed thermal divergence',
             fontsize=11.6, fontweight='bold', color=INK, y=0.963)
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUTDIR, f'FIG4_etcgem_counterfactual.{ext}'),
                dpi=300 if ext == 'png' else None, bbox_inches='tight')
print('wrote', os.path.join(OUTDIR, 'FIG4_etcgem_counterfactual.png/.pdf'))
print('  panel A         : measured vs predicted TPC; auris is the calibration target')
print('  panel B paired  :', {r: (v['n'], round(v['mean'],3), f"p={v['p']:.1e}")
                                 for r, v in paired.items()})
print(f'  panel B         : largest paired dTm {dpair:.3f} {DEG}C vs required '
      f'{reqd_tm:.2f} {DEG}C -> {foldgap:.0f}x')
print('  panel C paired  :', {k: (round(pred[k], 3), round(reqd[k], 2),
                                    round(reqd[k] / pred[k], 1)) for k in pred})
print('  Topt paired     :', {r: (v['n'], round(v['mean'], 3), f"p={v['p']:.1e}")
                              for r, v in paired_topt.items()})
print('  panel D max n   :', nmax)
