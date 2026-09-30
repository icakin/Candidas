#!/usr/bin/env python3
"""ribo_seq2tm.py -- ribosomal protein stability from Seq2Tm, treated exactly
like the metabolic enzymes (task 4).  Replaces the Zeldovich/IVYWREL route
with the machinery already used for u(T).

Inputs (repo root):
  ribo_tm_predictions.csv          Seq2Tm Tm of the RBH ribosomal proteins
                                   (124 per species for auris/hae/duo, 78 para)
  gem/tables/machinery_ribosome.faa  sequences -> lengths for the two-state model
  gem/tables/rbh/pid_*.tsv         auris -> ortholog id maps (hae, duo pairable;
                                   para predictions use XP_ ids, RBH uses CPAR2_,
                                   so para is species-level only)
  model_extra/tmr_uncertainty.csv  fitted TmR + CIs (optional, for panel d)

What it does
  1. species-level summaries of ribosomal Tm (mean, median, min, 5th pct)
  2. paired ortholog differences auris - hae, auris - duo (Wilcoxon, sign,
     bootstrap CI), i.e. the ribosomal analogue of Fig. 4b
  3. u_R(T): mean unfolded fraction across ribosomal proteins per species from
     the same two-state model as the enzymes, at raw Seq2Tm scale, plus the
     weakest-link version (fraction of proteins unfolded > 50 %)
  4. the amplification the fitted dTmR would require relative to the paired
     Seq2Tm shift, set against the S. cerevisiae / S. uvarum anchor used in
     Fig. 4c (8 degC growth-limit gap <-> 1.6 degC mean ortholog Tm gap, x5)

Output: ribo_seq2tm_summary.csv, ribo_seq2tm_pairs.csv, ribo_seq2tm.png
"""
import sys, os, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, binomtest, gaussian_kde

HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from model_v5_core import SPS,COLS,two_state,REPO
os.chdir(REPO)

def fasta(p):
    h,s=None,[]
    for l in open(p):
        if l.startswith('>'):
            if h: yield h,''.join(s)
            h,s=l[1:].strip(),[]
        else: s.append(l.strip())
    if h: yield h,''.join(s)

# ── load ──────────────────────────────────────────────────────────────────────
L={h.split()[0].split('|')[-1]:len(s) for h,s in fasta('gem/tables/machinery_ribosome.faa')}
d=pd.read_csv('ribo_tm_predictions.csv')
d['sp']=d.id.str.split('|').str[0]; d['g']=d.id.str.split('|').str[1]; d['L']=d.g.map(L)
d=d[d.sp.isin(SPS)]
print(f"ribosomal proteins with Seq2Tm Tm: "+", ".join(f"{s} {int((d.sp==s).sum())}" for s in SPS))
print(f"  lengths available: {int(d.L.notna().sum())}/{len(d)}  (missing -> species median length)")
for s in SPS:
    m=d[(d.sp==s)&d.L.notna()].L.median(); d.loc[(d.sp==s)&d.L.isna(),'L']=m

# ── 1. species summaries ──────────────────────────────────────────────────────
S=d.groupby('sp').pred_tm.agg(n='size',mean='mean',median='median',sd='std',min='min',
                              q05=lambda x:x.quantile(.05),q25=lambda x:x.quantile(.25)).reindex(SPS)
S['mean_len']=d.groupby('sp').L.mean().reindex(SPS)
print("\n══ RIBOSOMAL PROTEIN Tm (Seq2Tm, raw scale) ══")
print(S.round(2).to_string())
print("  note: parapsilosis has 78 (not 124) proteins and shorter mean length; not like-for-like")

# ── 2. paired ortholog differences ────────────────────────────────────────────
A=d[d.sp=='auris'].set_index('g').pred_tm
pairs=[]; PAIR={}
for sp in ['haemulonii','duobushaemulonii','parapsilosis']:
    t=pd.read_csv(f'gem/tables/rbh/pid_{sp}.tsv',sep='\t',header=None,names=['a','b','p'])
    m=dict(zip(t.a,t.b)); B=d[d.sp==sp].set_index('g').pred_tm
    x=[(g,m[g],A[g],B[m[g]]) for g in A.index if g in m and m[g] in B.index]
    if not x: print(f"\n  auris vs {sp}: 0 pairable proteins (id systems differ) -> species-level only"); continue
    P=pd.DataFrame(x,columns=['auris_id','other_id','tm_auris','tm_other']); P['sp']=sp; P['d']=P.tm_auris-P.tm_other
    pairs.append(P); PAIR[sp]=P
    rng=np.random.default_rng(0); bs=[rng.choice(P.d,len(P),replace=True).mean() for _ in range(5000)]
    lo,hi=np.percentile(bs,[2.5,97.5])
    print(f"\n  auris − {sp}: n={len(P)} pairs   mean Δ={P.d.mean():+.2f}°C [{lo:+.2f}, {hi:+.2f}]"
          f"   median={P.d.median():+.2f}   auris higher in {(P.d>0).mean()*100:.0f}%"
          f"   Wilcoxon p={wilcoxon(P.d).pvalue:.1e}   sign p={binomtest(int((P.d>0).sum()),int((P.d!=0).sum())).pvalue:.1e}")
if pairs: pd.concat(pairs).to_csv(os.path.join(HERE,'ribo_seq2tm_pairs.csv'),index=False)

# ── 3. u_R(T) from the two-state model, raw Seq2Tm scale ──────────────────────
Tg=np.arange(20,50.01,0.5)
UR={}; WL={}
for s in SPS:
    sub=d[d.sp==s]; F=two_state(sub.L.values,sub.pred_tm.values,Tg)     # (nT, nprot) folded fraction
    UR[s]=1-F.mean(axis=1); WL[s]=(F<0.5).mean(axis=1)
print("\n══ u_R(T): mean unfolded fraction of ribosomal proteins (two-state, raw Seq2Tm) ══")
print(f"{'T':>4} "+" ".join(f"{s[:6]:>8}" for s in SPS)+"    weakest-link fraction (Tm<T)")
for T in (36,38,40,42,44,46,48):
    i=int(np.where(np.isclose(Tg,T))[0][0])
    print(f"{T:>4} "+" ".join(f"{UR[s][i]:8.4f}" for s in SPS)+"    "+" ".join(f"{WL[s][i]:6.3f}" for s in SPS))
print("  -> at the raw Seq2Tm scale essentially nothing is unfolded at 44 °C in any species,")
print("     so sequence-predicted ribosomal stability, used literally, predicts no growth limit and")
print("     no species difference over the assay range (same conclusion as Fig. 4a for the enzymes).")

# ── 4. required amplification vs anchor ───────────────────────────────────────
fit={'auris':38.98,'haemulonii':34.71,'duobushaemulonii':34.93,'parapsilosis':33.43}; ci_lo={}; ci_hi={}; src='model_v5 point estimates'
tu=os.path.join(HERE,'tmr_uncertainty.csv')
if os.path.exists(tu):
    t=pd.read_csv(tu); t=t[t.quantity.str.startswith('TmR_')]
    fx=t[t.variant.str.startswith('Ea fixed')] if 'variant' in t else t
    if len(fx):
        for _,r in fx.iterrows():
            s=r.quantity[4:]
            if s in SPS: fit[s]=r.point; ci_lo[s]=r.boot_lo95; ci_hi[s]=r.boot_hi95
        src='tmr_uncertainty.csv (Ea fixed)'
print(f"\n══ REQUIRED vs PREDICTED (fitted TmR from {src}) ══")
ANCHOR=8.0/1.6
rows=[]
for sp,P in PAIR.items():
    need=fit['auris']-fit[sp]; got=P.d.mean()
    rows.append(dict(contrast=f'auris-{sp}',fitted_dTmR=need,seq2tm_paired_dTm=got,amplification_required=need/got if got else np.nan,
                     anchor_amplification=ANCHOR,explained_at_anchor=got*ANCHOR,frac_explained_at_anchor=got*ANCHOR/need))
    print(f"  {sp:17} fitted ΔTmR={need:+.2f}   Seq2Tm paired Δ={got:+.2f}   amplification needed ×{need/got:.1f}"
          f"   | with the S. cerevisiae/uvarum anchor (×{ANCHOR:.0f}): {got*ANCHOR:+.2f}°C = {got*ANCHOR/need*100:.0f}% of the fitted gap")
R=pd.DataFrame(rows); R.to_csv(os.path.join(HERE,'ribo_seq2tm_summary.csv'),index=False)
S.to_csv(os.path.join(HERE,'ribo_seq2tm_species.csv'))

# ── figure ────────────────────────────────────────────────────────────────────
fig,ax=plt.subplots(2,2,figsize=(13,9))
a=ax[0,0]; xs=np.linspace(44,72,300)
for s in SPS:
    v=d[d.sp==s].pred_tm; a.plot(xs,gaussian_kde(v)(xs),color=COLS[s],lw=2,label=f"C. {s} (n={len(v)}, mean {v.mean():.1f})")
    a.axvline(v.mean(),color=COLS[s],lw=1,ls='--',alpha=.7)
a.set_xlabel('Seq2Tm Tm of ribosomal proteins (°C)'); a.set_ylabel('density'); a.legend(fontsize=7.5)
a.set_title('a  Per-protein ribosomal Tm from sequence (raw Seq2Tm)',fontsize=9.5,loc='left')
a=ax[0,1]
for sp,P in PAIR.items():
    a.hist(P.d,bins=30,color=COLS[sp],alpha=.55,label=f"auris − {sp}: mean {P.d.mean():+.2f}, Wilcoxon p={wilcoxon(P.d).pvalue:.1e}")
a.axvline(0,color='k',lw=1); a.set_xlabel('ΔTm per ortholog pair, auris − relative (°C)'); a.legend(fontsize=7.5)
a.set_title('b  Paired ortholog differences (ribosomal analogue of Fig. 4b)',fontsize=9.5,loc='left')
a=ax[1,0]
for s in SPS: a.plot(Tg,UR[s],color=COLS[s],lw=2,label=f"C. {s}")
a.axvspan(22,44,color='grey',alpha=.08); a.text(33,a.get_ylim()[1]*0.9 if a.get_ylim()[1]>0 else 0.01,'assay range',ha='center',fontsize=8,color='grey')
a.set_xlabel('Temperature (°C)'); a.set_ylabel('u_R(T): mean unfolded fraction'); a.legend(fontsize=7.5)
a.set_title('c  Ribosomal u_R(T) from the same two-state model as the enzymes',fontsize=9.5,loc='left')
a=ax[1,1]
for s in SPS:
    x=S.loc[s,'mean']; y=fit[s]
    if s in ci_lo: a.errorbar(x,y,yerr=[[y-ci_lo[s]],[ci_hi[s]-y]],fmt='none',ecolor=COLS[s],lw=1.2,capsize=3)
    a.scatter(x,y,s=140,color=COLS[s],edgecolor='k',zorder=5,marker='s' if s=='parapsilosis' else 'o')
    a.annotate(f"C. {s}"+(' (78 shorter proteins)' if s=='parapsilosis' else ''),(x,y),xytext=(7,4),textcoords='offset points',fontsize=8)
a.set_xlabel('mean Seq2Tm ribosomal Tm (°C)'); a.set_ylabel('fitted TmR (°C)')
a.set_title('d  Fitted TmR against sequence-predicted ribosomal Tm (n = 4; para not like-for-like)',fontsize=9.5,loc='left')
fig.suptitle('Ribosomal protein stability from sequence (Seq2Tm) treated like the enzymes: consistent in sign, ~10× too small',fontsize=10.5)
plt.tight_layout(); plt.savefig(os.path.join(HERE,'ribo_seq2tm.png'),dpi=150,bbox_inches='tight')
print("\nSaved: ribo_seq2tm_species.csv, ribo_seq2tm_pairs.csv, ribo_seq2tm_summary.csv, ribo_seq2tm.png")
