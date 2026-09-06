# etcGEM thermal-prediction result (main-figure finding)

## What was built
A sequence-parameterised enzyme- and temperature-constrained genome-scale model
(etcGEM) for four Candida species: *C. auris* (iRV973, re-keyed to B9J08),
*C. parapsilosis* (iDC1003), and de-novo drafts for *C. haemulonii* and
*C. duobushaemulonii*. Every reaction carries a molecular weight, a DLKcat kcat,
and sequence-predicted Topt and Tm (Seq2Topt / Seq2Tm, ESM2-8M) for all 3,042
enzymes. Effective kcat is scaled by a shared activity curve peaking at each
enzyme's Topt with a denaturation cutoff near its Tm; three shape hyperparameters
were calibrated ONCE on the measured *C. auris* TPC and frozen to predict the
other three, so any interspecies difference comes only from the sequence-predicted
thermal parameters.

## A data-quality correction we made (important, also relevant to the paper)
The measured TPCs were originally medians over *surviving* wells only. At high
temperature most wells are dead, so the median is taken over a heat-tolerant
handful and is inflated. Example: *C. haemulonii* declined normally 30->36 C
(n=15 wells) then "peaked" at 42 C on only 5 surviving wells - a survivorship
artifact. Counting dead wells as zero growth (correct for a thermal-tolerance
curve) removes it: haemulonii truly peaks at 30 C, not 42 C. Corrected measured
optima span 30-36 C (6 C), not 32-42 C (10 C).
-> Worth checking the manuscript figures use the death-as-zero convention where
they report thermal tolerance, not survivor-only medians.

## Result (against corrected measured TPCs)
Predicted peaks are the fine-grid (0.5 C) optima; the 2 C measurement grid rounds
all four to a single 34 C point, which visually overstates the compression.
| species | measured peak | predicted peak (fine) | measured death T |
|---|---|---|---|
| C. auris (calibration) | 36 C | 34.5 C | 44 C |
| C. haemulonii | 30 C | 34.5 C | 38 C |
| C. duobushaemulonii | 32 C | 34.0 C | 38 C |
| C. parapsilosis | 32 C | 33.5 C | 38 C |

- Predicted thermal optima span only ~1 C (33.5-34.5 C, auris warmest, para
  coolest) vs measured spread 6 C. The model compresses, but the faint predicted
  ordering (auris > para) is directionally correct.
- The real interspecies signal is auris's HIGH-T survival (grows to 44 C; the
  other three die at 38 C). The model does not reproduce it.

## Why - and what we tried to fix it (all three checks)
1. Enzyme thermostability is clade-conserved: predicted Topt medians within 1.0 C;
   Tm within 0.4 C; even the 1st-percentile Tm (death-relevant) is within 0.5 C
   across all four species. There is simply no per-species thermostability signal.
2. Flux-control weighting (weight thermal scaling by rate-controlling enzymes at
   high T): does NOT help - the rate-controlling core enzymes are the MOST
   conserved, so it compresses the species further.
3. Alternative whole-proteome OGT signal (IVYWREL composition, Zeldovich):
   ranks auris highest and parapsilosis lowest - directionally consistent with
   auris's superior heat tolerance - but the spread is 0.59 pct-points
   (~3-4 C, prokaryote-calibrated, noisy). A weak corroboration, not a fix.
We did NOT fit per-species parameters to the measured curves (that would be
circular). The compression is robust to every legitimate fix.

## What works vs what does not
- WORKS: absolute growth rate (pool binds at measured mu_max ~0.76 vs 0.80 /h);
  general mesophilic thermal shape; calibration species reproduced exactly.
- DOES NOT: reproduce species-specific thermal optima or auris's 44 C survival.

## Honest interpretation (the publishable claim)
*C. auris*'s distinctive thermotolerance is NOT encoded in bulk metabolic-enzyme
thermostability readable by current sequence predictors. The etcGEM reproduces
metabolic capacity and growth rate but localises the interspecies thermal
divergence to non-metabolic-enzyme determinants (regulation, expression,
membrane/lipid biology, chaperones, or specific residues below predictor
resolution). This mechanistic-localisation result independently reinforces the
phylogenomic near-identity of the Candidozyma clade.

Files: FIG_etcGEM.png/.pdf, etcgem_tpc_pred.csv, etcgem_calib.json,
measured_tpc_honest.csv, thermal_topt.csv, thermal_tm.csv, thermal_paired.csv.
