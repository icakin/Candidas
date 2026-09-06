# The within-*C. auris* clade spread in CUE optimum: verdict

**Do not quote it.** Gab's instinct was right, and the between-species contrast should
carry the premise on its own.

## What was quoted

`phylo/phylogenomics_results.md` states that the four *C. auris* clades differ by 2.74 °C
in CUE optimum while their genome signature is flat. That number came from the **fitted**
(Sharpe–Schoolfield, per-cell) optima, which inherit the N0 = N_inoc·e^(r·δ)
back-projection. Within-*auris* fine structure from that route is exactly the class of
result the paper already cut for N0 fragility.

## What was tested

The **model-free** optimum does not inherit that problem: CUE = 1/(1 + a·K/(r·e^(r·δ)))
with *a* temperature-independent, so the optimum is argmin(K/r) and every conversion
constant cancels. So the fair test is whether the clade differences survive there.

Marginal intervals are not that test. Two overlapping intervals can still have a
difference excluding zero, and two non-overlapping ones always do. So the difference was
bootstrapped directly, resampling series within temperature within each clade, 2000
replicates, using the estimator of `02_scale_free.R` verbatim.

## Result

Model-free optima: Clade1 27.42 [26.86, 27.69], Clade2 28.23 [27.69, 33.72],
Clade3 30.06 [24.61, 31.60], Clade4 31.48 [27.97, 32.54]. Point-estimate spread 4.06 °C,
not 2.74 °C.

| contrast | diff (°C) | 95% CI | P(same sign) | excludes 0 | survives multiplicity |
|---|---|---|---|---|---|
| Clade1 − Clade2 | −0.81 | [−6.37, −0.19] | 0.998 | yes | yes |
| Clade1 − Clade3 | −2.64 | [−4.27, +2.62] | 0.963 | no | no |
| Clade1 − Clade4 | −4.06 | [−5.22, −0.47] | 0.983 | yes | no |
| Clade2 − Clade3 | −1.83 | [−3.52, +4.15] | 0.833 | no | no |
| Clade2 − Clade4 | −3.25 | [−4.37, +3.78] | 0.832 | no | no |
| Clade3 − Clade4 | −1.42 | [−5.21, +2.92] | 0.716 | no | no |

**Two of six exclude zero nominally. At most one of six survives correction for six comparisons**
(Bonferroni-style threshold P(same sign) ≥ 0.9958). The table above is the Python run of the
estimator. The R script `reports/C11_scale_free/R/05_clade_contrasts.R`, run on the Mac with the
same seed logic, gives identical point estimates and P(same sign) = 0.995 for Clade1 − Clade2,
which falls just below the 0.9958 threshold: **0 of 6 survive in R**. The contrast sits on the
line; whether it is counted or not depends on bootstrap noise, which is itself the point.

## Why even that one is not quotable

The single surviving contrast, Clade1 − Clade2, has a point estimate of **0.81 °C** and an
interval running to **−6.37 °C**: the interval is roughly eight times wider than the effect
it brackets, and it clears zero at −0.19. Its width is driven by Clade2's own bootstrap
instability (that clade's optimum interval spans 27.7 to 33.7 °C). Clades 2, 3 and 4 each
have optimum intervals 4–9 °C wide, so the fine ordering among them is simply not resolved
by these data.

The one comparison with a clean-looking separation, Clade1 versus Clade4 at −4.06 °C, does
not survive multiplicity, and picking it out of six after the fact is the error the
correction exists to prevent.

## Action taken

The 2.74/2.7 °C spread is removed from the phylogenomics write-up and from the grant
materials. The premise it was supporting — that thermal physiology varies without
corresponding genomic variation — rests on the between-species contrasts, which are large,
survive every N0 treatment, and need no fine structure.

Reproduce: `Rscript reports/C11_scale_free/R/05_clade_contrasts.R`
