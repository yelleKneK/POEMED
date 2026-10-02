# Simulate the Real-Data-Motivated Heterogeneous Mediation Setting

Generates a linear mediation data set from the article's
*real-data-motivated* heterogeneous setting, whose coefficients are
calibrated to the DNA-methylation case study of Guo et al. (2022).
Unlike the homogeneous and contrasting presets of
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md),
this setting has eleven active mediators with a genuine **mix of
positive and negative** effects that neither all agree in sign
(homogeneous) nor exactly cancel (contrasting): the messy middle ground
that real data tends to present. The calibration constants are the
[guo_calibration](https://yelleknek.github.io/POEMED/reference/guo_calibration.md)
object; the total indirect effect is `guo_calibration$beta_per_c1` times
\\c_1\\, about \\-1.597 c_1\\. (The article reports \\-1.5977 c_1\\,
from the full-precision Guo coefficients; the values shipped here are
those coefficients rounded to three decimals.)

## Usage

``` r
simulate_guo_mediation(
  c1,
  c2 = 0.5,
  confounders = FALSE,
  alpha_m = "setting0",
  n = NULL,
  seed = NULL
)
```

## Arguments

- c1:

  Signal-strength scale for the exposure-mediator coefficients
  (`Gamma_x = c1 * Gamma_x_G`), a single finite number. At `c1 = 0`
  there is no mediation (the Type I error setting); the article sweeps
  `c1` over `c(0, +/- 0.1, ..., +/- 1)`.

- c2:

  Direct effect (exposure-outcome coefficient), a single finite number.
  Default 0.5.

- confounders:

  `TRUE` or `FALSE`; if `TRUE`, include the calibrated confounders in
  the data-generating process (the article's "with confounders"
  scenario). Default `FALSE`.

- alpha_m:

  Outcome-mediator coefficients at the eleven active loci: either a
  numeric vector of 11 finite values, or one of the strings `"setting0"`
  (the default, the values 1, 0.9, 0.8, -0.9, -0.8, -0.7, 0.6, 0.5, 0.4,
  0.3, and 0.2), `"homogeneous_like"`, or `"contrasting_like"` (the two
  alternative sets the article also studies, which remain heterogeneous
  because they are paired with the scattered `Gamma_x`).

- n:

  Number of observations, a single whole number of at least 2. Default
  `NULL`, the case-study size, 85.

- seed:

  Optional integer seed. When supplied it is set for the duration of the
  call and the caller's random number generator state is restored on
  exit. `NULL` (the default) sets no seed: the draws come from the
  session's random number stream, which the call advances as any random
  function does.

## Value

A list with the same shape as
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md):
the exposure `X` (`n` by 1), mediators `M` (`n` by 1008), outcome `Y`,
confounders `Z` (`n` by 8, or `NULL`), the full coefficient vectors
`alpha_m` and `Gamma_x`, the total indirect effect `beta`, the
`active_mediators` (the loci nonzero in both paths: all eleven under the
default `alpha_m`, the first five under `"homogeneous_like"`, the first
four under `"contrasting_like"`, and none when `c1 = 0`), and the
dimensions.

## Details

The mediator dimension (\\p = 1008\\), the eleven active loci, and the
coefficient vectors come from Guo et al. (2022) and are fixed. The
exposure and confounders, which the article took from the actual
case-study data, are simulated here (standard normal) so the function is
self-contained. The mediator and outcome coefficients, the part of the
calibration that matters scientifically, are the values printed in the
supplement of Yu and Kelley (in press), where the Guo estimates appear
to three decimals.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

Guo, X., Li, R., Liu, J., & Zeng, M. (2022). High-dimensional mediation
analysis for selecting DNA methylation loci mediating childhood trauma
and cortisol stress reactivity. *Journal of the American Statistical
Association, 117*(539), 1110–1121.
[doi:10.1080/01621459.2022.2053136](https://doi.org/10.1080/01621459.2022.2053136)

## See also

[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md),
[guo_calibration](https://yelleknek.github.io/POEMED/reference/guo_calibration.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md).

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEMED/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEMED/reference/pe_simulation_study.md),
[`plot.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/plot.poemed_tbl.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md),
[`ss_power_pe_mediation()`](https://yelleknek.github.io/POEMED/reference/ss_power_pe_mediation.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_guo_mediation(c1 = 0.5)
dim(d$M)               # 85 x 1008
#> [1]   85 1008
d$active_mediators     # the 11 active loci
#>  [1]  1  2  3  4  5  6  7  8  9 10 11
round(d$beta, 4)       # beta_per_c1 * 0.5, i.e. -0.7985
#> [1] -0.7985
# At the case study's n = 85 and c1 = 0.5 the signal is faint for either
# test. With c1 = 1 and n = 150 both tests reject the global null, and the
# power-enhanced test also selects seven of the eleven active loci. The
# fit searches a 20-value grid to keep the example fast; the default grid
# selects the same seven loci and takes five times as long.
d2 <- simulate_guo_mediation(c1 = 1, n = 150, seed = 113)
fit <- pe_mediation(d2$X, d2$Y, d2$M, outcome = "continuous",
                    lambda_grid = seq(0.05, 2, length.out = 20))
fit[fit$term %in% c("pval_hdmm", "pval_pe", "n_selected_mediators"), ]
#>  term                 value   
#>  pval_hdmm            < 0.0001
#>  pval_pe              < 0.0001
#>  n_selected_mediators 7       
#> 
#> Outcome model: continuous (linear)
#> Selected mediators (7): 1, 4, 5, 6, 9, 10, 11
#> Tuning parameter (HBIC): lambda = 0.153 from 20 values in [0.05, 2]
```
