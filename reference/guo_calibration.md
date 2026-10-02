# Calibration Constants for the Real-Data-Motivated Heterogeneous Setting

The fixed constants for the real-data-motivated heterogeneous mediation
simulation of Yu and Kelley (in press), calibrated to the
DNA-methylation case study of Guo et al. (2022): a \\p = 1008\\-mediator
linear model with eleven active loci whose effects mix positive and
negative signs.
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md)
generates data from these constants; the object is exported so the exact
calibration is inspectable.

## Usage

``` r
guo_calibration
```

## Format

A list with components:

- p, s, n:

  The mediator count (1008), confounder count (9, an intercept plus
  eight covariates), and case-study sample size (85).

- locations:

  Indices of the eleven active mediators.

- alpha_m:

  Outcome-mediator coefficients at the active loci (the designed
  simulation values).

- Gamma_x:

  Exposure-mediator coefficients at the active loci.

- alpha_z:

  Confounder-outcome coefficients (length 9).

- Gamma_z:

  Confounder-mediator coefficients at the active loci (an 11 by 9
  matrix).

- beta_per_c1:

  The total indirect effect per unit `c1`, equal to
  `sum(Gamma_x * alpha_m)` \\\approx -1.597\\ (the article reports
  -1.5977 from the unrounded coefficients).

- alpha_m_estimated:

  The raw Guo et al. (2022) estimate of the outcome-mediator
  coefficients (for reference; the simulation uses the designed
  `alpha_m`).

- alpha_m_variants:

  Two alternative outcome-mediator coefficient sets the article also
  studies, `homogeneous_like` and `contrasting_like`, which remain
  heterogeneous when paired with `Gamma_x`.

## Source

Transcribed from the supplement of Yu and Kelley (in press) (its
parameter-configuration section), which in turn calibrates to the
`simulation_allS.Rdata` of Guo et al. (2022).

## Details

Pairing the exposure-mediator coefficients `Gamma_x` with the
outcome-mediator coefficients `alpha_m` gives a total indirect effect of
`beta_per_c1` \\\approx -1.597\\ per unit of the signal scale \\c_1\\.
(The article reports \\-1.5977\\, computed from the full-precision Guo
coefficients; the values shipped here are those coefficients rounded to
the three decimals printed in the supplement, which give \\-1.597\\.)

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

[`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md)

Other mediation simulation:
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEMED/reference/pe_simulation_study.md),
[`plot.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/plot.poemed_tbl.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md),
[`ss_power_pe_mediation()`](https://yelleknek.github.io/POEMED/reference/ss_power_pe_mediation.md)

## Examples

``` r
str(guo_calibration, max.level = 1)
#> List of 11
#>  $ p                : int 1008
#>  $ s                : int 9
#>  $ n                : int 85
#>  $ locations        : int [1:11] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ alpha_m          : num [1:11] 1 0.9 0.8 -0.9 -0.8 -0.7 0.6 0.5 0.4 0.3 ...
#>  $ Gamma_x          : num [1:11] -0.251 -0.221 -0.233 0.251 0.295 0.282 -0.332 -0.359 0.335 -0.345 ...
#>  $ alpha_z          : num [1:9] -0.336 -0.07 0.665 0.278 0.315 0.201 0.173 0.51 0.315
#>  $ Gamma_z          : num [1:11, 1:9] -0.045 -0.197 -0.076 -0.052 -0.033 -0.012 0.203 0.017 -0.364 -0.001 ...
#>  $ beta_per_c1      : num -1.6
#>  $ alpha_m_estimated: num [1:11] 0.166 0.243 0.248 -0.049 -0.294 -0.187 0.148 0.087 0.112 0.223 ...
#>  $ alpha_m_variants :List of 2
# The eleven individual indirect effects: nine negative and two positive
sign(guo_calibration$Gamma_x * guo_calibration$alpha_m)
#>  [1] -1 -1 -1 -1 -1 -1 -1 -1  1 -1  1
# Their sum is the total indirect effect per unit of c1
sum(guo_calibration$Gamma_x * guo_calibration$alpha_m)
#> [1] -1.597
```
