# Simulate High-Dimensional Mediation Data

Generates a data set from the mediation models studied in the article,
for a continuous, binary, or count outcome, under the homogeneous or
contrasting (heterogeneous) mediation patterns. The mediators carry an
AR(1) correlation structure, a handful are truly active, and the rest
are null, so the data exercise exactly the situation the power-enhanced
tests are designed for. This is the generator behind the package's
example data sets and behind
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md).

## Usage

``` r
simulate_mediation_data(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  pattern = c("homogeneous", "contrasting", "heterogeneous"),
  c1 = 1,
  c2 = 0.5,
  q = 1L,
  d = 0L,
  alpha_m = NULL,
  tau = NULL,
  rho = 0.5,
  sigma_y = 0.5,
  alpha_z = NULL,
  Gamma_z = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Number of observations.

- p:

  Number of candidate mediators.

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- pattern:

  Mediation pattern: `"homogeneous"` (active effects share a sign) or
  `"contrasting"` (active effects of opposite sign). `"heterogeneous"`
  is accepted as a synonym for `"contrasting"`. Ignored when `alpha_m`
  is supplied directly. The presets are the article's. With the default
  `tau`, the contrasting presets for the continuous and binary outcomes
  cancel to a total indirect effect of zero, while the count preset, (0,
  0, 0, 0.8, -0.7), gives a total indirect effect of -0.03 times `c1`.

- c1:

  Scale on the exposure-on-mediator coefficients \\\Gamma_x = c_1
  \tau\\. `c1 = 0` makes every mediator inactive (the global null).
  Default 1.

- c2:

  The direct effect \\\alpha_x\\ of each exposure on the outcome.
  Default 0.5.

- q:

  Number of exposures. Default 1. For `q > 1` you must supply `tau` as a
  `q` by `p` matrix.

- d:

  Number of confounders. Default 0 (none). When positive, confounders
  are drawn standard normal and enter both models through `alpha_z` and
  `Gamma_z`.

- alpha_m:

  Optional numeric vector of `p` finite mediator-on-outcome
  coefficients, overriding the `pattern` preset.

- tau:

  Optional base pattern for \\\Gamma_x\\, finite numbers. By default it
  is the article's (0.1, 0.2, 0.3, 0.4, 0.5, noise), with the noise
  loadings drawn from N(0, 0.5^2), for `q = 1`. Required when `q > 1`,
  as a `q` by `p` matrix (one row per exposure) or a vector of its
  values in column order; when `q = 1` a vector of length `p` serves.

- rho:

  AR(1) correlation parameter for the mediator noise. Default 0.5.

- sigma_y:

  Outcome noise standard deviation for a continuous outcome. Default
  0.5.

- alpha_z, Gamma_z:

  Optional confounder coefficients: a vector of `d` finite numbers and a
  `d` by `p` matrix (when `d = 1`, a vector of length `p` serves for
  `Gamma_z`). Each defaults to zeros when `d > 0`. They are used only
  when `d > 0`; supplying either with `d = 0`, where it would be
  ignored, is an error.

- seed:

  Optional integer seed. When supplied it is set for the duration of the
  call and the caller's random number generator state is restored on
  exit. `NULL` (the default) sets no seed: the draws come from the
  session's random number stream, which the call advances as any random
  function does.

## Value

A list with the simulated data and the truth used to generate it: `X`
(`n` by `q`), `M` (`n` by `p`), `Y` (length `n`), `Z` (`n` by `d`, or
`NULL`), the coefficient values `alpha_m`, `Gamma_x`, `alpha_x`, the
total indirect effect `beta` (`Gamma_x` times `alpha_m`), and
`active_mediators` (the indices of the truly active mediators, those
nonzero in both paths), together with `n`, `p`, `q`, `outcome`, and
`pattern`.

## Details

One exposure column (or \\q\\ of them) is drawn standard normal. The
mediators follow \\M = X \Gamma_x + Z \Gamma_z + \varepsilon_m\\ with
\\\varepsilon_m \sim N(0, \Sigma)\\, \\\Sigma\\ the AR(1) covariance
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md)
with parameter `rho`. The exposure-on-mediator coefficient is \\\Gamma_x
= c_1 \tau\\, where the base pattern \\\tau\\ has small increasing
loadings on the first five mediators and random noise loadings on the
rest, so `c1` scales the overall exposure-to-mediator signal (and
`c1 = 0` gives the global null of no mediation). The outcome is
generated from its model with mediator coefficients \\\alpha_m\\ (the
preset for the chosen `pattern`, or a user-supplied vector) and direct
effect \\\alpha_x = c_2\\:

- continuous: \\Y = M\alpha_m + X\alpha_x + Z\alpha_z + \varepsilon_y\\,
  \\\varepsilon_y \sim N(0, \sigma_y^2)\\;

- binary: \\Y \sim \mathrm{Bernoulli}(\mathrm{logit}^{-1}( M\alpha_m +
  X\alpha_x + Z\alpha_z))\\;

- count: \\Y \sim \mathrm{Poisson}(\exp(\eta))\\, with each value of the
  log-mean \\\eta\\ clamped to the interval from -5 to 5 (a value
  outside it is set to the nearer endpoint) to guard against overflow.

Every design argument is checked before any random number is drawn: a
coefficient vector or matrix of the wrong length or shape, or one
holding a missing or infinite value, stops with a message naming it.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
to analyze the data,
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
to sweep `c1`,
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md)
for the mediator covariance.

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEMED/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEMED/reference/pe_simulation_study.md),
[`plot.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/plot.poemed_tbl.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md),
[`ss_power_pe_mediation()`](https://yelleknek.github.io/POEMED/reference/ss_power_pe_mediation.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 100, p = 50, outcome = "continuous",
                             pattern = "contrasting")
dim(d$M)
#> [1] 100  50
d$active_mediators       # truly active mediators
#> [1] 1 2 3 4
round(d$beta, 10)        # total indirect effect, zero up to rounding error
#> [1] 0
```
