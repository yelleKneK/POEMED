# Monte Carlo Size and Power Curve for the PE Mediation Tests

Runs the article's simulation design: for each value of the
signal-strength scale \\c_1\\ on a grid, it simulates many data sets,
applies both the benchmark Wald test and the power-enhanced test, and
returns their empirical rejection rates. At \\c_1 = 0\\ (the global null
of no mediation) the rejection rate estimates the Type I error rate; at
nonzero \\c_1\\ it estimates power. The headline finding is visible
directly in the table: under a contrasting pattern the power-enhanced
test climbs toward one as \\\|c_1\|\\ grows, while the benchmark test,
whose target is a total indirect effect that cancels, gains power much
more slowly.

## Usage

``` r
pe_power_curve(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  pattern = c("homogeneous", "contrasting", "heterogeneous"),
  c1_grid = seq(0, 1, by = 0.25),
  c2 = 0.5,
  n_rep = 100,
  alpha_level = 0.05,
  method = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
  lambda_grid = seq(0.05, 10, length.out = 100),
  outcome_args = list(),
  cores = 1L,
  seed = NULL,
  progress = FALSE
)
```

## Arguments

- n, p:

  Number of observations and candidate mediators per data set.

- outcome:

  Outcome type passed to
  [`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md):
  `"continuous"`, `"binary"`, or `"count"`.

- pattern:

  Mediation pattern: `"homogeneous"` or `"contrasting"`
  (`"heterogeneous"` is a synonym for `"contrasting"`).

- c1_grid:

  Numeric vector of signal-strength scales to sweep. Default
  `seq(0, 1, by = 0.25)`. Include 0 to estimate the Type I error rate;
  the article also uses negative values.

- c2:

  Direct effect used in the simulation. Default 0.5, the value of the
  article's linear simulations; its logistic figures use 1 and its
  Poisson figures 0.4.

- n_rep:

  Number of Monte Carlo replications per grid point. Default 100. The
  article uses 1000.

- alpha_level:

  Significance level for the global test. A replication counts as a
  rejection when its p-value is at most `alpha_level`. Default 0.05.

- method, error_level:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md):
  the multiplicity method and target error rate for the
  mediator-screening step.

- lambda_grid:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md):
  the tuning grid for the penalized fit, by default
  `seq(0.05, 10, length.out = 100)` (see Details). A shorter or narrower
  grid speeds a study up.

- outcome_args:

  A named list of further arguments forwarded to
  [`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md)
  (for example `rho`, `q`, `d`, `tau`, `alpha_m`). Default empty. Each
  name must identify one argument of
  [`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md);
  an unknown or ambiguous name stops with an error rather than reaching
  the simulator. The design arguments this function sets itself (`n`,
  `p`, `outcome`, `pattern`, `c1`, `c2`) and `seed` may not appear here,
  whether spelled in full or abbreviated; a `seed` inside `outcome_args`
  would make every replication draw the same data.

- cores:

  Number of CPU cores for the Monte Carlo replications. Values above 1
  fork via the base parallel package (Unix only). With parallel cores
  the replications run on the `"L'Ecuyer-CMRG"` generator (the caller's
  generator kind is restored on exit), so a seeded parallel run is
  reproducible across runs at the same `cores` but need not match a
  serial run. The parallel streams are handed out afresh from the same
  generator state at every point of `c1_grid`, so with `cores` above 1
  each grid point replays the same `n_rep` random draws (common random
  numbers): the rates at neighboring grid points share their Monte Carlo
  noise, and the curve looks smoother than the separate errors of its
  points suggest. With `cores = 1` every grid point draws fresh data.
  Default 1.

- seed:

  Optional integer seed. When supplied it is set for the duration of the
  call and the caller's random number generator state is restored on
  exit. `NULL` (the default) sets no seed: the draws come from the
  session's random number stream, which the call advances as any random
  function does.

- progress:

  `TRUE` or `FALSE`; if `TRUE`, print a line per grid point as it
  completes. Default `FALSE`.

## Value

A tidy `data.frame` of class `poemed_tbl` with one row per grid point
and columns `c1`, `rejection_hdmm`, and `rejection_pe` (the empirical
rejection rates of the benchmark and power-enhanced tests), `n_valid`
(replications whose fit succeeded, the denominator of the rates), and
`n_empty` (replications whose penalized fit selected no mediator,
counted as non-rejections). The simulation settings are recorded in the
attributes `outcome`, `pattern`, `n`, `p`, `n_rep`, `alpha_level`,
`error_level`, and `lambda_grid`.

## Details

A rejection rate from `n_rep` replications carries a Monte Carlo
standard error of about \\\sqrt{r(1 - r) / n\_{rep}}\\; at 100
replications a rate near 0.05 is known to about 0.02 and a rate near 0.5
to about 0.05. The article uses 1000 replications. A replication whose
penalized fit selects no mediator counts as a non-rejection (both
p-values are 1), POEMED's reporting convention (the article states only
that \\J_m = 0\\ when the selection is empty); the number of such
replications is reported in `n_empty`. A replication whose fit fails
outright is dropped from the denominator and reported in `n_valid`, with
a warning that quotes the first error. A design that cannot be simulated
(for example a `p` too small for the pattern, or a misspelled name in
`outcome_args`) stops before any replication runs, with the same message
at any value of `cores`.

The tuning grid matters for what these curves show. Every replication
searches the same grid: the package default
`seq(0.05, 10, length.out = 100)`, or the grid passed as `lambda_grid`.
A study of many replications is where a grid tuned to the situation
pays: fewer values, or a range narrowed to where the HBIC minimum has
been seen to fall, make each fit faster at the price of a coarser
choice, and whatever the grid, each fit keeps the best value the grid
offers. The single-fit warning for a choice at an end of the grid (class
`poemed_grid_boundary`) is not raised inside a study, which searches a
fixed grid by design. The reproduction vignette states the grids of the
article's own scripts.

Under contrasting mediation the benchmark can have power although the
total indirect effect is zero: a penalized fit that keeps only part of
the canceling set leaves a total indirect effect, over the kept
mediators, that no longer cancels, and how often that happens depends on
the tuning grid. The article's logistic and Poisson figures fix \\c_2\\
at 1 and 0.4, not at this function's default of 0.5, so pass `c2` to
match them.

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md).

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEMED/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEMED/reference/pe_simulation_study.md),
[`plot.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/plot.poemed_tbl.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md),
[`ss_power_pe_mediation()`](https://yelleknek.github.io/POEMED/reference/ss_power_pe_mediation.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# n_rep = 5 keeps this example fast; a rate from five replications has a
# Monte Carlo standard error of up to 0.22, so read the shape, not the
# numbers. A reported study uses the article's design (n = 300, p = 500)
# and 1000 replications. A study is also where a shorter tuning grid
# pays: every replication searches the whole grid, so this one passes 20
# values from 0.05 to 2 instead of the default 100 values from 0.05 to 10
# (see the lambda_grid argument of ?pe_mediation).
set.seed(113)
pe_power_curve(n = 200, p = 60, outcome = "continuous",
               pattern = "contrasting", c1_grid = c(0, 0.5, 1),
               n_rep = 5, lambda_grid = seq(0.05, 2, length.out = 20))
#>  c1  rejection_hdmm rejection_pe n_valid n_empty
#>  0   0              0            5       0      
#>  0.5 0              0.8          5       0      
#>  1   0              1            5       0      
#> 
#> Outcome: continuous
```
