# Run the Power Study Under Several Mediation Patterns at Once

A convenience wrapper around
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
that runs the size and power study under more than one mediation pattern
in a single call and stacks the results, so a homogeneous and a
contrasting curve (the two panels of a figure in the article) come back
together. Each pattern is swept over the same signal grid.

## Usage

``` r
pe_simulation_study(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  patterns = c("homogeneous", "contrasting"),
  c1_grid = seq(0, 1, by = 0.25),
  c2 = 0.5,
  n_rep = 100,
  alpha_level = 0.05,
  method = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
  lambda_grid = seq(0.05, 10, length.out = 100),
  cores = 1L,
  seed = NULL
)
```

## Arguments

- n, p:

  Observations and candidate mediators per simulated data set.

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- patterns:

  Character vector of mediation patterns to run, each named once. The
  default runs `"homogeneous"` and `"contrasting"`.

- c1_grid, c2, n_rep, alpha_level, method, error_level, lambda_grid,
  cores:

  Passed to
  [`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md).

- seed:

  Optional integer seed. When supplied it is set for the duration of the
  call and the caller's random number generator state is restored on
  exit. `NULL` (the default) sets no seed: the draws come from the
  session's random number stream, which the call advances as any random
  function does.

## Value

A tidy `data.frame` (class `poemed_tbl`) with a leading `pattern` column
and the
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
columns (`c1`, `rejection_hdmm`, `rejection_pe`, `n_valid`, `n_empty`)
for each pattern. The settings are recorded in the attributes `outcome`,
`n`, `p`, `n_rep`, `alpha_level`, `error_level`, and `lambda_grid`.

## Details

Everything
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
says about Monte Carlo error, empty fits, and the tuning grid applies
here; the article's logistic and Poisson figures fix `c2` at 1 and 0.4.
The patterns take their data from one continuing random stream, so no
two patterns share their draws. With `cores` above 1 the points of
`c1_grid` within a pattern do share them (see `cores` on
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)).

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## See also

[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
for a single pattern,
[`plot.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/plot.poemed_tbl.md)
to draw one pattern's curves.

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEMED/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md),
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
pe_simulation_study(n = 120, p = 50, outcome = "continuous",
                    c1_grid = c(0, 1), n_rep = 5,
                    lambda_grid = seq(0.05, 2, length.out = 20))
#>  pattern     c1 rejection_hdmm rejection_pe n_valid n_empty
#>  homogeneous 0  0              0            5       0      
#>  homogeneous 1  0.6            1            5       0      
#>  contrasting 0  0              0            5       0      
#>  contrasting 1  0              1            5       0      
#> 
#> Outcome: continuous
```
