# Simulation-Based Sample Size Planning for the Power-Enhanced Mediation Test

Estimates, by Monte Carlo simulation, the sample size needed for the
power-enhanced test to reach a target power under a given mediation
pattern and signal strength. It evaluates the empirical power of the PE
test (and, for comparison, the benchmark Wald test) at each candidate
sample size on a grid, and reports the smallest grid value that reaches
the target. This is the design counterpart of the analysis function
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md):
use it to plan a study, then
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
to analyze it.

## Usage

``` r
ss_power_pe_mediation(
  outcome = c("continuous", "binary", "count"),
  pattern = c("homogeneous", "contrasting", "heterogeneous"),
  p,
  n_grid,
  c1 = 1,
  c2 = 0.5,
  target_power = 0.8,
  n_rep = 200,
  alpha_level = 0.05,
  method = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
  lambda_grid = seq(0.05, 10, length.out = 100),
  cores = 1L,
  seed = NULL
)
```

## Arguments

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- pattern:

  Mediation pattern: `"homogeneous"` or `"contrasting"`
  (`"heterogeneous"` is a synonym for `"contrasting"`).

- p:

  Number of candidate mediators.

- n_grid:

  Vector of candidate sample sizes to evaluate: whole numbers of at
  least 10 (a value listed twice is evaluated once). They are evaluated,
  and reported, in increasing order.

- c1:

  Signal strength (the exposure-on-mediator scale), a single finite
  nonzero number; at `c1 = 0` the rejection rate is the test's size,
  which no sample size raises to a power target. Default 1.

- c2:

  Direct effect. Default 0.5, the value of the article's linear
  simulations; its logistic figures use 1 and its Poisson figures 0.4.

- target_power:

  Desired power for the PE test, a single number between 0 and 1 (see
  Details for the replications it needs). Default 0.8.

- n_rep:

  Monte Carlo replications per candidate `n`. Default 200.

- alpha_level:

  Significance level. Default 0.05.

- method, error_level:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md).

- lambda_grid:

  Passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md):
  the tuning grid searched at each candidate `n`, by default
  `seq(0.05, 10, length.out = 100)`. A plan runs `n_rep` fits at every
  candidate `n`, so a shorter or narrower grid speeds it up.

- cores:

  Number of CPU cores (Unix forking). A seeded parallel run is
  reproducible across runs at the same `cores` but need not match a
  serial run. Default 1.

- seed:

  Optional integer seed. When supplied it is set for the duration of the
  call and the caller's random number generator state is restored on
  exit. `NULL` (the default) sets no seed: the draws come from the
  session's random number stream, which the call advances as any random
  function does.

## Value

A tidy `data.frame` (class `poemed_ss_plan`, a `poemed_tbl`) with one
row per candidate sample size and columns `n`, `power_pe`, `power_hdmm`,
`n_valid`, and `n_empty` (see
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)).
The smallest `n` reaching `target_power` for the PE test is stored in
the `"recommended_n"` attribute (`NA` if no grid value reaches it), with
`"target_power"`, `"outcome"`, `"pattern"`, `"c1"`, `"c2"`, `"n_rep"`,
and `"alpha_level"`. Printing the table adds a line with the target and
the recommendation.

## Details

This is a simulation-based procedure, not a closed-form power formula:
no analytic power expression exists for the PE test, so sample size is
instead determined by direct Monte Carlo simulation. Power is estimated
by simulating `n_rep` data sets at each candidate `n` with
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md)
and recording how often the test rejects at level `alpha_level`. Because
the estimate is a Monte Carlo proportion, it carries simulation error of
roughly \\\sqrt{power(1 - power) / n\\rep}\\; raise `n_rep` for a
smoother curve and a more stable recommendation. The grid approach
(rather than a root search) is deliberate: the power curve is monotone
but noisy, so a search can stop early on a lucky draw, whereas the grid
shows the whole trajectory.

An estimated power from `n_rep` replications moves in steps of
`1 / n_rep`, so a target above `1 - 1 / n_rep` can be reached only when
every replication rejects, which is weak evidence that the true power
reaches it. Such a target draws a warning of class
`poemed_target_unresolvable` naming the smallest `n_rep` that resolves
it. When no candidate sample size reaches the target, the recommendation
is `NA` and a warning of class `poemed_target_not_reached` gives the
largest estimated power and where it occurred, so the grid can be
extended. The printed table ends with the target and the recommendation
(or the statement that the grid did not reach the target).

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## See also

[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md),
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md).

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEMED/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md),
[`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md),
[`pe_simulation_study()`](https://yelleknek.github.io/POEMED/reference/pe_simulation_study.md),
[`plot.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/plot.poemed_tbl.md),
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md),
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# n_rep = 5 keeps this example fast; a power estimate from five
# replications has a Monte Carlo standard error of up to 0.22, so the
# recommended n here is a demonstration, not a plan. Planning a study
# deserves n_rep of several hundred. A plan is also where a shorter
# tuning grid pays: every replication searches the whole grid, so this
# one passes 20 values from 0.05 to 2 instead of the default 100 values
# from 0.05 to 10 (see the lambda_grid argument of ?pe_mediation).
set.seed(113)
plan <- ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                              p = 50, n_grid = c(150, 250), n_rep = 5,
                              lambda_grid = seq(0.05, 2, length.out = 20))
plan                        # the last line gives the recommendation
#>  n   power_pe power_hdmm n_valid n_empty
#>  150 1        0.2        5       0      
#>  250 1        0          5       0      
#> 
#> Outcome: continuous
#> Target power 0.8 for the power-enhanced test: recommended n = 150, the smallest n in the grid that reaches it (5 replications per n).
attr(plan, "recommended_n") # 150, the smallest grid n reaching 0.8
#> [1] 150
```
