# Monte Carlo Study of Individual-Mediator Identification (FWER / FDR)

Evaluates how well the power-enhanced screen recovers the *individual*
active mediators, the experiment reported in the article's supplement.
For each signal-strength scale \\c_1\\ on a grid it simulates many data
sets, selects a set of mediators under each multiplicity method, and
scores those selections against the known truth, returning the empirical
familywise error rate, false discovery rate, precision, and recall. This
complements
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md),
which evaluates the *global* test rather than which mediators are
flagged.

## Usage

``` r
pe_identification_study(
  n,
  p,
  outcome = c("continuous", "binary", "count"),
  pattern = c("homogeneous", "contrasting", "heterogeneous"),
  c1_grid = seq(0, 1, by = 0.25),
  c2 = 0.5,
  n_rep = 200,
  methods = c("Bonferroni", "BH", "BY"),
  error_level = 0.05,
  truth = c("outcome_effect", "mediation"),
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
  `seq(0, 1, by = 0.25)`. Include 0 to estimate the null familywise
  error rate.

- c2:

  Direct effect used in the simulation. Default 0.5.

- n_rep:

  Number of Monte Carlo replications per grid point. Default 200. The
  article uses 1000.

- methods:

  Multiplicity methods to evaluate, one or more of `"Bonferroni"`
  (familywise error rate), `"BH"`, and `"BY"` (false discovery rate),
  each named once. Default all three.

- error_level:

  Target error rate for the mediator-screening step (the FWER level for
  Bonferroni, the FDR level for BH and BY). Default 0.05.

- truth:

  Which mediators count as truly active: `"outcome_effect"` (nonzero
  outcome coefficient, the convention of the supplement's identification
  tables) or `"mediation"` (both paths nonzero, the article's definition
  of the true active set in its Theorem 3). See Details.

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
  an unknown or ambiguous name stops with an error. The design arguments
  this function sets itself (`n`, `p`, `outcome`, `pattern`, `c1`, `c2`)
  and `seed` may not appear here, whether spelled in full or
  abbreviated.

- cores:

  Number of CPU cores for the replications. Values above 1 fork via the
  base parallel package (Unix only) and run the replications on the
  `"L'Ecuyer-CMRG"` generator (the caller's generator kind is restored
  on exit); a seeded parallel run is reproducible across runs at the
  same `cores` but need not match a serial run. The parallel streams are
  handed out afresh from the same generator state at every point of
  `c1_grid`, so with `cores` above 1 each grid point replays the same
  `n_rep` random draws (common random numbers) and the rows at different
  grid points share their Monte Carlo noise. With `cores = 1` every grid
  point draws fresh data. Default 1.

- seed:

  Optional integer seed. When supplied it is set for the duration of the
  call and the caller's random number generator state is restored on
  exit. `NULL` (the default) sets no seed: the draws come from the
  session's random number stream, which the call advances as any random
  function does.

- progress:

  `TRUE` or `FALSE`; if `TRUE`, print a line per grid point. Default
  `FALSE`.

## Value

A tidy `data.frame` of class `poemed_tbl` with one row per (\\c_1\\,
method) combination and columns `c1`, `method`, `fwer`, `fdr`,
`precision`, `recall` (as defined in Details), `n_valid` (replications
whose fit succeeded), and `n_empty` (replications whose penalized fit
selected no mediator). The settings are recorded in the attributes
`outcome`, `pattern`, `n`, `p`, `n_rep`, `error_level`, `truth`, and
`lambda_grid`.

## Details

The supplement of Yu and Kelley (in press) reports four metrics,
empirical FWER, empirical FDR, precision, and recall, without defining
them; this function computes them as follows. In each replication the
selected set is compared with the truth set. The familywise error rate
is the proportion of replications selecting at least one mediator
outside the truth set. The false discovery rate is the mean false
discovery proportion (zero when nothing is selected). Precision is the
mean proportion of selected mediators that are true (scored zero when
nothing is selected), and recall is the mean proportion of true
mediators selected. The article uses 1000 replications.

The familywise error rate is a proportion of replications, so its Monte
Carlo standard error is \\\sqrt{r(1 - r) / n\_{valid}}\\, with \\r\\ the
reported rate. The other three metrics are means of per-replication
fractions between 0 and 1, for which that formula is only an upper
bound. Their standard error is the standard deviation of the
per-replication values divided by \\\sqrt{n\_{valid}}\\, and for recall
it can be as little as a third of the bound.

What counts as a true mediator is set by `truth`. Under
`truth = "outcome_effect"` (the default) a mediator is true when its
outcome coefficient \\\alpha\_{m,j}\\ is nonzero. This is the convention
of the supplement's identification tables: their \\c_1 = 0\\ rows report
positive precision and recall, which only a truth set that is nonempty
when the exposure-on-mediator paths are all zero allows. There the
familywise error rate is the probability of selecting a mediator with no
outcome effect. Under `truth = "mediation"` a mediator is true only when
both paths are nonzero. That is the definition of the true active set in
Theorem 3 of the article, against which its familywise error guarantee
is stated, and of the simulator's `active_mediators`. At \\c_1 = 0\\ no
mediator is then active, any selection is a false positive, and recall
is `NA`. For nonzero \\c_1\\ the two definitions agree under the
article's designs, whose exposure-on-mediator loadings are all nonzero.

The tuning grid governs these rates as much as the screen does. Every
replication searches the same grid: the package default
`seq(0.05, 10, length.out = 100)`, or the grid passed as `lambda_grid`,
and each fit keeps the best value the grid offers. A shorter or narrower
grid speeds a study up; the single-fit warning for a choice at an end of
the grid is not raised inside a study, which searches a fixed grid by
design. The reproduction vignette states the grids of the article's own
scripts.

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
for the global test,
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md)
for the per-method selected set of a single fit.

Other mediation simulation:
[`guo_calibration`](https://yelleknek.github.io/POEMED/reference/guo_calibration.md),
[`mediation_ar1_cov()`](https://yelleknek.github.io/POEMED/reference/mediation_ar1_cov.md),
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md),
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
pe_identification_study(n = 200, p = 60, outcome = "continuous",
                        pattern = "contrasting", c1_grid = c(0, 1),
                        n_rep = 5,
                        lambda_grid = seq(0.05, 2, length.out = 20))
#>  c1 method     fwer fdr  precision recall n_valid n_empty
#>  0  Bonferroni 0    0    0         0      5       0      
#>  0  BH         0    0    0         0      5       0      
#>  0  BY         0    0    0         0      5       0      
#>  1  Bonferroni 0    0    1         0.6    5       0      
#>  1  BH         0.2  0.08 0.92      0.7    5       0      
#>  1  BY         0    0    1         0.6    5       0      
#> 
#> Outcome: continuous
```
