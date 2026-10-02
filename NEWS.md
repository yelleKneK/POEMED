# POEMED 1.0.0

This is the first release of POEMED (Power Enhancement for High-Dimensional
Mediation Analysis), the companion package to Yu and Kelley (in press),
"Power Enhancement in High-Dimensional Heterogeneous Mediation Analysis," in
the *Journal of the American Statistical Association*. The name comes from
POwer-Enhanced MEDiation and is pronounced "POE-med".

POEMED asks whether any of a large set of candidate mediators carries the
effect of an exposure to an outcome and, when the answer is yes, which ones
do. The benchmark Wald test of Guo et al. (2022, 2023, 2024) tests the total
indirect effect, the sum of the individual indirect effects, so it loses its
power when active mediators act in opposite directions and their effects
cancel. The power-enhanced test of Yu and Kelley (in press), which follows
the power enhancement principle of Fan, Liao, and Yao (2015), adds a
component that accumulates the magnitudes of the marginal signals of the
mediators that pass a significance screen, so it detects active mediators
whether their effects are homogeneous, heterogeneous, or exactly
contrasting. The package works for continuous, binary, and count outcomes.

## Testing for Mediation

* `pe_mediation()` is the front end. Give it the exposure, the outcome, the
  candidate mediators, and any confounders, name the outcome type
  (`"continuous"`, `"binary"`, or `"count"`), and it reports the benchmark
  Wald test and the power-enhanced test side by side, the estimated total
  indirect effect, and the mediators the screen selected. One or several
  exposures may be analyzed.
* `pe_mediation_linear()`, `pe_mediation_logistic()`, and
  `pe_mediation_poisson()` are the outcome-specific workers that
  `pe_mediation()` dispatches to, and they can be called directly.
* `pe_mediate()` is a formula and data-frame front end for those who prefer
  naming columns to building matrices. It analyzes the complete cases and
  reports how many rows it dropped and which columns held the missing
  values; `pe_mediation()` and its workers take complete data.
* In the outcome model the mediators are penalized and the exposures and
  confounders are not, and the tuning parameter is chosen by the
  high-dimensional BIC, HBIC (Wang, Kim, and Li, 2013). For continuous and
  count outcomes the penalty is SCAD (Fan and Li, 2001), fit by the
  one-step local linear approximation (Zou and Li, 2008; Fan, Xue, and Zou,
  2014); for a binary outcome the fit is the partial penalized likelihood
  estimator of Guo et al. (2024), computed by coordinate descent with an
  adaptively rescaled SCAD penalty (see `?pe_mediation_logistic`).
* Recoding a mediator, or a continuous or binary outcome, in the opposite
  direction leaves both test statistics, both p-values, and the selected
  set unchanged; signed estimates, such as the total indirect effect when
  the outcome is reversed, change only in sign.
* The "When to Use POEMED" section of `?pe_mediation` states the method's
  assumptions, shows its practical limits on simulated designs (small
  samples, strongly correlated mediators, selection in finite samples, rare
  binary outcomes, and overdispersed counts), and points to other tools for
  a single mediator or a few. The same page states that the power-enhanced
  test controls its size asymptotically and that, in finite samples, it can
  reject somewhat more often than the benchmark, by up to about
  `error_level / log(log(n))`.
* `example_continuous`, `example_binary`, and `example_count` are small data
  sets simulated with `simulate_mediation_data()`, one for each outcome
  type, each stored with the active mediators that generated it and ready
  for a first fit (see `?example_data`).

## Identifying the Active Mediators

* The screen inside the power enhancement component selects individual
  mediators at the level set by `error_level`, with familywise error rate
  control (`method = "Bonferroni"`, the default) or false discovery rate
  control (`method = "BH"` for Benjamini-Hochberg or `method = "BY"` for
  Benjamini-Yekutieli). Like the size of the test, this control holds
  asymptotically.
* `pe_mediators()` returns the detail behind a fit for each mediator the
  penalized fit kept: its two path statistics (mediator to outcome and
  exposure to mediator), its screening p-value, and whether the screen
  selected it.
* `pe_selection()` reports the power-enhanced statistic, its enhancement
  component `J_m`, the global p-value, and the selected set under the
  method the fit used, or under all three multiplicity methods side by side
  when the fit was made with `report_all_methods = TRUE`.

## Tuning

* The default tuning grid, `seq(0.05, 10, length.out = 100)` for every
  outcome family, is the default value of `lambda_grid`, and any other grid
  can be passed in its place. HBIC keeps the best value the grid offers, so
  the grid can be suited to the task: shortened to save time when many fits
  are run, extended or refined when the chosen value sits at an end of the
  grid.
* Every fit records the grid it searched and the value HBIC chose in a
  `"tuning"` attribute and reports that value when it prints. A fit whose
  choice sits at either end of the grid warns (class
  `poemed_grid_boundary`) and suggests extending the grid and using a finer
  partition. A fit whose penalized model keeps no mediator warns (class
  `poemed_empty_fit`), says so when it prints, and is reported as a
  non-rejection (both tests 0, with p-values of 1), which is expected when
  no mediator is active but is not evidence that none is.

## Simulation, Power, and Planning

* `simulate_mediation_data()` generates data under the homogeneous and
  contrasting mediation patterns of the article for all three outcome
  types, and `mediation_ar1_cov()` builds the autoregressive mediator
  covariance the simulations use.
* `simulate_guo_mediation()`, with the `guo_calibration` constants, generates
  the article's real-data-motivated heterogeneous setting: a continuous
  outcome with eleven active mediators of mixed sign, calibrated to the
  DNA-methylation study of Guo et al. (2022).
* `pe_power_curve()` runs the Monte Carlo size and power study across a
  grid of signal strengths, and `pe_simulation_study()` runs it across
  mediation patterns.
* `pe_identification_study()` estimates how well the screen recovers the
  individual active mediators: the familywise error rate, false discovery
  rate, precision, and recall of the selected set under each multiplicity
  method.
* `ss_power_pe_mediation()` plans the sample size that reaches a target
  power by simulation.
* The Monte Carlo functions run in parallel through `cores` on Unix-like
  systems such as macOS and Linux (serially on Windows), report the number
  of empty fits at each design point, and take a `seed` that reproduces a
  run exactly at the same value of `cores` while leaving the caller's
  random-number state as it was.

## The WHO Health-Expenditure Analysis

* `WHO_health_mediation` holds the article's empirical data, a panel of 91
  World Health Organization members from 2000 to 2021 with economic growth
  as the exposure, 57 health-expenditure indicators as candidate mediators,
  and five health outcomes; `WHO_indicator_codebook` defines the
  indicators.
* `WHO_mediation_design()` assembles the exposure, outcome, mediator, and
  confounder matrices for any outcome and any region or income group.
* `WHO_mediation_analysis()` runs the article's subgroup analyses (global,
  by region, by income group, and by region and income group) with the
  preprocessing and the tuning grid, `seq(0.1, 2, length.out = 100)`,
  behind the article's Tables 1 and 2 and supplement Tables S.8 to S.12. It
  reproduces those tables apart from a few cells that
  `?WHO_mediation_analysis` lists and explains, and the reproduction
  vignette shows how to fit by hand the six cells that the function reports
  as `NA`. With `full_table = TRUE` it reports the power-enhanced p-value
  and the selected mediators under each multiplicity method.

## Results, Printing, and Plotting

* Results are tidy data frames. A fit from `pe_mediation()`, its workers, or
  `pe_mediate()` is a long table with a `term` column and a numeric `value`
  column; the tables from `pe_mediators()`, `pe_selection()`, the Monte
  Carlo and planning functions, and `WHO_mediation_analysis()` are wide,
  with one row per mediator, multiplicity method, design point, sample
  size, or group. Values are stored at full precision and rounded only for
  display: whole numbers print without decimals, p-values print to four
  decimals with `< 0.0001` below that, and other values print to four
  significant figures (`options(poemed.digits = 6)` shows more).
* Column names travel to the output: with a named `M`, the printed fit and
  the mediator tables name the mediators, and with several named exposures
  the per-exposure rows carry the exposures' names.
* `plot()` draws the power curves of a `pe_power_curve()` result, the
  benchmark and power-enhanced rejection rates against signal strength,
  and, for a fit, the screening evidence for each mediator the penalized
  fit kept, with the selected ones highlighted. `summary()` of a grouped
  table such as the output of `WHO_mediation_analysis()` counts the groups
  in which each test detects mediation, lists those in which the
  power-enhanced test detects mediation that the benchmark misses, and
  names the mediators selected most often. The broom verbs `tidy()` and
  `glance()` are re-exported so they work after `library(POEMED)` alone.

## Vignettes

* Three vignettes accompany the package: a guided tour,
  `vignette("POEMED")`; a head-to-head comparison of the power-enhanced test
  with the benchmark Wald test across mediation patterns, outcome types, and
  the WHO data, closing with where to find HILMA, GlobalTest, HDMT, and
  DACT, `vignette("poemed-vs-competitors")`; and a guide to rerunning the
  article's simulation designs and empirical analysis,
  `vignette("reproducing-the-article")`. `browseVignettes("POEMED")` lists
  them.
* An installation from CRAN includes the vignettes. An installation from
  GitHub includes them only when they are built during installation, with
  `remotes::install_github("yelleKneK/POEMED", build_vignettes = TRUE)`;
  without them, `vignette()` finds nothing to open. All three can also be
  read on the package website, <https://yelleknek.github.io/POEMED/>, where
  the guided tour is the "Get started" page and the other two are under
  Articles.

## Validation

* The test suite checks the fits of the three example data sets against
  stored reference values, checks the article's Tables 1, 2, and S.8 to
  S.12 cell by cell, and confirms that every function with a `seed`
  argument reproduces its results and leaves the caller's random-number
  state as it was.
