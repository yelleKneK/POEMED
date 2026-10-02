# POEMED: Power Enhancement for High-Dimensional Mediation Analysis

POEMED implements powerful global tests for high-dimensional mediation
analysis, for continuous, binary, and count outcomes. It answers the
question *“is there any active mediator among a large set of candidate
mediators?”* and, when the answer is yes, reports which individual
mediators the screen selects. The name comes from POwer-Enhanced
MEDiation and is pronounced “POE-med”.

## Why Power Enhancement

The conventional high-dimensional mediation test targets the **total
indirect effect** `beta = Gamma_x' alpha_m`, the sum of every mediator’s
individual indirect effect. When some indirect effects are positive and
others negative they can cancel, making `beta = 0` even though active
mediators exist, and a test built on `beta` is then powerless. The
power-enhanced (PE) test adds a component `J_m` that accumulates the
marginal signal from each individual mediator. Because `J_m` sums
magnitudes, opposite-signed effects reinforce rather than cancel, so the
test stays powerful under heterogeneous and contrasting mediation. Under
the global null the enhancement is built to vanish: as the sample size
and the number of candidate mediators grow, the test’s Type I error rate
approaches the nominal level, as the benchmark’s does, and in finite
samples it rejects slightly more often than the benchmark.

## When to Use POEMED

POEMED answers the global question, *“is there any active mediator among
many candidates?”*, and then names the mediators it selects. Reach for
it when the candidate mediators number in the dozens or more (`p` may
exceed `n`), the outcome is continuous, binary, or a count, and the
individual mediation effects may differ in sign. For a single mediator,
or a few mediators fit as one structural model (indirect effects with
confidence intervals, likelihood ratio tests of arbitrary indirect
effects, moderated mediation), use
[DMAR](https://CRAN.R-project.org/package=DMAR) instead. Among
high-dimensional global tests, POEMED implements the power-enhanced
tests of Yu and Kelley (in press) and the HDMM benchmark they extend;
HILMA, GlobalTest, HDMT, and DACT are separate methods in their own
packages. The vignette `poemed-vs-competitors` shows how to run HILMA
and GlobalTest beside a POEMED fit and where to find HDMT and DACT. The
guided-tour vignette and
[`?pe_mediation`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
state the method’s practical limits (small `n`, strongly correlated
mediators, rare binary outcomes, overdispersed counts, missing data).

## Installation

Once POEMED is on CRAN:

``` r

install.packages("POEMED")
```

The development version, from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("yelleKneK/POEMED", build_vignettes = TRUE)
```

POEMED imports `glmnet` and `ncvreg`; `knitr` and `rmarkdown` build the
vignettes. After installing, `browseVignettes("POEMED")` opens the
guides.

## Quick Start

``` r

library(POEMED)

# A contrasting setting: active mediators whose indirect effects cancel,
# so the total indirect effect is zero.
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                             pattern = "contrasting", c1 = 1)

# One call reports the benchmark Wald test and the power-enhanced test
# side by side, and names the mediators the screen selected.
pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
```

The benchmark test (`pval_hdmm`) does not reject the (zero) total
indirect effect; the power-enhanced test (`pval_pe`) rejects decisively
and selects three of the four active mediators. The footer also reports
the tuning parameter the high-dimensional BIC chose and the grid it was
chosen from, `seq(0.05, 10, length.out = 100)` by default for every
outcome type; the grid can be shortened to save time or extended when
the value chosen sits at an end of it, which the fit warns about (see
[`?pe_mediation`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)).
Run as printed, this fit chooses 0.05, the smallest value, and warns; a
grid extended below it and spaced more finely
(`seq(0.01, 10, length.out = 200)`) chooses an interior value and
selects the same three mediators.

## Reproducing the Article

All data needed to rerun the article’s analyses ship with the package;
nothing is downloaded at run time. The vignette
**`reproducing-the-article`** reruns the empirical analysis cell by cell
and gives the calls that run the simulation designs; the table below
maps each article result to the function that runs its analysis.

| Article result | Function |
|----|----|
| Section 4.1, Figure 3: linear-model size and power | `pe_power_curve(outcome = "continuous", pattern = "homogeneous" / "contrasting")` |
| Section 4.2, Figures 4 and 5: logistic and Poisson size and power | `pe_power_curve(outcome = "binary", c2 = 1)` and `pe_power_curve(outcome = "count", c2 = 0.4)` |
| Supplement Section S.4.1, Tables S.1 to S.3: FWER/FDR mediator-identification study | [`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md) |
| Supplement Section S.4.4: real-data-motivated heterogeneous setting | [`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md) with the shipped `guo_calibration` constants |
| Section 5, Tables 1 (IMR) and 2 (PoU); supplement Tables S.10 (LEB), S.11 (U5MR), and S.12 (LBW) | `WHO_mediation_analysis("imr")`, `"pou"`, `"leb"`, `"u5mr"`, `"lbw"` |
| Supplement Section S.5, Tables S.8 to S.12: PE, PE-BH, and PE-BY columns | `WHO_mediation_analysis(..., full_table = TRUE)` |
| Where PE detects mediation the benchmark misses | [`summary()`](https://rdrr.io/r/base/summary.html) of a [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md) result |

The Monte Carlo studies use the article’s `n = 300`, `p = 500`, and
1,000 replications (the supplement’s real-data-motivated setting its
case-study `n = 85` and `p = 1008`), and each panel takes hours, so the
vignette does not rerun them; it gives the calls, on the grids of the
article’s scripts where the article records one and on the package
default grid for the real-data-motivated setting, with `cores > 1` to
parallelize (near-linear). The real-data analyses are deterministic:
under its default tuning grid,
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
reproduces the selected mediators of Tables 1, 2, and S.8 to S.12 (every
selected set but one BH set) and their benchmark p-values to the printed
four decimals, apart from a few cells that differ in the fourth decimal
or below 0.0001 and six cells it does not fit: four in which two
indicators are exactly collinear and two in which the outcome is
constant. The reproduction vignette names those six cells and shows how
to fit them by hand. The package’s test suite pins every global, region,
and income row of Tables 1 and 2 (`tests/testthat/test-benchmark.R`) and
checks every populated cell of Tables 1, 2, and S.8 to S.12 against the
printed values (`tests/testthat/test-who_article_tables.R`).

``` r

# The global, region, and income models of the article's infant-mortality
# analysis (Table 1); add "region_income" to groupings for the
# region-by-income cells.
WHO_mediation_analysis("imr", groupings = c("global", "region", "income"))
```

## Vignettes

``` r

browseVignettes("POEMED")
```

- **`POEMED`**: a guided tour of the method and the package.
- **`reproducing-the-article`**: the empirical analysis rerun cell by
  cell, and the calls for the simulation designs.
- **`poemed-vs-competitors`**: a head-to-head against the benchmark
  test.

## Main Functions

| Function | Purpose |
|----|----|
| [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md) | Front end; dispatches on outcome type (continuous / binary / count) |
| [`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md) / `_logistic()` / `_poisson()` | Outcome-specific workers |
| [`pe_mediate()`](https://yelleknek.github.io/POEMED/reference/pe_mediate.md) | Data-frame / formula front end |
| [`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md) / [`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md) | Per-mediator screen detail; PE statistic and selected set per multiplicity method |
| [`plot()`](https://rdrr.io/r/graphics/plot.default.html) / [`tidy()`](https://yelleknek.github.io/POEMED/reference/poemed_broom.md) / [`glance()`](https://yelleknek.github.io/POEMED/reference/poemed_broom.md) / [`summary()`](https://rdrr.io/r/base/summary.html) | Plot, broom verbs, and a cross-group comparison digest |
| [`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md) | Generate data under homogeneous / contrasting patterns |
| [`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md) / `guo_calibration` | The real-data-motivated heterogeneous setting and its constants |
| [`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md) / [`pe_simulation_study()`](https://yelleknek.github.io/POEMED/reference/pe_simulation_study.md) | Monte Carlo size and power study |
| [`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md) | Monte Carlo FWER / FDR / precision / recall of mediator identification |
| [`ss_power_pe_mediation()`](https://yelleknek.github.io/POEMED/reference/ss_power_pe_mediation.md) | Sample-size planning for a target power |
| [`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md) / [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md) | Build / run the benchmark health-expenditure analysis |

Every test reports the benchmark Wald test and the power-enhanced test
side by side, plus the estimated total indirect effect. No confidence
interval is reported: the article reports none, and the Wald interval
built on the penalized estimate under-covers in finite samples.
`report_all_methods = TRUE` adds the selected set and PE p-value under
all three multiplicity methods, `drop_constant = TRUE` tolerates
zero-variance mediators, and `cores` parallelizes the Monte Carlo and
benchmark routines.

## Data

- `WHO_health_mediation`: the empirical panel analyzed in the article:
  91 World Health Organization member states over 2000–2021, with GDP
  per capita growth as the exposure, 57 WHO Global Health Expenditure
  Database indicators as candidate mediators, and five World Bank health
  outcomes (infant and under-five mortality, life expectancy, low
  birthweight, undernourishment).
- `WHO_indicator_codebook`: definitions of the 57 indicators.
- `guo_calibration`: constants for the real-data-motivated heterogeneous
  simulation.
- `example_continuous`, `example_binary`, `example_count`: small
  simulated data sets, one per outcome type.

## Reference

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

If you use POEMED in published work, please cite the article. Run
`citation("POEMED")` for the BibTeX entry.

## Authors

- **Xiufan Yu**, Department of Applied and Computational Mathematics and
  Statistics, University of Notre Dame.
- **Ken Kelley** (maintainer), Department of IT, Analytics, and
  Operations, University of Notre Dame (<kkelley@nd.edu>).

## License

GPL (\>= 3).
