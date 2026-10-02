# Reproducing the Article's Analyses With POEMED

This vignette accompanies the article “Power Enhancement in
High-Dimensional Heterogeneous Mediation Analysis” (Yu and Kelley, in
press at the *Journal of the American Statistical Association*). Its
second part reruns the article’s empirical data analysis (its Section 5)
with POEMED and reads the package’s results against the printed tables,
cell by cell. Its first part describes the article’s simulation designs
(its Section 4, and Sections S.4.1 and S.4.4 of its supplement) and
gives the calls that run each of them at the article’s scale, on the
grids of the article’s own scripts where the article records one. The
companion vignette *Power Enhancement for High-Dimensional Mediation: A
Guided Tour* introduces the method itself.

**Scope of the comparison.** The article benchmarks the power-enhanced
tests (PE-HDMM for linear outcomes, PE-HDGMM for generalized outcomes)
against several competitors: HDMM, HILMA, and GlobalTest. POEMED
implements the proposed PE tests and the **HDMM** benchmark (the
total-indirect-effect Wald test that the PE component augments). Every
POEMED result therefore reports the central contrast of the article’s
figures and tables, **PE versus HDMM**; HILMA and GlobalTest are
separate methods outside this package.

Load the package and fix the random-number seed before running any of
the code below:

``` r

library(POEMED)
set.seed(113)
```

## Part I. The Article’s Simulation Designs (Article Section 4)

The simulation study asks two questions of every test. Under the global
null of no active mediator (`c1 = 0`) the rejection rate estimates the
**Type I error rate**, which should sit near the nominal level,
`alpha_level = 0.05`. Away from the null (`c1 != 0`) the rejection rate
estimates **power**.
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
sweeps `c1` and returns both the benchmark (`rejection_hdmm`) and
power-enhanced (`rejection_pe`) rejection rates at each grid point, and
[`pe_simulation_study()`](https://yelleknek.github.io/POEMED/reference/pe_simulation_study.md)
runs several patterns at once and stacks the results.

The article’s designs use `n = 300`, `p = 500`, and **1,000
replications** per design point, and its scripts search a tuning grid
set by hand for each outcome family: 20 values from 0.2 to 0.39 for the
linear model, 15 values from 0.04 to 0.2 for the logistic model, and 15
values from 0.7 to 5 for the Poisson model. (The linear script searched
a second grid, 20 values from 0.27 to 0.46, for the reduced model that
the benchmark test needs; the package searches one grid for both fits,
so the calls below pass the full-model grid.) The calls below for
Sections 4.1, 4.2, and S.4.1 pass those grids as `lambda_grid`, so they
run the article’s designs as its scripts did, whatever the package’s
default grid (100 values from 0.05 to 10); the article records no grid
for the real-data-motivated setting of Section S.4.4, whose call uses
the package default. The chunks in this part are shown but not run: a
linear or logistic panel takes one to a few hours of computing at this
scale, and the Poisson homogeneous panel days. Set `cores` above 1 to
parallelize; a seeded run repeats exactly at the same number of cores.

### Section 4.1: Linear Mediation Models

The article studies a univariate exposure with `n = 300` continuous
outcomes and `p = 500` candidate mediators whose noise has an
autoregressive `0.5^|i-j|` covariance. The mediator-exposure coefficient
is `Gamma_x = c1 * tau`, where `tau` is 0.1, 0.2, 0.3, 0.4, 0.5 for the
first five mediators and an independent N(0, 0.5^2) draw for each of the
others, the direct effect is `c2 = 0.5`, and two outcome-mediator
patterns are contrasted:

- **Setting (i), homogeneous:**
  `alpha_m = (1, 0.8, 0.6, 0.4, 0.2, 0, ...)`, so all active mediators
  push the same way and the total indirect effect is `beta = 0.7 * c1`
  (nonzero). This is the regime where ordinary tests already work; the
  question is whether PE *adds* power.
- **Setting (ii), contrasting:**
  `alpha_m = (1, -0.5, 0.4, -0.3, 0, ...)`, four active mediators whose
  effects cancel so the total indirect effect is `beta = 0` for every
  `c1`. This is the regime that defeats a total-indirect-effect test,
  and where PE is designed to rescue power.

POEMED maps these to `pattern = "homogeneous"` and
`pattern = "contrasting"` in
[`simulate_mediation_data()`](https://yelleknek.github.io/POEMED/reference/simulate_mediation_data.md)
and
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md).
The article’s Figure 3 sweeps `c1` from -1 to 1 (the symmetric design
fills in the negative half), one panel per pattern:

``` r

c1_grid <- seq(-1, 1, by = 0.1)
lm_study <- pe_simulation_study(n = 300, p = 500, outcome = "continuous",
                                patterns = c("homogeneous", "contrasting"),
                                c1_grid = c1_grid, c2 = 0.5, n_rep = 1000,
                                lambda_grid = seq(0.2, 0.39, length.out = 20),
                                cores = 4, seed = 113)
```

The result has one row per pattern and `c1`, with the two rejection
rates and the number of replications whose penalized fit selected no
mediator (`n_empty`).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
[`pe_power_curve()`](https://yelleknek.github.io/POEMED/reference/pe_power_curve.md)
result draws the two curves against `c1`, as the guided-tour vignette
shows.

#### The Real-Data-Motivated Heterogeneous Setting (Supplement Section S.4.4)

The supplement adds a third, harder pattern in its Section S.4.4: a
setting calibrated to the DNA-methylation case study of Guo et al.
(2022), with **eleven active mediators of mixed sign** that neither all
agree (homogeneous) nor exactly cancel (contrasting).
[`simulate_guo_mediation()`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md)
generates it from the calibration constants shipped as
`guo_calibration`; the total indirect effect is
`guo_calibration$beta_per_c1 * c1`, about `-1.597 * c1`. (The article
reports `-1.5977`, computed from the full-precision Guo coefficients;
the package ships those coefficients rounded to the three decimals
printed in the supplement, which give `-1.597`.)

``` r

d <- simulate_guo_mediation(c1 = 1, seed = 113)
dim(d$M)                              # 85 x 1008, the case-study dimensions
#> [1]   85 1008
d$active_mediators                    # the 11 active loci
#>  [1]  1  2  3  4  5  6  7  8  9 10 11
round(d$beta, 4)                      # beta_per_c1 * c1, i.e. about -1.597
#> [1] -1.597
```

The supplement’s Figure S.14 sweeps `c1` over
`c(0, +/- 0.1, ..., +/- 1)` at the case-study size `n = 85`, without the
calibrated confounders in its panel (a) and with them in its panel (b),
averaged over 1,000 replications. The article records no tuning grid for
this setting, so the fits below use the package default (the example on
[`?simulate_guo_mediation`](https://yelleknek.github.io/POEMED/reference/simulate_guo_mediation.md)
passes a shorter grid to stay fast and selects the same loci as the
default):

``` r

rej <- function(c1, confounders) {
  ps <- replicate(1000, {
    d <- simulate_guo_mediation(c1 = c1, confounders = confounders)
    f <- pe_mediation(d$X, d$Y, d$M, Z = d$Z, outcome = "continuous")
    c(f$value[f$term == "pval_hdmm"], f$value[f$term == "pval_pe"])
  })
  rowMeans(ps <= 0.05)
}
guo_without <- sapply(seq(-1, 1, by = 0.1), rej, confounders = FALSE)
guo_with <- sapply(seq(-1, 1, by = 0.1), rej, confounders = TRUE)
```

### Section 4.2: Mediation Models With Generalized Linear Outcomes

The same construction is repeated for non-continuous outcomes through a
link function: a **logistic** model for binary outcomes (PE-HDGMM versus
HDGMM) and a **Poisson** model for counts. Again `n = 300`, `p = 500`,
with a homogeneous and a contrasting pattern apiece; the article fixes
the direct effect at `c2 = 1` for the logistic models and `c2 = 0.4` for
the Poisson models. The calls below run the four panels of the article’s
Figures 4 and 5. The Poisson homogeneous fits take tens of seconds each
at this design, so that panel alone needs well over 100 hours of
computing (days, even on 4 cores); each of the other three panels needs
a few hours.

``` r

c1_grid <- seq(-1, 1, by = 0.1)
# Logistic: c2 = 1
glm_logit <- pe_simulation_study(n = 300, p = 500, outcome = "binary",
                                 patterns = c("homogeneous", "contrasting"),
                                 c1_grid = c1_grid, c2 = 1, n_rep = 1000,
                                 lambda_grid = seq(0.04, 0.2, length.out = 15),
                                 cores = 4, seed = 113)
# Poisson: c2 = 0.4
glm_pois <- pe_simulation_study(n = 300, p = 500, outcome = "count",
                                patterns = c("homogeneous", "contrasting"),
                                c1_grid = c1_grid, c2 = 0.4, n_rep = 1000,
                                lambda_grid = seq(0.7, 5, length.out = 15),
                                cores = 4, seed = 113)
```

### Supplement Section S.4.1: Identifying Individual Mediators (FWER and FDR)

Beyond the global test, the article’s supplement studies how well the PE
machinery **identifies which** individual mediators are active, and
whether it controls the familywise error rate (FWER) or false discovery
rate (FDR). POEMED exposes this through the `method` argument
(`"Bonferroni"` for FWER, `"BH"` or `"BY"` for FDR) and the
`report_all_methods` switch, which reports all three side by side as the
supplement does; for a single fit,
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md)
shows the selected set and the power-enhanced p-value under each method
(the guided-tour vignette shows one).

[`pe_identification_study()`](https://yelleknek.github.io/POEMED/reference/pe_identification_study.md)
runs the supplement’s identification *study*: it scores each method’s
selections against the truth set across a grid of `c1` and reports the
empirical FWER, FDR, precision, and recall, with the supplement’s
definitions. A mediator counts as truly active when its outcome
coefficient is nonzero, at every `c1` including `c1 = 0` (five mediators
under the homogeneous pattern, four under the contrasting pattern); pass
`truth = "mediation"` for the both-paths definition, under which nothing
is active at `c1 = 0` and `recall` is `NA` there. The supplement’s
Tables S.1 to S.3 come from the calls below, one per pattern:

``` r

id_study <- lapply(c("homogeneous", "contrasting"), function(pattern)
  pe_identification_study(n = 300, p = 500, outcome = "continuous",
                          pattern = pattern, c1_grid = seq(0, 1, by = 0.1),
                          n_rep = 1000,
                          lambda_grid = seq(0.2, 0.39, length.out = 20),
                          cores = 4, seed = 113))
```

## Part II. Empirical Data Analysis (Article Section 5)

The article asks whether health-care expenditure mediates the
relationship between economic growth and population health across the
world’s countries. POEMED ships that exact data set,
`WHO_health_mediation`: a country-by-year panel (2000–2021) with

- **Exposure** `X`: annual growth in GDP per capita (World Bank),
- **Mediators** `M`: 57 health-expenditure indicators (WHO Global Health
  Expenditure Database),
- **Outcomes** `Y`: five health outcomes, infant mortality (`imr`),
  under-five mortality (`u5mr`), life expectancy (`leb`), low
  birthweight (`lbw`), and undernourishment (`pou`),

with WHO region and World Bank income group available as the grouping
and confounding variables.

``` r

data(WHO_health_mediation)
dim(WHO_health_mediation)
#> [1] 2002   68
table(WHO_health_mediation$region)
#> 
#>  AFR  AMR  EMR  EUR SEAR  WPR 
#>  704  616  176  198  110  198
```

### Section 5, Table 1: Infant Mortality (IMR)

[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
builds a table in the layout of the article’s per-outcome tables: it
fits the global model (region, income, and year as confounders), then a
separate model within each WHO region (income and year as confounders)
and each income group (region and year as confounders), and reports the
HDMM and PE p-values for the total indirect effect together with the
mediators the PE component selects.

The function searches one tuning grid, `seq(0.1, 2, length.out = 100)`,
for the global all-country model and for every subgroup model alike. It
is the grid behind the article’s real-data tables (the article’s text
gives no grid; it states only that the tuning parameter is chosen by the
HBIC criterion), and it is the default of both `lambda_grid`, for the
subgroup models, and `lambda_grid_global`, for the global model, so no
explicit grid arguments are needed.

``` r

imr <- WHO_mediation_analysis("imr", groupings = c("global", "region", "income"))
imr[, c("group", "n_countries", "pval_hdmm", "pval_pe", "selected_mediators")]
#>  group        n_countries pval_hdmm pval_pe  selected_mediators
#>  ALL          91          0.0644    < 0.0001 gge_gdp           
#>  AFR          32          0.4434    0.4434   none              
#>  AMR          28          0.4972    0.4972   none              
#>  EMR          8           0.0261    0.0261   none              
#>  EUR          9           < 0.0001  < 0.0001 none              
#>  SEAR         5           0.6483    < 0.0001 ext_usd2021_pc    
#>  WPR          9           0.8681    < 0.0001 chi_che           
#>  Low          16          0.3091    < 0.0001 pvtd_usd2021      
#>  Lower-middle 27          0.9543    0.9543   none              
#>  Upper-middle 26          0.4481    0.4481   none              
#>  High         22          0.0389    < 0.0001 oops_che, shi_che 
#> 
#> Outcome: imr
```

Read this against the article’s Table 1. The highlights:

- **Global (ALL):** general government expenditure as a share of GDP
  (`gge_gdp`) is the mediator the screen selects; HDMM is borderline
  (`pval_hdmm` = 0.064) while PE is essentially zero.
- **SEAR (South-East Asia):** HDMM sees nothing (`pval_hdmm` = 0.65)
  while PE detects external health expenditure per capita in constant
  2021 US dollars (`ext_usd2021_pc`).
- **WPR (Western Pacific):** HDMM sees nothing (`pval_hdmm` = 0.87)
  while PE detects compulsory health insurance as a share of current
  health expenditure (`chi_che`).
- **Low income:** HDMM sees nothing (`pval_hdmm` = 0.31) while PE
  detects domestic private health expenditure (`pvtd_usd2021`).
- **EMR (Eastern Mediterranean), EUR (Europe), and high income:** the
  benchmark rejects (`pval_hdmm` = 0.026, 0.000094, and 0.039; Table 1
  prints the EUR value as 0.0001). In the high-income group PE
  identifies out-of-pocket spending and social health insurance as
  shares of current health expenditure (`oops_che`, `shi_che`); in EMR
  and EUR no single indicator passes the screen, so the PE p-value
  equals the benchmark’s.
- **AFR, AMR, lower-middle income, and upper-middle income:** neither
  test detects mediation.

In the global model, SEAR, WPR, and the low-income group, heterogeneous
mediation hides the signal from the conventional test and the PE
component recovers it. Among the income groups the strongest signals
appear in the low- and high-income groups, as Yu and Kelley (in press)
report in Table 1.

The table above agrees with the article’s Table 1 row for row. Every
group’s identified set is the one the table lists, and every benchmark
p-value matches the printed value to four decimals except SEAR, which is
0.6483 here where the table prints 0.6482. Where the table prints a PE
p-value below `1e-12`, the package prints `< 0.0001`, its display floor;
the values behind those entries are below `1e-12` as well.

[`summary()`](https://rdrr.io/r/base/summary.html) digests the whole
table, highlighting the groups where the power-enhanced test detects
mediation the benchmark misses (the article’s Table 1 sets those four
rows in bold):

``` r

summary(imr)
#> POEMED comparison: IMR
#>   groups: 11 (11 with data), alpha_level = 0.05
#>   detected by benchmark (HDMM): 3   by power-enhanced (PE): 7
#>   PE detects mediation in 4 groups the benchmark misses:
#>  group n_countries pval_hdmm  pval_pe selected_mediators
#>    ALL          91    0.0644 < 0.0001            gge_gdp
#>   SEAR           5    0.6483 < 0.0001     ext_usd2021_pc
#>    WPR           9    0.8681 < 0.0001            chi_che
#>    Low          16    0.3091 < 0.0001       pvtd_usd2021
#>   most-flagged mediators:
#>     chi_che            1
#>     ext_usd2021_pc     1
#>     gge_gdp            1
#>     oops_che           1
#>     pvtd_usd2021       1
#>     shi_che            1
```

The supplement also gives an **extended** version of each table (its
Table S.8 for IMR) that breaks the power-enhanced result out by
multiplicity method. Pass `full_table = TRUE` to get that layout, with
the PE p-value and selected mediators reported separately under
Bonferroni (FWER), BH, and BY (FDR):

``` r

WHO_mediation_analysis("imr", groupings = "global", full_table = TRUE)
#>  grouping group n_countries pval_hdmm pval_pe_bonferroni selected_bonferroni
#>  global   ALL   91          0.0644    < 0.0001           gge_gdp            
#>  pval_pe_bh selected_bh      pval_pe_by selected_by
#>  < 0.0001   chi_che, gge_gdp < 0.0001   gge_gdp    
#> 
#> Outcome: imr
```

### The Other Outcomes, and the Region-by-Income Cells

The remaining article outcomes (undernourishment in the article’s Table
2, and life expectancy, under-five mortality, and low birthweight in the
supplement’s Tables S.10 to S.12) run identically; only the outcome code
changes:

``` r

for (y in c("u5mr", "leb", "lbw", "pou"))
  print(WHO_mediation_analysis(y, groupings = c("global", "region", "income")))
```

The article also drills into the **region-by-income** cells.
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
returns the region-by-income block of the table in one call, with `NA`
for the cells that have no countries or whose fit fails:

``` r

imr_cells <- WHO_mediation_analysis("imr", groupings = "region_income")
```

To look inside one cell, build it with
[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md),
which assembles the exposure, outcome, mediator, and confounder blocks
for a chosen subset, and fit it the way
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
does. A small cell can contain indicators that never vary inside it.
Such a column carries no signal and cannot be standardized, so the code
below drops it by hand before standardizing, as the function does, and
the print footer then names the selected mediators by their indicator
codes. The code also passes the function’s default grid, so this fit is
the AFR-Upper-middle row of the block above:

``` r

cell <- WHO_mediation_design("imr", region = "AFR", income = "Upper-middle")
c(n_countries = cell$n_countries, mediators = ncol(cell$M))
#> n_countries   mediators 
#>           4          57

keep <- apply(cell$M, 2L, sd) > 0          # FALSE for a constant indicator
cell$indicators[!keep]
#> [1] "chi_pvt_che"
fit <- pe_mediation(scale(cell$X), cell$Y - mean(cell$Y),
                    scale(cell$M[, keep, drop = FALSE]),
                    Z = cell$Z, outcome = "continuous", scale = FALSE,
                    lambda_grid = seq(0.1, 2, length.out = 100))
fit[fit$term %in% c("pval_hdmm", "pval_pe"), ]
#>  term      value   
#>  pval_hdmm 0.0148  
#>  pval_pe   < 0.0001
#> 
#> Outcome model: continuous (linear)
#> Selected mediators (1): pvtd_ncu2021_pc
#> Tuning parameter (HBIC): lambda = 0.138 from 100 values in [0.1, 2]
```

The selected mediator, domestic private health expenditure per capita in
constant 2021 national currency units (`pvtd_ncu2021_pc`), is the one
the article’s Table 1 lists for this cell, and the benchmark p-value,
0.0148, is the value the table prints.

### Six Cells Reported as NA

Across the region-by-income blocks of the article’s tables for its five
outcomes (Table 1, Table 2, and Tables S.8 to S.12), six cells that the
article prints with a value come back from
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
as `NA`, each with a warning that names the cell and quotes the reason.
They are of two kinds.

**Four rank-deficient cells.** In the EMR-Lower-middle cells of IMR and
U5MR and the AFR-Upper-middle cells of LEB and U5MR, two of the 57
indicators are exactly collinear: domestic and external health
expenditure as shares of current health expenditure (`dom_che`,
`ext_che`) sum to 100 in every country-year of those cells. The
AFR-Upper-middle cells add a second identity, since compulsory and
social health insurance as shares of current health expenditure
(`chi_che`, `shi_che`) coincide there (compulsory private insurance,
`chi_pvt_che`, is zero throughout, which is why it was dropped as
constant above). The benchmark test inverts the covariance matrix of the
mediators the penalized fit selected, and when that selection keeps both
members of a pair the matrix is singular; the fit stops with the error
that the system is computationally singular, and the function reports
the cell as `NA`. (The same two identities hold in the IMR cell above,
where the selection keeps `ext_che` and `chi_che` without their
partners, so the test goes through.) In three of the four hand fits
below HBIC chooses 0.1, the smallest value of the article’s grid, for
the full model, the reduced model, or both; the fits keep that grid as
the article’s analysis did, so the warning a single
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
fit gives at an end of its grid is not shown in this section. The
article’s values come back when the collinearity is removed before the
fit, as the article’s analysis did: it dropped `dom_che` and `vfa_che`
(near-exact partners of `ext_che` and `cfa_che`) in the EMR-Lower-middle
cells, and `shi_che` and `ext_che` in the AFR-Upper-middle cells:

``` r

cell <- WHO_mediation_design("imr", region = "EMR", income = "Lower-middle")
range(cell$M[, "dom_che"] + cell$M[, "ext_che"])   # 100 in every row
#> [1] 100 100

fit_cell <- function(outcome, region, income, drop) {
  cell <- WHO_mediation_design(outcome, region = region, income = income)
  keep <- apply(cell$M, 2L, sd) > 0 & !(cell$indicators %in% drop)
  fit <- pe_mediation(scale(cell$X), cell$Y - mean(cell$Y),
                      scale(cell$M[, keep, drop = FALSE]),
                      Z = cell$Z, outcome = "continuous", scale = FALSE,
                      lambda_grid = seq(0.1, 2, length.out = 100))
  selected <- cell$indicators[keep][attr(fit, "selected_mediators")]
  data.frame(outcome = outcome, cell = paste(region, income, sep = "-"),
             dropped = paste(drop, collapse = ", "),
             pval_hdmm = fit$value[fit$term == "pval_hdmm"],
             selected = if (length(selected)) paste(selected, collapse = ", ")
                        else "none")
}
rd <- rbind(fit_cell("imr",  "EMR", "Lower-middle",
                     drop = c("dom_che", "vfa_che")),
            fit_cell("u5mr", "EMR", "Lower-middle",
                     drop = c("dom_che", "vfa_che")),
            fit_cell("leb",  "AFR", "Upper-middle",
                     drop = c("shi_che", "ext_che")),
            fit_cell("u5mr", "AFR", "Upper-middle",
                     drop = c("shi_che", "ext_che")))
rd
#>   outcome             cell          dropped  pval_hdmm        selected
#> 1     imr EMR-Lower-middle dom_che, vfa_che 0.08372121            none
#> 2    u5mr EMR-Lower-middle dom_che, vfa_che 0.09144904            none
#> 3     leb AFR-Upper-middle shi_che, ext_che 0.13114741            none
#> 4    u5mr AFR-Upper-middle shi_che, ext_che 0.03894879 pvtd_ncu2021_pc
```

The four benchmark p-values, 0.0837, 0.0914, 0.1311, 0.0389, are the
values the article prints for these cells (Table 1 for IMR, Table S.10
for LEB, and Table S.11 for U5MR), and the U5MR AFR-Upper-middle cell
recovers the mediator Table S.11 lists, `pvtd_ncu2021_pc`.

**Two constant-outcome cells.** In the EUR-High and WPR-High cells of
the undernourishment analysis (PoU, the article’s Table 2), the outcome
is 2.5 in every country-year of the cell: the World Bank publishes a
prevalence of undernourishment below 2.5 percent as 2.5, its lower
reporting bound, and the high-income countries of those two regions are
all at that bound. There is no variation in the outcome to mediate, so
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
refuses the fit with the message that `Y` is constant, and the function
reports the two cells as `NA`. Table 2 prints a p-value of 1 for both.

``` r

sapply(c("EUR", "WPR"), function(r)
  unique(WHO_mediation_design("pou", region = r, income = "High")$Y))
#> EUR WPR 
#> 2.5 2.5
```

### Why a Mediator? Opening Up the Screening Step

For any fit,
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md)
shows the evidence the PE component acted on: the two path statistics,
the screening p-value (the larger of the two path p-values), and the
`selected` flag, for every candidate the penalized fit selected. This is
the per-mediator detail behind the last column of the article’s tables.
In this global model HBIC picks 0.1, the smallest value of the grid, as
it does in about half of the article’s cells; the warning a single fit
gives at an end of its grid is not shown here either.

``` r

des <- WHO_mediation_design("imr")                  # global IMR design
gfit <- pe_mediation(scale(des$X), des$Y - mean(des$Y), scale(des$M),
                     Z = des$Z, outcome = "continuous", scale = FALSE,
                     lambda_grid = seq(0.1, 2, length.out = 100))
pe_mediators(gfit)
#>  mediator name            t_outcome t_exposure screen_p selected
#>  2        che_pc_usd      8.43      -1.686     0.0918   FALSE   
#>  4        gghed           3.027     0.6476     0.5172   FALSE   
#>  6        ext             1.141     1.534      0.2540   FALSE   
#>  9        pvtd_che        1.058     0.537      0.5913   FALSE   
#>  10       oops_che        0.9397    0.9471     0.3474   FALSE   
#>  14       gghed_gge       -2.405    1.231      0.2183   FALSE   
#>  16       pvtd_pc_usd     2.309     -0.8475    0.3967   FALSE   
#>  17       oop_pc_usd      -8.1      -1.315     0.1886   FALSE   
#>  19       cfa_che         -2.355    0.6753     0.4995   FALSE   
#>  21       chi_che         -3.603    3.111      0.0019   FALSE   
#>  22       shi_che         2.082     3.703      0.0374   FALSE   
#>  25       vhi_che         -3.755    -2.494     0.0126   FALSE   
#>  26       gge_gdp         -5.025    -3.732     0.0002    TRUE   
#>  29       gghed_usd       -4.748    0.4069     0.6841   FALSE   
#>  31       ext_usd         2.95      3.451      0.0032   FALSE   
#>  35       ext_ncu_pc      -3.875    -1.373     0.1697   FALSE   
#>  38       pvtd_ppp_pc     -0.5555   -0.2293    0.8187   FALSE   
#>  39       ext_ppp_pc      -1.532    2.224      0.1256   FALSE   
#>  40       pvtd_gdp        2.492     -1.064     0.2875   FALSE   
#>  41       ext_gdp         -4.303    0.325      0.7452   FALSE   
#>  45       ext_ncu2021     -3.454    1.731      0.0834   FALSE   
#>  52       pvtd_ncu2021_pc -4.381    0.1229     0.9022   FALSE   
#>  53       ext_ncu2021_pc  5.762     -1.532     0.1255   FALSE   
#>  57       ext_usd2021_pc  7.531     1.361      0.1735   FALSE
```

## Reproducibility Notes

- The global, regional, and income analyses above are deterministic.
  Called with no explicit grid arguments,
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
  searches `seq(0.1, 2, length.out = 100)`, the grid behind the
  article’s real-data tables, for the all-country model and for every
  subgroup model, and returns the selected mediators and p-values shown
  here.
- Part I runs no simulation study. Its `eval = FALSE` blocks give the
  calls that run the article’s designs at their scale and with 1,000
  replications: `n = 300` and `p = 500` on the grids of the article’s
  scripts (passed as `lambda_grid`) for Sections 4.1, 4.2, and S.4.1,
  and the case-study `n = 85` and `p = 1008` on the package default grid
  for supplement Section S.4.4, for which the article records no grid.
  Set `cores` above 1 to parallelize and `seed` for exact
  reproducibility; a seeded run with `cores` above 1 repeats exactly at
  the same number of cores.
- Exact real-data p-values depend on the preprocessing. The article’s
  text says that covariates are standardized before fitting.
  [`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md)
  instead applies the preprocessing of the authors’ analysis code behind
  the article’s tables: it standardizes the exposure and mediators,
  centers the outcome, and leaves the confounders (the region and income
  indicators and the year) unscaled. The region and income factors enter
  through treatment contrasts with alphabetical reference levels. The
  guided-tour vignette shows how to do this by hand and what
  standardizing the confounders as well changes.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

Guo, X., Li, R., Liu, J., & Zeng, M. (2022). High-dimensional mediation
analysis for selecting DNA methylation loci mediating childhood trauma
and cortisol stress reactivity. *Journal of the American Statistical
Association, 117*(539), 1110–1121.
<https://doi.org/10.1080/01621459.2022.2053136>

## Session Information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] POEMED_1.0.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.6         knitr_1.52        rlang_1.3.0       xfun_0.61        
#>  [5] otel_0.2.0        generics_0.1.4    textshaping_1.0.5 jsonlite_2.0.0   
#>  [9] htmltools_0.5.9   ragg_1.5.2        sass_0.4.10       glmnet_5.1       
#> [13] rmarkdown_2.32    grid_4.6.1        evaluate_1.0.5    jquerylib_0.1.4  
#> [17] fastmap_1.2.0     foreach_1.5.2     yaml_2.3.12       lifecycle_1.0.5  
#> [21] compiler_4.6.1    codetools_0.2-20  fs_2.1.0          Rcpp_1.1.2       
#> [25] lattice_0.22-9    systemfonts_1.3.2 digest_0.6.39     R6_2.6.1         
#> [29] ncvreg_3.16.0     splines_4.6.1     shape_1.4.6.1     bslib_0.12.0     
#> [33] Matrix_1.7-5      withr_3.0.3       tools_4.6.1       iterators_1.0.14 
#> [37] survival_3.8-6    pkgdown_2.2.1     cachem_1.1.0      desc_1.4.3
```
