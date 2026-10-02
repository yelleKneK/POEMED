# Power-Enhanced Mediation Test From a Data Frame (Formula Interface)

A convenience front end to
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
for users who have a data frame rather than ready-made matrices. Give it
the outcome and exposure as a formula, name the mediator columns, and
(optionally) name confounder columns; it assembles the exposure,
outcome, mediator, and confounder matrices (turning factor confounders
into indicator variables) and runs the test. This is the gentler entry
point; the matrix interface
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
is the workhorse.

## Usage

``` r
pe_mediate(
  formula,
  data,
  mediators,
  confounders = NULL,
  outcome = c("continuous", "binary", "count"),
  ...
)
```

## Arguments

- formula:

  A two-sided formula `outcome ~ exposure` naming the outcome and one or
  more exposure columns in `data`. A `.` stands for every other column
  of `data`, so either name the exposures explicitly or subtract the
  mediator and confounder columns (`y ~ . - m1 - m2`); a `.` that brings
  a mediator or confounder in as an exposure stops.

- data:

  A `data.frame` containing the outcome, exposure, mediator, and any
  confounder columns.

- mediators:

  Character vector of mediator column names in `data` (the candidate
  mediators `M`). At least two columns, each listed once, and none of
  them the outcome, an exposure, or a confounder.

- confounders:

  Optional character vector of confounder column names. Factor and
  character columns are expanded into indicator variables; numeric
  columns enter as is. Default `NULL`; an empty vector, `character(0)`,
  also means no confounders.

- outcome:

  Outcome type: `"continuous"`, `"binary"`, or `"count"`.

- ...:

  Further arguments passed to
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
  (for example `method`, `error_level`, `report_all_methods`). Among
  them is `lambda_grid`, the grid of candidate tuning parameters that
  HBIC searches, `seq(0.05, 10, length.out = 100)` by default for every
  outcome family: a shorter or narrower grid makes the fit faster, and a
  wider or finer one is worth trying when the value chosen is the
  smallest or largest in the grid (a non-empty fit warns, and the print
  footer, which reports the full model's choice, says so). See the
  argument on
  [`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md).

## Value

The tidy `poemed_tbl` returned by
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md).

## Details

Each column plays one role. The outcome and the exposures come from
`formula`, the mediators from `mediators`, and the confounders from
`confounders`. A column listed in two roles, or a mediator listed twice,
stops with a message naming it: the outcome listed as a mediator, for
example, would be "identified" as its own mediator.

Missing values follow one complete-case rule across every column the
model uses (outcome, exposures, mediators, and confounders). A row with
a missing value in any of them is dropped, and a message reports how
many rows were dropped and which columns held the missing values. The
method has no missing-data mechanism, so the test is computed on the
complete cases; impute beforehand if dropping rows is not appropriate
for your data. The matrix interface
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
drops nothing: it stops on a missing value and leaves the choice of
complete cases to the user.

The outcome must be numeric (or logical, for a binary outcome). A factor
or character outcome stops with a message, because a factor's numeric
codes are not its values: a count stored as a factor would otherwise be
analyzed as a different count.

## How to Cite

If you use POEMED in published work, please cite Yu and Kelley (in
press), the article that introduces its methods and the package.
`citation("POEMED")` gives the full reference and a BibTeX entry.

## See also

[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md)
for the matrix interface,
[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md)
for the worked health-expenditure design.

Other mediation tests:
[`pe_mediation()`](https://yelleknek.github.io/POEMED/reference/pe_mediation.md),
[`pe_mediation_linear()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_linear.md),
[`pe_mediation_logistic()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_logistic.md),
[`pe_mediation_poisson()`](https://yelleknek.github.io/POEMED/reference/pe_mediation_poisson.md),
[`pe_mediators()`](https://yelleknek.github.io/POEMED/reference/pe_mediators.md),
[`pe_selection()`](https://yelleknek.github.io/POEMED/reference/pe_selection.md),
[`summary.poemed_tbl()`](https://yelleknek.github.io/POEMED/reference/summary.poemed_tbl.md)

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
# Build a small data frame and test through the formula interface.
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, outcome = "continuous",
                             pattern = "contrasting", c1 = 1)
df <- data.frame(y = d$Y, x = d$X[, 1], d$M)
meds <- grep("^X", names(df), value = TRUE)   # the mediator columns
# The mediator columns are named, so the footer names the selected ones.
pe_mediate(y ~ x, data = df, mediators = meds, outcome = "continuous")
#> Warning: HBIC selected lambda = 0.05, the smallest value of `lambda_grid`, for the full model and the reduced model of the benchmark test. The criterion is minimized over the grid alone, so its minimum may lie beyond that end or between the end and its neighbor. Consider extending the grid past it and using a finer partition (more values), and compare the selected mediators and p-values across grids.
#>  term                  value   
#>  stat_hdmm             1.302   
#>  pval_hdmm             0.2538  
#>  stat_pe               865.5   
#>  j_pe                  864.2   
#>  pval_pe               < 0.0001
#>  total_indirect_effect -0.08396
#>  n_selected_mediators  3       
#>  df                    1       
#>  n_candidate_mediators 60      
#>  n_observations        200     
#> 
#> Outcome model: continuous (linear)
#> Selected mediators (3): X2, X3, X4
#> Tuning parameter (HBIC): lambda = 0.05 from 100 values in [0.05, 10] (the grid's lower end)
# HBIC chose the smallest value of the default grid here, so the fit
# warned; see the lambda_grid argument of ?pe_mediation.
```
