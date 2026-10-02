# Tidy Printing for POEMED Result Tables

Every estimation and testing function in POEMED returns a tidy
`data.frame` with a `term` column and one or more numeric columns. A
single numeric column routinely holds quantities on very different
scales: a whole-number count of mediators next to a chi-square statistic
in the hundreds next to a p-value of \\10^{-12}\\. The base
[`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html) method
formats a whole column with one common format, which forces either a
wall of trailing zeros or a slide into scientific notation. The
`poemed_tbl` class supplies `print` and `format` methods that format
each value on its own terms: whole numbers (counts, sample sizes,
mediator indices) print without a decimal part, other values print to a
few significant figures, and p-values print to a fixed number of decimal
places with a “\< 0.0001” floor.

## Usage

``` r
# S3 method for class 'poemed_tbl'
format(x, digits = getOption("poemed.digits", 4L), digits_p = 4L, ...)

# S3 method for class 'poemed_tbl'
print(x, digits = getOption("poemed.digits", 4L), digits_p = 4L, ...)

# S3 method for class 'poemed_tbl'
x[...]

# S3 method for class 'poemed_tbl'
rbind(..., deparse.level = 1)
```

## Arguments

- x:

  A `poemed_tbl` object (a tidy `data.frame` returned by a POEMED
  function).

- digits:

  Number of significant figures for non-integer values, a single whole
  number from 1 to 22 (the range base R's
  [`format()`](https://rdrr.io/r/base/format.html) accepts). Defaults to
  `getOption("poemed.digits", 4L)`, so `options(poemed.digits = 6)`
  changes it for the session.

- digits_p:

  Number of decimal places for p-values, a single whole number from 1 to
  15 (about the number of decimal digits a double holds reliably).
  Defaults to 4. A p-value below `10^(-digits_p)` prints as “\< ”
  followed by that bound, so “\< 0.0001” at the default and “\< 0.01”
  with `digits_p = 2`.

- ...:

  For `print`, further arguments passed to
  [`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html);
  `row.names` and `right` default to `FALSE` for the tidy look. For
  `format`, ignored. For `[`, the row and column indices and `drop`, as
  for a data frame. For `rbind`, the tables to combine and any options
  of [`rbind.data.frame`](https://rdrr.io/r/base/cbind.html), such as
  `make.row.names`.

- deparse.level:

  Passed to [`rbind()`](https://rdrr.io/r/base/cbind.html).

## Value

`print.poemed_tbl` returns `x` invisibly. `format.poemed_tbl` returns a
`data.frame` whose numeric columns have been formatted to character for
display. The `[` method returns a `poemed_tbl` with the attributes of
`x` (or a vector, when a single column is extracted), and
`rbind.poemed_tbl` returns a `poemed_tbl` whose attributes follow the
rule under Subsetting and Combining.

## Details

The stored numeric values are never rounded. Only the display changes,
so any arithmetic you do on the returned object (a difference of two
statistics, a further calculation) uses full precision. To see more
digits, raise `digits` or read the column directly with `x$value`.

This mirrors the display convention of the DMAR package, whose house
style POEMED follows.

## Subsetting and Combining

A `poemed_tbl` carries its display information as attributes: which rows
of the `value` column hold p-values, the outcome model, the selected
mediators, and the tuning record that the print footer reports.
Subsetting with `[` or [`subset()`](https://rdrr.io/r/base/subset.html)
keeps all of them, whether you select rows, columns, or both, so a
subset still prints its p-values with the floor and names the fit it
came from. (Base R would drop them whenever columns are selected.)

Combining tables with [`rbind()`](https://rdrr.io/r/base/cbind.html)
keeps these attributes only when every argument is a `poemed_tbl` and
all of them carry the same ones, as when you split one table and put it
back together. Tables from different fits describe different models, so
their combination keeps only the list of p-value rows and prints no
single-fit footer. A row whose `term` begins with `pval_` prints as a
p-value even when that list is missing.

## Author

Xiufan Yu and Ken Kelley

## Examples

``` r
set.seed(113)
d <- simulate_mediation_data(n = 200, p = 60, pattern = "contrasting",
                             outcome = "continuous")
fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
#> Warning: HBIC selected lambda = 0.05, the smallest value of `lambda_grid`, for the full model and the reduced model of the benchmark test. The criterion is minimized over the grid alone, so its minimum may lie beyond that end or between the end and its neighbor. Consider extending the grid past it and using a finer partition (more values), and compare the selected mediators and p-values across grids.
# HBIC chose the smallest value of the default grid here, so the fit
# warned; see the lambda_grid argument of ?pe_mediation.
fit                       # rounded for reading, with the p-value floor
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
#> Selected mediators (3): 2, 3, 4
#> Tuning parameter (HBIC): lambda = 0.05 from 100 values in [0.05, 10] (the grid's lower end)
fit$value[fit$term == "pval_pe"]   # full precision underneath
#> [1] 3.135368e-190

# A subset keeps the p-value format and the footer of its fit.
subset(fit, term %in% c("stat_pe", "pval_pe"))
#>  term    value   
#>  stat_pe 865.5   
#>  pval_pe < 0.0001
#> 
#> Outcome model: continuous (linear)
#> Selected mediators (3): 2, 3, 4
#> Tuning parameter (HBIC): lambda = 0.05 from 100 values in [0.05, 10] (the grid's lower end)
```
