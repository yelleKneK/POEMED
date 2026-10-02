# Codebook for the WHO Health-Expenditure Indicators

The 57 health-expenditure indicators that serve as candidate mediators
in
[WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md),
with their short codes and descriptions. The codes and their order are
those of Table S.7 in the supplement to Yu and Kelley (in press); the
descriptions are adapted from that table rather than copied,
abbreviating some spelled-out denominators (for example "as % of GDP"
for "as % of Gross Domestic Product (GDP)" and "as % of CHE" for "as %
of Current Health Expenditure (CHE)"). Many indicators are alternative
normalizations of the same underlying spending (per capita, percent of
GDP, percent of current health expenditure, current versus constant 2021
currency, national currency versus US dollars versus purchasing power
parity), which is why they are strongly intercorrelated.

## Usage

``` r
WHO_indicator_codebook
```

## Format

A `data.frame` with 57 rows and 2 columns:

- indicator:

  The short indicator code, matching a column name in
  [WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md)
  (for example `gge_gdp`, `shi_che`).

- description:

  The full indicator name. Common abbreviations: CHE (current health
  expenditure), GGHE-D (domestic general government health expenditure),
  PVT-D (domestic private health expenditure), EXT (external health
  expenditure), OOPS (out-of-pocket spending), SHI (social health
  insurance), GGE (general government expenditure), GDP (gross domestic
  product), NCU (national currency unit), PPP (purchasing power parity).

## Source

Adapted from Table S.7 ("A List of Health Spending Indicators and the
Associated Indicator Code used in Real Data Analysis") of the supplement
to Yu and Kelley (in press); underlying definitions from the WHO Global
Health Expenditure Database.

## References

Yu, X., & Kelley, K. (in press). Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis. *Journal of the American Statistical
Association*.

## See also

[WHO_health_mediation](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md),
[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md).

Other WHO mediation data:
[`WHO_health_mediation`](https://yelleknek.github.io/POEMED/reference/WHO_health_mediation.md),
[`WHO_mediation_analysis()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_analysis.md),
[`WHO_mediation_design()`](https://yelleknek.github.io/POEMED/reference/WHO_mediation_design.md)

## Examples

``` r
data(WHO_indicator_codebook)
head(WHO_indicator_codebook)
#>    indicator
#> 1    che_gdp
#> 2 che_pc_usd
#> 3        che
#> 4      gghed
#> 5       pvtd
#> 6        ext
#>                                                             description
#> 1 Current Health Expenditure (CHE) as % of Gross Domestic Product (GDP)
#> 2                    Current Health Expenditure (CHE) per capita in US$
#> 3                                      Current Health Expenditure (CHE)
#> 4               Domestic General Government Health Expenditure (GGHE-D)
#> 5                           Domestic Private Health Expenditure (PVT-D)
#> 6                                     External Health Expenditure (EXT)
# Look up a mediator flagged in the article's analysis.
WHO_indicator_codebook[WHO_indicator_codebook$indicator == "gge_gdp", ]
#>    indicator                                      description
#> 26   gge_gdp General Government Expenditure (GGE) as % of GDP
```
