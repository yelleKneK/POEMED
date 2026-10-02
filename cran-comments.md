# cran-comments

## Submission

This is a new submission.

POEMED implements the methods of the accompanying peer-reviewed article,
Yu, X., & Kelley, K. (in press), "Power Enhancement in High-Dimensional
Heterogeneous Mediation Analysis," *Journal of the American Statistical
Association*.

## Test environments

* local macOS (Apple Silicon), R 4.6.1: full `R CMD check --as-cran`
  with both manuals built, from the release tarball, built from a clean
  archive of the repository by the maintainer's release gate
  (2026-10-01).
* win-builder R-release and R-devel: [[FILL: dates and status]].
* GitHub Actions R-CMD-check matrix (Ubuntu release, devel, and oldrel,
  macOS, Windows): [[FILL: result after the first push]].

## R CMD check results

0 errors, 0 warnings, 1 NOTE:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Ken Kelley <kkelley@nd.edu>'

New submission
```

The local log counts two NOTEs, the second from the HTML manual check,
"Skipping checking HTML validation: 'tidy' doesn't look like recent
enough HTML Tidy", which is environmental (the CRAN check farm has a
current HTML Tidy) and so is not listed above.

## Check time

The package contains no `\donttest{}` and no `\dontrun{}`, so the
examples run in a single pass: about 12 s across the 24 help
pages with examples, the slowest page about 1 s. The tests run in
about 39 s on CRAN's path (the Monte Carlo blocks and the full
real-data reproduction carry `skip_on_cran()` and run in the maintainer's
release gate, where the whole suite of 1,440 expectations passes) and the
three vignettes build in about 47 s. The whole check took
about 134 s locally. On win-builder the examples ran in
[[FILL: s]], the tests in [[FILL: s]], and the vignette rebuild in
[[FILL: s]].

## Possibly misspelled words in DESCRIPTION

The incoming check flags four words; all are intended.

* **POwer** and **MEDiation** are deliberate: they show where the package
  name's letters come from (POwer-Enhanced MEDiation).
* **Yu** is the first author's surname (Xiufan Yu), cited in the
  Description as "Yu and Kelley (in press)".
* **familywise** is the standard spelling of the error rate (familywise
  error rate), as in the American Psychological Association's style and
  the multiple-testing literature.

## Downstream dependencies

This is a new package; no reverse dependencies on CRAN.
