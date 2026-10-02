# Live oracle checks for POEMED: the package against the reference
# implementation that produced the numbers in the article.
#
# POEMED's estimation internals are ported from the reference `PEmediation`
# package (Yu and Kelley's review-era implementation, functions
# PEHDGMM_linear / _logistic / _poisson). One deliberate difference besides
# the p-value's upper tail: POEMED evaluates the SCAD derivative in the LLA
# step at |coefficient|, where the reference evaluated it at the signed
# coefficient. The two agree whenever no Step 1 coefficient falls below
# -lambda, which holds for the example data this harness compares, so the
# comparison below still applies bit for bit. The shipped test suite pins the
# article's values as constants, which is the operative guard on CRAN and
# for users. This script re-runs the LIVE comparison on the maintainer's
# machine, where the reference source exists, so drift between POEMED and the
# code behind the article is caught at release time without the reference
# ever being a dependency.
#
# For each outcome family it fits the shipped example data with both
# implementations at matching settings and requires: the power-enhanced
# statistic equal to 1e-8 (relative), the total indirect effect equal to
# 1e-8, and the same selected mediator set. The p-value is compared only
# where the reference's `1 - pchisq()` has not cancelled to zero; POEMED's
# `pchisq(lower.tail = FALSE)` is the one intentional numerical difference.
#
# Usage, from the package root, with POEMED installed (the release gate
# installs the tarball it certifies into a scratch library first):
#   Rscript tools/oracle_checks.R
# Set POEMED_REFERENCE to override the reference source location. A missing
# reference skips with a message rather than failing, so the script
# degrades gracefully on machines without it.
#
# This directory is excluded from the tarball by .Rbuildignore (^tools$).

ref_src <- Sys.getenv("POEMED_REFERENCE", file.path(
  "/Users/kkelley/Dropbox/Research",
  "Mediation, Power Enhancement High Dimensional/Archive/PriorPackage/PEmediation"))
if (!dir.exists(ref_src)) {
  message("oracle checks: 0 ok, 3 skipped (reference source not found at ", ref_src, ")")
  quit(status = 0L)
}

ref_lib <- file.path(tempdir(), "poemed_oracle_lib")
dir.create(ref_lib, showWarnings = FALSE)
st <- system2("R", c("CMD", "INSTALL", "--no-docs", "--no-test-load",
                     paste0("--library=", shQuote(ref_lib)), shQuote(ref_src)),
              stdout = FALSE, stderr = FALSE)
if (st != 0L) stop("could not install the reference PEmediation from ", ref_src, call. = FALSE)
suppressPackageStartupMessages(library(PEmediation, lib.loc = ref_lib))
suppressPackageStartupMessages(library(POEMED))

cases <- list(
  continuous = list(data = POEMED::example_continuous, family = "gaussian"),
  binary     = list(data = POEMED::example_binary,     family = "binomial"),
  count      = list(data = POEMED::example_count,      family = "poisson"))

rel <- function(a, b) abs(a - b) / max(1, abs(b))
ok <- 0L; failed <- character()
for (outcome in names(cases)) {
  d <- cases[[outcome]]$data
  ref <- PEmediation::PEHDGMM(X = d$X, Y = d$Y, M = d$M, Z = d$Z,
                              Y_family = cases[[outcome]]$family,
                              method = "Bonferroni", scale = TRUE,
                              FWER_level = 0.05,
                              lamb_grid = seq(0.05, 1, length.out = 20))
  fit <- POEMED::pe_mediation(d$X, d$Y, d$M, Z = d$Z, outcome = outcome,
                            method = "Bonferroni", scale = TRUE,
                            error_level = 0.05,
                            lambda_grid = seq(0.05, 1, length.out = 20))
  val <- function(term) fit$value[fit$term == term]
  poemed_selected <- sort(as.integer(attr(fit, "selected_mediators")))
  ref_selected  <- sort(as.integer(ref$selected_mediators))
  problems <- c(
    if (rel(val("stat_pe"), ref$stat_PEHDGMM) > 1e-8)
      sprintf("stat_pe %.10g vs reference %.10g", val("stat_pe"), ref$stat_PEHDGMM),
    if (rel(val("total_indirect_effect"), sum(ref$beta_hat)) > 1e-8)
      sprintf("total_indirect_effect %.10g vs reference %.10g",
              val("total_indirect_effect"), sum(ref$beta_hat)),
    if (!identical(poemed_selected, ref_selected))
      sprintf("selected set {%s} vs reference {%s}",
              paste(poemed_selected, collapse = ","), paste(ref_selected, collapse = ",")),
    if (ref$pval_PEHDGMM > 1e-12 && rel(val("pval_pe"), ref$pval_PEHDGMM) > 1e-6)
      sprintf("pval_pe %.10g vs reference %.10g", val("pval_pe"), ref$pval_PEHDGMM))
  if (length(problems)) {
    failed <- c(failed, paste0(outcome, ": ", paste(problems, collapse = "; ")))
  } else {
    ok <- ok + 1L
    message(sprintf("== %-10s stat_pe %.6g, total indirect %.6g, %d selected: matches the reference",
                    outcome, val("stat_pe"), val("total_indirect_effect"), length(poemed_selected)))
  }
}
message(sprintf("oracle checks: %d ok, 0 skipped, %d failed", ok, length(failed)))
if (length(failed)) stop(paste(failed, collapse = "\n"), call. = FALSE)
