# Tests for the display layer outside print.poemed_tbl itself: the plot
# method (both views, the cap for a screening p-value stored as 0, and the
# graphical parameters passed through its dots), the broom verbs, and the
# empty-fit warning and footer. Plots draw to a pdf file in tempdir(), and
# the drawn coordinates are read back from the device's display list.

# The points drawn while `code` runs on a pdf device in tempdir(): one entry
# per call that drew points (plot, lines, legend), each with its x and y and
# the plotting symbols used. Returns NULL when the display list does not have
# the layout these tests read, so a caller can skip rather than fail.
drawn_points <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  grDevices::dev.control("enable")
  force(code)
  pts <- tryCatch({
    dl <- grDevices::recordPlot()[[1L]]
    lapply(dl, function(e) {
      a <- as.list(e[[2L]])
      if (length(a) >= 4L && is.list(a[[1L]]) &&
          identical(a[[1L]]$name, "C_plotXY"))
        list(x = a[[2L]]$x, y = a[[2L]]$y, pch = a[[4L]])
      else NULL
    })
  }, error = function(e) NULL)
  pts <- Filter(Negate(is.null), pts)
  if (length(pts)) pts else NULL
}

# Draw `code` to a pdf file in tempdir() and discard the file.
to_pdf <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  force(code)
}

# A fit-shaped table carrying a hand-made mediator table, so the plot's
# handling of a screening p-value stored as exactly 0 can be tested without
# a strong-signal fit.
fake_fit <- function(screen_p, selected) {
  tab <- POEMED:::.as_poemed_tbl(data.frame(term = "pval_pe", value = 0.01),
                                 p_terms = "pval_pe")
  attr(tab, "mediator_table") <- data.frame(
    mediator = seq_along(screen_p), t_outcome = 1, t_exposure = 1,
    screen_p = screen_p, selected = selected)
  tab
}

test_that("a screening p-value stored as 0 is drawn at a marked cap", {
  ev <- POEMED:::.screen_evidence(c(0, 1e-3, 0.5, NA))
  cap <- -log10(.Machine$double.eps)
  expect_equal(ev$cap, cap)
  expect_identical(ev$capped, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(ev$height, c(cap, 3, -log10(0.5), NA))
  # A finite height above the default cap raises the cap to it.
  ev2 <- POEMED:::.screen_evidence(c(0, 1e-40))
  expect_equal(ev2$cap, 40)
  expect_equal(ev2$height, c(40, 40))
  # No zero, no cap in use.
  expect_false(any(POEMED:::.screen_evidence(c(0.2, 0.01))$capped))

  fit <- fake_fit(c(0, 0, 1e-6, 0.4), c(TRUE, TRUE, TRUE, FALSE))
  pts <- drawn_points(expect_invisible(plot(fit)))
  skip_if(is.null(pts), "the display list layout is not the one these tests read")
  expect_equal(pts[[1L]]$y, c(cap, cap, 6, -log10(0.4)))
  expect_equal(pts[[1L]]$pch, c(17, 17, 19, 19))       # capped points are triangles
  # The legend gains an entry (a third point) for the capped points.
  expect_equal(pts[[length(pts)]]$pch, c(19, 19, 17))
})

test_that("a strong-signal fit is drawn with finite heights and no 307-decade point", quiet_grid({
  set.seed(113)
  tau1 <- matrix(0, 1, 40)
  tau1[1, 1:4] <- 0.8
  d <- simulate_mediation_data(n = 200, p = 40, q = 1, outcome = "continuous",
                               pattern = "homogeneous", c1 = 1, tau = tau1)
  fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous")
  tab <- pe_mediators(fit)
  stored <- tab$screen_p
  pts <- drawn_points(plot(fit))
  expect_identical(pe_mediators(fit)$screen_p, stored)   # plotting changes nothing
  skip_if(is.null(pts), "the display list layout is not the one these tests read")
  y <- pts[[1L]]$y
  expect_true(all(is.finite(y)))
  expect_true(all(y < 300))
  ev <- POEMED:::.screen_evidence(stored)
  expect_identical(ev$capped, stored == 0)
  expect_equal(y, ev$height)
}))

test_that("the per-mediator view accepts the common graphical parameters", {
  fit <- fake_fit(c(0, 1e-6, 0.4), c(TRUE, TRUE, FALSE))
  expect_no_error(to_pdf(plot(fit, xlab = "X", ylab = "Y", pch = 2, col = "red",
                              main = "M", xaxt = "s", ylim = c(0, 20))))
  pts <- drawn_points(plot(fit, pch = 2))
  skip_if(is.null(pts), "the display list layout is not the one these tests read")
  expect_true(all(pts[[1L]]$pch == 2))                   # the user's pch wins
  expect_equal(pts[[length(pts)]]$pch, c(2, 2, 2))       # and reaches the legend
})

test_that("the power-curve view draws both curves and accepts graphical parameters", {
  pc <- pe_power_curve(n = 80, p = 30, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(0, 1), n_rep = 2,
                       seed = 113)
  pts <- drawn_points(expect_invisible(plot(pc)))
  expect_no_error(to_pdf(plot(pc, ylim = c(0, 0.5), xlab = "c", ylab = "rate",
                              pch = 2, type = "l", col = "blue", lty = 3)))
  skip_if(is.null(pts), "the display list layout is not the one these tests read")
  expect_equal(pts[[1L]]$x, pc$c1)
  expect_equal(pts[[1L]]$y, pc$rejection_pe)
  expect_equal(pts[[2L]]$y, pc$rejection_hdmm)
})

test_that("a table with no plot view says so, and the broom verbs pass it through", {
  plain <- POEMED:::.as_poemed_tbl(data.frame(a = 1))
  expect_message(to_pdf(plot(plain)), "No POEMED plot")
  pc <- pe_power_curve(n = 80, p = 30, outcome = "continuous",
                       pattern = "contrasting", c1_grid = c(0, 1), n_rep = 2,
                       seed = 113)
  expect_identical(tidy(pc), as.data.frame(unclass(pc)))
  expect_identical(glance(pc), as.data.frame(unclass(pc)))
})

test_that("tidy() numbers its rows 1..k whatever the mediator column holds", quiet_grid({
  db <- simulate_mediation_data(n = 100, p = 30, outcome = "binary",
                                pattern = "homogeneous", c1 = 1, seed = 113)
  fb <- pe_mediation(db$X, db$Y, db$M, outcome = "binary")
  td <- tidy(fb)
  expect_identical(rownames(td), as.character(seq_len(nrow(td))))
  # The values are those of the per-mediator table.
  expect_equal(td$mediator, pe_mediators(fb)$mediator)
  expect_equal(td$screen_p, pe_mediators(fb)$screen_p)
}))

test_that("the empty-fit warning and footer describe the HBIC-selected fit", {
  set.seed(113)
  d <- simulate_mediation_data(n = 100, p = 30, outcome = "continuous",
                               pattern = "homogeneous", c1 = 0)
  w <- expect_warning(
    fit <- pe_mediation(d$X, d$Y, d$M, outcome = "continuous",
                        lambda_grid = c(5, 10)),
    class = "poemed_empty_fit")
  msg <- conditionMessage(w)
  expect_match(msg, "The penalized fit that HBIC selected (lambda = 10, from 2 values",
               fixed = TRUE)
  expect_no_match(msg, "at any of the", fixed = TRUE)
  expect_no_match(msg, "smaller values may select", fixed = TRUE)
  # The p = 1 convention is the package's reporting convention; the article
  # states only that J_m = 0 when the selection is empty.
  expect_match(msg, "POEMED's reporting convention for an empty selection",
               fixed = TRUE)
  expect_no_match(msg, "Yu and Kelley", fixed = TRUE)
  out <- capture.output(print(fit))
  expect_true(any(grepl("empty at the HBIC-selected lambda", out, fixed = TRUE)))
  expect_false(any(grepl("empty at every lambda", out, fixed = TRUE)))
})
