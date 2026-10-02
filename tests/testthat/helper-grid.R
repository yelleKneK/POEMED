# The fits in the test suite run at small designs on the default tuning grid,
# where HBIC usually picks the grid's smallest lambda and the fit warns (class
# poemed_grid_boundary) that the grid could be extended and refined. The
# warning is tested on its own in test-default_grid.R; in the other tests it
# is noise, so the tests that fit on the default grid run under this handler.
# Every other warning still surfaces.
quiet_grid <- function(expr) {
  withCallingHandlers(expr,
                      poemed_grid_boundary = function(w)
                        invokeRestart("muffleWarning"))
}

# A 20-value grid for tests whose subject is not the grid: about five times
# faster per fit than the default 100 values, which keeps the suite short on
# CRAN's check machines. Tests of the default itself are in test-default_grid.R.
fast_grid <- seq(0.05, 2, length.out = 20)
