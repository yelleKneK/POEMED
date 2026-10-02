# The seed contract of every exported function with a `seed` argument: a
# supplied seed reproduces the run bit for bit, the caller's `.Random.seed` is
# byte-identical afterwards whether or not one existed before the call, the
# generator kind is restored (with no saved state before the call included),
# and `seed = NULL` leaves the stream advancing. Verified at runtime, never by
# reading the source.

seeded_calls <- list(
  simulate_mediation_data = function(seed)
    simulate_mediation_data(n = 40, p = 12, outcome = "continuous",
                            pattern = "contrasting", seed = seed),
  simulate_guo_mediation = function(seed)
    simulate_guo_mediation(c1 = 0.5, n = 30, seed = seed),
  pe_power_curve = function(seed)
    pe_power_curve(n = 60, p = 15, outcome = "continuous",
                   pattern = "contrasting", c1_grid = 1, n_rep = 2,
                   lambda_grid = fast_grid, seed = seed),
  pe_simulation_study = function(seed)
    pe_simulation_study(n = 60, p = 15, outcome = "continuous", c1_grid = 1,
                        n_rep = 2, lambda_grid = fast_grid, seed = seed),
  pe_identification_study = function(seed)
    pe_identification_study(n = 60, p = 15, outcome = "continuous",
                            pattern = "contrasting", c1_grid = 1, n_rep = 2,
                            methods = "Bonferroni", lambda_grid = fast_grid,
                            seed = seed),
  # Two replications cannot resolve the default 0.8 target, so the
  # planner's two classed warnings are expected and muffled here;
  # test-study_guards.R makes each one fire.
  ss_power_pe_mediation = function(seed)
    withCallingHandlers(
      ss_power_pe_mediation(outcome = "continuous", pattern = "contrasting",
                            p = 15, n_grid = c(40, 60), n_rep = 2,
                            lambda_grid = fast_grid, seed = seed),
      poemed_target_unresolvable = function(w) invokeRestart("muffleWarning"),
      poemed_target_not_reached = function(w) invokeRestart("muffleWarning")))

strip <- function(x) { attributes(x) <- NULL; x }

for (nm in names(seeded_calls)) {
  f <- seeded_calls[[nm]]
  test_that(paste0(nm, "() reproduces a seeded run and restores the RNG state"), {
    # A seed reproduces the run bit for bit.
    a <- f(113); b <- f(113)
    expect_identical(a, b)
    # A different seed gives different data (the seed is not ignored); the
    # Monte Carlo tables are too coarse at two replications to assert this.
    if (nm %in% c("simulate_mediation_data", "simulate_guo_mediation"))
      expect_false(identical(strip(unclass(a)), strip(unclass(f(114)))))
    # Existing state: byte-identical afterwards, kind untouched.
    set.seed(999)
    before <- get(".Random.seed", envir = globalenv())
    kind_before <- RNGkind()
    invisible(f(113))
    expect_identical(get(".Random.seed", envir = globalenv()), before)
    expect_identical(RNGkind(), kind_before)
    # Absent state: still absent afterwards, and the kind is unchanged.
    if (exists(".Random.seed", envir = globalenv()))
      rm(".Random.seed", envir = globalenv())
    kind_before <- RNGkind()
    invisible(f(113))
    expect_false(exists(".Random.seed", envir = globalenv()))
    expect_identical(RNGkind(), kind_before)
    # seed = NULL draws from the user's stream and advances it.
    set.seed(999); invisible(f(NULL))
    after <- get(".Random.seed", envir = globalenv())
    set.seed(999)
    expect_false(identical(after, get(".Random.seed", envir = globalenv())))
  })
}

test_that("a seeded call restores a non-default generator kind", {
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(5)
  before <- get(".Random.seed", envir = globalenv())
  invisible(seeded_calls$pe_power_curve(113))
  expect_identical(RNGkind()[1L], "L'Ecuyer-CMRG")
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("`seed` must be a single whole number in the integer range, or NULL", {
  expect_error(simulate_mediation_data(n = 40, p = 12, seed = "a"), "`seed` must be")
  expect_error(simulate_mediation_data(n = 40, p = 12, seed = c(1, 2)), "`seed` must be")
  expect_error(pe_power_curve(n = 40, p = 12, n_rep = 2, seed = NA), "`seed` must be")
  # A fractional seed is rejected, not truncated to the integer below it.
  expect_error(simulate_mediation_data(n = 40, p = 12, seed = 1.7),
               "`seed` must be a single whole number")
  # A seed outside the integer range is rejected by name, before set.seed().
  expect_error(simulate_mediation_data(n = 40, p = 12, seed = 1e10),
               "`seed` must be a single whole number")
  expect_error(simulate_mediation_data(n = 40, p = 12, seed = -2^31),
               "`seed` must be a single whole number")
  # The integer range's ends are valid seeds.
  expect_no_error(simulate_mediation_data(n = 40, p = 12,
                                          seed = .Machine$integer.max))
})

# RNGkind() returns three strings through R 4.6 (generator, normal kind,
# sample kind) and four from R 4.7, which adds the binomial sampler's kind
# (`binom.kind`, with the fix of rbinom()'s BTPE algorithm). The package
# switches and restores the first three, so the tests below that name the
# kinds compare RNGkind()[1:3]; the tests that compare RNGkind() with its
# own earlier value cover whatever further component the running R reports.
#
# A seeded parallel call switches the generator to L'Ecuyer-CMRG for the call.
# When the caller had no saved state (a fresh session), withr removes the
# state it created but does not put the kind back, so the package restores
# the kind itself. The switch does not depend on forking, so this runs on
# every platform and on CRAN, at two cores.
test_that("a seeded parallel call with no saved state restores the generator kind", {
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  calls <- list(
    function() pe_power_curve(n = 60, p = 15, outcome = "continuous",
                              pattern = "contrasting", c1_grid = 1, n_rep = 2,
                              lambda_grid = fast_grid, cores = 2L, seed = 113),
    function() pe_identification_study(n = 60, p = 15, outcome = "continuous",
                                       pattern = "contrasting", c1_grid = 1,
                                       n_rep = 2, methods = "Bonferroni",
                                       lambda_grid = fast_grid,
                                       cores = 2L, seed = 113))
  for (f in calls) {
    # RNGkind() creates a saved state, so set the kind first, then remove it.
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    if (exists(".Random.seed", envir = globalenv()))
      rm(".Random.seed", envir = globalenv())
    invisible(f())
    expect_false(exists(".Random.seed", envir = globalenv()))
    expect_identical(RNGkind()[1:3], c("Mersenne-Twister", "Inversion", "Rejection"))
    # The user's own seeded draws afterwards are the Mersenne-Twister ones.
    set.seed(42)
    expect_equal(rnorm(1), 1.37095844714667, tolerance = 1e-12)
  }
})

# The restore puts back the generator and normal kinds and sets the sample
# kind again only if it changed, so a caller who chose the "Rounding"
# sampler (whose setting warns) sees no warning from a parallel call, seeded
# or not, and keeps that sampler.
test_that("a parallel call keeps a caller's Rounding sampler without a warning", {
  old <- RNGkind()
  on.exit(suppressWarnings(do.call(RNGkind, as.list(old))), add = TRUE)
  rounding <- c("Mersenne-Twister", "Inversion", "Rounding")
  suppressWarnings(do.call(RNGkind, as.list(rounding)))
  par_call <- function(seed)
    pe_power_curve(n = 60, p = 15, outcome = "continuous",
                   pattern = "contrasting", c1_grid = 1, n_rep = 2,
                   lambda_grid = fast_grid, cores = 2L, seed = seed)
  # Seeded, with a saved state: the state and all three kinds come back.
  set.seed(999)
  before <- get(".Random.seed", envir = globalenv())
  expect_no_warning(invisible(par_call(113)))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(RNGkind()[1:3], rounding)
  # Seeded, with no saved state.
  rm(".Random.seed", envir = globalenv())
  expect_no_warning(invisible(par_call(113)))
  expect_false(exists(".Random.seed", envir = globalenv()))
  expect_identical(RNGkind()[1:3], rounding)
  # Unseeded: the kind is switched for the call only.
  set.seed(999)
  expect_no_warning(invisible(par_call(NULL)))
  expect_identical(RNGkind()[1:3], rounding)
  # The identification study takes the same path.
  expect_no_warning(invisible(
    pe_identification_study(n = 60, p = 15, outcome = "continuous",
                            pattern = "contrasting", c1_grid = 1, n_rep = 2,
                            methods = "Bonferroni", lambda_grid = fast_grid,
                            cores = 2L, seed = 113)))
  expect_identical(RNGkind()[1:3], rounding)
})

test_that("a sample kind changed during a parallel call is put back", {
  old <- RNGkind()
  on.exit(suppressWarnings(do.call(RNGkind, as.list(old))), add = TRUE)
  inner <- function(seed) {
    POEMED:::.poemed_local_seed(seed, parallel = TRUE)
    suppressWarnings(RNGkind(sample.kind = "Rounding"))
    invisible(NULL)
  }
  for (seed in list(NULL, 113)) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(1)
    inner(seed)
    expect_identical(RNGkind()[1:3], c("Mersenne-Twister", "Inversion", "Rejection"))
  }
})

test_that("a seeded parallel run reproduces itself and restores the RNG state", {
  skip_on_cran()
  skip_on_os("windows")
  par_call <- function(seed, cores)
    pe_power_curve(n = 60, p = 15, outcome = "continuous", pattern = "contrasting",
                   c1_grid = 1, n_rep = 4, cores = cores, seed = seed)
  a <- par_call(113, 2L); b <- par_call(113, 2L)
  expect_identical(a, b)
  set.seed(999)
  before <- get(".Random.seed", envir = globalenv()); kind_before <- RNGkind()
  invisible(par_call(113, 2L))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(RNGkind(), kind_before)
  if (exists(".Random.seed", envir = globalenv()))
    rm(".Random.seed", envir = globalenv())
  kind_before <- RNGkind()
  invisible(par_call(113, 2L))
  expect_false(exists(".Random.seed", envir = globalenv()))
  expect_identical(RNGkind(), kind_before)
  # pe_identification_study takes the same path
  id1 <- pe_identification_study(n = 60, p = 15, outcome = "continuous",
                                 pattern = "contrasting", c1_grid = 1, n_rep = 4,
                                 methods = "Bonferroni", cores = 2L, seed = 113)
  id2 <- pe_identification_study(n = 60, p = 15, outcome = "continuous",
                                 pattern = "contrasting", c1_grid = 1, n_rep = 4,
                                 methods = "Bonferroni", cores = 2L, seed = 113)
  expect_identical(id1, id2)
  # an unseeded parallel call switches the kind for the call only (the kind
  # is not reset first, so a kind left behind by the calls above would show)
  expect_identical(RNGkind()[1L], "Mersenne-Twister")
  kind_before <- RNGkind()
  invisible(par_call(NULL, 2L))
  expect_identical(RNGkind(), kind_before)
})
