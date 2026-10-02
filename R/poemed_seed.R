# Seed discipline for every function that takes a `seed` argument.
#
# A supplied seed is set for the duration of the calling function only. The
# caller's random number generator state, including the generator kind, is
# restored when that function exits, so a call made reproducible by `seed`
# does not perturb the random stream of the script around it. The saving and
# restoring is delegated to withr::local_seed(), the accepted way to do this
# on CRAN (the CRAN policies do not allow a package to modify the global
# environment, which a hand-rolled assignment of `.Random.seed` there does);
# the package itself never reads or writes `.Random.seed`. A NULL seed is a
# no-op: the draws come from the user's current generator state, exactly as
# if the function had no `seed` argument.
#
# A seed must be a single whole number that fits in R's integer range, the
# values set.seed() accepts. A fractional seed is rejected rather than
# truncated (1.7 would otherwise run as seed 1), and an out-of-range seed is
# rejected here rather than by set.seed() after a coercion warning.
#
# `parallel = TRUE` runs the call on the L'Ecuyer-CMRG generator, whose
# streams parallel::mclapply() hands out to forked replications so a seeded
# parallel run is reproducible at the same core count. With a seed, withr
# sets the kind and the state together. On exit withr puts back the caller's
# saved state, which carries its kind, when the caller had one. When the
# caller had no saved state (a fresh session, or a script whose first random
# draw is this call), withr removes the state it created but leaves the
# L'Ecuyer-CMRG kind in force, so every later draw in the session, the
# user's own set.seed() included, would come from the other generator. The
# helper therefore records the caller's kind first and registers its own
# restore after withr's. Deferred actions run last in, first out, so the
# kind is put back before withr removes or restores the state, and the
# caller ends with the kind and the state (or its absence) it started with.
# With no seed only the kind is switched, and switched back on exit, so the
# user's stream keeps advancing through the call exactly as it did before (an
# unseeded parallel call in a seeded wrapper such as pe_simulation_study()
# must not hand every pattern the same draws).
#
# The restore (.poemed_restore_kind()) puts back what the switch changed:
# the generator and the normal kind. The sample kind is set again only if
# it changed during the call, because setting it to "Rounding" warns
# ("non-uniform 'Rounding' sampler used"), and a caller who chose that
# sampler would otherwise see the warning after every parallel call.
#
# Call it from the body of the function whose exit should restore the state
# (the default `envir` is that function's frame), never from a helper.
#' @keywords internal
#' @noRd
.poemed_local_seed <- function(seed, parallel = FALSE, envir = parent.frame()) {
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != round(seed) || abs(seed) > .Machine$integer.max)
      stop("`seed` must be a single whole number between -",
           .Machine$integer.max, " and ", .Machine$integer.max,
           ", or NULL for no seed.", call. = FALSE)
    if (isTRUE(parallel)) {
      # RNGkind() with no arguments reports the kind without creating a
      # saved state, so it is safe to call when the caller has none.
      old_kind <- RNGkind()
      withr::local_seed(as.integer(seed), .local_envir = envir,
                        .rng_kind = "L'Ecuyer-CMRG")
      withr::defer(.poemed_restore_kind(old_kind), envir = envir)
    } else {
      withr::local_seed(as.integer(seed), .local_envir = envir)
    }
  } else if (isTRUE(parallel)) {
    old_kind <- RNGkind("L'Ecuyer-CMRG")
    withr::defer(.poemed_restore_kind(old_kind), envir = envir)
  }
  invisible(seed)
}

# Put back the generator kind and normal kind recorded in `old_kind` (the
# strings RNGkind() returns), and the sample kind only if it differs from
# the one in force. Only the change of generator draws from the stream (to
# seed the new generator); the normal and sample kinds are set without a
# draw, so leaving an unchanged sample kind alone gives the same stream as
# setting all three. From R 4.7 RNGkind() also reports the binomial
# sampler's kind as a fourth string; the package never changes it, and it
# is put back only if something inside the call did. Not exported.
#' @keywords internal
#' @noRd
.poemed_restore_kind <- function(old_kind) {
  RNGkind(kind = old_kind[1L], normal.kind = old_kind[2L])
  if (!identical(RNGkind()[3L], old_kind[3L]))
    suppressWarnings(RNGkind(sample.kind = old_kind[3L]))
  if (length(old_kind) >= 4L && !identical(RNGkind()[4L], old_kind[4L]))
    do.call(RNGkind, list(binom.kind = old_kind[4L]))
  invisible(NULL)
}
