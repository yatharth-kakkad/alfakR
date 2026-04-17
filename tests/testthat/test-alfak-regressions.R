make_counts <- function(values, rownames_vec, colnames_vec) {
  x <- matrix(values, nrow = length(rownames_vec), byrow = TRUE)
  rownames(x) <- rownames_vec
  colnames(x) <- colnames_vec
  x
}

reference_project_forward_log <- function(x0, f, timepoints) {
  out <- matrix(NA_real_, nrow = length(x0), ncol = length(timepoints))
  log_x0 <- log(x0)
  for (i in seq_along(timepoints)) {
    lv <- log_x0 + f * timepoints[i]
    denom <- alfakR:::logSumExp(lv)
    out[, i] <- exp(lv - denom)
  }
  out
}

reference_neg_log_lik <- function(param, counts, timepoints) {
  K <- nrow(counts)
  free_idx <- if (K > 1) seq_len(K - 1) else integer(0)
  f_free <- param[free_idx]
  f_full <- c(f_free, -sum(f_free))
  log_x0 <- param[K:(2 * K - 1)]
  nll <- 0
  for (i in seq_len(ncol(counts))) {
    lv <- log_x0 + f_full * timepoints[i]
    denom <- alfakR:::logSumExp(lv)
    for (k in seq_len(K)) {
      if (counts[k, i] > 0) {
        nll <- nll - counts[k, i] * (lv[k] - denom)
      }
    }
  }
  nll
}

reference_neighbor_objective <- function(fc_param, parent_fitness, pij_values,
                                         parent_birth_times, timepoints, parent_xfit,
                                         child_obs, ntot, parent_fitness_mean,
                                         prior_mean, prior_sd, do_prior, tol) {
  xc_est <- colSums(do.call(rbind, lapply(seq_along(parent_fitness), function(i) {
    tt <- pmax(0, timepoints - parent_birth_times[i])
    alfakR:::fExp_stable(fc_param, parent_fitness[i], pij_values[i], tt, tol = tol) * parent_xfit[i, ]
  })))
  xc_est <- pmax(0, pmin(1, xc_est))
  res <- stats::dbinom(child_obs, ntot, prob = xc_est, log = TRUE)
  if (do_prior && is.finite(parent_fitness_mean)) {
    res <- c(res, stats::dnorm(fc_param - parent_fitness_mean, mean = prior_mean, sd = prior_sd, log = TRUE))
  }
  res[!is.finite(res)] <- -(10^9)
  -sum(res)
}

make_simple_yi <- function(x, dt = 1) {
  list(x = x, dt = dt)
}

reference_qr_accum <- function(x_trim, dx_dt) {
  K <- nrow(x_trim)
  Q_accum <- matrix(0, nrow = K, ncol = K)
  r_accum <- numeric(K)
  for (t_idx in seq_len(ncol(x_trim))) {
    xt <- x_trim[, t_idx]
    M_t <- diag(as.numeric(xt), nrow = length(xt), ncol = length(xt)) - outer(xt, xt)
    Q_accum <- Q_accum + M_t %*% M_t
    r_accum <- as.numeric(r_accum + M_t %*% dx_dt[, t_idx])
  }
  list(Q_accum = Q_accum, r_accum = r_accum)
}

test_that("minobs includes karyotypes exactly at the threshold", {
  yi <- list(
    x = make_counts(
      c(10, 10,
        15, 15),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  res <- alfakR:::solve_fitness_bootstrap(
    yi,
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4
  )

  expect_setequal(colnames(res$final_fitness), c("2.2.2", "2.2.1"))
})

test_that("matrix-like count inputs are accepted and coerced at entry", {
  yi_df <- list(
    x = data.frame(
      "0" = c(10, 15),
      "1" = c(10, 15),
      row.names = c("2.2.2", "2.2.1"),
      check.names = FALSE
    ),
    dt = 1
  )

  res_df <- alfakR:::solve_fitness_bootstrap(
    yi_df,
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4
  )
  expect_setequal(colnames(res_df$final_fitness), c("2.2.2", "2.2.1"))

  skip_if_not_installed("Matrix")
  x_sparse <- Matrix::Matrix(c(10, 15, 10, 15), nrow = 2, sparse = TRUE)
  rownames(x_sparse) <- c("2.2.2", "2.2.1")
  colnames(x_sparse) <- c("0", "1")

  res_sparse <- alfakR:::solve_fitness_bootstrap(
    list(x = x_sparse, dt = 1),
    minobs = 20,
    nboot = 1,
    n0 = 1e4,
    nb = 1e6,
    pm = 1e-4
  )
  expect_setequal(colnames(res_sparse$final_fitness), c("2.2.2", "2.2.1"))
})

test_that("strict karyotype parsing rejects malformed IDs and mixed dimensions", {
  parsed <- alfakR:::parse_karyotype_ids(c("2.2", "2.3"))
  expect_identical(dim(parsed), c(2L, 2L))
  expect_true(is.integer(parsed))
  expect_equal(unname(parsed[1, ]), c(2L, 2L))

  expect_error(alfakR:::parse_karyotype_ids(c("2a.2", "2.3")), "Invalid karyotype ID")
  expect_error(alfakR:::parse_karyotype_ids("2..2"), "Invalid karyotype ID")
  expect_error(alfakR:::parse_karyotype_ids(c("2.0.2", "2.3.2")), "Invalid karyotype ID")
  expect_error(alfakR:::parse_karyotype_ids(c("2.-1.2", "2.3.2")), "Invalid karyotype ID")
  expect_error(alfakR:::parse_karyotype_ids(c("2.2", "2.2.2")), "same number of dot-separated components")
})

test_that("build_W_rcpp validates p, Nmax, and karyotype strings safely", {
  expect_silent(alfakR::build_W_rcpp(c("2.2", "2.3"), p = 0.01, Nmax = Inf))
  expect_error(alfakR::build_W_rcpp(c("2.2", "2a.3"), p = 0.01), "Invalid karyotype ID")
  expect_error(alfakR::build_W_rcpp(c("2.2", "2.2.2"), p = 0.01), "same number of dot-separated components")
  expect_error(alfakR::build_W_rcpp(c("2.2", "2.3"), p = NA_real_), "`p`")
  expect_error(alfakR::build_W_rcpp(c("2.2", "2.3"), p = 1.5), "`p`")
})

test_that("solve_fitness_bootstrap rejects malformed karyotype rownames before bootstrapping", {
  yi_bad <- list(
    x = make_counts(
      c(10, 10,
        20, 20),
      rownames_vec = c("2a.2", "2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  expect_error(
    alfakR:::solve_fitness_bootstrap(
      yi_bad,
      minobs = 1,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = 1e-4
    ),
    "Invalid karyotype ID"
  )
})

test_that("count-matrix validation rejects invalid values and rounds non-integers once", {
  x_inf <- make_counts(
    c(10, Inf,
      15, 15),
    rownames_vec = c("2.2.2", "2.2.1"),
    colnames_vec = c("0", "1")
  )
  expect_error(
    alfakR:::coerce_count_matrix(x_inf),
    "finite count values"
  )

  x_neg <- make_counts(
    c(10, -1,
      15, 15),
    rownames_vec = c("2.2.2", "2.2.1"),
    colnames_vec = c("0", "1")
  )
  expect_error(
    alfakR:::coerce_count_matrix(x_neg),
    "non-negative count values"
  )

  x_non_integer <- data.frame(
    "0" = c(10.2, 15.7),
    "1" = c(9.8, 14.3),
    row.names = c("2.2.2", "2.2.1"),
    check.names = FALSE
  )
  expect_warning(
    rounded <- alfakR:::coerce_count_matrix(x_non_integer),
    "rounding to the nearest integer once at entry"
  )
  expect_equal(rounded, round(as.matrix(x_non_integer)))
})

test_that("zero-depth timepoints are rejected before normalization or bootstrap", {
  x_zero <- make_counts(
    c(10, 0,
      15, 0),
    rownames_vec = c("2.2.2", "2.2.1"),
    colnames_vec = c("0", "1")
  )

  expect_error(
    alfakR:::solve_fitness_bootstrap(
      make_simple_yi(x_zero),
      minobs = 1,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = 1e-4
    ),
    "zero-depth column"
  )

  expect_error(
    alfakR::alfak(
      yi = make_simple_yi(x_zero),
      outdir = tempfile("alfak_zero_depth_"),
      minobs = 1,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = 1e-4
    ),
    "zero-depth column"
  )
})

test_that("alfak validates optional arguments before running heavy work", {
  yi <- make_simple_yi(
    make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    )
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      returned <- invisible(alfakR::alfak(
        yi = yi,
        outdir = tempfile("alfak_default_passage_"),
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4
      ))
      expect_identical(returned, 0.1)

      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 0, n0 = 1e4, nb = 1e6, pm = 1e-4), "`nboot`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = -1, n0 = 1e4, nb = 1e6, pm = 1e-4), "`nboot`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1.5, n0 = 1e4, nb = 1e6, pm = 1e-4), "`nboot`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 0, nb = 1e6, pm = 1e-4), "`n0`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = Inf, nb = 1e6, pm = 1e-4), "`n0`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = NA_real_, pm = 1e-4), "`nb`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = NA), "`correct_efflux`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = "TRUE"), "`correct_efflux`")
      expect_error(alfakR::alfak(yi = yi, outdir = tempfile(), minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = 1), "`correct_efflux`")
    },
    solve_fitness_bootstrap = function(...) {
      seen$solve_called <- TRUE
      list(
        initial_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
        final_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
        initial_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
        final_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
        nn_fitness = matrix(numeric(0), nrow = 1, ncol = 0)
      )
    },
    fitKrig = function(...) {
      list(
        summary_stats = data.frame(k = "2.2.2", mean = 0, median = 0, sd = 0, fq = TRUE, nn = FALSE),
        posterior_samples = matrix(0, nrow = 1, ncol = 1),
        krig_stable_mean = NULL,
        krig_stable_median = NULL
      )
    },
    xval = function(...) list(R2R = 0.1),
    .package = "alfakR"
  )

  expect_true(isTRUE(seen$solve_called))
})

test_that("solve_fitness_bootstrap validates bootstrap controls and pm before neighbour generation", {
  yi <- make_simple_yi(
    make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    )
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 0, n0 = 1e4, nb = 1e6, pm = 1e-4), "`nboot`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 0, nb = 1e6, pm = 1e-4), "`n0`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = NA_real_, pm = 1e-4), "`nb`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = NA_real_), "`pm`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = NA), "`correct_efflux`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = "TRUE"), "`correct_efflux`")
      expect_error(alfakR:::solve_fitness_bootstrap(yi, minobs = 1, nboot = 1, n0 = 1e4, nb = 1e6, pm = 1e-4, correct_efflux = 1), "`correct_efflux`")
    },
    gen_nn_info = function(...) {
      seen$gen_nn_called <- TRUE
      list()
    },
    .package = "alfakR"
  )

  expect_false(isTRUE(seen$gen_nn_called))
})

test_that("correct_efflux stops before bootstrap when viability is non-positive", {
  yi <- list(
    x = make_counts(
      c(5, 5),
      rownames_vec = "2.2.2",
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  expect_error(
    alfakR:::solve_fitness_bootstrap(
      yi,
      minobs = 1,
      nboot = 2,
      n0 = 1e4,
      nb = 1e6,
      pm = 0.2,
      correct_efflux = TRUE
    ),
    "correct_efflux viability pre-check failed before bootstrap: pm=0.2"
  )
})

test_that("correct_efflux warns once when viability is positive but tiny", {
  yi <- list(
    x = make_counts(
      c(5, 5,
        6, 6),
      rownames_vec = c("4.4.4.4.4", "2.2.2.2.2"),
      colnames_vec = c("0", "2")
    ),
    dt = 99
  )
  pm <- 1 - ((1 + 5e-7) / 2)^(1 / 20)

  expect_warning(
    alfakR:::solve_fitness_bootstrap(
      yi,
      minobs = 1,
      nboot = 1,
      n0 = 1e4,
      nb = 1e6,
      pm = pm,
      correct_efflux = TRUE,
      passage_times = c(0, 2.5)
    ),
    "0 < viability <"
  )
})

test_that("fExp_stable stays finite and matches the analytic limit near fc == fp", {
  tt <- c(0, 1, 2.5, 5)
  limit_val <- 0.2 * 0.3 * tt

  exact <- alfakR:::fExp_stable(fc_arg = 0.3, fp_arg = 0.3, pij_val = 0.2, tt_arg = tt)
  near <- alfakR:::fExp_stable(fc_arg = 0.3 + 1e-12, fp_arg = 0.3, pij_val = 0.2, tt_arg = tt)

  expect_true(all(is.finite(exact)))
  expect_true(all(is.finite(near)))
  expect_equal(exact, limit_val, tolerance = 1e-12)
  expect_equal(near, limit_val, tolerance = 1e-12)
})

test_that("C++ numerical kernels match the previous R reference calculations", {
  x0 <- c(0.3, 0.7)
  f <- c(0.1, -0.1)
  timepoints <- c(0, 1.5, 4)
  project_ref <- reference_project_forward_log(x0, f, timepoints)
  project_cpp <- alfak_project_forward_log_cpp(x0, f, timepoints)
  expect_equal(project_cpp, project_ref, tolerance = 1e-12)

  counts <- matrix(c(10, 5, 7,
                     4, 11, 9), nrow = 2, byrow = TRUE)
  param <- c(0.2, log(0.4), log(0.6))
  expect_equal(
    alfak_neg_log_lik_cpp(param, counts, timepoints),
    reference_neg_log_lik(param, counts, timepoints),
    tolerance = 1e-12
  )

  parent_fitness <- c(0.15, 0.05)
  pij_values <- c(0.2, 0.1)
  parent_birth_times <- c(-1, 0.5)
  parent_xfit <- matrix(c(0.4, 0.35, 0.3,
                          0.2, 0.25, 0.3), nrow = 2, byrow = TRUE)
  child_obs <- c(0, 2, 3)
  ntot <- c(10, 12, 14)
  expect_equal(
    alfak_neighbor_objective_cpp(
      fc_param = 0.11,
      parent_fitness = parent_fitness,
      pij_values = pij_values,
      parent_birth_times = parent_birth_times,
      timepoints = timepoints,
      parent_xfit = parent_xfit,
      child_obs = child_obs,
      ntot = ntot,
      parent_fitness_mean = 0.12,
      prior_mean = 0.01,
      prior_sd = 0.2,
      do_prior = TRUE,
      tol = alfakR:::ALFAK_FEXP_DELTA_TOL
    ),
    reference_neighbor_objective(
      fc_param = 0.11,
      parent_fitness = parent_fitness,
      pij_values = pij_values,
      parent_birth_times = parent_birth_times,
      timepoints = timepoints,
      parent_xfit = parent_xfit,
      child_obs = child_obs,
      ntot = ntot,
      parent_fitness_mean = 0.12,
      prior_mean = 0.01,
      prior_sd = 0.2,
      do_prior = TRUE,
      tol = alfakR:::ALFAK_FEXP_DELTA_TOL
    ),
    tolerance = 1e-12
  )

  x_trim <- matrix(c(0.3, 0.4,
                     0.7, 0.6), nrow = 2, byrow = TRUE)
  dx_dt <- matrix(c(0.1, -0.05,
                    -0.1, 0.05), nrow = 2, byrow = TRUE)
  qr_ref <- reference_qr_accum(x_trim, dx_dt)
  qr_cpp <- alfak_qr_accum_cpp(x_trim, dx_dt)
  expect_equal(qr_cpp$Q_accum, qr_ref$Q_accum, tolerance = 1e-12)
  expect_equal(qr_cpp$r_accum, qr_ref$r_accum, tolerance = 1e-12)
})

test_that("C++ numerical kernels validate dimensions and non-finite inputs", {
  expect_error(
    alfak_project_forward_log_cpp(c(0.5, 0.5), c(0.1), c(0, 1)),
    "same length"
  )
  expect_error(
    alfak_project_forward_log_cpp(c(0, 0), c(0.1, -0.1), c(0, 1)),
    "positive finite value"
  )

  counts_bad <- matrix(c(1, NA_real_, 2, 3), nrow = 2)
  expect_error(
    alfak_neg_log_lik_cpp(c(0.1, log(0.5), log(0.5)), counts_bad, c(0, 1)),
    "finite non-negative values"
  )

  expect_error(
    alfak_neighbor_objective_cpp(
      fc_param = 0.1,
      parent_fitness = c(0.2, 0.1),
      pij_values = 0.1,
      parent_birth_times = c(0, 1),
      timepoints = c(0, 1),
      parent_xfit = matrix(0.5, nrow = 2, ncol = 2),
      child_obs = c(0, 1),
      ntot = c(10, 10),
      parent_fitness_mean = 0.15,
      prior_mean = 0,
      prior_sd = 0.1,
      do_prior = FALSE,
      tol = 1e-8
    ),
    "matching lengths/rows"
  )
})

test_that("negative parent exposure times are clamped at zero in neighbour objective", {
  full_obj <- alfak_neighbor_objective_cpp(
    fc_param = 0.1,
    parent_fitness = c(0.2, 0.2),
    pij_values = c(0.3, 0.4),
    parent_birth_times = c(10, 0),
    timepoints = c(0, 1),
    parent_xfit = matrix(0.5, nrow = 2, ncol = 2),
    child_obs = c(0, 0),
    ntot = c(10, 10),
    parent_fitness_mean = 0.2,
    prior_mean = 0,
    prior_sd = 0.1,
    do_prior = FALSE,
    tol = 1e-8
  )

  clamped_reference <- alfak_neighbor_objective_cpp(
    fc_param = 0.1,
    parent_fitness = 0.2,
    pij_values = 0.4,
    parent_birth_times = 0,
    timepoints = c(0, 1),
    parent_xfit = matrix(0.5, nrow = 1, ncol = 2),
    child_obs = c(0, 0),
    ntot = c(10, 10),
    parent_fitness_mean = 0.2,
    prior_mean = 0,
    prior_sd = 0.1,
    do_prior = FALSE,
    tol = 1e-8
  )

  expect_equal(full_obj, clamped_reference, tolerance = 1e-12)
})

test_that("passage_times defines one validated internal time axis everywhere", {
  yi <- list(
    x = make_counts(
      c(10, 11, 12,
        20, 21, 22),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1", "2")
    ),
    dt = 99
  )
  passage_times <- c(0, 2.5, 7)
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      res <- alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        passage_times = passage_times
      )

      expect_equal(seen$compute_dx_dt, passage_times)
      expect_equal(seen$optimize_initial_frequencies, passage_times)
      expect_equal(seen$joint_optimize, passage_times)
      expect_equal(seen$project_forward_log, passage_times)
      expect_equal(seen$find_birth_times, c(-1000, max(passage_times)))
      expect_equal(
        unname(res$final_fitness[1, ]),
        rep(log(1e6 / 1e4) / diff(passage_times)[1], 2),
        tolerance = 1e-12
      )
    },
    compute_dx_dt = function(x, timepoints) {
      seen$compute_dx_dt <- timepoints
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    optimize_initial_frequencies = function(x_obs, f, timepoints) {
      seen$optimize_initial_frequencies <- timepoints
      rep(1 / nrow(x_obs), nrow(x_obs))
    },
    joint_optimize = function(counts, timepoints, f_init, x0_init) {
      seen$joint_optimize <- timepoints
      list(f = rep(0, nrow(counts)), x0 = rep(1 / nrow(counts), nrow(counts)))
    },
    project_forward_log = function(x0, f, timepoints) {
      seen$project_forward_log <- timepoints
      matrix(rep(x0, length(timepoints)), nrow = length(x0), ncol = length(timepoints))
    },
    find_birth_times = function(opt_res, time_range, minF) {
      seen$find_birth_times <- time_range
      rep(0, length(opt_res$f))
    },
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
      list(solution = rep(0, nrow(Dmat)))
    },
    gen_nn_info = function(fq, pm) {
      list()
    },
    .package = "alfakR"
  )
})

test_that("resolve_time_axis rejects non-increasing supplied passage_times", {
  yi <- list(
    x = make_counts(
      c(10, 11, 12,
        20, 21, 22),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1", "2")
    ),
    dt = 1
  )

  expect_error(
    alfakR:::resolve_time_axis(yi, passage_times = c(0, 2, 2)),
    "strictly increasing"
  )
})

test_that("birth-time fallback keeps neighbour estimation finite when roots are all missing", {
  yi <- list(
    x = make_counts(
      c(10, 12, 11),
      rownames_vec = "2.2.2",
      colnames_vec = c("0", "1", "2")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  expect_warning(
    res <- testthat::with_mocked_bindings(
      {
        alfakR:::solve_fitness_bootstrap(
          yi,
          minobs = 1,
          nboot = 1,
          n0 = 1e4,
          nb = 1e6,
          pm = 1e-4
        )
      },
      compute_dx_dt = function(x, timepoints) {
        matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
      },
      optimize_initial_frequencies = function(x_obs, f, timepoints) {
        1
      },
      joint_optimize = function(counts, timepoints, f_init, x0_init) {
        list(f = 0, x0 = 1)
      },
      project_forward_log = function(x0, f, timepoints) {
        matrix(1, nrow = 1, ncol = length(timepoints),
               dimnames = list("2.2.2", NULL))
      },
      find_birth_times = function(opt_res, time_range, minF) {
        rep(NA_real_, length(opt_res$f))
      },
      run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
        list(solution = 0)
      },
      gen_nn_info = function(fq, pm) {
        nn <- list(list(ni = "2.2.3", nj = "2.2.2", pij = 0.1))
        names(nn) <- "2.2.3"
        nn
      },
      run_optimise_checked = function(f, interval, ..., context) {
        objective <- f(mean(interval), ...)
        seen$objective <- objective
        list(minimum = mean(interval), objective = objective)
      },
      .package = "alfakR"
    ),
    "finite fallback birth times"
  )

  expect_true(is.finite(seen$objective))
  expect_true(all(is.finite(res$nn_fitness)))
})

test_that("xval samples one real bootstrap replicate for all folds", {
  fq_boot <- list(
    final_fitness = matrix(
      c(1, 10, 100,
        2, 20, 200,
        3, 30, 300),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(NULL, c("1.1", "3.3", "5.5"))
    ),
    nn_fitness = matrix(numeric(0), nrow = 3, ncol = 0)
  )

  testthat::with_mocked_bindings(
    {
      set.seed(123)
      res <- alfakR:::xval(fq_boot)
      observed <- as.numeric(res$tmp[, "test_f"])
      is_real_row <- vapply(
        seq_len(nrow(fq_boot$final_fitness)),
        function(i) identical(observed, as.numeric(fq_boot$final_fitness[i, ])),
        logical(1)
      )
      expect_true(any(is_real_row))
    },
    build_cached_krig_fit = function(ktrain, y, kpred = NULL, give_warnings = TRUE) {
      list(fit = structure(list(train_f = y), class = "mock_krig"), pred_dist = NULL)
    },
    predict_cached_krig = function(object, x, dist_mat, ...) {
      rep(mean(object$train_f), nrow(x))
    },
    .package = "alfakR"
  )
})

test_that("cached Krig refits match fresh fields::Krig fits", {
  set.seed(1)
  ktrain <- matrix(runif(16), ncol = 2)
  y_initial <- c(0.2, -0.1, 0.3, 0.05, -0.2, 0.4, 0.1, -0.15)
  y_updated <- c(0.4, -0.25, 0.15, 0.1, -0.1, 0.35, 0.05, -0.05)
  kpred <- matrix(
    c(0.25, 0.5,
      0.75, 0.5,
      0.5, 0.25),
    ncol = 2,
    byrow = TRUE
  )

  cache <- suppressWarnings(
    alfakR:::build_cached_krig_fit(ktrain, y_initial, kpred = kpred, give_warnings = FALSE)
  )
  refit <- suppressWarnings(
    alfakR:::refit_cached_krig(
      cache,
      y = y_updated,
      x_pred = kpred,
      pred_dist = cache$pred_dist,
      give_warnings = FALSE
    )
  )

  fresh_fit <- suppressWarnings(
    fields::Krig(
      ktrain,
      y_updated,
      cov.function = "stationary.cov",
      cov.args = list(Covariance = "Matern", smoothness = 1.5),
      nstep.cv = alfakR:::ALFAK_KRIG_NSTEP_CV,
      give.warnings = FALSE
    )
  )
  fresh_preds <- as.numeric(stats::predict(fresh_fit, kpred))

  expect_equal(refit$fit$lambda, fresh_fit$lambda, tolerance = 1e-10)
  expect_equal(refit$fit$eff.df, fresh_fit$eff.df, tolerance = 1e-10)
  expect_equal(as.numeric(refit$preds), fresh_preds, tolerance = 1e-8)
  expect_equal(as.numeric(stats::predict(refit$fit, kpred)), fresh_preds, tolerance = 1e-8)
})

test_that("fitKrig returns structured NA bootstrap outputs when a bootstrap fit is not trainable", {
  fq_boot <- list(
    final_fitness = matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "2.2.2")),
    nn_fitness = matrix(numeric(0), nrow = 3, ncol = 0)
  )

  res <- suppressWarnings(alfakR:::fitKrig(fq_boot, nboot = 2))

  expect_identical(dim(res$posterior_samples), c(1L, 2L))
  expect_true(all(is.na(res$posterior_samples)))
  expect_length(res$boot_results, 2)
  expect_true(all(vapply(res$boot_results, is.list, logical(1))))
  expect_true(all(vapply(res$fit_boot_list, is.null, logical(1))))
  expect_true(all(vapply(res$boot_results, function(x) all(is.na(x$preds)), logical(1))))
  expect_true(all(is.na(res$summary_stats$mean)))
  expect_true(all(is.na(res$summary_stats$median)))
  expect_true(all(is.na(res$summary_stats$sd)))
})

test_that("nn_prior = 'none' disables latent-neighbour prior contribution", {
  yi <- list(
    x = make_counts(
      c(10, 11,
        5, 4),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  seen <- new.env(parent = emptyenv())

  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior = "none"
      )
    },
    compute_dx_dt = function(x, timepoints) {
      matrix(0, nrow = nrow(x), ncol = ncol(x) - 1)
    },
    optimize_initial_frequencies = function(x_obs, f, timepoints) {
      1
    },
    joint_optimize = function(counts, timepoints, f_init, x0_init) {
      list(f = 0, x0 = 1)
    },
    project_forward_log = function(x0, f, timepoints) {
      matrix(1, nrow = 1, ncol = length(timepoints),
             dimnames = list("2.2.2", NULL))
    },
    find_birth_times = function(opt_res, time_range, minF) {
      0
    },
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) {
      list(solution = 0)
    },
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    run_optimise_checked = function(f, interval, ..., context) {
      dots <- list(...)
      if (grepl("latent child 2.2.1", context, fixed = TRUE)) {
        seen$latent_do_prior <- isTRUE(dots$do_prior_param)
      }
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )

  expect_false(isTRUE(seen$latent_do_prior))
})

test_that("empirical prior SD uses floor and user-supplied nn_prior_sd is respected", {
  yi <- list(
    x = make_counts(
      c(10, 11,
        5, 4),
      rownames_vec = c("2.2.2", "2.2.3"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )

  floor_capture <- new.env(parent = emptyenv())
  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior_sd_floor = 0.123
      )
    },
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    optimize_initial_frequencies = function(x_obs, f, timepoints) 1,
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0, x0 = 1),
    project_forward_log = function(x0, f, timepoints) matrix(1, nrow = 1, ncol = length(timepoints), dimnames = list("2.2.2", NULL)),
    find_birth_times = function(opt_res, time_range, minF) 0,
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = 0),
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior) && is.finite(prior_sd) && is.null(floor_capture$prior_sd)) {
        floor_capture$prior_sd <- prior_sd
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )
  expect_equal(floor_capture$prior_sd, 0.123, tolerance = 1e-12)

  user_capture <- new.env(parent = emptyenv())
  testthat::with_mocked_bindings(
    {
      alfakR:::solve_fitness_bootstrap(
        yi,
        minobs = 20,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        nn_prior_sd = 0.456,
        nn_prior_sd_floor = 0.123
      )
    },
    compute_dx_dt = function(x, timepoints) matrix(0, nrow = nrow(x), ncol = ncol(x) - 1),
    optimize_initial_frequencies = function(x_obs, f, timepoints) 1,
    joint_optimize = function(counts, timepoints, f_init, x0_init) list(f = 0, x0 = 1),
    project_forward_log = function(x0, f, timepoints) matrix(1, nrow = 1, ncol = length(timepoints), dimnames = list("2.2.2", NULL)),
    find_birth_times = function(opt_res, time_range, minF) 0,
    run_solve_qp_checked = function(Dmat, dvec, Amat, bvec, meq, context) list(solution = 0),
    gen_nn_info = function(fq, pm) {
      nn <- list(
        list(ni = "2.2.3", nj = "2.2.2", pij = 0.2),
        list(ni = "2.2.1", nj = "2.2.2", pij = 0.2)
      )
      names(nn) <- c("2.2.3", "2.2.1")
      nn
    },
    alfak_neighbor_objective_cpp = function(fc_param, parent_fitness, pij_values,
                                            parent_birth_times, timepoints, parent_xfit,
                                            child_obs, ntot, parent_fitness_mean,
                                            prior_mean, prior_sd, do_prior, tol) {
      if (isTRUE(do_prior) && is.finite(prior_sd) && is.null(user_capture$prior_sd)) {
        user_capture$prior_sd <- prior_sd
      }
      0
    },
    run_optimise_checked = function(f, interval, ..., context) {
      f(mean(interval))
      list(minimum = mean(interval), objective = 0)
    },
    .package = "alfakR"
  )
  expect_equal(user_capture$prior_sd, 0.456, tolerance = 1e-12)
})

test_that("find_steady_state matches long-time row-vector ODE dynamics", {
  skip_if_not_installed("deSolve")
  skip_if_not_installed("Matrix")

  lscape <- data.frame(
    k = c("2.2", "3.1"),
    mean = c(1, 1),
    stringsAsFactors = FALSE
  )
  A <- Matrix::Matrix(
    matrix(c(0.1, 0.9,
             0.8, 0.2), nrow = 2, byrow = TRUE),
    sparse = TRUE,
    dimnames = list(lscape$k, lscape$k)
  )
  ode_out <- deSolve::ode(
    y = c(0.8, 0.2),
    times = seq(0, 100, by = 0.1),
    func = alfakR:::chrmod_rel,
    parms = list(A = A)
  )
  terminal <- as.numeric(ode_out[nrow(ode_out), -1])
  terminal <- terminal / sum(terminal)

  ss <- testthat::with_mocked_bindings(
    {
      alfakR::find_steady_state(lscape, p = 0.01)
    },
    build_W_rcpp = function(karyotype_strings, p, Nmax = Inf) {
      A
    },
    .package = "alfakR"
  )

  expect_equal(unname(ss), terminal, tolerance = 1e-4)
})

test_that("ABM parent accounting consumes each dividing parent exactly once", {
  res <- alfakR:::run_karyotype_abm(
    initial_population_r = stats::setNames(list(10), "1"),
    fitness_map_r = stats::setNames(list(1, 1), c("1", "2")),
    p_missegregation = 1,
    dt = 1,
    n_steps = 1L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 123L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )

  expect_true("1" %in% names(res))
  step1_counts <- res[["1"]]
  expect_equal(sum(as.numeric(step1_counts)), 10)
})

test_that("landscape_data_output controls whether landscape_data.Rds is written", {
  yi <- list(
    x = make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  fq_boot_stub <- list(
    initial_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    initial_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )
  landscape_stub <- list(
    summary_stats = data.frame(
      k = "2.2.2",
      mean = 0,
      median = 0,
      sd = 0,
      fq = TRUE,
      nn = FALSE
    ),
    posterior_samples = matrix(0, nrow = 1, ncol = 1),
    krig_stable_mean = list(tag = "mean"),
    krig_stable_median = list(tag = "median")
  )
  xval_stub <- list(
    tmp = matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("test_f", "est_f"))),
    R2R = 0
  )
  outdir_false <- file.path(tempdir(), "alfak_landscape_data_false")
  outdir_true <- file.path(tempdir(), "alfak_landscape_data_true")
  unlink(outdir_false, recursive = TRUE)
  unlink(outdir_true, recursive = TRUE)

  testthat::with_mocked_bindings(
    {
      invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir_false,
        passage_times = NULL,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        landscape_data_output = FALSE
      ))
      expect_false(file.exists(file.path(outdir_false, "landscape_data.Rds")))

      invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir_true,
        passage_times = NULL,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4,
        landscape_data_output = TRUE
      ))
      expect_true(file.exists(file.path(outdir_true, "landscape_data.Rds")))
    },
    solve_fitness_bootstrap = function(...) fq_boot_stub,
    fitKrig = function(...) landscape_stub,
    xval = function(...) xval_stub,
    .package = "alfakR"
  )
})

test_that("softmax is stable for large logits and rejects non-finite values", {
  expect_equal(alfakR:::softmax(c(1000, 1000)), c(0.5, 0.5), tolerance = 1e-12)
  expect_error(alfakR:::softmax(c(0, Inf)), "non-finite logits")
})

test_that("alfak saves xval.Rds as a scalar R2R for downstream compatibility", {
  yi <- list(
    x = make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  fq_boot_stub <- list(
    initial_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    initial_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )
  landscape_stub <- list(
    summary_stats = data.frame(
      k = "2.2.2",
      mean = 0,
      median = 0,
      sd = 0,
      fq = TRUE,
      nn = FALSE
    ),
    posterior_samples = matrix(0, nrow = 1, ncol = 1),
    krig_stable_mean = NULL,
    krig_stable_median = NULL
  )
  xval_stub <- list(
    tmp = matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("test_f", "est_f"))),
    R2R = 0.42
  )
  outdir <- file.path(tempdir(), "alfak_xval_scalar")
  unlink(outdir, recursive = TRUE)

  testthat::with_mocked_bindings(
    {
      returned <- invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir,
        passage_times = NULL,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4
      ))
      saved <- readRDS(file.path(outdir, "xval.Rds"))
      expect_identical(returned, 0.42)
      expect_identical(saved, 0.42)
    },
    solve_fitness_bootstrap = function(...) fq_boot_stub,
    fitKrig = function(...) landscape_stub,
    xval = function(...) xval_stub,
    .package = "alfakR"
  )
})

test_that("alfak accepts scalar NA_real_ cross-validation outputs and still writes core files", {
  yi <- list(
    x = make_counts(
      c(10, 12,
        20, 18),
      rownames_vec = c("2.2.2", "2.2.1"),
      colnames_vec = c("0", "1")
    ),
    dt = 1
  )
  fq_boot_stub <- list(
    initial_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_fitness = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    initial_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    final_frequencies = matrix(1, nrow = 1, ncol = 1, dimnames = list(NULL, "2.2.2")),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )
  landscape_stub <- list(
    summary_stats = data.frame(
      k = "2.2.2",
      mean = 0,
      median = 0,
      sd = 0,
      fq = TRUE,
      nn = FALSE
    ),
    posterior_samples = matrix(0, nrow = 1, ncol = 1),
    krig_stable_mean = NULL,
    krig_stable_median = NULL
  )
  outdir <- file.path(tempdir(), "alfak_xval_na")
  unlink(outdir, recursive = TRUE)

  testthat::with_mocked_bindings(
    {
      returned <- invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir,
        passage_times = NULL,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4
      ))
      saved <- readRDS(file.path(outdir, "xval.Rds"))
      expect_true(is.numeric(returned) && length(returned) == 1 && is.na(returned))
      expect_true(is.numeric(saved) && length(saved) == 1 && is.na(saved))
      expect_true(file.exists(file.path(outdir, "bootstrap_res.Rds")))
      expect_true(file.exists(file.path(outdir, "landscape.Rds")))
      expect_true(file.exists(file.path(outdir, "landscape_posterior_samples.Rds")))
    },
    solve_fitness_bootstrap = function(...) fq_boot_stub,
    fitKrig = function(...) landscape_stub,
    xval = function(...) NA_real_,
    .package = "alfakR"
  )
})

test_that("constant-response cross-validation returns NA_real_ instead of NaN", {
  fq_boot <- list(
    final_fitness = matrix(
      c(1, 1, 1,
        1, 1, 1),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("2.2", "3.3", "4.4"))
    ),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )

  expect_warning(
    res <- alfakR:::xval(fq_boot),
    "Not enough valid observations"
  )
  expect_true(is.numeric(res$R2R) && length(res$R2R) == 1 && is.na(res$R2R))
  expect_true(is.na(alfakR:::extract_xval_r2r(res)))
})

test_that("alfak saves scalar NA_real_ when Krig fitting fails during xval", {
  yi <- make_simple_yi(
    make_counts(
      c(10, 12,
        20, 18,
        5, 6),
      rownames_vec = c("2.2.2", "2.2.1", "2.2.3"),
      colnames_vec = c("0", "1")
    )
  )
  fq_boot_stub <- list(
    initial_fitness = matrix(c(0.1, 0.2, 0.3,
                               0.15, 0.25, 0.35), nrow = 2, byrow = TRUE,
                             dimnames = list(NULL, c("2.2.2", "2.2.1", "2.2.3"))),
    final_fitness = matrix(c(0.1, 0.2, 0.3,
                             0.15, 0.25, 0.35), nrow = 2, byrow = TRUE,
                           dimnames = list(NULL, c("2.2.2", "2.2.1", "2.2.3"))),
    initial_frequencies = matrix(c(0.3, 0.4, 0.3,
                                   0.25, 0.45, 0.3), nrow = 2, byrow = TRUE,
                                 dimnames = list(NULL, c("2.2.2", "2.2.1", "2.2.3"))),
    final_frequencies = matrix(c(0.3, 0.4, 0.3,
                                 0.25, 0.45, 0.3), nrow = 2, byrow = TRUE,
                               dimnames = list(NULL, c("2.2.2", "2.2.1", "2.2.3"))),
    nn_fitness = matrix(numeric(0), nrow = 2, ncol = 0)
  )
  landscape_stub <- list(
    summary_stats = data.frame(k = c("2.2.2", "2.2.1", "2.2.3"), mean = 0, median = 0, sd = 0, fq = TRUE, nn = FALSE),
    posterior_samples = matrix(0, nrow = 3, ncol = 1),
    krig_stable_mean = NULL,
    krig_stable_median = NULL
  )
  outdir <- file.path(tempdir(), "alfak_xval_krig_failure")
  unlink(outdir, recursive = TRUE)

  testthat::with_mocked_bindings(
    {
      returned <- invisible(alfakR::alfak(
        yi = yi,
        outdir = outdir,
        minobs = 1,
        nboot = 1,
        n0 = 1e4,
        nb = 1e6,
        pm = 1e-4
      ))
      saved <- readRDS(file.path(outdir, "xval.Rds"))
      expect_true(is.numeric(returned) && length(returned) == 1 && is.na(returned))
      expect_true(is.numeric(saved) && length(saved) == 1 && is.na(saved))
    },
    solve_fitness_bootstrap = function(...) fq_boot_stub,
    fitKrig = function(...) landscape_stub,
    build_cached_krig_fit = function(...) stop("mock Krig failure"),
    .package = "alfakR"
  )
})

test_that("ABM treats zero and negative fitness as no-division and validates p_missegregation", {
  zero_res <- alfakR:::run_karyotype_abm(
    initial_population_r = stats::setNames(list(10), "1"),
    fitness_map_r = stats::setNames(list(0), "1"),
    p_missegregation = 0,
    dt = 1,
    n_steps = 1L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 123L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )
  expect_equal(as.numeric(zero_res[["1"]]), 10)

  negative_res <- alfakR:::run_karyotype_abm(
    initial_population_r = stats::setNames(list(10), "1"),
    fitness_map_r = stats::setNames(list(-1), "1"),
    p_missegregation = 0,
    dt = 1,
    n_steps = 1L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 123L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )
  expect_equal(as.numeric(negative_res[["1"]]), 10)

  expect_error(
    alfakR:::run_karyotype_abm(
      initial_population_r = stats::setNames(list(10), "1"),
      fitness_map_r = stats::setNames(list(1), "1"),
      p_missegregation = 1.2,
      dt = 1,
      n_steps = 1L,
      max_population_size = 0,
      culling_survival_fraction = 0.1,
      record_interval = 1L,
      seed = 123L,
      grf_centroids = matrix(0, 0, 0),
      grf_lambda = NA_real_
    ),
    "p_missegregation"
  )
})

test_that("ABM supports large population counts without int-sized random-distribution limits", {
  res <- alfakR:::run_karyotype_abm(
    initial_population_r = stats::setNames(list(3e9), "1"),
    fitness_map_r = stats::setNames(list(0), "1"),
    p_missegregation = 0,
    dt = 1,
    n_steps = 1L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 123L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )
  expect_equal(as.numeric(res[["1"]]), 3e9)
})

test_that("prediction helpers reject invalid probabilities, times, and x0 inputs", {
  lscape <- data.frame(k = c("2.2", "3.1"), mean = c(0.1, 0.2), stringsAsFactors = FALSE)
  x0 <- c("2.2" = 0.5, "3.1" = 0.5)

  expect_error(alfakR::predict_evo(lscape, p = NA_real_, times = c(0, 1), x0 = x0), "`p`")
  expect_error(alfakR::predict_evo(lscape, p = NaN, times = c(0, 1), x0 = x0), "`p`")
  expect_error(alfakR::predict_evo(lscape, p = 1.2, times = c(0, 1), x0 = x0), "`p`")
  expect_error(alfakR::predict_evo(lscape, p = 0.1, times = c(0, Inf), x0 = x0), "`times`")
  expect_error(alfakR::predict_evo(lscape, p = 0.1, times = c(0, 1), x0 = c("2.2" = NA_real_, "3.1" = 1)), "`x0`")

  expect_error(
    alfakR::run_abm_simulation_grf(
      centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
      lambda = 1,
      p = 1.2,
      times = c(0, 1),
      x0 = c("2.2" = 1)
    ),
    "`p`"
  )
  expect_error(
    alfakR::run_abm_simulation_grf(
      centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
      lambda = 1,
      p = 0.1,
      times = c(0, Inf),
      x0 = c("2.2" = 1)
    ),
    "`times`"
  )
})
