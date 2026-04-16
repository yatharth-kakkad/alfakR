make_counts <- function(values, rownames_vec, colnames_vec) {
  x <- matrix(values, nrow = length(rownames_vec), byrow = TRUE)
  rownames(x) <- rownames_vec
  colnames(x) <- colnames_vec
  x
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
        predict = function(object, newdata, ...) {
          rep(mean(object$train_f), nrow(newdata))
        },
        .package = "stats"
      )
    },
    Krig = function(x, Y, ...) {
      structure(list(train_f = Y), class = "mock_krig")
    },
    .package = "fields"
  )
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
    nn_fitness = matrix(numeric(0), nrow = 1, ncol = 0)
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
