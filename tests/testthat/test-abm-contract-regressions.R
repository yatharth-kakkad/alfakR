test_that("largest remainder allocation preserves exact large integer-valued totals", {
  total_size <- .Machine$integer.max + 25
  alloc <- alfakR:::largest_remainder_allocate(c(a = 0.2, b = 0.3, c = 0.5), total_size)

  expect_identical(names(alloc), c("a", "b", "c"))
  expect_true(is.double(alloc))
  expect_identical(length(alloc), 3L)
  expect_true(all(alloc >= 0))
  expect_true(all(alloc == floor(alloc)))
  expect_identical(sum(alloc), total_size)
})

test_that("ABM wrappers keep large initial counts above the 32-bit range integer-valued", {
  total_size <- .Machine$integer.max + 25
  lscape <- data.frame(k = c("2.2", "3.1"), mean = c(0.1, 0.2), stringsAsFactors = FALSE)
  x0 <- c("2.2" = 0.4, "3.1" = 0.6)
  captured_regular <- NULL
  captured_grf <- NULL

  testthat::with_mocked_bindings(
    {
      suppressMessages(
        alfakR::predict_evo(
          lscape = lscape,
          p = 0.01,
          times = 0,
          x0 = x0,
          prediction_type = "ABM",
          abm_pop_size = total_size,
          abm_delta_t = 0.1,
          abm_record_interval = 1,
          abm_seed = 1
        )
      )
    },
    run_karyotype_abm = function(initial_population_r, fitness_map_r, p_missegregation, dt,
                                 n_steps, max_population_size, culling_survival_fraction,
                                 record_interval, seed, grf_centroids, grf_lambda) {
      captured_regular <<- initial_population_r
      list("0" = stats::setNames(as.numeric(unlist(initial_population_r, use.names = FALSE)),
                                 names(initial_population_r)))
    },
    .package = "alfakR"
  )

  expect_true(all(as.numeric(unlist(captured_regular, use.names = FALSE)) == floor(as.numeric(unlist(captured_regular, use.names = FALSE)))))
  expect_identical(sum(as.numeric(unlist(captured_regular, use.names = FALSE))), total_size)

  testthat::with_mocked_bindings(
    {
      alfakR::run_abm_simulation_grf(
        centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
        lambda = 1,
        p = 0.01,
        times = 0,
        x0 = c("2.2" = 0.4, "3.1" = 0.6),
        abm_pop_size = total_size,
        abm_delta_t = 0.1,
        abm_record_interval = 1,
        abm_seed = 1
      )
    },
    run_karyotype_abm = function(initial_population_r, fitness_map_r, p_missegregation, dt,
                                 n_steps, max_population_size, culling_survival_fraction,
                                 record_interval, seed, grf_centroids, grf_lambda) {
      captured_grf <<- initial_population_r
      list("0" = stats::setNames(as.numeric(unlist(initial_population_r, use.names = FALSE)),
                                 names(initial_population_r)))
    },
    .package = "alfakR"
  )

  expect_true(all(as.numeric(unlist(captured_grf, use.names = FALSE)) == floor(as.numeric(unlist(captured_grf, use.names = FALSE)))))
  expect_identical(sum(as.numeric(unlist(captured_grf, use.names = FALSE))), total_size)
})

test_that("run_karyotype_abm rejects invalid initial population counts", {
  common_args <- list(
    fitness_map_r = stats::setNames(list(0), "1"),
    p_missegregation = 0,
    dt = 1,
    n_steps = 0L,
    max_population_size = 0,
    culling_survival_fraction = 0.1,
    record_interval = 1L,
    seed = 1L,
    grf_centroids = matrix(0, 0, 0),
    grf_lambda = NA_real_
  )

  expect_error(
    do.call(alfakR:::run_karyotype_abm, c(list(initial_population_r = stats::setNames(list(1.5), "1")), common_args)),
    "integer-valued"
  )
  expect_error(
    do.call(alfakR:::run_karyotype_abm, c(list(initial_population_r = stats::setNames(list(Inf), "1")), common_args)),
    "finite"
  )
  expect_error(
    do.call(alfakR:::run_karyotype_abm, c(list(initial_population_r = stats::setNames(list(-1), "1")), common_args)),
    "non-negative"
  )
})

test_that("abm_max_pop accepts unlimited and capped semantics in regular and GRF wrappers", {
  lscape <- data.frame(k = c("2.2", "3.1"), mean = c(0.1, 0.2), stringsAsFactors = FALSE)
  x0 <- c("2.2" = 0.5, "3.1" = 0.5)
  regular_caps <- numeric(0)
  grf_caps <- numeric(0)
  counts_out <- stats::setNames(c(50, 50), c("2.2", "3.1"))

  testthat::with_mocked_bindings(
    {
      for (cap in c(-1, 0, 5)) {
        suppressMessages(
          alfakR::predict_evo(
            lscape = lscape,
            p = 0.01,
            times = 0,
            x0 = x0,
            prediction_type = "ABM",
            abm_pop_size = 100,
            abm_delta_t = 0.1,
            abm_max_pop = cap,
            abm_record_interval = 1,
            abm_seed = 1
          )
        )
      }
    },
    run_karyotype_abm = function(initial_population_r, fitness_map_r, p_missegregation, dt,
                                 n_steps, max_population_size, culling_survival_fraction,
                                 record_interval, seed, grf_centroids, grf_lambda) {
      regular_caps <<- c(regular_caps, max_population_size)
      list("0" = counts_out)
    },
    .package = "alfakR"
  )

  expect_identical(regular_caps, c(-1, 0, 5))

  testthat::with_mocked_bindings(
    {
      for (cap in c(-1, 0, 5)) {
        alfakR::run_abm_simulation_grf(
          centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
          lambda = 1,
          p = 0.01,
          times = 0,
          x0 = c("2.2" = 0.5, "3.1" = 0.5),
          abm_pop_size = 100,
          abm_delta_t = 0.1,
          abm_max_pop = cap,
          abm_record_interval = 1,
          abm_seed = 1
        )
      }
    },
    run_karyotype_abm = function(initial_population_r, fitness_map_r, p_missegregation, dt,
                                 n_steps, max_population_size, culling_survival_fraction,
                                 record_interval, seed, grf_centroids, grf_lambda) {
      grf_caps <<- c(grf_caps, max_population_size)
      list("0" = counts_out)
    },
    .package = "alfakR"
  )

  expect_identical(grf_caps, c(-1, 0, 5))
})

test_that("ABM wrappers return time columns exactly as requested and preserve duplicates", {
  lscape <- data.frame(k = c("2.2", "3.1"), mean = c(0.1, 0.2), stringsAsFactors = FALSE)
  x0 <- c("2.2" = 0.25, "3.1" = 0.75)
  requested_times <- c(0, 0, 0.1)
  regular_interval <- NULL
  grf_interval <- NULL

  regular_res <- testthat::with_mocked_bindings(
    {
      suppressMessages(
        alfakR::predict_evo(
          lscape = lscape,
          p = 0.01,
          times = requested_times,
          x0 = x0,
          prediction_type = "ABM",
          abm_pop_size = 100,
          abm_delta_t = 0.1,
          abm_record_interval = -1,
          abm_seed = 1
        )
      )
    },
    run_karyotype_abm = function(initial_population_r, fitness_map_r, p_missegregation, dt,
                                 n_steps, max_population_size, culling_survival_fraction,
                                 record_interval, seed, grf_centroids, grf_lambda) {
      regular_interval <<- record_interval
      list(
        "0" = stats::setNames(c(25, 75), c("2.2", "3.1")),
        "1" = stats::setNames(c(40, 60), c("2.2", "3.1"))
      )
    },
    .package = "alfakR"
  )

  expect_identical(regular_interval, 1L)
  expect_identical(regular_res$time, requested_times)
  expect_identical(nrow(regular_res), length(requested_times))
  expect_identical(as.list(regular_res[1, c("2.2", "3.1")]), as.list(regular_res[2, c("2.2", "3.1")]))

  grf_res <- testthat::with_mocked_bindings(
    {
      alfakR::run_abm_simulation_grf(
        centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
        lambda = 1,
        p = 0.01,
        times = requested_times,
        x0 = c("2.2" = 1),
        abm_pop_size = 100,
        abm_delta_t = 0.1,
        abm_record_interval = -1,
        abm_seed = 1
      )
    },
    run_karyotype_abm = function(initial_population_r, fitness_map_r, p_missegregation, dt,
                                 n_steps, max_population_size, culling_survival_fraction,
                                 record_interval, seed, grf_centroids, grf_lambda) {
      grf_interval <<- record_interval
      list(
        "0" = stats::setNames(c(100), "2.2"),
        "1" = stats::setNames(c(100), "2.2")
      )
    },
    .package = "alfakR"
  )

  expect_identical(grf_interval, 1L)
  expect_identical(grf_res$time, requested_times)
  expect_identical(nrow(grf_res), length(requested_times))
  expect_identical(as.list(grf_res[1, "2.2", drop = FALSE]), as.list(grf_res[2, "2.2", drop = FALSE]))
})

test_that("ABM wrappers reject times that do not align with the ABM step grid", {
  lscape <- data.frame(k = c("2.2", "3.1"), mean = c(0.1, 0.2), stringsAsFactors = FALSE)
  x0 <- c("2.2" = 0.5, "3.1" = 0.5)

  expect_error(
    suppressMessages(
      alfakR::predict_evo(
        lscape = lscape,
        p = 0.01,
        times = c(0, 0.15),
        x0 = x0,
        prediction_type = "ABM",
        abm_pop_size = 100,
        abm_delta_t = 0.1,
        abm_record_interval = 1,
        abm_seed = 1
      )
    ),
    "align exactly with the ABM step grid"
  )

  expect_error(
    alfakR::run_abm_simulation_grf(
      centroids = matrix(c(2, 2, 3, 1), ncol = 2, byrow = TRUE),
      lambda = 1,
      p = 0.01,
      times = c(0, 0.15),
      x0 = c("2.2" = 1),
      abm_pop_size = 100,
      abm_delta_t = 0.1,
      abm_record_interval = 1,
      abm_seed = 1
    ),
    "align exactly with the ABM step grid"
  )
})

test_that("transition treatment ABM returns tau2 and endpoint summaries", {
  res <- alfakR:::run_transition_treatment_abm(
    initial_population_r = stats::setNames(list(50, 50), c("2.2", "3.1")),
    untreated_fitness_map_r = stats::setNames(as.list(c(0.2, 0.1, 0.05)), c("2.2", "3.1", "2.3")),
    treated_fitness_map_r = stats::setNames(as.list(c(0.05, 0.2, 0.3)), c("2.2", "3.1", "2.3")),
    adjacency_r = list("2.2" = c("3.1", "2.3"), "3.1" = c("2.2", "2.3"), "2.3" = c("2.2", "3.1")),
    transition_karyotypes_r = "2.3",
    p_missegregation = 0.02,
    base_death_rate = 0.01,
    base_birth_rate = 0.02,
    fitness_birth_scale = 0.1,
    second_treatment_strength = 0.9,
    tau1_step = 31L,
    dt = 1,
    n_steps = 60L,
    record_interval = 10L,
    seed = 1L
  )

  expect_named(res, c(
    "tau2_step", "tau2_time", "tau2_transition_fraction",
    "condition1", "condition2", "metrics", "endpoint_summary"
  ))
  expect_true(res$tau2_step >= 31L)
  expect_true(is.data.frame(res$metrics))
  expect_true(is.data.frame(res$endpoint_summary))
  expect_true(all(c("total_population", "diversity", "transition_fraction") %in% names(res$metrics)))
  expect_true(is.finite(res$endpoint_summary$population_reduction_pct))
})

test_that("transition treatment wrapper initializes on untreated peaks and targets top ranked transitions", {
  node_metadata <- data.frame(
    karyotype = c("2.2", "3.1", "2.3", "1.3"),
    untreated_fitness = c(0.2, 0.1, 0.05, 0.03),
    treated_fitness = c(0.05, 0.2, 0.3, 0.25),
    is_untreated_peak = c(TRUE, TRUE, FALSE, FALSE),
    transition_rank = c(NA, NA, 2, 1),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = c("2.2", "3.1", "2.3"),
    to = c("2.3", "2.3", "1.3"),
    stringsAsFactors = FALSE
  )

  captured <- NULL
  testthat::with_mocked_bindings(
    {
      res <- alfakR::run_transition_karyotype_abm(
        node_metadata = node_metadata,
        edges = edges,
        times = c(0, 1),
        tau1 = 0,
        transition_top_n = 1,
        abm_delta_t = 1,
        abm_seed = 1
      )
      expect_identical(res$inputs$transition_karyotypes, "1.3")
      expect_identical(res$inputs$initial_untreated_peaks, c("2.2", "3.1"))
    },
    run_transition_treatment_abm = function(initial_population_r, untreated_fitness_map_r,
                                            treated_fitness_map_r, adjacency_r,
                                            transition_karyotypes_r, ...) {
      captured <<- list(
        initial_population_r = initial_population_r,
        transition_karyotypes_r = transition_karyotypes_r,
        adjacency_r = adjacency_r
      )
      list(
        tau2_step = 1L,
        tau2_time = 1,
        tau2_transition_fraction = 0,
        condition1 = list("0" = stats::setNames(c(50, 50), c("2.2", "3.1"))),
        condition2 = list("0" = stats::setNames(c(50, 50), c("2.2", "3.1"))),
        metrics = data.frame(),
        endpoint_summary = data.frame()
      )
    },
    .package = "alfakR"
  )

  expect_identical(names(captured$initial_population_r), c("2.2", "3.1"))
  expect_equal(as.numeric(unlist(captured$initial_population_r)), c(1, 1))
  expect_identical(captured$transition_karyotypes_r, "1.3")
  expect_true("2.2" %in% names(captured$adjacency_r))
  expect_true("2.3" %in% captured$adjacency_r[["2.2"]])
})

test_that("transition treatment wrapper passes tau1 at the report treatment step", {
  node_metadata <- data.frame(
    karyotype = c("2.2", "3.1", "2.3"),
    untreated_fitness = c(0.2, 0.1, 0.05),
    treated_fitness = c(0.05, 0.2, 0.3),
    is_untreated_peak = c(TRUE, TRUE, FALSE),
    transition_rank = c(NA, NA, 1),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = c("2.2", "3.1"),
    to = c("2.3", "2.3"),
    stringsAsFactors = FALSE
  )
  captured_tau1_step <- NULL

  testthat::with_mocked_bindings(
    {
      res <- alfakR::run_transition_karyotype_abm(
        node_metadata = node_metadata,
        edges = edges,
        times = c(0, 30),
        tau1 = 30,
        transition_top_n = 1,
        abm_delta_t = 1,
        abm_seed = 1
      )
      expect_identical(res$inputs$tau1_step, 30L)
      expect_identical(res$inputs$abm_pop_size, 2L)
      expect_null(res$inputs$abm_pop_size_requested)
    },
    run_transition_treatment_abm = function(initial_population_r, untreated_fitness_map_r,
                                            treated_fitness_map_r, adjacency_r,
                                            transition_karyotypes_r, tau1_step, ...) {
      captured_tau1_step <<- tau1_step
      list(
        tau2_step = 30L,
        tau2_time = 30,
        tau2_transition_fraction = 0,
        condition1 = list("0" = stats::setNames(c(1, 1), c("2.2", "3.1"))),
        condition2 = list("0" = stats::setNames(c(1, 1), c("2.2", "3.1"))),
        metrics = data.frame(),
        endpoint_summary = data.frame()
      )
    },
    .package = "alfakR"
  )

  expect_identical(captured_tau1_step, 30L)
})

test_that("chrmod kernels validate transition-matrix shape and still accept legal inputs", {
  expect_error(
    alfakR:::chrmod_cpp(0, c(0.5, 0.5), list()),
    "parms\\$A"
  )
  expect_error(
    alfakR:::chrmod_cpp(0, c(0.5, 0.5), list(A = matrix(1, nrow = 2, ncol = 3))),
    "square matrix"
  )
  expect_error(
    alfakR:::chrmod_rel_cpp(0, c(0.5, 0.5), list(A = diag(3))),
    "length\\(x\\)"
  )
  expect_no_error(alfakR:::chrmod_cpp(0, c(0.5, 0.5), list(A = diag(2))))
  expect_no_error(alfakR:::chrmod_rel_cpp(0, c(0.5, 0.5), list(A = diag(2))))
})

test_that("alfak_neighbor_objective_cpp rejects non-integer or impossible observed counts", {
  base_args <- list(
    fc_param = 0.1,
    parent_fitness = 0.2,
    pij_values = 0.1,
    parent_birth_times = 0,
    timepoints = c(0, 1),
    parent_xfit = matrix(c(0.5, 0.5), nrow = 1),
    child_obs = c(0, 1),
    ntot = c(10, 10),
    parent_fitness_mean = 0.2,
    prior_mean = 0,
    prior_sd = 0.1,
    do_prior = FALSE,
    tol = 1e-8
  )

  expect_error(
    do.call(alfakR:::alfak_neighbor_objective_cpp, modifyList(base_args, list(child_obs = c(0, 1.5)))),
    "integer-valued"
  )
  expect_error(
    do.call(alfakR:::alfak_neighbor_objective_cpp, modifyList(base_args, list(ntot = c(10, 10.5)))),
    "integer-valued"
  )
  expect_error(
    do.call(alfakR:::alfak_neighbor_objective_cpp, modifyList(base_args, list(child_obs = c(0, 11)))),
    "must not exceed"
  )
})

test_that("parse_karyotype_ids rejects components beyond the supported integer range", {
  expect_error(
    alfakR:::parse_karyotype_ids(c("2147483648.2", "2.2")),
    "exceeds supported integer range"
  )
})
