# ======================================================================
# prediction_functions.R
#
# Combines ODE and ABM simulation methods for karyotype evolution.
# Provides a master function `predict_evo` to dispatch simulations.
# ======================================================================
# -------------------------------------------------------------
# 1. Karyotype String Parsing and Utility Functions
# -------------------------------------------------------------

#' Parse Karyotype Strings to Integer Matrix (Internal)
#'
#' Converts a character vector of karyotype strings (e.g., "2.2.1...")
#' into an integer matrix. All strings must represent karyotypes with the
#' same number of chromosome types. The number of types is inferred from
#' the first string.
#'
#' @param kvec Character vector of karyotype strings.
#' @return An integer matrix (N_karyotypes x N_chromosome_types).
#' @keywords internal
#' @noRd
parse_karyotypes <- function(kvec) {
  parse_karyotype_ids(kvec)
}

#' Convert Integer Vector Karyotype to String Tag (Internal)
#' @param v Integer vector representing a karyotype.
#' @return A single string representation (e.g., "2.2.1...").
#' @keywords internal
#' @noRd
vec_to_tag <- function(v) {
  paste(v, collapse = ".")
}

validate_finite_numeric_vector <- function(x, name, allow_empty = FALSE) {
  if (!is.numeric(x) || (!allow_empty && length(x) == 0) || any(!is.finite(x))) {
    qualifier <- if (allow_empty) "a numeric vector containing only finite values" else "a non-empty numeric vector containing only finite values"
    stop(sprintf("`%s` must be %s.", name, qualifier), call. = FALSE)
  }
  invisible(NULL)
}

validate_scalar_finite_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
    stop(sprintf("`%s` must be a single finite numeric value.", name), call. = FALSE)
  }
  invisible(NULL)
}

validate_probability_closed <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("`%s` must be a single finite numeric value in [0, 1].", name), call. = FALSE)
  }
  invisible(NULL)
}

validate_times_vector <- function(times, name = "times", non_negative = FALSE) {
  validate_finite_numeric_vector(times, name)
  if (is.unsorted(times, strictly = FALSE)) {
    stop(sprintf("`%s` must be sorted in non-decreasing order.", name), call. = FALSE)
  }
  if (non_negative && any(times < 0)) {
    stop(sprintf("`%s` must contain only non-negative timepoints for ABM simulation.", name), call. = FALSE)
  }
  invisible(NULL)
}

validate_integerish_scalar <- function(x, name, min_value = NULL, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != floor(x)) {
    stop(sprintf("`%s` must be a single finite integer-like value.", name), call. = FALSE)
  }
  if (!is.null(min_value)) {
    if (allow_zero) {
      if (x < min_value) {
        stop(sprintf("`%s` must be >= %s.", name, format(min_value, trim = TRUE)), call. = FALSE)
      }
    } else if (x <= min_value) {
      stop(sprintf("`%s` must be > %s.", name, format(min_value, trim = TRUE)), call. = FALSE)
    }
  }
  invisible(NULL)
}

validate_cpp_integerish_scalar <- function(x, name, min_value = NULL, allow_negative_one = FALSE,
                                           max_value = ALFAK_MAX_EXACT_INTEGER, target = c("long long", "int")) {
  target <- match.arg(target)
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != floor(x)) {
    stop(sprintf("`%s` must be a single finite integer-valued scalar.", name), call. = FALSE)
  }
  if (abs(x) > max_value) {
    stop(sprintf("`%s` must be exactly representable in R and no larger than %.0f in magnitude.", name, max_value), call. = FALSE)
  }
  if (!is.null(min_value) && x < min_value && !(allow_negative_one && x == -1)) {
    stop(sprintf("`%s` must be >= %s.", name, format(min_value, trim = TRUE)), call. = FALSE)
  }
  if (target == "int" && (x < -.Machine$integer.max - 1 || x > .Machine$integer.max)) {
    stop(sprintf("`%s` exceeds the supported C++ int range.", name), call. = FALSE)
  }
  invisible(NULL)
}

largest_remainder_allocate <- function(prob, total_size) {
  validate_finite_numeric_vector(prob, "prob")
  original_names <- names(prob)
  if (any(prob < 0)) {
    stop("`prob` must contain only non-negative values.", call. = FALSE)
  }
  validate_cpp_integerish_scalar(total_size, "total_size", min_value = 0, target = "long long")
  if (!length(prob)) {
    return(integer(0))
  }
  if (sum(prob) <= 0) {
    stop("`prob` must sum to a positive value.", call. = FALSE)
  }
  prob <- prob / sum(prob)
  raw <- prob * total_size
  counts <- floor(raw)
  remaining <- as.integer(total_size - sum(counts))
  if (remaining > 0) {
    fractional <- raw - counts
    order_idx <- order(-fractional, seq_along(fractional))
    counts[order_idx[seq_len(remaining)]] <- counts[order_idx[seq_len(remaining)]] + 1
  }
  counts <- as.integer(counts)
  if (!is.null(original_names)) {
    names(counts) <- original_names
  }
  counts
}

validate_named_frequency_vector <- function(x, expected_names = NULL, expected_dim = NULL, name = "x0") {
  validate_finite_numeric_vector(x, name)
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    stop(sprintf("`%s` must be a named numeric vector.", name), call. = FALSE)
  }
  if (anyDuplicated(names(x))) {
    stop(sprintf("`%s` must not contain duplicate names.", name), call. = FALSE)
  }
  if (any(x < 0)) {
    stop(sprintf("`%s` must contain only non-negative values.", name), call. = FALSE)
  }
  if (abs(sum(x) - 1) > 1e-6) {
    stop(sprintf("`%s` must sum to 1 (within tolerance).", name), call. = FALSE)
  }
  parsed <- parse_karyotype_ids(names(x))
  if (!is.null(expected_dim) && ncol(parsed) != expected_dim) {
    stop(sprintf("Karyotype names in `%s` do not have %d chromosome counts.", name, expected_dim), call. = FALSE)
  }
  if (!is.null(expected_names) && !setequal(names(x), expected_names)) {
    stop(sprintf("Names of `%s` must exactly match the karyotypes in `lscape$k`.", name), call. = FALSE)
  }
  invisible(parsed)
}

# -------------------------------------------------------------
# ODE Simulation Components
# -------------------------------------------------------------

#' Build Sparse Transition Matrix W using Rcpp Helpers
#' @param karyotype_strings Character vector of unique karyotype strings.
#' @param p Single chromosome missegregation probability.
#' @param Nmax Optional limit to the number of missegregations allowable
#' @return A sparse matrix (dgCMatrix) representing W (Parent x Daughter transitions).
#' @importFrom Matrix sparseMatrix
#' @export
#' @examples
#' \dontrun{
#'   k_str =c("2.2.1.2","2.2.2.3","2.2.3.2")
#'   A <- get_A(k_str,0.01)
#'   A <- get_A(k_str,0.01,Nmax=1) ## approximation of 1 missegregation case
#' }
#' 
build_W_rcpp <- function(karyotype_strings, p,Nmax=Inf) {
  parse_karyotype_ids(karyotype_strings)
  validate_probability_closed(p, "p")
  if (!(is.infinite(Nmax) || (is.numeric(Nmax) && length(Nmax) == 1 && is.finite(Nmax) && Nmax >= 0))) {
    stop("`Nmax` must be Inf or a non-negative finite number.", call. = FALSE)
  }
  w_structure <- tryCatch(
    get_A_inputs(karyotype_strings,p,Nmax), # C++ function
    error = function(e) stop("Error in get_A_inputs(C++ call): ", e$message, call. = FALSE)
  )
  W <- tryCatch(
    Matrix::sparseMatrix(i = w_structure$i, j = w_structure$j, x = w_structure$x,
                         dims = w_structure$dims, dimnames = w_structure$dimnames,
                         repr = "C"), 
    error = function(e) stop("Error creating sparseMatrix from Rcpp results: ", e$message, call. = FALSE)
  )
  W
}

#' A function suitable for use with `deSolve::ode`.
#' @param time 
#' @param x state vector.
#' @param parms A list containing the transition matrix, named A.
#' @return derivatives
#' @keywords internal
#' @noRd
chrmod_rel <- function(time, x, parms) {
  A <- parms$A
  g <- as.numeric(x %*% A)        # per‐type contributions
  phi <- sum(g)                   # average growth
  list(g - x * phi)               # dx/dt
}

#' Run ODE Simulation (Internal)
#' @param lscape Data frame with 'k' (karyotype string) and 'mean' (fitness).
#' @param p Missegregation probability.
#' @param times Numeric vector of time points for output.
#' @param x0 Named numeric vector of initial frequencies (must sum to 1).
#' @param ode_method Solver method for `deSolve::ode`.
#' @param Nmax Optional limit to the number of missegregations allowable.
#' @return Data frame with 'time' column and frequency columns for each karyotype.
#' @importFrom deSolve ode
#' @keywords internal
#' @noRd
run_ode_simulation <- function(lscape, p, times, x0, ode_method,Nmax=Inf) {
  message("Setting up ODE simulation...")
  k_strings <- lscape$k
  r_values <- lscape$mean 
  if(any(r_values < 0, na.rm = TRUE)) warning("Negative fitness values found in lscape$mean; ODE behavior might be unexpected.", call. = FALSE)
  
  message("Building sparse W matrix via Rcpp...")
  W <- build_W_rcpp(k_strings, p, Nmax) 
  W <- W*r_values
  ode_output <- tryCatch( 
    deSolve::ode(y = x0, times = times, func = chrmod_rel, parms = list(A=W), method = ode_method),
    error = function(e) stop(sprintf("ODE integration failed with method '%s': %s", ode_method, e$message), call. = FALSE)
  )
  
  message("ODE simulation finished.")
  res_df <- as.data.frame(ode_output)
  
  if(anyNA(res_df[,-1, drop=FALSE])) {
    warning("NA values detected in ODE results. This might indicate numerical instability or issues with the ODE system/parameters.", call. = FALSE)
  }
  final_sums <- rowSums(res_df[,-1, drop=FALSE], na.rm = TRUE) # Use na.rm=TRUE for sum if NAs are present
  if(any(abs(final_sums - 1) > 1e-5, na.rm = TRUE)) { 
    warning("ODE row sums deviate significantly from 1. Check model stability and parameters.", call. = FALSE)
  }
  
  res_df
}

# -------------------------------------------------------------
# ABM Simulation Components
# -------------------------------------------------------------

#' Run Agent-Based Model (ABM) Simulation (Internal)
#' @param lscape Data frame with 'k' (karyotype string) and 'mean' (fitness).
#' @param p Missegregation probability.
#' @param times Numeric vector of target time points for output.
#' @param x0 Named numeric vector of initial frequencies (must sum to 1).
#' @param abm_pop_size Initial total population size.
#' @param abm_delta_t Time duration of a single ABM step.
#' @param abm_max_pop Maximum population size (carrying capacity), <= 0 for unlimited.
#' @param abm_culling_survival Survival fraction if max_pop is exceeded.
#' @param abm_record_interval Record state every N steps.
#' @param abm_seed RNG seed (-1 for random).
#' @return Data frame with 'time' column and frequency columns for each karyotype.
#' @importFrom stats setNames
#' @importFrom tidyr pivot_wider
#' @keywords internal
#' @noRd
run_abm_simulation <- function(lscape, p, times, x0, abm_pop_size, abm_delta_t,
                               abm_max_pop, abm_culling_survival,
                               abm_record_interval, abm_seed) {
  
  message("Setting up ABM simulation...")
  if (!is.data.frame(lscape) || !all(c("k", "mean") %in% names(lscape))) {
    stop("'lscape' must be a data.frame with columns 'k' and 'mean'.", call. = FALSE)
  }
  if (!is.character(lscape$k) || length(lscape$k) == 0) {
    stop("'lscape$k' must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(lscape$k)) {
    stop("'lscape$k' must contain unique karyotype strings.", call. = FALSE)
  }
  validate_finite_numeric_vector(lscape$mean, "lscape$mean")
  if (length(lscape$mean) != length(lscape$k)) {
    stop("'lscape$mean' must have the same length as 'lscape$k'.", call. = FALSE)
  }
  parse_karyotypes(lscape$k)
  validate_probability_closed(p, "p")
  validate_times_vector(times, non_negative = TRUE)
  validate_named_frequency_vector(x0, expected_names = lscape$k)
  validate_cpp_integerish_scalar(abm_pop_size, "abm_pop_size", min_value = 1, target = "long long")
  validate_positive_finite(abm_delta_t, "abm_delta_t")
  validate_cpp_integerish_scalar(abm_max_pop, "abm_max_pop", min_value = 0, target = "long long")
  validate_probability_closed(abm_culling_survival, "abm_culling_survival")
  validate_cpp_integerish_scalar(abm_record_interval, "abm_record_interval", min_value = 1, target = "int")
  if (abm_record_interval == 0) {
    stop("`abm_record_interval` must not be zero.", call. = FALSE)
  }
  validate_cpp_integerish_scalar(abm_seed, "abm_seed", min_value = 0, allow_negative_one = TRUE, target = "int")
  
  initial_counts <- largest_remainder_allocate(x0, abm_pop_size)
  if(any(initial_counts < 0)) {
    warning("Negative counts generated for ABM initial population after rounding; treating as 0.", call.=FALSE)
    initial_counts[initial_counts < 0] <- 0 
  }
  initial_pop_list_unfiltered <- as.list(initial_counts)
  initial_pop_list <- initial_pop_list_unfiltered[initial_counts > 0] 
  
  if(length(initial_pop_list) == 0) stop("Initial population for ABM is zero after filtering zero counts.", call. = FALSE)
  
  fitness_map_list <- stats::setNames(as.list(lscape$mean), lscape$k)
  
  max_time <- max(times, na.rm = TRUE) 
  num_steps <- ceiling(max_time / abm_delta_t)
  if (!is.finite(num_steps) || num_steps < 0) stop("Number of ABM steps is invalid (max_time / abm_delta_t). Check 'times' and 'abm_delta_t'.", call. = FALSE)
  if (num_steps > 1e7) {
    warning("ABM simulation requires a very large number of steps; check `times` and `abm_delta_t`.", call. = FALSE)
  }
  if (num_steps > .Machine$integer.max) {
    stop("Computed number of ABM steps exceeds the supported integer range.", call. = FALSE)
  }
  
  message(sprintf("Starting ABM simulation for %d steps (up to time %.2f)...", num_steps, max_time))
  sim_results_list_cpp <- tryCatch(
    run_karyotype_abm(
      initial_population_r = initial_pop_list,
      fitness_map_r        = fitness_map_list,
      p_missegregation     = p,
      dt                   = abm_delta_t,
      n_steps              = as.integer(num_steps),
      max_population_size  = abm_max_pop,
      culling_survival_fraction = abm_culling_survival,
      record_interval      = as.integer(abm_record_interval),
      seed                 = as.integer(abm_seed),
      grf_centroids        = matrix(numeric(0), nrow = 0, ncol = 0), # CORRECT R equivalent
      grf_lambda           = NA_real_                                # R equivalent
    ), error = function(e) {
      # It's useful to see the full error structure from Rcpp if possible
      print(e) # print the full error object
      stop("Error during run_karyotype_abm (C++ call) execution: ", e$message, call. = FALSE)
    }
  )
  message("ABM simulation finished.")
  
  message("Processing ABM results...")
  if (length(sim_results_list_cpp) == 0) {
    warning("ABM simulation returned no results from C++.", call. = FALSE)
    empty_df_res <- data.frame(time = numeric(0))
    ktypes_all <- names(x0)
    for (kt_name in ktypes_all) empty_df_res[[kt_name]] <- numeric(0)
    return(empty_df_res)
  }
  
  # Convert list of named vectors (from C++) to a long data frame
  results_df_list <- lapply(names(sim_results_list_cpp), function(step_name_str) {
    counts_vec <- sim_results_list_cpp[[step_name_str]]
    step_num <- as.integer(step_name_str) 
    time_point <- step_num * abm_delta_t
    
    if (length(counts_vec) > 0 && sum(counts_vec, na.rm=TRUE) > 0) { 
      total_count <- sum(counts_vec, na.rm=TRUE)
      freq_vec <- counts_vec / total_count
      karyo_names <- names(freq_vec)
      if(is.null(karyo_names) && length(freq_vec) > 0) { # Should have names from C++
        warning(paste0("Step ", step_name_str, ": ABM counts vector missing names. Using V1, V2..."), call.=FALSE)
        karyo_names <- paste0("V", seq_along(freq_vec))
      }
      
      data.frame(time = time_point, 
                 Karyotype = karyo_names, 
                 Frequency = as.numeric(freq_vec),
                 stringsAsFactors = FALSE)
    } else { 
      data.frame(time = time_point, Karyotype = character(0), Frequency = numeric(0),
                 stringsAsFactors = FALSE)
    }
  })
  results_long_df <- do.call(rbind, results_df_list)
  
  all_karyotypes_initial <- names(x0) 
  
  if (nrow(results_long_df) > 0 && "Karyotype" %in% names(results_long_df)) { # Check Karyotype col exists
    results_wide_df <- tidyr::pivot_wider(results_long_df,
                                          names_from = "Karyotype",
                                          values_from = "Frequency",
                                          values_fill = 0.0) 
    
    missing_cols <- setdiff(all_karyotypes_initial, names(results_wide_df))
    if (length(missing_cols) > 0) {
      for(col_name in missing_cols) results_wide_df[[col_name]] <- 0.0
    }
    
    # Ensure "time" column is first, then others
    time_col_present <- "time" %in% names(results_wide_df)
    if(!time_col_present && nrow(results_wide_df) > 0) stop("Internal error: 'time' column lost during pivot_wider in ABM processing.", call. = FALSE)
    
    # Select and order columns
    final_col_order <- intersect(c("time", all_karyotypes_initial), names(results_wide_df))
    # Add back any karyotypes from all_karyotypes_initial that might have been completely absent in results
    # (already handled by missing_cols loop mostly)
    
    results_final_df <- results_wide_df[, final_col_order, drop = FALSE]
    
  } else {
    message("ABM processing resulted in empty data frame or no 'Karyotype' column; returning structure based on initial times and karyotypes.")
    results_final_df <- data.frame(time = if(length(times) > 0) unique(times) else numeric(0))
    for (kt_name in all_karyotypes_initial) results_final_df[[kt_name]] <- 0.0
    if (nrow(results_final_df) == 0 && length(times) == 0) { # Truly empty case
      # Construct an empty df with correct column names if all_karyotypes_initial is also empty
      col_names_for_empty <- "time"
      if(length(all_karyotypes_initial) > 0) col_names_for_empty <- c("time", all_karyotypes_initial)
      results_final_df <- data.frame(matrix(ncol = length(col_names_for_empty), nrow = 0,
                                            dimnames=list(NULL, col_names_for_empty)))
    }
  }
  results_final_df
}
# -------------------------------------------------------------
# Master Prediction Function
# -------------------------------------------------------------

#' Predict Karyotype Evolution using ODE or ABM
#'
#' Master function to run karyotype evolution simulations. Assumes Rcpp functions
#' are compiled and available from the package.
#'
#' @param lscape data.frame. Must contain columns 'k' (character, karyotype string,
#'   e.g., "2.2.1...") and 'mean' (numeric, fitness r_k). All karyotype strings
#'   must represent karyotypes with the same number of chromosome types.
#' @param p numeric. Probability of single chromosome missegregation per division (0 <= p <= 1).
#' @param times numeric vector. Time points at which to report frequencies. Must be sorted ascending and non-empty.
#' @param x0 numeric vector. Initial frequencies of karyotypes. Must be named with
#'   karyotype strings matching all unique karyotypes present in `lscape$k`.
#'   The vector will be reordered to match `lscape$k` if names are a setequal match. Must sum to 1.
#' @param prediction_type character. Either "ODE" or "ABM". Default is "ODE".
#' @param ode_method character. Method for `deSolve::ode` (e.g., "lsoda", "lsodes").
#'   Used only if `prediction_type` is "ODE". Default is "lsoda".
#' @param abm_pop_size integer. Initial total population size for ABM.
#'   Used only if `prediction_type` is "ABM". Default is 1e4.
#' @param abm_delta_t numeric. Duration of one time step in ABM.
#'   Used only if `prediction_type` is "ABM". Default is 0.1.
#' @param abm_max_pop numeric. Carrying capacity for ABM (values <= 0 typically mean unlimited population growth,
#'   as handled by the C++ function). Used only if `prediction_type` is "ABM". Default is 1e7.
#' @param abm_culling_survival numeric. Fraction of population surviving when `abm_max_pop`
#'   is exceeded (0 <= x <= 1). Used only if `prediction_type` is "ABM". Default is 0.1.
#' @param abm_record_interval integer. Record ABM state every N steps.
#'   Used only if `prediction_type` is "ABM". Default is 10.
#' @param Nmax Optional limit to the number of missegregations allowable (ODE model only).
#' @param abm_seed integer. Seed for ABM's random number generator (-1 for a random seed based on device,
#'   any other integer for a fixed seed). Used only if `prediction_type` is "ABM". Default is -1.
#'
#' @return A data.frame with the first column 'time' and subsequent columns named by
#'   karyotype strings from `lscape$k`, containing their frequencies at the requested `times` (for ODE)
#'   or at recorded ABM steps that cover the requested time range.
#'
#' @export
#' @examples
#' \dontrun{
#' # This example requires the C++ functions to be compiled and available
#' # within the alfakR package.
#'
#' # 1. Define a landscape (fitness for different karyotypes)
#' # Ensuring all karyotype strings have the same number of elements.
#' landscape_df <- data.frame(
#'   k = c("2.2", "3.1", "1.3"), # Example for 2 chromosome types
#'   mean = c(0.1, 0.12, 0.08)  # Fitness values (r_k)
#' )
#'
#' # 2. Define initial frequencies
#' initial_freq <- c(0.9, 0.05, 0.05)
#' names(initial_freq) <- landscape_df$k # Ensure names match
#'
#' # 3. Define simulation parameters
#' prob_missegregation <- 0.001
#' simulation_times <- seq(0, 50, by = 1)
#'
#' # Ensure alfakR is loaded (e.g. after devtools::load_all() or library(alfakR))
#' # If running this example outside of a context where alfakR is loaded,
#' # the C++ functions (run_karyotype_abm etc.) will not be found.
#'
#'
#' ode_results <- predict_evo(lscape = landscape_df,
#'                            p = prob_missegregation,
#'                            times = simulation_times,
#'                            x0 = initial_freq,
#'                            prediction_type = "ODE")
#' print(head(ode_results))
#' abm_results <- predict_evo(lscape = landscape_df,
#'                            p = prob_missegregation,
#'                            times = simulation_times,
#'                            x0 = initial_freq,
#'                            prediction_type = "ABM",
#'                            abm_pop_size = 1000, # Smaller for quick example
#'                            abm_delta_t = 0.1,
#'                            abm_max_pop = 10000,
#'                            abm_record_interval = 10)
#' print(head(abm_results))
#' }
predict_evo <- function(lscape, p, times, x0, prediction_type = "ODE",
                        ode_method = "lsoda",
                        Nmax=Inf,
                        abm_pop_size = 1e4, abm_delta_t = 0.1, abm_max_pop = 1e7,
                        abm_culling_survival = 0.1, abm_record_interval = 10,
                        abm_seed = -1) {
  
  # --- Input Validation ---
  if (!is.data.frame(lscape) || !all(c("k", "mean") %in% names(lscape))) {
    stop("'lscape' must be a data.frame with columns 'k' (karyotype strings) and 'mean' (fitness).", call. = FALSE)
  }
  if (!is.character(lscape$k) || length(lscape$k) == 0) {
    stop("'lscape$k' must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(lscape$k)) {
    stop("'lscape$k' contains duplicate karyotype strings. Please provide unique karyotypes in lscape.", call. = FALSE)
  }
  validate_finite_numeric_vector(lscape$mean, "lscape$mean")
  if (length(x0) != nrow(lscape)) {
    stop("Length of 'x0' must match the number of rows (karyotypes) in 'lscape'.", call. = FALSE)
  }
  validate_probability_closed(p, "p")
  validate_times_vector(times)
  validate_named_frequency_vector(x0, expected_names = lscape$k)
  parse_karyotypes(lscape$k)
  # Reorder x0 to match lscape$k order for consistency downstream
  x0 <- x0[lscape$k] 
  if(anyNA(x0)){ # Check after reordering if any names in lscape$k were not in original x0 names
    stop("Mismatch between names(x0) and lscape$k led to NAs after reordering x0. Ensure all lscape$k are in names(x0).", call. = FALSE)
  }
  
  if (!(prediction_type %in% c("ODE", "ABM"))) stop("'prediction_type' must be either 'ODE' or 'ABM'.", call. = FALSE)
  
  result_df <- NULL 
  if (prediction_type == "ODE") {
    result_df <- run_ode_simulation(lscape = lscape, p = p, times = times, x0 = x0, ode_method = ode_method,Nmax=Nmax)
  } else if (prediction_type == "ABM") {
    # ABM specific parameter validation
    validate_times_vector(times, non_negative = TRUE)
    validate_cpp_integerish_scalar(abm_pop_size, "abm_pop_size", min_value = 1, target = "long long")
    validate_positive_finite(abm_delta_t, "abm_delta_t")
    validate_cpp_integerish_scalar(abm_max_pop, "abm_max_pop", min_value = 0, target = "long long")
    validate_probability_closed(abm_culling_survival, "abm_culling_survival")
    validate_cpp_integerish_scalar(abm_record_interval, "abm_record_interval", min_value = 1, target = "int")
    if (abm_record_interval == 0) {
      stop("'abm_record_interval' must not be zero.", call. = FALSE)
    }
    validate_cpp_integerish_scalar(abm_seed, "abm_seed", min_value = 0, allow_negative_one = TRUE, target = "int")
    
    result_df <- run_abm_simulation(lscape = lscape, p = p, times = times, x0 = x0,
                                    abm_pop_size = abm_pop_size, abm_delta_t = abm_delta_t,
                                    abm_max_pop = abm_max_pop, abm_culling_survival = abm_culling_survival,
                                    abm_record_interval = abm_record_interval, abm_seed = abm_seed)
  } else {
    # This case should be caught by the earlier check of prediction_type
    stop("Internal error: Unknown prediction_type somehow passed validation.", call. = FALSE)
  }
  
  message("Simulation process complete.")
  return(result_df)
}

#' Find Steady State Karyotype Distribution
#'
#' Calculates the steady-state frequency distribution of karyotypes by finding the
#' dominant eigenvector of the transition matrix M = t(W) %*% diag(r).
#' Assumes Rcpp functions for W matrix construction are compiled and available from the package.
#'
#' @param lscape data.frame. Must contain columns 'k' (character, karyotype string,
#'   e.g., "2.2.1...") and 'mean' (numeric, fitness r_k). All karyotype strings
#'   in `lscape$k` must be unique and represent karyotypes with the same number of chromosome types.
#' @param p numeric. Probability of single chromosome missegregation per division (0 <= p <= 1).
#' @param Nmax Optional limit to the number of missegregations allowable.
#' @return A named numeric vector of steady-state frequencies, ordered according to
#'   `lscape$k`. Returns `NULL` if calculation fails (e.g., invalid inputs,
#'   eigen decomposition issues).
#' @export
#' @importFrom Matrix Diagonal sparseMatrix 
#' @importFrom RSpectra eigs
#' @importFrom methods is
#' @examples
#' \dontrun{
#' # This example requires the C++ functions to be compiled and available
#' # within the alfakR package.
#'
#' landscape_df <- data.frame(
#'   k = c("2.2", "3.1", "1.3"), # Example for 2 chromosome types
#'   mean = c(0.1, 0.12, 0.08)
#' )
#' prob_missegregation <- 0.001
#'
#' # if (requireNamespace("alfakR", quietly = TRUE) &&
#' #     exists("rcpp_prepare_W_structure", where = "package:alfakR", mode="function")) {
#' #
#' #   ss_dist <- find_steady_state(lscape = landscape_df, p = prob_missegregation)
#' #   if (!is.null(ss_dist)) {
#' #     print("Steady-state distribution:")
#' #     print(ss_dist)
#' #     print(paste("Sum of frequencies:", sum(ss_dist)))
#' #   }
#' # } else {
#' #  message("Skipping find_steady_state example: alfakR not loaded or Rcpp fns not found.")
#' #  message("Run devtools::load_all() in your package project first, or install the package.")
#' # }
#' }
find_steady_state <- function(lscape, p, Nmax=Inf) {
  
  # --- Input Validation ---
  if (!is.data.frame(lscape) || !all(c("k", "mean") %in% names(lscape))) {
    stop("Invalid 'lscape' input: must be a data.frame with 'k' and 'mean' columns.", call. = FALSE)
  }
  if (!is.character(lscape$k) || length(lscape$k) == 0) {
    stop("'lscape$k' must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(lscape$k)) { # Check for duplicates that would cause issues
    stop("'lscape$k' contains duplicate karyotype strings. Steady state calculation requires unique karyotypes in lscape.", call. = FALSE)
  }
  validate_finite_numeric_vector(lscape$mean, "lscape$mean")
  if (length(lscape$mean) != length(lscape$k)) {
    stop("'lscape$mean' must be numeric and have the same length as 'lscape$k'.", call. = FALSE)
  }
  validate_probability_closed(p, "p")
  
  # Validate karyotype strings before passing to C++
  # If parse_karyotypes fails, it will stop.
  parsed_k_validation <- parse_karyotypes(lscape$k) 
  
  k_strings <- lscape$k
  r_values <- lscape$mean 
  names(r_values) <- k_strings 
  N_k_states <- length(k_strings) 
  
  W <- tryCatch(build_W_rcpp(k_strings, p, Nmax), error = function(e) {
    # Pass error message from build_W_rcpp which might come from C++
    stop(sprintf("Failed to build W matrix in find_steady_state: %s", e$message), call. = FALSE);
  })
  # build_W_rcpp should stop on C++ error, so W should not be NULL here unless tryCatch changed.
  # But good practice to check.
  if (is.null(W) || !methods::is(W, "sparseMatrix")) { # Using methods::is for S4 objects
    stop("Internal error: build_W_rcpp did not return a valid sparse Matrix.", call. = FALSE)
  }
  
  # The dimnames of W from rcpp_prepare_W_structure should match k_strings
  # because k_strings are validated, duplicates would have stopped it,
  # and rcpp_prepare_W_structure uses the (now unique) input for dimnames.
  # Thus, r_values (named by original k_strings) should align.
  # A robust check:
  if (!identical(colnames(W), k_strings) || !identical(rownames(W), k_strings)) {
    stop("Internal error: Dimension names of W matrix do not match input lscape$k. This indicates an issue in C++ W matrix construction or R-side k_string handling.", call. = FALSE)
  }
  
  M_matrix <- W*r_values
  # chrmod_rel() evolves row vectors via x %*% M_matrix, so the steady state is
  # the dominant left eigenvector of M_matrix, equivalently the dominant right
  # eigenvector of t(M_matrix).
  M_eigs <- Matrix::t(M_matrix)
  
  eig_result <- NULL
  if (nrow(M_eigs) < 3) {
    dense_eigs <- eigen(as.matrix(M_eigs))
    dominant_idx <- which.max(Re(dense_eigs$values))
    eig_result <- list(
      values = dense_eigs$values[dominant_idx],
      vectors = dense_eigs$vectors[, dominant_idx, drop = FALSE]
    )
  } else {
    eig_result <- tryCatch(RSpectra::eigs(M_eigs, k = 1, which = "LR"), error = function(e_lr) {
      warning("RSpectra::eigs with 'LR' failed; falling back to dense eigen decomposition: ", e_lr$message, call. = FALSE)
      dense_eigs <- eigen(as.matrix(M_eigs))
      dominant_idx <- which.max(Re(dense_eigs$values))
      list(
        values = dense_eigs$values[dominant_idx],
        vectors = dense_eigs$vectors[, dominant_idx, drop = FALSE]
      )
    })
  }
  
  if (is.null(eig_result) || length(eig_result$vectors) == 0) { # Check vectors specifically
    stop("Eigen decomposition (RSpectra::eigs) failed to return eigenvectors.", call. = FALSE)
  }
  
  ss_vector <- Re(eig_result$vectors[, 1]) 
  
  # Normalize to be a probability distribution
  max_abs_val_idx <- which.max(abs(ss_vector))
  if(length(max_abs_val_idx) == 0) { # all ss_vector elements might be NA/NaN or vector is empty
    stop("Internal error: Could not determine sign for eigenvector normalization (all elements problematic or empty).", call.=FALSE)
  }
  ss_vector <- ss_vector * sign(ss_vector[max_abs_val_idx[1]]) 
  ss_vector[ss_vector < 1e-10] <- 0 
  
  sum_ss <- sum(ss_vector)
  if (is.na(sum_ss) || sum_ss <= 1e-9) { 
    stop(sprintf("Steady state vector sums to %s (close to zero or NA); cannot normalize. Eigenvalue was ~%f. Check input fitness values and transition rates.", 
                 as.character(sum_ss), Re(eig_result$values[1])), call. = FALSE)
  }
  ss_vector <- ss_vector / sum_ss
  
  # Ensure names are from lscape$k (which should match W's dimnames at this point)
  names(ss_vector) <- k_strings 
  
  return(ss_vector)
}

# -------------------------------------------------------------
# Public wrapper: ABM simulation under a GRF fitness landscape
# -------------------------------------------------------------

#' Simulate karyotype dynamics in a Gaussian‑Random‑Field (GRF) fitness landscape
#'
#' This function runs the C++ agent‑based model (`run_karyotype_abm`) in
#' **GRF mode**: fitness values are generated on‑the‑fly from a set of
#' centroid points and a wavelength λ, instead of using a lookup table.
#'
#' @param centroids Numeric matrix with one centroid per **row**  
#'   (dimensions *n_centroids × K*), where *K* is the number of
#'   chromosome types.
#' @param lambda Positive scalar. The GRF wavelength λ; smaller values give
#'   a more rugged landscape.
#' @param p Missegregation probability (0 ≤ \code{p} ≤ 1).
#' @param times Numeric vector of time points to sample.
#' @param x0 **Named** numeric vector of initial karyotype frequencies
#'   (must sum to 1).  The names must be karyotype strings of length *K*
#'   matching the dimension of \code{centroids}.
#' @param abm_pop_size Initial total population size.
#' @param abm_delta_t Duration of one ABM step.
#' @param abm_max_pop Carrying capacity. Use \code{<= 0} for unlimited.
#' @param abm_culling_survival Fraction of cells retained when the population
#'   exceeds \code{abm_max_pop}.
#' @param abm_record_interval Record population state every N steps. If a negative value 
#' is provided, then population is recorded every passage.
#' @param abm_seed RNG seed.  Use \code{-1} for a random seed.
#' @param normalize_freq Should ABM counts be normalized to frequencies?
#' @return A **wide data‑frame**: first column \code{time}, remaining columns
#'   one per karyotype, giving relative frequencies at each sampled time.
#'
#' @examples
#' # Two‑chromosome example with 4 centroids
#' cents <- matrix(c(2,2,
#'                   3,1,
#'                   1,3,
#'                   4,4), ncol = 2, byrow = TRUE)
#' init  <- c("2.2" = 1)        # start at (2,2)
#' times <- seq(0, 5, by = 1)
#'
#' sim <- run_abm_simulation_grf(
#'   centroids = cents, lambda = 2,
#'   p = 0.02, times = times, x0 = init,
#'   abm_pop_size = 5e3, abm_delta_t = 0.05,
#'   abm_record_interval = 20, abm_seed = 42
#' )
#' head(sim)
#'
#' @export
run_abm_simulation_grf <- function(centroids, lambda, p, times, x0,
                                   abm_pop_size        = 1e4,
                                   abm_delta_t         = 0.1,
                                   abm_max_pop         = 1e7,
                                   abm_culling_survival = 0.1,
                                   abm_record_interval  = 10,
                                   abm_seed             = -1,
                                   normalize_freq=T) {
  
  ## -- validation (same as internal draft, trimmed for brevity) -------------
  if(!is.matrix(centroids) || !is.numeric(centroids) || nrow(centroids) == 0)
    stop("'centroids' must be a non‑empty numeric matrix.", call. = FALSE)
  if (any(!is.finite(centroids))) {
    stop("'centroids' must contain only finite numeric values.", call. = FALSE)
  }
  validate_positive_finite(lambda, "lambda")
  validate_probability_closed(p, "p")
  validate_times_vector(times, non_negative = TRUE)
  validate_named_frequency_vector(x0, expected_dim = ncol(centroids))
  validate_cpp_integerish_scalar(abm_pop_size, "abm_pop_size", min_value = 1, target = "long long")
  validate_positive_finite(abm_delta_t, "abm_delta_t")
  validate_cpp_integerish_scalar(abm_max_pop, "abm_max_pop", min_value = 0, target = "long long")
  validate_probability_closed(abm_culling_survival, "abm_culling_survival")
  validate_cpp_integerish_scalar(abm_record_interval, "abm_record_interval", min_value = 1, target = "int")
  if (abm_record_interval == 0) {
    stop("'abm_record_interval' must not be zero.", call. = FALSE)
  }
  validate_cpp_integerish_scalar(abm_seed, "abm_seed", min_value = 0, allow_negative_one = TRUE, target = "int")
  validate_scalar_logical(normalize_freq, "normalize_freq")
  K <- ncol(centroids)
  
  ## -- initial population ----------------------------------------------------
  init_counts <- largest_remainder_allocate(x0, abm_pop_size)
  init_counts[init_counts < 0] <- 0
  init_list   <- as.list(init_counts)[init_counts > 0]  
  if(!length(init_list)) stop("Initial population is zero.", call. = FALSE)
  
  ## -- run C++ ---------------------------------------------------------------
  steps <- ceiling(max(times) / abm_delta_t)
  if (steps > 1e7) {
    warning("ABM simulation requires a very large number of steps; check `times` and `abm_delta_t`.", call. = FALSE)
  }
  if (!is.finite(steps) || steps < 0 || steps > .Machine$integer.max) {
    stop("Computed number of ABM steps exceeds the supported integer range.", call. = FALSE)
  }
  elapsed <- system.time({
    cpp_res <- run_karyotype_abm(
      initial_population_r      = init_list,
      fitness_map_r             = setNames(list(), character(0)),
      p_missegregation          = p,
      dt                        = abm_delta_t,
      n_steps                   = as.integer(steps),
      max_population_size       = abm_max_pop,
      culling_survival_fraction = abm_culling_survival,
      record_interval           = as.integer(abm_record_interval),
      seed                      = as.integer(abm_seed),
      grf_centroids             = centroids,
      grf_lambda                = lambda
    )
  })
  
  cat(sprintf("Simulation completed in %.2f seconds.\n", elapsed["elapsed"]))
  ## -- convert to wide data‑frame (unchanged) --------------------------------
  if(!length(cpp_res)) {
    warning("C++ returned no results.")
    out <- data.frame(time = numeric(0)); for(nm in names(x0)) out[[nm]] <- numeric(0)
    return(out)
  }
  long <- lapply(names(cpp_res), function(s) {
    t <- as.numeric(s) * abm_delta_t
    cnt <- cpp_res[[s]]
    if(length(cnt) && sum(cnt) > 0) {
      freq <- cnt 
      if(normalize_freq) freq <- cnt / sum(cnt)
      data.frame(time = t, Karyotype = names(freq), Frequency = as.numeric(freq))
    } else data.frame(time = t, Karyotype = character(0), Frequency = numeric(0))
  })
  long <- do.call(rbind, long)
  
  if(nrow(long)) {
    wide <- tidyr::pivot_wider(long, names_from = "Karyotype",
                               values_from = "Frequency", values_fill = 0)
    miss <- setdiff(names(x0), names(wide))
    for(m in miss) wide[[m]] <- 0
    karyo_cols <- setdiff(names(wide), "time")
    wide[, c("time", karyo_cols), drop = FALSE]
  } else {
    data.frame(time = unique(times),
               t(matrix(0, nrow = length(times), ncol = length(x0),
                        dimnames = list(NULL, names(x0)))))
  }
}


# ======================================================================
# End of prediction_functions.R
# ======================================================================
