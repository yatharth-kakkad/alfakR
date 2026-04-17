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
  if(!is.character(kvec)) stop("Input 'kvec' must be a character vector.", call. = FALSE)
  if(length(kvec) == 0) return(matrix(integer(0), ncol = 0, nrow = 0)) # Handle empty input
  
  split_chr <- strsplit(kvec, "[.]", perl = TRUE)
  
  # Determine expected length from the first non-empty, parsable string
  expected_len <- -1
  first_valid_k_for_len <- NULL
  for(i in seq_along(split_chr)){
    if(length(split_chr[[i]]) > 0 && !all(split_chr[[i]] == "")){
      expected_len <- length(split_chr[[i]])
      first_valid_k_for_len <- kvec[i]
      break
    }
  }
  
  if (expected_len == -1) { # All strings were empty or only dots
    if(all(sapply(kvec, function(s) s == "" || grepl("^\\.*$", s)))){
      # If all strings are truly empty or just dots, return empty matrix with 0 cols (or error)
      warning("All karyotype strings are empty or malformed (e.g. '...'). Cannot determine chromosome count.", call. = FALSE)
      return(matrix(integer(0), ncol = 0, nrow = length(kvec)))
    }
    stop("Cannot determine number of chromosome types from input 'kvec'. First element was problematic.", call. = FALSE)
  }
  if (expected_len == 0 && !is.null(first_valid_k_for_len) && nchar(first_valid_k_for_len)>0 ) {
    # This case implies the first string was e.g. "." which split into 0 length useful segments.
    stop(sprintf("First parsable karyotype string '%s' resulted in 0 chromosome types. Check format.", first_valid_k_for_len), call. = FALSE)
  }
  
  
  len_ok <- vapply(split_chr, length, integer(1)) == expected_len
  if (!all(len_ok)) {
    stop(sprintf("All karyotypes must have the same number of integers (expected %d, derived from first valid karyotype). Problem at input indices: %s",
                 expected_len, paste(which(!len_ok), collapse=", ")), call. = FALSE)
  }
  
  # Suppressing potential warnings from as.numeric on non-numbers, error handles it
  num_values <- unlist(split_chr, use.names = FALSE)
  if (length(num_values) == 0 && length(kvec) > 0) { # e.g. kvec was list("", "") and expected_len became 0
    # This state means expected_len was probably 0 from something like kvec = c("", "")
    # which should have been caught by earlier checks for expected_len.
    # If expected_len is > 0, then num_values cannot be empty if all len_ok.
    if(expected_len > 0) stop("Internal error: num_values became empty despite passing length checks.", call. = FALSE)
    # If expected_len was 0 (e.g. kvec=c("")), create 0-col matrix
    return(matrix(integer(0), nrow=length(kvec), ncol=0))
  }
  
  
  mat_num <- tryCatch(
    matrix(as.numeric(num_values), ncol = expected_len, byrow = TRUE),
    warning = function(w) {
      if (grepl("NAs introduced by coercion", w$message, fixed = TRUE)) {
        stop("Non-numeric values in karyotype strings led to NAs during numeric conversion.", call. = FALSE)
      }
      warning(w) 
      matrix(NA_real_, ncol = expected_len, nrow = length(split_chr)) 
    }
  )
  
  if (anyNA(mat_num)) stop("Parsing failed: Non-numeric values found or NAs introduced.", call. = FALSE)
  
  # Check if numbers are actual integers after as.numeric
  # Compare with rounded version. If not equal, then it wasn't an integer.
  if (!all(mat_num == round(mat_num), na.rm = TRUE)) { # Check for non-integer numbers
    stop("Parsing failed: Karyotype components must be whole numbers.", call. = FALSE)
  }
  
  mat <- apply(mat_num, 2, as.integer) 
  if (anyNA(mat) && !anyNA(mat_num)) { # Should not happen if previous check passes
    stop("Internal error: Karyotype components could not be coerced to integer without NA, despite being whole numbers.", call. = FALSE)
  }
  
  if (any(mat <= 0, na.rm = TRUE)) stop("Parsing failed: Karyotype components must be positive integers.", call. = FALSE)
  storage.mode(mat) <- "integer"
  mat
}

#' Convert Integer Vector Karyotype to String Tag (Internal)
#' @param v Integer vector representing a karyotype.
#' @return A single string representation (e.g., "2.2.1...").
#' @keywords internal
#' @noRd
vec_to_tag <- function(v) {
  paste(v, collapse = ".")
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
  
  initial_counts_raw <- x0 * abm_pop_size
  initial_counts <- round(initial_counts_raw) 
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
  if (num_steps <= 0) stop("Number of ABM steps is non-positive (max_time / abm_delta_t). Check 'times' and 'abm_delta_t'.", call. = FALSE)
  
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
                                          names_from = .data$Karyotype, 
                                          values_from = .data$Frequency,
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
  if (!is.numeric(lscape$mean) || length(lscape$mean) != length(lscape$k)) {
    stop("'lscape$mean' must be numeric and have the same length as 'lscape$k'.", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1 || p < 0 || p > 1) {
    stop("'p' (missegregation probability) must be a single numeric value between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(times) || anyNA(times) || is.unsorted(times, strictly = FALSE)) {
    stop("'times' must be a sorted numeric vector without NAs.", call. = FALSE)
  }
  if (length(times) == 0) stop("'times' must not be empty.", call. = FALSE)
  
  if (!is.numeric(x0) || is.null(names(x0))) {
    stop("'x0' (initial frequencies) must be a named numeric vector.", call. = FALSE)
  }
  if (length(x0) != nrow(lscape)) {
    stop("Length of 'x0' must match the number of rows (karyotypes) in 'lscape'.", call. = FALSE)
  }
  if(!setequal(names(x0), lscape$k)) {
    stop("Names of 'x0' must exactly match all unique karyotypes in 'lscape$k'. Order does not matter initially, but all names must be present in both.", call. = FALSE)
  }
  # Reorder x0 to match lscape$k order for consistency downstream
  x0 <- x0[lscape$k] 
  if(anyNA(x0)){ # Check after reordering if any names in lscape$k were not in original x0 names
    stop("Mismatch between names(x0) and lscape$k led to NAs after reordering x0. Ensure all lscape$k are in names(x0).", call. = FALSE)
  }
  
  if (abs(sum(x0) - 1.0) > 1e-6) stop("'x0' frequencies must sum to 1 (within tolerance).", call. = FALSE)
  if (any(x0 < -1e-9)) stop("'x0' must contain non-negative frequencies (allowing for small numerical errors near zero).", call. = FALSE) # Allow tiny negatives
  x0[x0<0] <- 0 # Correct tiny negatives before use
  
  if (!(prediction_type %in% c("ODE", "ABM"))) stop("'prediction_type' must be either 'ODE' or 'ABM'.", call. = FALSE)
  
  # Further validation (moved from internal functions to be user-facing checks)
  # This implicitly checks karyotype string format and consistency.
  # If parse_karyotypes fails, it will stop here.
  parsed_k_for_validation <- parse_karyotypes(lscape$k)
  
  result_df <- NULL 
  if (prediction_type == "ODE") {
    result_df <- run_ode_simulation(lscape = lscape, p = p, times = times, x0 = x0, ode_method = ode_method,Nmax=Nmax)
  } else if (prediction_type == "ABM") {
    # ABM specific parameter validation
    if (!is.numeric(abm_pop_size) || length(abm_pop_size) != 1 || abm_pop_size <= 0 || floor(abm_pop_size) != abm_pop_size) {
      stop("'abm_pop_size' must be a single positive integer.", call. = FALSE)
    }
    if (!is.numeric(abm_delta_t) || length(abm_delta_t) != 1 || abm_delta_t <= 0) {
      stop("'abm_delta_t' must be a single positive number.", call. = FALSE)
    }
    if (!is.numeric(abm_max_pop) || length(abm_max_pop) != 1) {
      stop("'abm_max_pop' must be a single numeric value.", call. = FALSE)
    }
    if (!is.numeric(abm_culling_survival) || length(abm_culling_survival) != 1 || abm_culling_survival < 0 || abm_culling_survival > 1) {
      stop("'abm_culling_survival' must be a single number between 0 and 1.", call. = FALSE)
    }
    if (!is.numeric(abm_record_interval) || length(abm_record_interval) != 1 || abm_record_interval <= 0 || floor(abm_record_interval) != abm_record_interval) {
      stop("'abm_record_interval' must be a single positive integer.", call. = FALSE)
    }
    if (!is.numeric(abm_seed) || length(abm_seed) != 1 || floor(abm_seed) != abm_seed ) { # Ensure abm_seed is integer-like
      stop("'abm_seed' must be a single integer value.", call. = FALSE)
    }
    
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
    warning("Invalid 'lscape' input: must be a data.frame with 'k' and 'mean' columns.", call. = FALSE); return(NULL)
  }
  if (!is.character(lscape$k) || length(lscape$k) == 0) {
    warning("'lscape$k' must be a non-empty character vector.", call. = FALSE); return(NULL)
  }
  if (anyDuplicated(lscape$k)) { # Check for duplicates that would cause issues
    stop("'lscape$k' contains duplicate karyotype strings. Steady state calculation requires unique karyotypes in lscape.", call. = FALSE)
  }
  if (!is.numeric(lscape$mean) || length(lscape$mean) != length(lscape$k)) {
    warning("'lscape$mean' must be numeric and have the same length as 'lscape$k'.", call. = FALSE); return(NULL)
  }
  if (!is.numeric(p) || length(p) != 1 || p < 0 || p > 1) {
    warning("'p' (missegregation probability) must be a single numeric value between 0 and 1.", call. = FALSE); return(NULL)
  }
  
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
  try_LR_on_LM_fail <- TRUE # Control if we try LR after LM

  if (nrow(M_eigs) < 3) {
    dense_eigs <- eigen(as.matrix(M_eigs))
    dominant_idx <- which.max(Re(dense_eigs$values))
    eig_result <- list(
      values = dense_eigs$values[dominant_idx],
      vectors = dense_eigs$vectors[, dominant_idx, drop = FALSE]
    )
  } else {
    eig_result <- tryCatch(RSpectra::eigs(M_eigs, k = 1, which = "LM"), error = function(e_lm) {
      warning("RSpectra::eigs with 'LM' failed: ", e_lm$message, ". ", 
              if(try_LR_on_LM_fail) "Trying 'LR'." else "Not trying 'LR'.", call. = FALSE)
      if(try_LR_on_LM_fail){
        return(tryCatch(RSpectra::eigs(M_eigs, k = 1, which = "LR"), error = function(e_lr) {
          stop(sprintf("RSpectra::eigs also failed with 'LR': %s", e_lr$message), call. = FALSE)
        }))
      } else {
        stop(sprintf("RSpectra::eigs with 'LM' failed: %s. Aborting.", e_lm$message), call. = FALSE)
      }
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
  if(!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
    stop("'lambda' must be a single positive numeric value.", call. = FALSE)
  K <- ncol(centroids)
  
  parse_len <- function(s) length(strsplit(s, "\\.")[[1]])
  bad <- names(x0)[vapply(names(x0), parse_len, 0L) != K]
  if(length(bad))
    stop("Karyotype names ", paste(bad, collapse = ", "),
         " do not have ", K, " chromosome counts.", call. = FALSE)
  if(abs(sum(x0) - 1) > 1e-6)
    stop("'x0' must sum to 1.", call. = FALSE)
  
  ## -- initial population ----------------------------------------------------
  init_counts <- round(x0 * abm_pop_size)
  init_counts[init_counts < 0] <- 0
  init_list   <- as.list(init_counts)[init_counts > 0]  
  if(!length(init_list)) stop("Initial population is zero.", call. = FALSE)
  
  ## -- run C++ ---------------------------------------------------------------
  steps <- ceiling(max(times) / abm_delta_t)
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
    wide <- tidyr::pivot_wider(long, names_from = .data$Karyotype,
                               values_from = .data$Frequency, values_fill = 0)
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
