#' ALFAK: Adaptive Landscape Fitness Inference from Karyotype Dynamics
#'
#' Performs fitness landscape inference using Allele Frequency Karyotype dynamics.
#' This function estimates fitness values for observed karyotypes and their
#' neighbors, performs Kriging to create a fitness landscape, and evaluates
#' the landscape using cross-validation.
#'
#' @param yi A list containing the input data. Expected elements:
#'   \itemize{
#'     \item `x`: A matrix of karyotype counts (rows are karyotypes as strings
#'       like "2.2.1", columns are timepoints). Rownames should be present.
#'       Colnames should represent time values that can be multiplied by `dt`.
#'     \item `dt`: A numeric value representing the scaling factor for timepoints
#'       (e.g., generation time if colnames of `x` are generations).
#'   }
#' @param outdir A string specifying the directory path where result files
#'   (RDS format) will be saved. The directory will be created if it doesn't exist.
#' @param passage_times An optional numeric vector giving the final internal time
#'   axis used by the fitter. If `NULL`, the time axis is calculated as
#'   `as.numeric(colnames(yi$x)) * yi$dt`. When supplied, `passage_times` is used
#'   as-is and must be numeric, finite, strictly increasing, and have length
#'   `ncol(yi$x)`.
#' @param minobs An integer, the minimum total number of observations (reads/counts)
#'   for a karyotype across all timepoints to be considered "frequent" and included
#'   in the analysis. A karyotype is considered frequent when
#'   `rowSums(yi$x) >= minobs`. Default is 20.
#' @param nboot An integer, the number of bootstrap iterations for fitness
#'   estimation and for the Kriging process in `fitKrig`. Default is 45.
#' @param n0 A numeric value, the initial effective population size at the
#'   start of a passage or growth phase, used for g0 calculation. Default is 1e5.
#' @param nb A numeric value, the bottleneck effective population size (population
#'   size after transfer), used for g0 calculation. Default is 1e7.
#' @param pm A numeric value, the per-locus mutation/error rate used in `pij`
#'   calculations. Default is 0.00005.
#' @param allow_noninteger_counts Logical; if `FALSE` (default), non-integer
#'   counts in `yi$x` are rejected. If `TRUE`, they are rounded once at entry
#'   with a warning.
#' @param correct_efflux Logical; if `TRUE`, apply the efflux correction after a
#'   one-time viability pre-check on the frequent karyotypes. The viability term
#'   currently depends only on total copy number through
#'   `2 * (1 - pm)^total_copy_number - 1`, so this is a ploidy-level approximation
#'   rather than a chromosome-specific model.
#' @param landscape_data_output Logical; if `TRUE`, also save the optional
#'   `landscape_data.Rds` file containing the stable Kriging mean and median
#'   model objects. Default is `FALSE`, so only the documented core outputs are
#'   written.
#' @param nn_prior Character; nearest-neighbour prior mode for latent children.
#'   `"empirical"` uses the current empirical child-minus-parent prior estimated
#'   from observed neighbours, while `"none"` disables that prior contribution.
#'   Default is `"empirical"`.
#' @param nn_prior_sd Optional numeric scalar. If supplied, this overrides the
#'   empirically estimated prior standard deviation for latent-neighbour fitting.
#' @param nn_prior_sd_floor Numeric scalar giving the minimum standard deviation
#'   used when the empirical prior variance is zero or too small. Default is
#'   `1e-3`.
#' @param krig_bootstrap_mode Character; `"joint"` (default) samples one full
#'   bootstrap row at a time when building Kriging posterior samples, preserving
#'   within-bootstrap dependence across karyotypes. `"marginal"` retains the old
#'   column-wise synthetic resampling behavior.
#'
#' @return Returns the cross-validation R-squared value (`Rxv`) invisibly.
#'   The function primarily saves its results to RDS files in the `outdir`:
#'   \itemize{
#'     \item `bootstrap_res.Rds`: Results from `solve_fitness_bootstrap`.
#'     \item `landscape.Rds`: Summary statistics (mean, median, sd) of the
#'       Kriging-inferred fitness landscape from `fitKrig`.
#'     \item `landscape_posterior_samples.Rds`: The full matrix of posterior
#'       samples from the Kriging bootstraps in `fitKrig`.
#'     \item `xval.Rds`: The cross-validation R-squared value (`Rxv`).
#'     \item `landscape_data.Rds`: Optional stable Kriging mean/median model
#'       objects, written only when `landscape_data_output = TRUE`.
#'   }
#'
#' @export
#' @importFrom quadprog solve.QP
#' @importFrom fields Krig
#' @importFrom stats rmultinom optim uniroot dbinom dnorm optimise dist predict median sd setNames complete.cases rnorm
#'
#' @examples
#' \dontrun{
#' # Create dummy data for yi
#' karyotypes_str <- c("2.2.2", "2.2.1", "2.1.2", "1.2.2", "2.2.3")
#' timepoints_num <- 5
#' counts_data <- matrix(
#'   abs(rpois(length(karyotypes_str) * timepoints_num, lambda = 20)), # Use pois for counts
#'   nrow = length(karyotypes_str),
#'   ncol = timepoints_num
#' )
#' rownames(counts_data) <- karyotypes_str
#' colnames(counts_data) <- 1:timepoints_num # Generations
#'
#' dummy_yi_data <- list(
#'   x = counts_data,
#'   dt = 1 # dt = 1 if colnames are generations
#' )
#'
#' temp_output_dir <- tempfile("alfak_example_")
#'
#' # Run alfak
#' result_r_squared <- alfak(
#'   yi = dummy_yi_data,
#'   outdir = temp_output_dir,
#'   passage_times = NULL,
#'   minobs = 5,      # Lowered for dummy data
#'   nboot = 10,      # Lowered for quick example
#'   n0 = 1e4,
#'   nb = 1e6,
#'   pm = 0.0001
#' )
#' print(paste("Cross-validation R-squared:", result_r_squared))
#'
#' # Check for created files
#' list.files(temp_output_dir)
#'
#' # Clean up
#' unlink(temp_output_dir, recursive = TRUE)
#' }
alfak <- function(yi, outdir, passage_times = NULL, minobs = 20,
                  nboot = 45,
                  n0 = 1e5,
                  nb = 1e7,
                  pm = 0.00005,
                  allow_noninteger_counts = FALSE,
                  correct_efflux=FALSE,
                  landscape_data_output = FALSE,
                  nn_prior = c("empirical", "none"),
                  nn_prior_sd = NULL,
                  nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR,
                  krig_bootstrap_mode = c("joint", "marginal")) {

  # Note: library calls removed, dependencies handled by @importFrom or DESCRIPTION

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  validate_positive_integer(nboot, "nboot")
  validate_positive_integer(n0, "n0")
  validate_positive_integer(nb, "nb")
  validate_probability(pm, "pm")
  validate_scalar_logical(allow_noninteger_counts, "allow_noninteger_counts")
  validate_scalar_logical(correct_efflux, "correct_efflux")
  validate_scalar_logical(landscape_data_output, "landscape_data_output")
  nn_prior <- validate_nn_prior_mode(nn_prior)
  krig_bootstrap_mode <- validate_krig_bootstrap_mode(krig_bootstrap_mode)
  validate_nn_prior_controls(nn_prior_sd = nn_prior_sd, nn_prior_sd_floor = nn_prior_sd_floor)
  yi$x <- coerce_count_matrix(yi$x, allow_noninteger_counts = allow_noninteger_counts)
  validate_positive_depth(yi$x)

  get_frequent_karyotypes(yi$x, minobs)
  resolve_time_axis(yi, passage_times)

  # Parallelism and cl related code removed

  fq_boot <- solve_fitness_bootstrap(yi, minobs = minobs, nboot = nboot,
                                     n0 = n0, nb = nb, pm = pm,
                                     allow_noninteger_counts = allow_noninteger_counts,
                                     passage_times = passage_times,correct_efflux=correct_efflux,
                                     nn_prior = nn_prior,
                                     nn_prior_sd = nn_prior_sd,
                                     nn_prior_sd_floor = nn_prior_sd_floor)
  saveRDS(fq_boot, file = file.path(outdir, "bootstrap_res.Rds"))

  landscape_data <- fitKrig(fq_boot, nboot, krig_bootstrap_mode = krig_bootstrap_mode)
  saveRDS(landscape_data$summary_stats, file = file.path(outdir, "landscape.Rds"))
  saveRDS(landscape_data$posterior_samples, file = file.path(outdir, "landscape_posterior_samples.Rds"))

  if (isTRUE(landscape_data_output)) {
    Krig_stable <- list(landscape_data$krig_stable_mean, landscape_data$krig_stable_median)
    names(Krig_stable) <- c("mean", "median")
    saveRDS(Krig_stable, file = file.path(outdir, "landscape_data.Rds"))
  }
  xval_res <- xval(fq_boot)
  Rxv <- extract_xval_r2r(xval_res)
  saveRDS(Rxv, file = file.path(outdir, "xval.Rds"))

  ##END HERE.

  invisible(Rxv) # Return Rxv invisibly as side-effect saving is primary
}

##########################################
# Helper functions (internal)
##########################################

#' Calculate the single-chromosome transition probability p_ij
#'
#' Validates `i`, `j`, and `beta` before calling the C++ transition kernel.
#'
#' @param i Positive integer parent copy number.
#' @param j Positive integer daughter copy number.
#' @param beta Finite mis-segregation probability in `[0, 1]`.
#' @return A transition probability in `[0, 1]`.
#' @export
pij <- function(i, j, beta) {
  validate_positive_integer(i, "i")
  validate_positive_integer(j, "j")
  validate_probability(beta, "beta", upper_inclusive = TRUE)
  pij_cpp(i, j, beta)
}

#' Convert string like "1.2.3" to numeric vector
#' @keywords internal
#' @noRd
parse_karyotype_ids <- function(ids) {
  if (!is.character(ids) || length(ids) == 0 || any(!nzchar(ids))) {
    stop("Karyotype IDs must be non-empty character strings.", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop("Karyotype IDs must be unique.", call. = FALSE)
  }
  pieces <- strsplit(ids, "\\.", perl = TRUE)
  lens <- lengths(pieces)
  if (any(lens == 0) || length(unique(lens)) != 1) {
    stop(
      "All karyotype IDs must have the same number of dot-separated components.",
      call. = FALSE
    )
  }
  bad <- vapply(pieces, function(x) any(!grepl("^[1-9][0-9]*$", x)), logical(1))
  if (any(bad)) {
    stop(
      sprintf("Invalid karyotype ID(s): %s", paste(ids[bad], collapse = ", ")),
      call. = FALSE
    )
  }
  mat <- matrix(as.integer(unlist(pieces, use.names = FALSE)), nrow = length(ids), byrow = TRUE)
  storage.mode(mat) <- "integer"
  rownames(mat) <- ids
  mat
}

#' Convert string like "1.2.3" to numeric vector
#' @keywords internal
#' @noRd
s2v <- function(s) {
  parsed <- parse_karyotype_ids(as.character(s))
  if (nrow(parsed) == 1) {
    return(as.numeric(parsed[1, ]))
  }
  parsed
}

#' Calculate R-squared
#' @keywords internal
#' @noRd
R2R <- function(obs, pred) {
  valid <- is.finite(obs) & is.finite(pred)
  obs <- obs[valid]
  pred <- pred[valid]
  if (length(obs) < 2) {
    return(NA_real_)
  }
  denom <- sum((obs - mean(obs))^2)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  1 - sum((pred - obs)^2) / denom
}

#' Extract the scalar R2R value from xval() output
#' @keywords internal
#' @noRd
extract_xval_r2r <- function(xval_result) {
  if (is.list(xval_result) && !is.null(xval_result$R2R)) {
    r2r_val <- xval_result$R2R
  } else {
    r2r_val <- xval_result
  }
  if (!is.numeric(r2r_val) || length(r2r_val) != 1 || is.nan(r2r_val) || is.infinite(r2r_val)) {
    stop("`xval()` must return a single numeric R2R value, optionally `NA_real_`, or a list containing scalar `R2R`.")
  }
  as.numeric(r2r_val)
}

ALFAK_FEXP_DELTA_TOL <- 1e-8
ALFAK_EFFLUX_VIABILITY_TOL <- 1e-6
ALFAK_NN_PRIOR_SD_FLOOR <- 1e-3
ALFAK_COUNT_INTEGER_TOL <- sqrt(.Machine$double.eps)
ALFAK_KRIG_NSTEP_CV <- 200L
ALFAK_MAX_EXACT_INTEGER <- 2^53 - 1

#' Validate scalar positive integer input
#' @keywords internal
#' @noRd
validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 1 || x != floor(x)) {
    stop(sprintf("`%s` must be a single positive integer.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate scalar positive finite numeric input
#' @keywords internal
#' @noRd
validate_positive_finite <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive finite numeric value.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate scalar logical input
#' @keywords internal
#' @noRd
validate_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("`%s` must be a single TRUE/FALSE value.", name), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate scalar probability input
#' @keywords internal
#' @noRd
validate_probability <- function(x, name, upper_inclusive = FALSE) {
  upper_ok <- if (upper_inclusive) isTRUE(x <= 1) else isTRUE(x < 1)
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < 0 || !upper_ok) {
    bound <- if (upper_inclusive) "[0, 1]" else "[0, 1)"
    stop(sprintf("`%s` must be a single finite numeric value in %s.", name, bound), call. = FALSE)
  }
  invisible(NULL)
}

#' Validate Kriging bootstrap sampling mode
#' @keywords internal
#' @noRd
validate_krig_bootstrap_mode <- function(mode) {
  match.arg(mode, c("joint", "marginal"))
}

#' Validate that every timepoint has positive sequencing depth
#' @keywords internal
#' @noRd
validate_positive_depth <- function(x) {
  zero_depth <- colSums(x) == 0
  if (any(zero_depth)) {
    depth_names <- colnames(x)
    if (is.null(depth_names)) {
      depth_names <- as.character(seq_len(ncol(x)))
    }
    stop(
      sprintf(
        paste0(
          "Each timepoint must have positive total counts; zero-depth column(s): %s. ",
          "Remove or explicitly impute missing timepoints."
        ),
        paste(depth_names[zero_depth], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Numerically stable softmax
#' @keywords internal
#' @noRd
softmax <- function(z) {
  if (!all(is.finite(z))) {
    stop("Cannot softmax non-finite logits.", call. = FALSE)
  }
  z <- z - max(z)
  ex <- exp(z)
  ex / sum(ex)
}

#' Safe mean that returns NA for all-missing vectors
#' @keywords internal
#' @noRd
mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  mean(x)
}

#' Safe median that returns NA for all-missing vectors
#' @keywords internal
#' @noRd
median_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NA_real_)
  }
  stats::median(x)
}

#' Safe standard deviation that returns NA when fewer than 2 finite values exist
#' @keywords internal
#' @noRd
sd_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(NA_real_)
  }
  stats::sd(x)
}

#' Validate nearest-neighbour prior mode
#' @keywords internal
#' @noRd
validate_nn_prior_mode <- function(nn_prior) {
  match.arg(nn_prior, c("empirical", "none"))
}

#' Validate nearest-neighbour prior controls
#' @keywords internal
#' @noRd
validate_nn_prior_controls <- function(nn_prior_sd = NULL, nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR) {
  if (!is.null(nn_prior_sd)) {
    validate_positive_finite(nn_prior_sd, "nn_prior_sd")
  }
  validate_positive_finite(nn_prior_sd_floor, "nn_prior_sd_floor")
  invisible(NULL)
}

#' Coerce supported count containers to a base numeric matrix
#' @keywords internal
#' @noRd
coerce_count_matrix <- function(x, allow_noninteger_counts = FALSE) {
  if (is.null(x)) {
    stop("`yi$x`/`data$x` must be a two-dimensional numeric count object.")
  }
  x_dim <- dim(x)
  if (length(x_dim) != 2) {
    stop("`yi$x`/`data$x` must be a two-dimensional numeric count object.")
  }

  x_mat <- try(as.matrix(x), silent = TRUE)
  if (inherits(x_mat, "try-error") || !is.matrix(x_mat)) {
    stop("`yi$x`/`data$x` must be coercible to a numeric matrix of karyotype counts.")
  }
  if (!is.numeric(x_mat)) {
    stop("`yi$x`/`data$x` must contain numeric karyotype counts.")
  }
  if (any(!is.finite(x_mat))) {
    stop("`yi$x`/`data$x` must contain only finite count values.")
  }
  if (any(x_mat < 0)) {
    stop("`yi$x`/`data$x` must contain non-negative count values.")
  }
  non_integer <- abs(x_mat - round(x_mat)) > ALFAK_COUNT_INTEGER_TOL
  if (any(non_integer)) {
    if (!isTRUE(allow_noninteger_counts)) {
      stop("Non-integer values detected in `yi$x`/`data$x`; set `allow_noninteger_counts = TRUE` to round once at entry.", call. = FALSE)
    }
    warning("Non-integer values detected in `yi$x`/`data$x`; rounding to the nearest integer once at entry.")
    x_mat <- round(x_mat)
  }
  x_mat
}

#' Convert a probability vector to K-1 free softmax logits
#' @keywords internal
#' @noRd
free_softmax_logits <- function(prob) {
  if (!is.numeric(prob) || length(prob) == 0 || any(!is.finite(prob)) || any(prob < 0)) {
    stop("`prob` must be a finite non-negative probability vector.", call. = FALSE)
  }
  if (length(prob) == 1) {
    return(numeric(0))
  }
  prob <- prob / sum(prob)
  ref <- prob[length(prob)]
  if (!is.finite(ref) || ref <= 0) {
    stop("The reference probability for free-softmax logits must be positive.", call. = FALSE)
  }
  log(prob[-length(prob)]) - log(ref)
}

#' Expand K-1 free logits into a K-vector of probabilities
#' @keywords internal
#' @noRd
softmax_from_free_logits <- function(logits_free) {
  softmax(c(logits_free, 0))
}

#' Ensure birth-time estimates are finite before neighbour estimation
#' @keywords internal
#' @noRd
sanitize_birth_times <- function(birth_times_est, peak_times, timepoints) {
  mean_risetime <- mean(peak_times - birth_times_est, na.rm = TRUE)
  fallback_used <- FALSE
  if (!is.finite(mean_risetime)) {
    mean_risetime <- 0
    fallback_used <- TRUE
  }

  missing_birth <- !is.finite(birth_times_est)
  if (any(missing_birth)) {
    birth_times_est[missing_birth] <- peak_times[missing_birth] - mean_risetime
    fallback_used <- TRUE
  }

  unresolved <- !is.finite(birth_times_est)
  if (any(unresolved)) {
    safe_fallback <- peak_times
    safe_fallback[!is.finite(safe_fallback)] <- min(timepoints)
    birth_times_est[unresolved] <- safe_fallback[unresolved]
    fallback_used <- TRUE
  }

  if (fallback_used) {
    warning("Using finite fallback birth times for nearest-neighbour estimation because root-finding did not return enough finite birth times.")
  }

  birth_times_est
}

#' Format a short karyotype preview for diagnostics
#' @keywords internal
#' @noRd
format_karyotype_preview <- function(karyotypes, max_show = 5) {
  if (length(karyotypes) == 0) {
    return("<none>")
  }
  shown <- utils::head(karyotypes, max_show)
  preview <- paste(shown, collapse = ", ")
  if (length(karyotypes) > max_show) {
    preview <- paste0(preview, ", ...")
  }
  preview
}

#' Select frequent karyotypes using the documented minobs rule
#' @keywords internal
#' @noRd
get_frequent_karyotypes <- function(x, minobs) {
  x <- coerce_count_matrix(x)
  if (is.null(rownames(x)) || any(!nzchar(rownames(x)))) {
    stop("`yi$x`/`data$x` must have non-empty rownames for karyotype IDs.")
  }
  parse_karyotype_ids(rownames(x))
  if (!is.numeric(minobs) || length(minobs) != 1 || !is.finite(minobs) || minobs < 0) {
    stop("`minobs` must be a single non-negative finite numeric value.")
  }
  minobs <- as.numeric(minobs)
  fq <- rownames(x)[rowSums(x) >= minobs]
  if (length(fq) == 0) {
    stop(sprintf(
      "no frequent karyotypes detected for minobs = %s; frequent karyotypes require rowSums(x) >= minobs",
      format(minobs, trim = TRUE)
    ))
  }
  fq
}

#' Resolve and validate the internal time axis
#' @keywords internal
#' @noRd
resolve_time_axis <- function(data, passage_times = NULL) {
  if (!is.list(data) || is.null(data$x)) {
    stop("`yi`/`data` must be a list containing a matrix-like `x` of karyotype counts.")
  }
  data$x <- coerce_count_matrix(data$x)
  if (ncol(data$x) < 2) {
    stop("At least two timepoints are required in `yi$x`/`data$x`.")
  }

  if (is.null(passage_times)) {
    if (is.null(data$dt) || !is.numeric(data$dt) || length(data$dt) != 1 || !is.finite(data$dt)) {
      stop("`yi$dt` must be a single finite numeric value when `passage_times` is NULL.")
    }
    raw_axis <- suppressWarnings(as.numeric(colnames(data$x)))
    if (length(raw_axis) != ncol(data$x) || any(!is.finite(raw_axis))) {
      stop("When `passage_times` is NULL, `colnames(yi$x)` must be numeric and finite so `colnames(yi$x) * yi$dt` defines the time axis.")
    }
    time_axis <- raw_axis * data$dt
  } else {
    if (!is.numeric(passage_times)) {
      stop("`passage_times` must be a numeric vector when supplied.")
    }
    if (length(passage_times) != ncol(data$x)) {
      stop(sprintf(
        "`passage_times` must have length %d to match ncol(yi$x); got %d.",
        ncol(data$x), length(passage_times)
      ))
    }
    if (any(!is.finite(passage_times))) {
      stop("`passage_times` must contain only finite values.")
    }
    time_axis <- as.numeric(passage_times)
  }

  delta_t <- diff(time_axis)
  if (any(delta_t <= 0)) {
    stop("The internal time axis must be strictly increasing.")
  }

  time_axis
}

#' Normalize counts to frequencies by column
#' @keywords internal
#' @noRd
normalize_columns <- function(count_matrix) {
  validate_positive_depth(count_matrix)
  totals <- colSums(count_matrix)
  sweep(count_matrix, 2, totals, "/")
}

#' Validate viability for efflux correction once before bootstrapping
#' @keywords internal
#' @noRd
prepare_efflux_viability <- function(fq_vec, pm, correct_efflux,
                                     viability_tol = ALFAK_EFFLUX_VIABILITY_TOL) {
  if (!is.numeric(pm) || length(pm) != 1 || !is.finite(pm) || pm < 0 || pm >= 1) {
    stop("`pm` must be a single finite numeric value in [0, 1).")
  }
  viability <- setNames(rep(1, nrow(fq_vec)), rownames(fq_vec))
  if (!isTRUE(correct_efflux)) {
    return(viability)
  }

  viability <- setNames(2 * (1 - pm)^rowSums(fq_vec) - 1, rownames(fq_vec))

  non_positive <- viability <= 0
  if (any(non_positive)) {
    affected <- names(viability)[non_positive]
    affected_vals <- viability[non_positive]
    stop(sprintf(
      paste0(
        "correct_efflux viability pre-check failed before bootstrap: pm=%s yields ",
        "%d frequent karyotype(s) with viability <= 0; affected=%s; viability range among affected=[%.6g, %.6g]."
      ),
      format(pm, trim = TRUE),
      length(affected),
      format_karyotype_preview(affected),
      min(affected_vals),
      max(affected_vals)
    ))
  }

  near_zero <- viability < viability_tol
  if (any(near_zero)) {
    affected <- names(viability)[near_zero]
    affected_vals <- viability[near_zero]
    warning(sprintf(
      paste0(
        "correct_efflux viability pre-check before bootstrap: pm=%s yields ",
        "%d frequent karyotype(s) with 0 < viability < %.1e; affected=%s; viability range among affected=[%.6g, %.6g]."
      ),
      format(pm, trim = TRUE),
      length(affected),
      viability_tol,
      format_karyotype_preview(affected),
      min(affected_vals),
      max(affected_vals)
    ))
  }

  viability
}

#' Run optim and surface non-convergence explicitly
#' @keywords internal
#' @noRd
run_optim_checked <- function(par, fn, ..., method = "BFGS", control = NULL, context) {
  opt <- try(stats::optim(par = par, fn = fn, ..., method = method, control = control), silent = TRUE)
  if (inherits(opt, "try-error")) {
    stop(sprintf("%s failed: %s", context, as.character(opt)))
  }
  if (!all(is.finite(opt$par)) || !is.finite(opt$value)) {
    stop(sprintf("%s returned non-finite parameters or objective values.", context))
  }
  if (!is.null(opt$convergence) && opt$convergence != 0) {
    warning(sprintf(
      "%s returned convergence code %d%s",
      context,
      opt$convergence,
      if (!is.null(opt$message) && nzchar(opt$message)) paste0(": ", opt$message) else "."
    ))
  }
  opt
}

#' Run optimise and warn on invalid scalar optima
#' @keywords internal
#' @noRd
run_optimise_checked <- function(f, interval, ..., context) {
  opt <- try(stats::optimise(f, interval = interval, ...), silent = TRUE)
  if (inherits(opt, "try-error")) {
    warning(sprintf("%s failed: %s", context, as.character(opt)))
    return(NULL)
  }
  if (!is.finite(opt$minimum) || !is.finite(opt$objective)) {
    warning(sprintf("%s returned a non-finite optimum.", context))
    return(NULL)
  }
  opt
}

#' Run solve.QP and fail loudly on invalid optimizer state
#' @keywords internal
#' @noRd
run_solve_qp_checked <- function(Dmat, dvec, Amat, bvec, meq, context) {
  qp_sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq), silent = TRUE)
  if (inherits(qp_sol, "try-error")) {
    stop(sprintf("%s failed: %s", context, as.character(qp_sol)))
  }
  if (is.null(qp_sol$solution) || any(!is.finite(qp_sol$solution))) {
    stop(sprintf("%s returned a non-finite solution.", context))
  }
  qp_sol
}

#' Weighted parent fitness for nearest-neighbour prior
#' @keywords internal
#' @noRd
weighted_parent_fitness <- function(nni_item, fpar) {
  valid_parents <- nni_item$nj[nni_item$nj %in% names(fpar)]
  if (length(valid_parents) == 0) {
    return(NA_real_)
  }
  parent_fitness <- fpar[valid_parents]
  parent_weights <- nni_item$pij[match(valid_parents, nni_item$nj)]
  if (length(parent_weights) == length(parent_fitness) &&
      all(is.finite(parent_weights)) &&
      all(parent_weights >= 0) &&
      sum(parent_weights) > 0) {
    return(stats::weighted.mean(parent_fitness, w = parent_weights))
  }
  mean(parent_fitness, na.rm = TRUE)
}

#' Numerically stable exposure term for neighbour estimation
#' @keywords internal
#' @noRd
fExp_stable <- function(fc_arg, fp_arg, pij_val, tt_arg, tol = ALFAK_FEXP_DELTA_TOL) {
  delta <- fc_arg - fp_arg
  if (abs(delta) < tol) {
    # This is the analytic limit of fp * (exp(tt * delta) - 1) / delta as delta -> 0.
    return(pij_val * fp_arg * tt_arg)
  }
  pij_val * fp_arg * expm1(tt_arg * delta) / delta
}

#' Generate all single-step neighbours for karyotype IDs
#' @keywords internal
#' @noRd
gen_all_neighbours <- function(ids, as.strings = TRUE, remove_nullisomes = TRUE) {
  if (as.strings) {
    parsed_ids <- parse_karyotype_ids(as.character(ids))
    ids <- lapply(seq_len(nrow(parsed_ids)), function(i) as.numeric(parsed_ids[i, ]))
  }
  nkern <- do.call(rbind, lapply(1:length(ids[[1]]), function(i) {
    x0 <- rep(0, length(ids[[1]]))
    x1 <- x0
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0, x1)
  }))
  n <- do.call(rbind, lapply(ids, function(ii) t(apply(nkern, 1, function(i) i + ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind, ids), n) # Intentionally add originals back
  n <- unique(n) # Keep unique set
  n <- n[-(1:nids), , drop=FALSE] # Remove the nids original ids that were just added
  # drop=FALSE ensures matrix structure even if 1 row left
  if (remove_nullisomes && nrow(n) > 0) # Check nrow > 0 before apply
    n <- n[apply(n, 1, function(ni) sum(ni < 1) == 0), , drop = FALSE]
  n
}

fields_krig_find_lambda <- local({
  krig_find_lambda <- get("KrigFindLambda", envir = asNamespace("fields"))
  function(...) krig_find_lambda(...)
})

fields_krig_coef <- local({
  krig_coef <- get("Krig.coef", envir = asNamespace("fields"))
  function(...) krig_coef(...)
})

fields_krig_parameters <- local({
  krig_parameters <- get("Krig.parameters", envir = asNamespace("fields"))
  function(...) krig_parameters(...)
})

fields_krig_ynew <- local({
  krig_ynew <- get("Krig.ynew", envir = asNamespace("fields"))
  function(...) krig_ynew(...)
})

krig_covariance_args <- function() {
  list(
    Covariance = "Matern",
    smoothness = 1.5
  )
}

predict_cached_krig <- function(object, x, dist_mat, Z = NULL, drop.Z = FALSE, just.fixed = FALSE) {
  if (is.null(x)) {
    x <- object$x
  } else {
    x <- as.matrix(x)
  }
  if (is.null(Z)) {
    Z <- object$Z
  } else {
    Z <- as.matrix(Z)
  }

  x_scaled <- scale(x, object$transform$x.center, object$transform$x.scale)
  null_fun <- get(object$null.function.name, envir = asNamespace("fields"))
  Tmatrix <- do.call(null_fun, c(object$null.args, list(x = x_scaled, Z = Z, drop.Z = drop.Z)))

  if (drop.Z) {
    pred <- Tmatrix %*% object$d[object$ind.drift]
  } else {
    pred <- Tmatrix %*% object$d
  }

  if (!just.fixed) {
    cov_args <- object$args
    cov_args$distMat <- dist_mat
    cov_fun <- get(object$cov.function.name, envir = asNamespace("fields"))
    pred <- pred + do.call(cov_fun, c(cov_args, list(x1 = x_scaled, x2 = object$knots, C = object$c)))
  }

  as.numeric(pred)
}

build_cached_krig_fit <- function(ktrain, y, kpred = NULL, give_warnings = TRUE) {
  train_dist <- fields::rdist(ktrain, ktrain)
  pred_dist <- NULL
  if (!is.null(kpred)) {
    pred_dist <- fields::rdist(kpred, ktrain)
  }
  fit <- fields::Krig(
    ktrain,
    y,
    cov.function = "stationary.cov",
    cov.args = krig_covariance_args(),
    nstep.cv = ALFAK_KRIG_NSTEP_CV,
    give.warnings = give_warnings
  )
  list(fit = fit, train_dist = train_dist, pred_dist = pred_dist)
}

refit_cached_krig <- function(cache, y, x_pred = NULL, pred_dist = cache$pred_dist, give_warnings = TRUE) {
  fit <- cache$fit
  gcv_out <- fields_krig_find_lambda(
    fit,
    nstep.cv = ALFAK_KRIG_NSTEP_CV,
    cost = fit$cost,
    offset = fit$offset,
    y = y,
    give.warnings = give_warnings
  )

  if (fit$method != "user") {
    fit$lambda <- gcv_out$lambda.est[fit$method, 1]
    fit$eff.df <- gcv_out$lambda.est[fit$method, 2]
  }

  y_info <- fields_krig_ynew(fit, y = y)
  coef_out <- fields_krig_coef(fit, lambda = fit$lambda, y = y)

  fit$gcv.grid <- gcv_out$gcv.grid
  fit$lambda.est <- gcv_out$lambda.est
  fit$warningTable <- gcv_out$warningTable
  fit$y <- as.numeric(y)
  fit$yM <- y_info$yM
  fit$c <- coef_out$c
  fit$d <- coef_out$d
  fit$tauHat.rep <- coef_out$tauHat.rep
  fit$tauHat.pure.error <- coef_out$tauHat.pure.error
  fit$pure.ss <- coef_out$pure.ss
  fit$fitted.values <- predict_cached_krig(fit, fit$x, dist_mat = cache$train_dist, Z = fit$Z)
  fit$residuals <- fit$y - fit$fitted.values

  fit_params <- fields_krig_parameters(fit)
  fit[names(fit_params)] <- fit_params

  if (fit$method == "user" && !is.na(fit$tau2)) {
    fit$best.model <- c(fit$lambda, fit$tau2, fit$sigma)
  } else {
    fit$best.model <- c(fit$lambda, fit$tauHat.MLE^2, fit$sigmahat)
  }
  fit$rhohat <- fit$sigmahat

  preds <- NULL
  if (!is.null(x_pred)) {
    preds <- predict_cached_krig(fit, x_pred, dist_mat = pred_dist)
  }

  list(fit = fit, preds = preds)
}

#' Bootstrap counts from data matrix
#' @keywords internal
#' @noRd
bootstrap_counts <- function(data) {
  num_species <- nrow(data)
  num_timepoints <- ncol(data)
  boot_data <- matrix(NA, nrow = num_species, ncol = num_timepoints)
  rownames(boot_data) <- rownames(data)
  for (i in seq_len(num_timepoints)) {
    total_counts <- sum(data[, i])
    if (total_counts == 0) {
      boot_data[, i] <- rep(0, num_species)
    } else {
      boot_data[, i] <- stats::rmultinom(1, size = total_counts, prob = data[, i] / total_counts)
    }
  }
  boot_data
}

#' Compute dx/dt (rate of change)
#' @keywords internal
#' @noRd
compute_dx_dt <- function(x, timepoints) {
  if (length(timepoints) != ncol(x)) {
    stop("`timepoints` must have length ncol(x).")
  }
  delta_t <- diff(timepoints)
  if (any(!is.finite(delta_t)) || any(delta_t <= 0)) {
    stop("`timepoints` must be finite and strictly increasing for `compute_dx_dt()`.")
  }
  sweep(x[, -1, drop = FALSE] - x[, -ncol(x), drop = FALSE], 2, delta_t, "/")
}

#' Calculate log-sum-exp for numerical stability
#' @keywords internal
#' @noRd
logSumExp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

#' Generate nearest neighbour information
#' @keywords internal
#' @noRd
gen_nn_info <- function(fq, pm = 0.00005) {
  validate_probability(pm, "pm")
  # fq is a character vector of karyotype strings
  nn_matrix <- gen_all_neighbours(fq) # Expects list of strings or char vector
  if(nrow(nn_matrix) == 0) return(list())

  nn_str <- as.character(apply(nn_matrix, 1, paste, collapse = "."))

  n_info <- lapply(nn_str, function(ni_string) { # Renamed ni to ni_string
    # gen_all_neighbours for a single string (as input `ids`)
    nj_matrix_inner <- gen_all_neighbours(ni_string) # Pass single string directly
    nj_strings_inner <- character(0)
    if(nrow(nj_matrix_inner) > 0) {
      nj_strings_inner <- as.character(apply(nj_matrix_inner, 1, paste, collapse = "."))
    }
    nj_filtered <- nj_strings_inner[nj_strings_inner %in% fq] # Renamed nj

    nivec <- as.numeric(parse_karyotype_ids(ni_string)[1, ])
    pij_vals <- sapply(nj_filtered, function(si_string) { # Renamed si to si_string
      si_vec <- as.numeric(parse_karyotype_ids(si_string)[1, ])
      prod(sapply(1:length(si_vec), function(k) pij(si_vec[k], nivec[k], pm)))
    })
    list(ni = ni_string, nj = nj_filtered, pij = pij_vals)
  })
  # Names are set in solve_fitness_bootstrap as per original: names(nn) <- sapply(nn,...)
  n_info
}

#' Negative log-likelihood calculation
#' @keywords internal
#' @noRd
neg_log_lik <- function(param, counts, timepoints) {
  alfak_neg_log_lik_cpp(param, counts, timepoints)
}

#' Jointly optimize fitness and initial frequencies
#' @keywords internal
#' @noRd
joint_optimize <- function(counts, timepoints, f_init, x0_init) {
  K <- length(f_init)
  if (K == 1) {
    return(list(f = 0, x0 = 1))
  }
  f_free_init <- f_init[seq_len(K - 1)]
  x0_init_log <- free_softmax_logits(x0_init)
  param_init <- c(f_free_init, x0_init_log)
  obj_fun <- function(par) neg_log_lik(par, counts, timepoints)
  opt <- run_optim_checked(par = param_init, fn = obj_fun,
                           method = "BFGS",
                           control = list(maxit = 200, reltol = 1e-8),
                           context = "joint_optimize")
  f_free_opt <- opt$par[seq_len(K - 1)]
  f_opt <- c(f_free_opt, -sum(f_free_opt))
  log_x0_opt <- opt$par[K:(2 * K - 2)]
  x0_opt <- softmax_from_free_logits(log_x0_opt)
  list(f = f_opt, x0 = x0_opt)
}

#' Project frequencies forward in time using log-space calculations
#' @keywords internal
#' @noRd
project_forward_log <- function(x0, f, timepoints) {
  alfak_project_forward_log_cpp(x0, f, timepoints)
}

#' Optimize initial frequencies given observed data and fitness values
#' @keywords internal
#' @noRd
optimize_initial_frequencies <- function(x_obs, f, timepoints) {
  K <- nrow(x_obs)
  if (K == 1) {
    return(1)
  }
  loss_function <- function(log_x0_free) {
    x0 <- softmax_from_free_logits(log_x0_free)
    x_pred <- project_forward_log(x0, f, timepoints)
    sum((x_pred - x_obs)^2)
  }
  x_ini <- x_obs[, 1] + 1e-6 # Original had this epsilon
  x_ini <- x_ini / sum(x_ini)
  opt_result <- run_optim_checked(par = free_softmax_logits(x_ini), fn = loss_function,
                                  method = "BFGS",
                                  control = list(maxit = 500, reltol = 1e-8),
                                  context = "optimize_initial_frequencies")
  softmax_from_free_logits(opt_result$par)
}

#' Find "birth times" for species based on reaching a minimum frequency
#' @keywords internal
#' @noRd
find_birth_times <- function(opt_res, time_range, minF) {
  f_est <- opt_res$f
  x0_est <- opt_res$x0
  num_species <- length(f_est)
  birth_times <- rep(NA, num_species)
  for (i in seq_len(num_species)) {
    if (f_est[i] <= min(f_est)) next
    birth_fn <- function(t) {
      log_x_t <- log(x0_est[i]) + f_est[i] * t # Original log(x0_est[i])
      denom <- logSumExp(log(x0_est) + f_est * t) # Original log(x0_est)
      exp(log_x_t - denom) - minF
    }
    root <- try(stats::uniroot(birth_fn, range(time_range), tol = 1e-6)$root, silent = TRUE)
    if (!inherits(root, "try-error")) {
      birth_times[i] <- root
    }
  }
  birth_times
}

##########################################
# Core functions (now serial)
##########################################

#' Solve fitness using bootstrap approach (Internal function)
#' @keywords internal
#' @noRd
solve_fitness_bootstrap <- function(data, minobs, nboot = 1000, epsilon = 1e-6, pm = 0.00005,
                                    n0, nb, passage_times = NULL, allow_noninteger_counts = FALSE, correct_efflux=FALSE,
                                    nn_prior = c("empirical", "none"),
                                    nn_prior_sd = NULL,
                                    nn_prior_sd_floor = ALFAK_NN_PRIOR_SD_FLOOR) {
  data$x <- coerce_count_matrix(data$x, allow_noninteger_counts = allow_noninteger_counts)
  validate_positive_depth(data$x)
  validate_positive_integer(nboot, "nboot")
  validate_positive_integer(n0, "n0")
  validate_positive_integer(nb, "nb")
  validate_probability(pm, "pm")
  validate_scalar_logical(allow_noninteger_counts, "allow_noninteger_counts")
  validate_scalar_logical(correct_efflux, "correct_efflux")
  nn_prior <- validate_nn_prior_mode(nn_prior)
  validate_nn_prior_controls(nn_prior_sd = nn_prior_sd, nn_prior_sd_floor = nn_prior_sd_floor)
  fq <- get_frequent_karyotypes(data$x, minobs)
  nn_info_list <- gen_nn_info(fq, pm) # Renamed 'nn' to 'nn_info_list' for clarity
  if (length(nn_info_list) > 0 && !is.null(nn_info_list[[1]]$ni)) { # Check if naming is needed
    names(nn_info_list) <- sapply(nn_info_list, function(nni) nni$ni)
  } else if (length(nn_info_list) == 0) {
    # names(nn_info_list) would be NULL, which is fine for later checks
  } else {
    warning("nn_info_list structure unexpected for naming in solve_fitness_bootstrap")
  }

  fq_vec <- do.call(rbind, lapply(fq, s2v))
  rownames(fq_vec) <- fq
  viability <- prepare_efflux_viability(fq_vec, pm = pm, correct_efflux = correct_efflux)

  # fq_nn <- which(as.matrix(stats::dist(fq_vec)) == 1) # fq_nn was not used
  timepoints <- resolve_time_axis(data, passage_times)
  num_species <- length(fq)
  num_timepoints <- ncol(data$x)

  bootstrap_iter <- function(b_iter_idx, current_data, current_fq, current_timepoints,
                             current_num_species, current_num_timepoints,
                             current_epsilon, current_n0, current_nb,
                             current_viability, current_nn_info) { # Renamed arguments
    boot_data <- bootstrap_counts(current_data$x) # Bootstrap from original full data
    x <- normalize_columns(boot_data[current_fq, , drop = FALSE])
    dx_dt <- compute_dx_dt(x, current_timepoints)
    x_trim <- x[, -1, drop = FALSE] # Use x(t+1) for M_t, as per original logic

    qr_terms <- alfak_qr_accum_cpp(x_trim, dx_dt)
    Q_accum <- qr_terms$Q_accum
    r_accum <- qr_terms$r_accum
    Dmat_boot <- 2 * Q_accum + diag(current_epsilon, current_num_species) # Original epsilon
    dvec_boot <- 2 * r_accum
    A_mat <- matrix(1, nrow = current_num_species, ncol = 1) # Renamed A
    bvec_val <- 0 # Renamed bvec

    qp_sol <- run_solve_qp_checked(Dmat_boot, dvec_boot, A_mat, bvec_val, meq = 1,
                                   context = sprintf("solve.QP bootstrap replicate %d", b_iter_idx))
    f_qp <- qp_sol$solution

    x0_init <- optimize_initial_frequencies(x, f_qp, current_timepoints)
    opt_res <- joint_optimize(boot_data[current_fq, , drop = FALSE], current_timepoints, f_qp, x0_init)

    # g0 still uses the first interval because the passage-level scaling step assumes a single
    # growth interval, but it now uses the same validated time_axis as every other time-dependent step.
    g0_val <- log(current_nb / current_n0) / diff(current_timepoints)[1] # Renamed g0

    if (correct_efflux) {
      viability_vec <- current_viability[current_fq]

      # Term 1: sum(x0 * f_rel / viability)
      sum_weighted_frel <- sum((opt_res$x0 * opt_res$f) / viability_vec)

      # Term 2: sum(x0 / viability)
      sum_weights <- sum(opt_res$x0 / viability_vec)

      # Solve for constant k
      k_const <- (sum_weighted_frel - g0_val) / sum_weights

      # Calculate absolute intrinsic division rates: (f_rel - k) / viability
      opt_res$f <- (opt_res$f - k_const) / viability_vec

    } else {
      # Original scaling: shifts mean to match g0
      opt_res$f <- opt_res$f + g0_val - sum(opt_res$x0 * opt_res$f)
    }

    birth_times_est <- find_birth_times(opt_res, time_range = c(-1000, max(current_timepoints)), minF = 1 / current_n0) # Renamed
    peak_times <- current_timepoints[apply(x, 1, which.max)]
    # Neighbour estimation relies on finite parent birth times; fall back to peak-aligned
    # values when root-finding cannot recover them.
    birth_times_est <- sanitize_birth_times(birth_times_est, peak_times = peak_times, timepoints = current_timepoints)

    x0par <- opt_res$x0
    names(x0par) <- current_fq
    fpar <- opt_res$f
    names(fpar) <- current_fq
    names(birth_times_est) <- current_fq

    # dfb <- as.matrix(stats::dist(as.numeric(fpar))) # dfb was not used
    # dfb[upper.tri(dfb)] <- dfb[upper.tri(dfb)] * (-1)
    f_final <- opt_res$f
    x0_final <- opt_res$x0

    xfit <- project_forward_log(x0par, fpar, current_timepoints)
    rownames(xfit) <- current_fq
    ntot <- colSums(boot_data) # Total counts from bootstrapped data
    fq_totals <- colSums(boot_data[current_fq, , drop = FALSE])
    fq_mass <- fq_totals / ntot
    parent_xfit_full <- sweep(xfit, 2, fq_mass, "*")
    ntot_rounded <- round(ntot)
    build_opt_fc <- function(nni_param, prior_mean_param = NaN, prior_sd_param = NaN, do_prior_param = FALSE) {
      if (length(nni_param$nj) == 0) {
        return(function(fc_param) 10^9)
      }
      child <- nni_param$ni
      parent_fitness_mean <- weighted_parent_fitness(nni_param, fpar)
      parent_fitness <- unname(fpar[nni_param$nj])
      parent_birth_times <- unname(birth_times_est[nni_param$nj])
      parent_xfit <- parent_xfit_full[nni_param$nj, , drop = FALSE]
      child_obs <- rep(0, length(current_timepoints))
      if (child %in% rownames(boot_data)) {
        child_obs <- as.numeric(boot_data[child, ])
      }

      function(fc_param) {
        alfak_neighbor_objective_cpp(
          fc_param = fc_param,
          parent_fitness = parent_fitness,
          pij_values = nni_param$pij,
          parent_birth_times = parent_birth_times,
          timepoints = current_timepoints,
          parent_xfit = parent_xfit,
          child_obs = child_obs,
          ntot = ntot_rounded,
          parent_fitness_mean = parent_fitness_mean,
          prior_mean = prior_mean_param,
          prior_sd = prior_sd_param,
          do_prior = do_prior_param,
          tol = ALFAK_FEXP_DELTA_TOL
        )
      }
    }
    search_interval <- range(fpar, na.rm=TRUE) # Added na.rm=TRUE
    interval_range <- diff(search_interval)
    if(length(search_interval) == 1 || interval_range == 0) { # Handle if all fpar are same or only one fpar
      interval_range <- abs(search_interval[1] * 0.5) + 1 # Create a sensible range
    }
    search_interval[1] <- search_interval[1] - interval_range
    search_interval[2] <- search_interval[2] + interval_range
    if(search_interval[1] == search_interval[2]){ # Final check for interval width
      search_interval[1] <- search_interval[1] - 1
      search_interval[2] <- search_interval[2] + 1
    }


    nn_present <- names(current_nn_info) %in% rownames(boot_data)
    if(length(nn_present) > 0 && any(nn_present)){ # Check if nn_present is not empty
      nn_present_indices <- which(nn_present)
      nn_present[nn_present_indices] <- nn_present[nn_present_indices] & sapply(names(current_nn_info)[nn_present_indices], function(ni_sapply) { #Renamed ni
        sum(boot_data[ni_sapply, ]) > 0
      })
    }

    fc <- rep(NaN, length(current_nn_info))
    names(fc) <- names(current_nn_info) # Pre-name fc

    if (any(nn_present)) {
      # Use names for sapply for robustness if current_nn_info can be sparse/differently ordered
      sapply_names <- names(current_nn_info)[nn_present]
      if(length(sapply_names) > 0) { # Ensure there are names to iterate over
        for (child_name in sapply_names) {
          objective_fn <- build_opt_fc(current_nn_info[[child_name]], do_prior_param = FALSE)
          res <- run_optimise_checked(objective_fn, interval = search_interval,
                                      context = sprintf("optimise nearest-neighbour fitness for observed child %s", child_name))
          if (!is.null(res)) {
            fc[child_name] <- res$minimum
          }
        }
      }
    }

    fc_prior_vals <- numeric(0) # Child-minus-parent-mean deltas for the neighbour prior
    if(any(nn_present)){
      # Calculate differences based on nn_info items that were present AND had their fc computed
      nn_present_and_fc_computed <- nn_present & !is.na(fc)
      if(any(nn_present_and_fc_computed)){
        fc_prior_vals <- unlist(lapply(current_nn_info[nn_present_and_fc_computed], function(nni_item) {
          parent_f_mean <- weighted_parent_fitness(nni_item, fpar)
          if(is.finite(parent_f_mean)){
            fc[nni_item$ni] - parent_f_mean
          } else {
            numeric(0) # No valid parents to compute difference from
          }
        }))
      }
    }

    use_empirical_prior <- nn_prior == "empirical"

    if (use_empirical_prior && any(!nn_present) && length(fc_prior_vals) > 0 && !all(is.na(fc_prior_vals))) {
      mean_fc_prior_val <- mean(fc_prior_vals, na.rm = TRUE) # Renamed mean_fc_prior
      if (is.null(nn_prior_sd)) {
        sd_fc_prior_val <- stats::sd(fc_prior_vals, na.rm = TRUE)
        if (!is.finite(sd_fc_prior_val)) {
          sd_fc_prior_val <- nn_prior_sd_floor
        }
        sd_fc_prior_val <- max(sd_fc_prior_val, nn_prior_sd_floor)
      } else {
        sd_fc_prior_val <- nn_prior_sd
      }

      sapply_names_not_present <- names(current_nn_info)[!nn_present]
      if(length(sapply_names_not_present) > 0) {
        for (child_name in sapply_names_not_present) {
          objective_fn <- build_opt_fc(current_nn_info[[child_name]],
                                       prior_mean_param = mean_fc_prior_val,
                                       prior_sd_param = sd_fc_prior_val,
                                       do_prior_param = TRUE)
          res <- run_optimise_checked(objective_fn, interval = search_interval,
                                      context = sprintf("optimise nearest-neighbour fitness with prior for latent child %s", child_name))
          if (!is.null(res)) {
            fc[child_name] <- res$minimum
          }
        }
      }
    } else if (any(!nn_present)) { # No prior to use
      sapply_names_not_present <- names(current_nn_info)[!nn_present]
      if(length(sapply_names_not_present) > 0) {
        for (child_name in sapply_names_not_present) {
          objective_fn <- build_opt_fc(current_nn_info[[child_name]], do_prior_param = FALSE)
          res <- run_optimise_checked(objective_fn, interval = search_interval,
                                      context = sprintf("optimise nearest-neighbour fitness without prior for latent child %s", child_name))
          if (!is.null(res)) {
            fc[child_name] <- res$minimum
          }
        }
      }
    }

    list(f_initial = f_qp,
         f_final = f_final,
         x0_initial = x0_init,
         x0_final = x0_final,
         f_nn = fc)
  }

  # Run bootstrap iterations serially using lapply
  boot_list <- lapply(seq_len(nboot), bootstrap_iter,
                      current_data = data, current_fq = fq, current_timepoints = timepoints,
                      current_num_species = num_species, current_num_timepoints = num_timepoints,
                      current_epsilon = epsilon, current_n0 = n0, current_nb = nb,
                      current_viability = viability, current_nn_info = nn_info_list)

  # Consolidate results
  f_initial_mat <- do.call(rbind, lapply(boot_list, function(x) x$f_initial))
  f_final_mat   <- do.call(rbind, lapply(boot_list, function(x) x$f_final))
  x0_initial_mat <- do.call(rbind, lapply(boot_list, function(x) x$x0_initial))
  x0_final_mat  <- do.call(rbind, lapply(boot_list, function(x) x$x0_final))
  f_nn_mat <- do.call(rbind, lapply(boot_list, function(x) x$f_nn))
  if (is.null(f_nn_mat)) {
    f_nn_mat <- matrix(numeric(0), nrow = length(boot_list), ncol = 0)
  }

  # Set column names if matrices are not empty and fq/names(nn_info_list) are not empty
  if(length(fq) > 0) {
    if(nrow(f_initial_mat) > 0) colnames(f_initial_mat) <- fq
    if(nrow(f_final_mat) > 0) colnames(f_final_mat) <- fq
    if(nrow(x0_initial_mat) > 0) colnames(x0_initial_mat) <- fq
    if(nrow(x0_final_mat) > 0) colnames(x0_final_mat) <- fq
  }
  if(length(names(nn_info_list)) > 0 && nrow(f_nn_mat) > 0) {
    colnames(f_nn_mat) <- names(nn_info_list)
  }

  list(initial_fitness = f_initial_mat,
       final_fitness = f_final_mat,
       initial_frequencies = x0_initial_mat,
       final_frequencies = x0_final_mat,
       nn_fitness = f_nn_mat)
}

#' Fit Kriging model to fitness data (Internal function)
#' @keywords internal
#' @noRd
fitKrig <- function(fq_boot, nboot, krig_bootstrap_mode = c("joint", "marginal")) {
  validate_positive_integer(nboot, "nboot")
  krig_bootstrap_mode <- validate_krig_bootstrap_mode(krig_bootstrap_mode)
  fboot <- cbind(fq_boot$final_fitness, fq_boot$nn_fitness)
  fq_str <- colnames(fq_boot$final_fitness)
  nn_str <- colnames(fq_boot$nn_fitness) # Will be NULL if nn_fitness is NULL or has no colnames

  # Handle cases where fq_str or nn_str might be NULL (e.g., if fq_boot$final_fitness is NULL)
  valid_fq_str <- if(is.null(fq_str)) character(0) else fq_str
  valid_nn_str <- if(is.null(nn_str)) character(0) else nn_str

  combined_strs <- c(valid_fq_str, valid_nn_str)
  if(length(combined_strs) == 0 || ncol(fboot) == 0) { # No data to train on
    warning("fitKrig: No valid fitness data (fq_str or nn_str) to train Kriging model.")
    empty_df <- data.frame(k = character(0), mean = numeric(0), median = numeric(0), sd = numeric(0),
                           fq = logical(0), nn = logical(0))
    return(list(summary_stats = empty_df,
                posterior_samples = matrix(numeric(0), ncol=0, nrow=0),
                boot_results = list(),
                fit_boot_list = list(),
                krig_stable_mean = NULL,
                krig_stable_median = NULL))
  }
  ktrain <- unname(parse_karyotype_ids(combined_strs))

  # Ensure nn_str is not NULL before passing to gen_all_neighbours
  ktest_neighbours_matrix <- matrix(numeric(0), ncol=ncol(ktrain)) # empty matrix with correct cols
  if (length(valid_nn_str) > 0 && !is.null(valid_nn_str)) {
    ktest_neighbours_matrix <- gen_all_neighbours(valid_nn_str)
  }

  ktest <- unique(rbind(ktest = ktrain, ktest_neighbours_matrix))
  ktest_str <- apply(ktest, 1, paste, collapse = ".")
  fq_ids <- ktest_str %in% valid_fq_str
  nn_ids <- ktest_str %in% valid_nn_str

  fboot_mean <- colMeans(fboot, na.rm = TRUE)
  fboot_median <- apply(fboot, 2, stats::median, na.rm = TRUE)
  if (!is.null(names(fboot_mean))) {
    fboot_mean <- fboot_mean[combined_strs]
    fboot_median <- fboot_median[combined_strs]
  }
  valid_mean <- is.finite(fboot_mean)
  valid_median <- is.finite(fboot_median)
  krig_stable_mean <- NULL
  krig_stable_median <- NULL
  safe_build_stable_krig <- function(label, train_x, train_y) {
    tryCatch(
      build_cached_krig_fit(
        train_x,
        train_y,
        give_warnings = TRUE
      )$fit,
      error = function(e) {
        warning(sprintf("fitKrig: Stable %s Kriging fit failed and will be omitted: %s", label, e$message))
        NULL
      }
    )
  }
  if (sum(valid_mean) >= 2 && length(unique(fboot_mean[valid_mean])) >= 2) {
    krig_stable_mean <- safe_build_stable_krig(
      "mean",
      ktrain[valid_mean, , drop = FALSE],
      fboot_mean[valid_mean]
    )
  } else {
    warning("fitKrig: Insufficient data for stable mean Kriging fit.")
  }
  if (sum(valid_median) >= 2 && length(unique(fboot_median[valid_median])) >= 2) {
    krig_stable_median <- safe_build_stable_krig(
      "median",
      ktrain[valid_median, , drop = FALSE],
      fboot_median[valid_median]
    )
  } else {
    warning("fitKrig: Insufficient data for stable median Kriging fit.")
  }

  boot_fit_cache <- new.env(parent = emptyenv())

  # Use lapply directly, as cl is removed
  boot_predictions_list <- lapply(seq_len(nboot), function(b) {
    boot_f <- if (krig_bootstrap_mode == "joint") {
      as.numeric(fboot[sample(seq_len(nrow(fboot)), 1), ])
    } else {
      boot_f_indices <- cbind(sample(seq_len(nrow(fboot)), ncol(fboot), replace = TRUE), seq_len(ncol(fboot)))
      as.vector(fboot[boot_f_indices])
    }

    valid_boot <- is.finite(boot_f)
    ktrain_boot <- ktrain[valid_boot, , drop = FALSE]
    boot_f_valid <- boot_f[valid_boot]

    # Every bootstrap iteration must return the same shape so failed Kriging fits
    # degrade to NA predictions without breaking downstream aggregation.
    if(nrow(ktrain_boot) < 2 ||
       length(boot_f_valid) < 2 ||
       nrow(unique(ktrain_boot)) < 2 ||
       length(unique(boot_f_valid)) < 2) {
      warning("fitKrig: Insufficient or incompatible data for Kriging in bootstrap iteration. Returning NAs.")
      return(list(
        fit_boot = NULL,
        preds = rep(NA_real_, nrow(ktest))
      ))
    }

    tryCatch({
      cache_key <- paste(which(valid_boot), collapse = ",")
      if (!nzchar(cache_key)) {
        cache_key <- "empty"
      }

      if (!exists(cache_key, envir = boot_fit_cache, inherits = FALSE)) {
        cache_entry <- build_cached_krig_fit(
          ktrain_boot,
          boot_f_valid,
          kpred = ktest,
          give_warnings = TRUE
        )
        assign(cache_key, cache_entry, envir = boot_fit_cache)
        fit_boot <- cache_entry$fit
        preds <- predict_cached_krig(fit_boot, ktest, dist_mat = cache_entry$pred_dist)
      } else {
        cache_entry <- get(cache_key, envir = boot_fit_cache, inherits = FALSE)
        refit <- refit_cached_krig(
          cache_entry,
          y = boot_f_valid,
          x_pred = ktest,
          pred_dist = cache_entry$pred_dist,
          give_warnings = TRUE
        )
        fit_boot <- refit$fit
        preds <- refit$preds
      }

      list(fit_boot = fit_boot, preds = preds)
    }, error = function(e) {
      warning(sprintf("fitKrig: Kriging bootstrap iteration failed and will contribute NA predictions: %s", e$message))
      list(
        fit_boot = NULL,
        preds = rep(NA_real_, nrow(ktest))
      )
    })
  })

  #boot_predictions <- do.call(cbind, boot_predictions_list)
  boot_predictions <- do.call(cbind, lapply(boot_predictions_list, `[[`, "preds"))
  fit_boot_list   <- lapply(boot_predictions_list, `[[`, "fit_boot")

  if(is.null(boot_predictions) || ncol(boot_predictions) == 0) { # Check if boot_predictions is empty
    pred_means <- rep(NA_real_, length(ktest_str))
    pred_medians <- rep(NA_real_, length(ktest_str))
    pred_sd <- rep(NA_real_, length(ktest_str))
  } else {
    pred_means <- apply(boot_predictions, 1, mean_or_na)
    pred_medians <- apply(boot_predictions, 1, median_or_na)
    pred_sd <- apply(boot_predictions, 1, sd_or_na)
  }

  summary_df <- data.frame(k = ktest_str, mean = pred_means, median = pred_medians, sd = pred_sd,
                           fq = fq_ids, nn = nn_ids)

  #list(summary_stats = summary_df, posterior_samples = boot_predictions)
  return(list(
    summary_stats     = summary_df,
    posterior_samples = boot_predictions,
    boot_results      = boot_predictions_list,
    fit_boot_list     = fit_boot_list,
    krig_stable_mean  = krig_stable_mean,
    krig_stable_median = krig_stable_median
  ))
}

#' Cross-validation for Kriging model (Internal function)
#' @keywords internal
#' @noRd
xval <- function(fq_boot) {
  fboot <- cbind(fq_boot$final_fitness, fq_boot$nn_fitness)
  fq_str <- colnames(fq_boot$final_fitness)
  nn_str <- colnames(fq_boot$nn_fitness) # Can be NULL

  valid_fq_str <- if(is.null(fq_str)) character(0) else fq_str
  valid_nn_str <- if(is.null(nn_str)) character(0) else nn_str

  combined_strs <- c(valid_fq_str, valid_nn_str)
  if(length(combined_strs) == 0 || ncol(fboot) == 0 || nrow(fboot) == 0) {
    warning("xval: No valid fitness data to perform cross-validation.")
    return(NA_real_)
  }
  ktrain <- unname(parse_karyotype_ids(combined_strs))

  # Original ids logic for xval
  if(length(valid_fq_str) == 0) {
    warning("xval: No fq_str defined, cannot perform original xval logic based on fq_str.")
    return(NA_real_)
  }

  ids <- unlist(lapply(seq_along(valid_fq_str), function(i_xval) {
    ki_neighbours_matrix <- gen_all_neighbours(valid_fq_str[i_xval])
    ki_neighbours_str <- character(0)
    if(nrow(ki_neighbours_matrix) > 0) {
      ki_neighbours_str <- as.character(apply(ki_neighbours_matrix, 1, paste, collapse = "."))
    }
    ki <- c(valid_fq_str[i_xval], ki_neighbours_str)
    idi <- rep(i_xval, length(ki))
    names(idi) <- ki
    idi
  }))
  # Fold ownership keeps the first assignment for overlapping neighbourhoods, which preserves
  # the long-standing heuristic while making the order-dependence explicit.
  ids <- ids[!duplicated(names(ids))] # Ensure unique names for ids
  uids <- unique(ids) # These are fold identifiers
  b <- sample(seq_len(nrow(fboot)), 1)
  fi <- fboot[b, ]

  # Use lapply directly
  tmp_list <- lapply(uids, function(id_fold) { # Renamed id to id_fold
    # Original logic for train/test split based on 'ids' and current 'id_fold'
    train_indices <- !(ids == id_fold)
    test_indices <- (ids == id_fold)

    # Check if ktrain (all possible karyotypes from combined_strs) aligns with names(ids)
    # names(ids) are the karyotype strings that were assigned a fold id
    # We need to map these back to rows in ktrain if ktrain corresponds to combined_strs

    # Create a mapping from karyotype string to row index in ktrain
    ktrain_map <- setNames(1:nrow(ktrain), combined_strs)

    train_k_names <- names(ids)[train_indices]
    test_k_names <- names(ids)[test_indices]

    # Ensure these names exist in ktrain_map
    train_k_names_valid <- train_k_names[train_k_names %in% names(ktrain_map)]
    test_k_names_valid <- test_k_names[test_k_names %in% names(ktrain_map)]

    if(length(train_k_names_valid) == 0 || length(test_k_names_valid) == 0) {
      return(matrix(NA, ncol=2, dimnames=list(NULL, c("test_f", "est_f")))) # Skip fold
    }

    train_k_rows <- ktrain_map[train_k_names_valid]
    test_k_rows <- ktrain_map[test_k_names_valid]

    train_k <- ktrain[train_k_rows, , drop = FALSE]
    # fi corresponds to combined_strs (colnames of fboot)
    # So, we need to subset fi based on names corresponding to train_k_names_valid
    train_f <- fi[train_k_names_valid]

    test_k <- ktrain[test_k_rows, , drop = FALSE]
    test_f <- fi[test_k_names_valid]

    # Filter NAs from training data, which could arise if some fi values were NA
    valid_train_points <- is.finite(train_f)
    train_k <- train_k[valid_train_points, , drop=FALSE]
    train_f <- train_f[valid_train_points]

    if(nrow(train_k) < 2 || nrow(unique(train_k)) < 2 || length(unique(train_f)) < 2) {
      warning(paste("Skipping xval fold for id_fold:", id_fold, "due to insufficient/non-unique training points after NA removal."))
      return(cbind(test_f = test_f, est_f = rep(NA, length(test_f))))
    }

    est_f <- tryCatch({
      fit_cache <- build_cached_krig_fit(train_k, train_f, kpred = test_k, give_warnings = TRUE)
      fit <- fit_cache$fit
      predict_cached_krig(fit, test_k, dist_mat = fit_cache$pred_dist)
    }, error = function(e) {
      warning(sprintf("xval fold fit failed for fold %s and will return NA predictions: %s", id_fold, e$message))
      rep(NA_real_, length(test_f))
    })
    cbind(test_f, est_f)
  })

  tmp <- do.call(rbind, tmp_list)
  tmp <- tmp[stats::complete.cases(tmp), , drop = FALSE] # Use stats::complete.cases

  r2r_val <- NA_real_
  if(nrow(tmp) < 2) { # R2R needs at least 2 points
    warning("Not enough valid observations after cross-validation to compute R2R.")
  }else{
    r2r_val <- R2R(tmp[, 1], tmp[, 2])
  }
  return(list(
    tmp  = tmp,
   R2R  = r2r_val
  ))
  #R2R(tmp[, 1], tmp[, 2])
}
