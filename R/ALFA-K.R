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
#' @param passage_times An optional numeric vector of passage times. If `NULL`,
#'   calculated from `colnames(yi$x) * yi$dt`.
#' @param minobs An integer, the minimum total number of observations (reads/counts)
#'   for a karyotype across all timepoints to be considered "frequent" and included
#'   in the analysis. Default is 20.
#' @param nboot An integer, the number of bootstrap iterations for fitness
#'   estimation and for the Kriging process in `fitKrig`. Default is 45.
#' @param n0 A numeric value, the initial effective population size at the
#'   start of a passage or growth phase, used for g0 calculation. Default is 1e5.
#' @param nb A numeric value, the bottleneck effective population size (population
#'   size after transfer), used for g0 calculation. Default is 1e7.
#' @param pm A numeric value, the per-locus mutation/error rate used in `pij`
#'   calculations. Default is 0.00005.
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
alfak <- function(yi, outdir, passage_times, minobs = 20,
                  nboot = 45,
                  n0 = 1e5,
                  nb = 1e7,
                  pm = 0.00005,
                  correct_efflux=FALSE) {
  
  # Note: library calls removed, dependencies handled by @importFrom or DESCRIPTION
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (max(rowSums(yi$x)) < minobs)
    stop(paste("no frequent karyotypes detected for minobs", minobs))
  
  # Parallelism and cl related code removed
  
  fq_boot <- solve_fitness_bootstrap(yi, minobs = minobs, nboot = nboot,
                                     n0 = n0, nb = nb, pm = pm,
                                     passage_times = passage_times,correct_efflux=correct_efflux)
  saveRDS(fq_boot, file = file.path(outdir, "bootstrap_res.Rds"))
  
  landscape_data <- fitKrig(fq_boot, nboot) # nboot is passed for Kriging iterations
  saveRDS(landscape_data$summary_stats, file = file.path(outdir, "landscape.Rds"))
  saveRDS(landscape_data$posterior_samples, file = file.path(outdir, "landscape_posterior_samples.Rds"))
  #saveRDS(landscape_data, file = file.path(outdir, "landscape_data.Rds"))
  Rxv <- xval(fq_boot)
  saveRDS(Rxv, file = file.path(outdir, "xval.Rds"))
  
  ##END HERE.
  
  invisible(Rxv) # Return Rxv invisibly as side-effect saving is primary
}

##########################################
# Helper functions (internal)
##########################################

#' Calculate transition probability p_ij
#' @keywords internal
#' @noRd
pij <- function(i, j, beta) {
  qij <- 0
  if (abs(i - j) > i) { ## not enough copies for i->j
    return(qij)
  }
  if (j == 0) j <- 2 * i
  s <- seq(abs(i - j), i, by = 2)
  for (z in s) {
    qij <- qij + choose(i, z) * beta^z * (1 - beta)^(i - z) *
      0.5^z * choose(z, (z + i - j) / 2)
  }
  return(qij)
}

#' Convert string like "1.2.3" to numeric vector
#' @keywords internal
#' @noRd
s2v <- function(s) as.numeric(unlist(strsplit(s, split = "[.]")))

#' Calculate R-squared
#' @keywords internal
#' @noRd
R2R <- function(obs, pred) {
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
}

#' Generate all single-step neighbours for karyotype IDs
#' @keywords internal
#' @noRd
gen_all_neighbours <- function(ids, as.strings = TRUE, remove_nullisomes = TRUE) {
  if (as.strings) 
    ids <- lapply(ids, function(ii) as.numeric(unlist(strsplit(ii, split = "[.]"))))
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
  num_species <- nrow(x)
  dx_dt <- matrix(NA, nrow = num_species, ncol = ncol(x) - 1)
  for (i in seq_len(num_species)) {
    dx_dt[i, ] <- diff(x[i, ]) / diff(timepoints)
  }
  dx_dt
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
    
    nivec <- s2v(ni_string)
    pij_vals <- sapply(nj_filtered, function(si_string) { # Renamed si to si_string
      si_vec <- as.numeric(unlist(strsplit(si_string, split = "[.]")))
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
  K <- nrow(counts)
  Tt <- ncol(counts)
  f_free <- param[1:(K - 1)]
  f_full <- c(f_free, -sum(f_free))
  log_x0 <- param[K:(2 * K - 1)]
  nll <- 0
  for (i in seq_len(Tt)) {
    lv <- log_x0 + f_full * timepoints[i]
    denom <- logSumExp(lv)
    for (k in seq_len(K)) {
      if (counts[k, i] > 0) {
        nll <- nll - counts[k, i] * (lv[k] - denom)
      }
    }
  }
  nll
}

#' Jointly optimize fitness and initial frequencies
#' @keywords internal
#' @noRd
joint_optimize <- function(counts, timepoints, f_init, x0_init) {
  K <- length(f_init)
  f_free_init <- f_init[1:(K - 1)]
  x0_init_log <- log(x0_init + 1e-12) # Original had this epsilon
  param_init <- c(f_free_init, x0_init_log)
  obj_fun <- function(par) neg_log_lik(par, counts, timepoints)
  opt <- stats::optim(par = param_init, fn = obj_fun, method = "BFGS",
                      control = list(maxit = 200, reltol = 1e-8))
  f_free_opt <- opt$par[1:(K - 1)]
  f_opt <- c(f_free_opt, -sum(f_free_opt))
  log_x0_opt <- opt$par[K:(2 * K - 1)]
  x0_opt <- exp(log_x0_opt)
  x0_opt <- x0_opt / sum(x0_opt)
  list(f = f_opt, x0 = x0_opt)
}

#' Project frequencies forward in time using log-space calculations
#' @keywords internal
#' @noRd
project_forward_log <- function(x0, f, timepoints) {
  K <- length(x0)
  out <- matrix(NA, nrow = K, ncol = length(timepoints))
  log_x0 <- log(x0) # Original did not have epsilon here
  for (i in seq_along(timepoints)) {
    lv <- log_x0 + f * timepoints[i]
    denom <- logSumExp(lv)
    out[, i] <- exp(lv - denom)
  }
  out
}

#' Optimize initial frequencies given observed data and fitness values
#' @keywords internal
#' @noRd
optimize_initial_frequencies <- function(x_obs, f, timepoints) {
  loss_function <- function(log_x0) {
    x0 <- exp(log_x0)
    x0 <- x0 / sum(x0)
    x_pred <- project_forward_log(x0, f, timepoints)
    sum((x_pred - x_obs)^2)
  }
  x_ini <- x_obs[, 1] + 1e-6 # Original had this epsilon
  x_ini <- x_ini / sum(x_ini)
  opt_result <- stats::optim(log(x_ini), loss_function, method = "BFGS")
  x0_opt <- exp(opt_result$par)
  x0_opt / sum(x0_opt)
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
                                    n0, nb, passage_times = NULL,correct_efflux=FALSE) {
  fq <- rownames(data$x)[rowSums(data$x) > minobs]
  nn_info_list <- gen_nn_info(fq, pm) # Renamed 'nn' to 'nn_info_list' for clarity
  if (length(nn_info_list) > 0 && !is.null(nn_info_list[[1]]$ni)) { # Check if naming is needed
    names(nn_info_list) <- sapply(nn_info_list, function(nni) nni$ni)
  } else if (length(nn_info_list) == 0) {
    # names(nn_info_list) would be NULL, which is fine for later checks
  } else {
    warning("nn_info_list structure unexpected for naming in solve_fitness_bootstrap")
  }
  
  fq_vec <- do.call(rbind, lapply(fq, s2v))
  
  P <- (1-pm)^rowSums(fq_vec)
  
  # fq_nn <- which(as.matrix(stats::dist(fq_vec)) == 1) # fq_nn was not used
  timepoints <- as.numeric(colnames(data$x)) * data$dt
  num_species <- length(fq)
  num_timepoints <- ncol(data$x)
  
  bootstrap_iter <- function(b_iter_idx, current_data, current_fq, current_timepoints, 
                             current_num_species, current_num_timepoints,
                             current_epsilon, current_n0, current_nb, 
                             current_passage_times, current_nn_info) { # Renamed arguments
    boot_data <- bootstrap_counts(current_data$x) # Bootstrap from original full data
    x <- apply(boot_data[current_fq, , drop = FALSE], 2, function(col) { # Subset fq after bootstrap
      s <- sum(col)
      if (s == 0) rep(0, length(col)) else col / s # Normalize, handle sum=0
    })
    dx_dt <- compute_dx_dt(x, current_timepoints)
    x_trim <- x[, -1, drop = FALSE] # Use x(t+1) for M_t, as per original logic
    
    Q_accum <- matrix(0, nrow = current_num_species, ncol = current_num_species)
    r_accum <- rep(0, current_num_species) # Original was rep(0, num_species)
    
    for (t_idx in 1:(current_num_timepoints - 1)) { # Renamed t to t_idx
      xt <- x_trim[, t_idx]
      M_t <- diag(xt) - outer(xt, xt)
      Q_accum <- Q_accum + M_t %*% M_t
      r_accum <- r_accum + M_t %*% dx_dt[, t_idx]
    }
    Dmat_boot <- 2 * Q_accum + diag(current_epsilon, current_num_species) # Original epsilon
    dvec_boot <- 2 * r_accum
    A_mat <- matrix(1, nrow = current_num_species, ncol = 1) # Renamed A
    bvec_val <- 0 # Renamed bvec
    
    # Using quadprog::solve.QP as it will be imported
    qp_sol <- quadprog::solve.QP(Dmat_boot, dvec_boot, A_mat, bvec_val, meq = 1)
    f_qp <- qp_sol$solution
    
    x0_init <- optimize_initial_frequencies(x, f_qp, current_timepoints)
    opt_res <- joint_optimize(boot_data[current_fq, , drop = FALSE], current_timepoints, f_qp, x0_init)
    
    g0_val <- log(current_nb / current_n0) / diff(current_timepoints)[1] # Renamed g0
    if (!is.null(current_passage_times))
      g0_val <- log(current_nb / current_n0) / diff(current_passage_times * current_data$dt)[1]
    
    if (correct_efflux) {
      viability <- 2 * P - 1
      
      # Term 1: sum(x0 * f_rel / viability)
      sum_weighted_frel <- sum((opt_res$x0 * opt_res$f) / viability)
      
      # Term 2: sum(x0 / viability)
      sum_weights <- sum(opt_res$x0 / viability)
      
      # Solve for constant k
      k_const <- (sum_weighted_frel - g0_val) / sum_weights
      
      # Calculate absolute intrinsic division rates: (f_rel - k) / viability
      opt_res$f <- (opt_res$f - k_const) / viability
      
    } else {
      # Original scaling: shifts mean to match g0
      opt_res$f <- opt_res$f + g0_val - sum(opt_res$x0 * opt_res$f)
    }
    
    birth_times_est <- find_birth_times(opt_res, time_range = c(-1000, max(current_timepoints)), minF = 1 / current_n0) # Renamed
    peak_times <- current_timepoints[apply(x, 1, which.max)]
    mean_risetime <- mean(peak_times - birth_times_est, na.rm = TRUE)
    birth_times_est[is.na(birth_times_est)] <- peak_times[is.na(birth_times_est)] - mean_risetime
    
    x0par <- opt_res$x0
    names(x0par) <- current_fq
    fpar <- opt_res$f
    names(fpar) <- current_fq
    names(birth_times_est) <- current_fq
    
    # dfb <- as.matrix(stats::dist(as.numeric(fpar))) # dfb was not used
    # dfb[upper.tri(dfb)] <- dfb[upper.tri(dfb)] * (-1)
    f_final <- opt_res$f
    x0_final <- opt_res$x0
    
    fExp <- function(fc_arg, fp_arg, pij_val, tt_arg) { # Renamed args
      pij_val * fp_arg / (fc_arg - fp_arg) * (exp(tt_arg * (fc_arg - fp_arg)) - 1)
    }
    xfit <- project_forward_log(x0par, fpar, current_timepoints)
    rownames(xfit) <- current_fq
    ntot <- colSums(boot_data) # Total counts from bootstrapped data
    
    opt_fc <- function(fc_param, nni_param, prior_mean_param = NULL, prior_sd_param = NULL, do_prior_param = FALSE) { #Renamed args
      child <- nni_param$ni
      xc_est <- colSums(do.call(rbind, lapply(1:length(nni_param$nj), function(i_loop) { #Renamed i
        parent_karyo <- nni_param$nj[i_loop] #Renamed par
        tt_val <- current_timepoints - birth_times_est[parent_karyo] #Renamed tt
        fExp(fc_param, fpar[parent_karyo], nni_param$pij[i_loop], tt_val) * xfit[parent_karyo, ]
      })))
      xc_est <- pmax(0, pmin(1, xc_est)) # Ensure probabilities are in [0,1]
      xc_obs <- rep(0, length(current_timepoints))
      if (child %in% rownames(boot_data)) # Check against full bootstrapped data
        xc_obs <- boot_data[nni_param$ni, ]
      
      res <- stats::dbinom(xc_obs, round(ntot), prob = xc_est, log = TRUE) # round(ntot) if ntot can be non-integer
      if (do_prior_param)
        res <- c(res, stats::dnorm(fc_param - fpar[nni_param$nj], mean = prior_mean_param, sd = prior_sd_param, log = TRUE))
      res[!is.finite(res)] <- -(10^9) # Penalize non-finite values
      -sum(res)
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
        results_sapply <- sapply(current_nn_info[sapply_names], function(nni_sapply2) { #Renamed nni
          res <- stats::optimise(opt_fc, interval = search_interval, nni_param = nni_sapply2, do_prior_param = FALSE)
          res$minimum
        })
        fc[sapply_names] <- results_sapply
      }
    }
    
    fc_prior_vals <- numeric(0) # Renamed fc_prior
    if(any(nn_present)){
      # Calculate differences based on nn_info items that were present AND had their fc computed
      nn_present_and_fc_computed <- nn_present & !is.na(fc)
      if(any(nn_present_and_fc_computed)){
        fc_prior_vals <- unlist(lapply(current_nn_info[nn_present_and_fc_computed], function(nni_item) {
          # Ensure fpar has names and nni_item$nj are valid names in fpar
          valid_parents <- nni_item$nj[nni_item$nj %in% names(fpar)]
          if(length(valid_parents)>0){
            # Original: fpar[nni_item$nj] - fc[nni_item$ni]
            # This implies multiple differences if nni_item$nj is a vector.
            # A common pattern is difference to the mean of parents, or each parent.
            # Let's assume difference to the mean parent fitness for now
            mean_parent_f <- mean(fpar[valid_parents], na.rm=TRUE)
            mean_parent_f - fc[nni_item$ni] 
          } else {
            numeric(0) # No valid parents to compute difference from
          }
        }))
      }
    }
    
    if (any(!nn_present) && length(fc_prior_vals) > 0 && !all(is.na(fc_prior_vals))) {
      mean_fc_prior_val <- mean(fc_prior_vals, na.rm = TRUE) # Renamed mean_fc_prior
      sd_fc_prior_val <- sd(fc_prior_vals, na.rm = TRUE)   # Renamed sd_fc_prior
      if(is.na(sd_fc_prior_val) || sd_fc_prior_val == 0) sd_fc_prior_val <- 1e-3 # Default small SD
      
      sapply_names_not_present <- names(current_nn_info)[!nn_present]
      if(length(sapply_names_not_present) > 0) {
        results_sapply_not_present <- sapply(current_nn_info[sapply_names_not_present], function(nni_sapply3) { #Renamed nni
          res <- stats::optimise(opt_fc, interval = search_interval, nni_param = nni_sapply3,
                                 prior_mean_param = mean_fc_prior_val,
                                 prior_sd_param = sd_fc_prior_val,
                                 do_prior_param = TRUE)
          res$minimum
        })
        fc[sapply_names_not_present] <- results_sapply_not_present
      }
    } else if (any(!nn_present)) { # No prior to use
      sapply_names_not_present <- names(current_nn_info)[!nn_present]
      if(length(sapply_names_not_present) > 0) {
        results_sapply_no_prior <- sapply(current_nn_info[sapply_names_not_present], function(nni_sapply3) { 
          res <- stats::optimise(opt_fc, interval = search_interval, nni_param = nni_sapply3, do_prior_param = FALSE)
          res$minimum
        })
        fc[sapply_names_not_present] <- results_sapply_no_prior
      }
    }
    
    list(f_initial = f_qp,
         f_final = f_final,
         x0_initial = x0_init,
         x0_final = x0_final,
         f_nn = fc)
  }
  
  # Run bootstrap iterations serially using lapply
  boot_list <- lapply(1:nboot, bootstrap_iter, 
                      current_data = data, current_fq = fq, current_timepoints = timepoints,
                      current_num_species = num_species, current_num_timepoints = num_timepoints,
                      current_epsilon = epsilon, current_n0 = n0, current_nb = nb,
                      current_passage_times = passage_times, current_nn_info = nn_info_list)
  
  # Consolidate results
  f_initial_mat <- do.call(rbind, lapply(boot_list, function(x) x$f_initial))
  f_final_mat   <- do.call(rbind, lapply(boot_list, function(x) x$f_final))
  x0_initial_mat <- do.call(rbind, lapply(boot_list, function(x) x$x0_initial))
  x0_final_mat  <- do.call(rbind, lapply(boot_list, function(x) x$x0_final))
  f_nn_mat <- do.call(rbind, lapply(boot_list, function(x) x$f_nn))
  
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
fitKrig <- function(fq_boot, nboot) {
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
    return(list(summary_stats = empty_df, posterior_samples = matrix(numeric(0), ncol=0, nrow=0)))
  }
  ktrain <- do.call(rbind, lapply(combined_strs, s2v))
  
  # Ensure nn_str is not NULL before passing to gen_all_neighbours
  ktest_neighbours_matrix <- matrix(numeric(0), ncol=ncol(ktrain)) # empty matrix with correct cols
  if (length(valid_nn_str) > 0 && !is.null(valid_nn_str)) {
    ktest_neighbours_matrix <- gen_all_neighbours(valid_nn_str)
  }
  
  ktest <- unique(rbind(ktest = ktrain, ktest_neighbours_matrix))
  ktest_str <- apply(ktest, 1, paste, collapse = ".")
  fq_ids <- ktest_str %in% valid_fq_str
  nn_ids <- ktest_str %in% valid_nn_str
  
  # Use lapply directly, as cl is removed
  boot_predictions_list <- lapply(1:nboot, function(b) {
    # Original sampling strategy to avoid spatially correlated errors
    boot_f_indices <- cbind(sample(1:nrow(fboot), ncol(fboot), replace = TRUE), 1:ncol(fboot))
    boot_f <- fboot[boot_f_indices, drop=FALSE] # Use drop=FALSE just in case
    if(ncol(fboot) == 1 && nrow(fboot) > 0) boot_f <- as.vector(boot_f[,1]) else boot_f <- as.vector(boot_f)
    
    
    # Ensure ktrain and boot_f have compatible dimensions and enough data for Krig
    if(nrow(ktrain) < 2 || length(boot_f) < 2 || length(unique(boot_f)) < 2 || nrow(ktrain) != length(boot_f)) {
      warning("fitKrig: Insufficient or incompatible data for Kriging in bootstrap iteration. Returning NAs.")
      return(rep(NA_real_, nrow(ktest)))
    }
    
    fit_boot <- fields::Krig(ktrain, boot_f,
                             cov.function = "stationary.cov",
                             cov.args = list(Covariance = "Matern", smoothness = 1.5))
    preds <- stats::predict(fit_boot, ktest)
    list(fit_boot = fit_boot, preds = preds)
  })
  
  #boot_predictions <- do.call(cbind, boot_predictions_list)
  boot_predictions <- do.call(cbind, lapply(boot_predictions_list, `[[`, "preds"))
  fit_boot_list   <- lapply(boot_predictions_list, `[[`, "fit_boot")

  if(is.null(boot_predictions) || ncol(boot_predictions) == 0) { # Check if boot_predictions is empty
    pred_means <- rep(NA_real_, length(ktest_str))
    pred_medians <- rep(NA_real_, length(ktest_str))
    pred_sd <- rep(NA_real_, length(ktest_str))
  } else {
    pred_means <- apply(boot_predictions, 1, mean, na.rm = TRUE) # Add na.rm=TRUE
    pred_medians <- apply(boot_predictions, 1, stats::median, na.rm = TRUE) # Add na.rm=TRUE
    pred_sd <- apply(boot_predictions, 1, stats::sd, na.rm = TRUE) # Add na.rm=TRUE
  }
  
  summary_df <- data.frame(k = ktest_str, mean = pred_means, median = pred_medians, sd = pred_sd,
                           fq = fq_ids, nn = nn_ids)
  
  #list(summary_stats = summary_df, posterior_samples = boot_predictions)
  return(list(
  summary_stats     = summary_df,
  posterior_samples = boot_predictions,
  boot_results     = boot_predictions_list,
  fit_boot_list = fit_boot_list
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
  ktrain <- do.call(rbind, lapply(combined_strs, s2v))
  
  # Original ids logic for xval
  if(length(valid_fq_str) == 0) {
    warning("xval: No fq_str defined, cannot perform original xval logic based on fq_str.")
    return(NA_real_)
  }
  
  ids <- unlist(lapply(1:length(valid_fq_str), function(i_xval) {
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
  ids <- ids[!duplicated(names(ids))] # Ensure unique names for ids
  uids <- unique(ids) # These are fold identifiers
  
  # Use lapply directly
  tmp_list <- lapply(uids, function(id_fold) { # Renamed id to id_fold
    fboot_shuffled <- apply(fboot, 2, sample) # Original shuffling strategy
    fi <- fboot_shuffled[1, ]
    
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
    valid_train_points <- !is.na(train_f)
    train_k <- train_k[valid_train_points, , drop=FALSE]
    train_f <- train_f[valid_train_points]
    
    if(nrow(train_k) < 2 || nrow(unique(train_k)) < 2 || length(unique(train_f))<1) { # Krig needs some variation
      warning(paste("Skipping xval fold for id_fold:", id_fold, "due to insufficient/non-unique training points after NA removal."))
      return(cbind(test_f = test_f, est_f = rep(NA, length(test_f))))
    }
    
    fit <- fields::Krig(train_k, train_f,
                        cov.function = "stationary.cov",
                        cov.args = list(Covariance = "Matern", smoothness = 1.5))
    est_f <- stats::predict(fit, test_k)
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