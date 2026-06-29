# Compare transition-karyotype targeting against a random-panel null model,
# across replicate seeds and a longer horizon so cells have time to spread
# across the graph (otherwise a random panel rarely overlaps occupied nodes).
#
# Run from the package root after installing/loading alfakR:
# Rscript inst/examples/run_transition_random_control_example.R

if (requireNamespace("alfakR", quietly = TRUE)) {
  library(alfakR)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  stop("Install alfakR or pkgload before running this example.", call. = FALSE)
}

node_metadata_path <- system.file("extdata", "Node_Metadata.csv", package = "alfakR")
edges_path <- system.file("extdata", "Edge_DF.csv", package = "alfakR")
transition_scores_path <- system.file("extdata", "transition_karyotype_scores.csv", package = "alfakR")

if (!nzchar(node_metadata_path) || !nzchar(edges_path) || !nzchar(transition_scores_path)) {
  node_metadata_path <- "inst/extdata/Node_Metadata.csv"
  edges_path <- "inst/extdata/Edge_DF.csv"
  transition_scores_path <- "inst/extdata/transition_karyotype_scores.csv"
}

common_args <- list(
  node_metadata = node_metadata_path,
  edges = edges_path,
  transition_scores = transition_scores_path,
  horizon_timestep = 180,
  transition_top_n = 10,
  tau1 = 30,
  p_missegregation = 0.02,
  base_death_rate = 0.01,
  base_birth_rate = 0.02,
  fitness_birth_scale = 0.1,
  second_treatment_strength = 0.9,
  abm_delta_t = 1,
  abm_record_interval = 18
)

seeds <- 1:20
rows <- lapply(seeds, function(s) {
  ranked_result <- do.call(run_transition_karyotype_abm, c(common_args, list(
    target_selection = "ranked", abm_seed = s
  )))
  random_result <- do.call(run_transition_karyotype_abm, c(common_args, list(
    target_selection = "random", random_target_seed = s, abm_seed = s
  )))
  data.frame(
    seed = s,
    condition1_final_population = ranked_result$endpoint_summary$condition1_final_population,
    condition1_final_diversity = ranked_result$endpoint_summary$condition1_final_diversity,
    ranked_final_population = ranked_result$endpoint_summary$condition2_final_population,
    ranked_final_diversity = ranked_result$endpoint_summary$condition2_final_diversity,
    ranked_population_reduction_pct = ranked_result$endpoint_summary$population_reduction_pct,
    ranked_diversity_reduction_pct = ranked_result$endpoint_summary$diversity_reduction_pct,
    random_final_population = random_result$endpoint_summary$condition2_final_population,
    random_final_diversity = random_result$endpoint_summary$condition2_final_diversity,
    random_population_reduction_pct = random_result$endpoint_summary$population_reduction_pct,
    random_diversity_reduction_pct = random_result$endpoint_summary$diversity_reduction_pct
  )
})
per_seed <- do.call(rbind, rows)

summarize <- function(x) c(mean = mean(x), sd = stats::sd(x))
panel_summary <- data.frame(
  panel = c("condition1_no_second_treatment", "condition2_ranked_transition_targets", "condition3_random_panel_null"),
  mean_final_population = c(
    mean(per_seed$condition1_final_population),
    mean(per_seed$ranked_final_population),
    mean(per_seed$random_final_population)
  ),
  sd_final_population = c(
    stats::sd(per_seed$condition1_final_population),
    stats::sd(per_seed$ranked_final_population),
    stats::sd(per_seed$random_final_population)
  ),
  mean_final_diversity = c(
    mean(per_seed$condition1_final_diversity),
    mean(per_seed$ranked_final_diversity),
    mean(per_seed$random_final_diversity)
  ),
  sd_final_diversity = c(
    stats::sd(per_seed$condition1_final_diversity),
    stats::sd(per_seed$ranked_final_diversity),
    stats::sd(per_seed$random_final_diversity)
  ),
  mean_population_reduction_pct_vs_c1 = c(NA, mean(per_seed$ranked_population_reduction_pct), mean(per_seed$random_population_reduction_pct)),
  sd_population_reduction_pct_vs_c1 = c(NA, stats::sd(per_seed$ranked_population_reduction_pct), stats::sd(per_seed$random_population_reduction_pct)),
  mean_diversity_reduction_pct_vs_c1 = c(NA, mean(per_seed$ranked_diversity_reduction_pct), mean(per_seed$random_diversity_reduction_pct)),
  sd_diversity_reduction_pct_vs_c1 = c(NA, stats::sd(per_seed$ranked_diversity_reduction_pct), stats::sd(per_seed$random_diversity_reduction_pct))
)

print(panel_summary)

output_dir <- file.path("inst", "extdata")
utils::write.csv(per_seed, file.path(output_dir, "transition_random_control_example_per_seed.csv"), row.names = FALSE)
utils::write.csv(panel_summary, file.path(output_dir, "transition_random_control_example_summary.csv"), row.names = FALSE)
