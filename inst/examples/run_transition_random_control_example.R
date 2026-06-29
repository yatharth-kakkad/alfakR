# Compare transition-karyotype targeting against a random-panel null model.
# null_model = TRUE runs the null-model condition automatically inside
# run_transition_karyotype_abm() alongside the no-treatment and ranked
# conditions, inheriting all the same inputs -- no separate call needed.
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
  null_model = TRUE,
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
  result <- do.call(run_transition_karyotype_abm, c(common_args, list(abm_seed = s)))
  es <- result$endpoint_summary
  data.frame(
    seed = s,
    condition1_final_population = es$condition1_final_population,
    condition1_final_diversity = es$condition1_final_diversity,
    ranked_final_population = es$condition2_final_population,
    ranked_final_diversity = es$condition2_final_diversity,
    ranked_population_reduction_pct = es$population_reduction_pct,
    ranked_diversity_reduction_pct = es$diversity_reduction_pct,
    random_final_population = es$condition3_final_population,
    random_final_diversity = es$condition3_final_diversity,
    random_population_reduction_pct = es$null_model_population_reduction_pct,
    random_diversity_reduction_pct = es$null_model_diversity_reduction_pct
  )
})
per_seed <- do.call(rbind, rows)

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
