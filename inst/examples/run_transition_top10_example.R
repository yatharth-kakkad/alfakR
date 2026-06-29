# Reproduce the top-10 transition-karyotype treatment ABM example.
#
# Run from the package root after installing/loading alfakR:
# Rscript inst/examples/run_transition_top10_example.R

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

result <- run_transition_karyotype_abm(
  node_metadata = node_metadata_path,
  edges = edges_path,
  transition_scores = transition_scores_path,
  horizon_timestep = 60,
  transition_top_n = 10,
  tau1 = 30,
  p_missegregation = 0.02,
  base_death_rate = 0.01,
  base_birth_rate = 0.02,
  fitness_birth_scale = 0.1,
  second_treatment_strength = 0.9,
  abm_delta_t = 1,
  abm_record_interval = 10,
  abm_seed = 123
)

print(result$inputs[c(
  "tau1",
  "tau1_step",
  "transition_top_n",
  "abm_pop_size",
  "second_treatment_strength"
)])
print(result$endpoint_summary)

output_dir <- file.path("inst", "extdata")
utils::write.csv(
  result$endpoint_summary,
  file.path(output_dir, "transition_top10_example_endpoint_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  result$metrics,
  file.path(output_dir, "transition_top10_example_metrics.csv"),
  row.names = FALSE
)
