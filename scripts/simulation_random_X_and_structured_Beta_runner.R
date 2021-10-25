library(IBMR)

source("scripts/simulation_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "final_results/simulations_random_X_and_structured_Beta_proxGamma_fixed"
dir.create(RESULT_PATH, recursive = TRUE)

methods = c("IBMR_int", "IBMR_no_Gamma", "subset", "relabel")
methods = c(methods, paste0(methods, "_ORC_fine"), paste0(methods, "_ORC_clean"), paste0(methods, "_ORC_fine_clean"))

defaults = list(
  category_mappings = simulate_category_mappings(3, c(6, 2, 2), list(rep(1, 6), rep(1, 6), rep(1, 6), rep(2, 6), c(rep(3, 3), rep(2, 3)), c(rep(2, 3), rep(3, 3)))),
  N = 4800,
  p = 500,
  nonsparsity = 0.2,
  pct_de = 0.1,
  b = 2,
  sigma = 1,
  rank = "int",
  batch_effect = 0.1
)

considered_values = list(
  batch_effect = c(0.025, 0.05, 0.1, 0.2, 0.4),
  nonsparsity = c(0.05, 0.1, 0.2, 0.4, 0.8),
  N = c(2400, 4800, 9600, 19200)
)

parameters = expand_parameters("random_X_and_structured_Beta", considered_values, defaults, 50, methods)

chunk_size = 8

for (i in 1:chunk_size) {

  PARAMETER_ID = (ARRAY_ID - 1) * chunk_size + i
  current_parameters = parameters[[PARAMETER_ID]]
  system.time({result = evaluate_parameters(current_parameters, generate_data_random_X_and_structured_Beta)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}
