library(IBMR)

source("scripts/simulation_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "results/simulations_1_new"
dir.create(RESULT_PATH, recursive = TRUE)

methods = c("IBMR", "IBMR_int", "IBMR_common_Gamma", "IBMR_no_Gamma", "glmnet_subset", "glmnet_split", "glmnet_relabel")
methods = c(methods, paste0(methods, "_ORC_clean"))

defaults = list(
  q = 16,
  N = 2000,
  K = 4,
  p = 500,
  nonzero = 75,
  b = 2,
  rank = 20,
  batch_effect = 0.1
)

considered_values = list(
  batch_effect = c(0.025, 0.05, 0.1, 0.2, 0.4),
  rank = c("int", 1, 5, 10, 20, 40, 80)
)

parameters = expand_parameters("fine_low_rank_batch_effect", considered_values, defaults, 50, methods)

defaults$rank = "int"

considered_values = list(
  batch_effect = c(0.025, 0.05, 0.1, 0.2, 0.4)
)

parameters = c(parameters, expand_parameters("fine_intercept_batch_effect", considered_values, defaults, 50, methods))


chunk_size = 6

for (i in 1:chunk_size) {

  PARAMETER_ID = (ARRAY_ID - 1) * chunk_size + i
  current_parameters = parameters[[PARAMETER_ID]]
  system.time({result = evaluate_parameters(current_parameters, generate_simulation_data_fine)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}
