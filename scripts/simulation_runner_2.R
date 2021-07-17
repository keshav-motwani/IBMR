library(IBMR)

source("scripts/simulation_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "results/simulations_2"
dir.create(RESULT_PATH, recursive = TRUE)

defaults = list(
  q = 16,
  N = 10000,
  K = 4,
  p = 20,
  nonzero = 20,
  b = 2,
  rank = 1,
  batch_effect = 0
)

considered_values = list(
  K = c(1, 2, 4, 8, 16)
)

methods = c("IBMR", "IBMR_common_Gamma", "IBMR_no_Gamma", "glmnet_subset", "glmnet_split", "glmnet_relabel")

parameters = expand_parameters("fine_clean_low_dim_simulations", considered_values, defaults, 50, methods)

current_parameters = parameters[[ARRAY_ID]]
system.time({result = evaluate_parameters(current_parameters, generate_simulation_data_fine)})
saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))
