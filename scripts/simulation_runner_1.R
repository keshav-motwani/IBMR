library(IBMR)

source("scripts/simulation_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "results/simulations_1"
dir.create(RESULT_PATH, recursive = TRUE)

defaults = list(
  q = 16,
  N = 2000,
  K = 4,
  p = 500,
  nonzero = 75,
  b = 2,
  r = 50
)

considered_values = list(
  # q = c(2, 4, 8, 16, 32),
  # K = c(1, 2, 4, 8, 16),
  nonzero = c(20, 50, 75, 100, 125) #,
  # b = c(0.25, 0.5, 1, 1.5, 2)
)

methods = c("IBMR", "IBMR_common_Gamma", "IBMR_no_Gamma", "elastic_net")

parameters = expand_parameters("fine_clean_simulations", considered_values, defaults, 50, methods)

current_parameters = parameters[[ARRAY_ID]]
system.time({result = evaluate_parameters(current_parameters, generate_simulation_data_fine_clean)})
saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))
