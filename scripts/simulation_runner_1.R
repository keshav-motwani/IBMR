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
  r = 50,
  rank = 20,
  batch_effect = 0.1
)

considered_values = list(
  batch_effect = c(0.025, 0.05, 0.1, 0.2, 0.4),
  rank = c(5, 10, 20, 40, 80),
  K = c(1, 2, 4, 8, 16),
  nonzero = c(20, 50, 75, 100, 125)
)

methods = c("IBMR", "IBMR_common_Gamma", "IBMR_int", "IBMR_no_Gamma",
            "IBMR_ORC_clean", "IBMR_common_Gamma_ORC_clean", "IBMR_int_ORC_clean", "IBMR_no_Gamma_ORC_clean")

parameters = expand_parameters("fine_simulations", considered_values, defaults, 50, methods)

for (i in 1:length(methods)) {

  PARAMETER_ID = (ARRAY_ID - 1) * length(methods) + i

  current_parameters = parameters[[PARAMETER_ID]]
  system.time({result = evaluate_parameters(current_parameters, generate_simulation_data_fine)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}
