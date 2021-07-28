library(IBMR)

source("scripts/simulation_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
DATA_PATH = "data/simulations"
RESULT_PATH = "results/simulations_3"
dir.create(RESULT_PATH, recursive = TRUE)

methods = c("IBMR", "IBMR_int", "IBMR_common_Gamma", "IBMR_no_Gamma", "glmnet_subset", "glmnet_relabel")
methods = c(methods, paste0(methods, "_ORC_clean"), paste0(methods, "_ORC_fine"), paste0(methods, "_ORC_fine_clean"))

defaults = list(
  sparsity = 0.84,
  N = 2400,
  rank = 20,
  batch_effect = 0.1
)

considered_values = list(
  batch_effect = c(0.025, 0.05, 0.1, 0.2, 0.4),
  sparsity = c(0.2, 0.5, 0.84, 0.87),
  N = c(600, 1200, 2400, 4800),
  rank = c("int", 1, 5, 20)
)

parameters = expand_parameters("low_rank_batch_effect", considered_values, defaults, 50, methods)

defaults$rank = "int"
considered_values$rank = NULL

parameters = c(parameters, expand_parameters("intercept_batch_effect", considered_values, defaults, 50, methods))

data = list(
  category_mappings = readRDS(file.path(DATA_PATH, "hao_category_mappings.rds")),
  X_star = readRDS(file.path(DATA_PATH, "hao_X_500.rds")),
  glmnet_fit = readRDS(file.path(DATA_PATH, "hao_glmnet_fit_500.rds"))
)

chunk_size = 24

for (i in 1:chunk_size) {

  PARAMETER_ID = (ARRAY_ID - 1) * chunk_size + i
  current_parameters = c(parameters[[PARAMETER_ID]], data)
  system.time({result = evaluate_parameters(current_parameters, generate_simulation_data_from_real)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}
