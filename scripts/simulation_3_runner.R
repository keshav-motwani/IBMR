library(IBMR)

source("scripts/simulation_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
DATA_PATH = "data/simulations"
RESULT_PATH = "results/simulations_3"
dir.create(RESULT_PATH, recursive = TRUE)

methods = c("IBMR_no_Gamma", "glmnet_subset", "glmnet_relabel")
methods = c(methods, paste0(methods, "_ORC_fine"))

defaults = list(
  sparsity = 0.2,
  N = 4800,
  category_mappings = 1,
  rank = 1,
  batch_effect = 0
)

considered_values = list(
  sparsity = c(0.05, 0.1, 0.2, 0.4, 0.8),
  N = c(2400, 4800, 9600, 19200)
)

parameters = expand_parameters("category_mappings_1", considered_values, defaults, 50, methods)

defaults$category_mappings = 2

parameters = c(parameters, expand_parameters("category_mappings_2", considered_values, defaults, 50, methods))

chunk_size = 3

data = list(
  X_star = readRDS(file.path(DATA_PATH, "hao_X_500.rds")),
  glmnet_fit = readRDS(file.path(DATA_PATH, "hao_glmnet_fit_500.rds"))
)

for (i in 1:chunk_size) {

  PARAMETER_ID = (ARRAY_ID - 1) * chunk_size + i

  current_parameters = c(parameters[[PARAMETER_ID]], data)
  current_parameters$category_mappings = readRDS(file.path(DATA_PATH, paste0("hao_category_mappings_", current_parameters$category_mappings, ".rds")))

  system.time({result = evaluate_parameters(current_parameters, generate_simulation_data_from_real)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}
