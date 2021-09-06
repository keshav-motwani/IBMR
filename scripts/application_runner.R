library(IBMR)

source("scripts/simulation_setup.R")
source("scripts/application_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "results/applications"
dir.create(RESULT_PATH, recursive = TRUE)

methods = c("IBMR_no_Gamma", "glmnet_subset", "glmnet_relabel")

defaults = list(
  cache_path = "../AnnotatedPBMC/data",
  split_index = 1,
  n_sample = 5000
)

considered_values = list(
  split_index = 1:20
)

parameters = expand_parameters("application", considered_values, defaults, 5, methods)

chunk_size = 1

for (i in 1:chunk_size) {

  PARAMETER_ID = (ARRAY_ID - 1) * chunk_size + i
  current_parameters = parameters[[PARAMETER_ID]]
  system.time({result = evaluate_parameters(current_parameters, prepare_real_data_application)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}
