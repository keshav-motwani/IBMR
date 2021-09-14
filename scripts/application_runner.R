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
  n_genes = 1000,
  n_sample = 5000
)

considered_values = list(
  split_index = 1:72
)

parameters = list()

for (n_genes in c(1000, 250, 500, 2000)) {
  defaults$n_genes = n_genes
  parameters = c(parameters, expand_parameters("n_genes; n_sample = 5000", considered_values, defaults, 5, methods))
}

defaults$n_genes = 1000

for (n_sample in c(5000, 1250, 2500, 10000)) {
  defaults$n_sample = n_sample
  parameters = c(parameters, expand_parameters("n_sample; n_genes = 1000 ", considered_values, defaults, 5, methods))
}

chunk_size = 3

for (i in 1:chunk_size) {

  PARAMETER_ID = (ARRAY_ID - 1) * chunk_size + i
  current_parameters = parameters[[PARAMETER_ID]]
  system.time({result = evaluate_parameters(current_parameters, prepare_real_data_application)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}
