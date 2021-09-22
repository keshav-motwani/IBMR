library(IBMR)

source("scripts/simulation_setup.R")
source("scripts/application_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "results/applications"
dir.create(RESULT_PATH, recursive = TRUE)

methods = c("IBMR_no_Gamma", "subset", "relabel")

defaults = list(
  cache_path = "../AnnotatedPBMC/data",
  split_index = 1,
  n_genes = 1000,
  n_sample = 5000
)

considered_values = list(
  split_index = 1:72
)

n_sample_sequence = c(5000, 1250, 2500, 10000)
n_genes_sequence = c(1000, 250, 500, 2000)

parameters = list()

for (n_sample in n_sample_sequence) {
  defaults$n_sample = n_sample
  parameters = c(parameters, expand_parameters(paste0("n_sample = ", defaults$n_sample, "; n_genes = ", defaults$n_genes), considered_values, defaults, 5, methods))
}

defaults$n_sample = 5000
for (n_genes in n_genes_sequence[-1]) {
  defaults$n_genes = n_genes
  parameters = c(parameters, expand_parameters(paste0("n_sample = ", defaults$n_sample, "; n_genes = ", defaults$n_genes), considered_values, defaults, 5, methods))
}

chunk_size = 3

for (i in 1:chunk_size) {

  PARAMETER_ID = (ARRAY_ID - 1) * chunk_size + i
  current_parameters = parameters[[PARAMETER_ID]]
  system.time({result = evaluate_parameters(current_parameters, prepare_real_data_application)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

  rm(result)
  print(gc())

}
