print(Sys.getenv('SLURM_JOB_ID'))

library(IBMR)

source("scripts/simulation_setup.R")
source("scripts/application_setup.R")

PARAMETER_ID = as.numeric(commandArgs(trailingOnly=TRUE)[1])
METHOD = commandArgs(trailingOnly=TRUE)[2]

RESULT_PATH = "results/application_R1"
dir.create(RESULT_PATH, recursive = TRUE)

if (METHOD == "Seurat") {
  methods = "Seurat"
} else {
  methods = c("IBMR_no_Gamma", "subset", "relabel", "IBMR_int", "SingleR")
}

defaults = list(
  cache_path = "../AnnotatedPBMC/data",
  split_index = 1,
  n_genes = 1000,
  n_sample = 10000
)

considered_values = list(
  split_index = 1:72
)

n_sample_sequence = c(10000, 5000, 15000, 20000)
n_genes_sequence = c(1000, 250, 500, 2000)

parameters = list()

for (n_sample in n_sample_sequence) {
  defaults$n_sample = n_sample
  parameters = c(parameters, expand_parameters(paste0("n_sample = ", defaults$n_sample, "; n_genes = ", defaults$n_genes), considered_values, defaults, 5, methods))
}

defaults$n_sample = 10000
for (n_genes in n_genes_sequence[-1]) {
  defaults$n_genes = n_genes
  parameters = c(parameters, expand_parameters(paste0("n_sample = ", defaults$n_sample, "; n_genes = ", defaults$n_genes), considered_values, defaults, 5, methods))
}

parameters = parameters

current_parameters = parameters[[PARAMETER_ID]]

print(PARAMETER_ID)
print(current_parameters)

path = file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds"))

if (!file.exists(path)) {

  system.time({result = evaluate_parameters(current_parameters, prepare_real_data_application)})
  saveRDS(result, path)

} else {

  print("Already computed.")

}
