library(IBMR)

source("scripts/simulation_setup.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "final_results/simulations_random_X_and_structured_Beta_updated_again_more_fine_more_nonzero_more_fine_sigma2_catmap62"
dir.create(RESULT_PATH, recursive = TRUE)

methods = c("IBMR_int", "IBMR_no_Gamma", "subset", "relabel")
methods = c(methods, paste0(methods, "_ORC_fine")) # , paste0(methods, "_ORC_clean"), paste0(methods, "_ORC_fine_clean"))

defaults = list(
  # 6 coarse categories, each split into 2 subcategories. 6 datasets, 3 observed at coarse resolution, 3 at fine resolution for all categories.
  category_mappings = simulate_category_mappings(2, c(6, 2), c(replicate(4, c(rep(1, 5), 2), simplify = FALSE), replicate(2, rep(2, 6), simplify = FALSE))),
  N = 1200,
  p = 500,
  nonzero = 100,
  s = 40, # number of genes which have shared coefficients within coarse categories
  sigma = 2, # sd of normal distribution for sampling coefficients
  rank = "int",
  batch_effect = 0.1 # norm of batch effect over norm of true predictors
)

considered_values = list(
  p = c(100, 200, 500, 1000, 2000),
  s = c(0, 20, 40, 60, 80),
  N = c(300, 600, 1200, 2400, 4800),
  batch_effect = c(0.025, 0.05, 0.1, 0.2, 0.4)
)

parameters = expand_parameters("random_X_and_structured_Beta", considered_values, defaults, 50, methods)

chunk_size = 4

for (i in 1:chunk_size) {

  PARAMETER_ID = (ARRAY_ID - 1) * chunk_size + i
  current_parameters = parameters[[PARAMETER_ID]]
  system.time({result = evaluate_parameters(current_parameters, generate_data_random_X_and_structured_Beta)})
  saveRDS(result, file.path(RESULT_PATH, paste0(gsub("___|__", "_", gsub(" |;|=|,", "_", current_parameters$run)), "_", current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}
