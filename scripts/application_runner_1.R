library(IBMR)

source("scripts/application_setup_1.R")

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
DATA_PATH = "data/application_1/"
RESULT_PATH = "results/application_1"
dir.create(RESULT_PATH, recursive = TRUE)

data = readRDS(file.path(DATA_PATH, "data.rds"))

methods = c("IBMR", "IBMR_int", "IBMR_common_Gamma", "IBMR_no_Gamma",
            "IBMR_ORC_FINE", "IBMR_int_ORC_FINE", "IBMR_common_Gamma_ORC_FINE", "IBMR_no_Gamma_ORC_FINE")

system.time({result = evaluate_parameters(data, methods[ARRAY_ID])})
saveRDS(result, file.path(RESULT_PATH, paste0(method, ".rds")))
