library(IBMR)

source("scripts/application_setup_1.R")

ARRAY_ID = 4 # as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
DATA_PATH = "data/application_1/"
RESULT_PATH = "results/application_1"
dir.create(RESULT_PATH, recursive = TRUE)

data = readRDS(file.path(DATA_PATH, "data.rds"))

for (set in c("train", "validation", "test")) {

for (k in 1:length(data[[set]]$X_list)) {
  indices = sample(1:nrow(data[[set]]$X_list[[k]]), 100)

    for (l in grep("list", names(data[[set]]), value = TRUE)) {
      if (grepl("Y", l)) {
        data[[set]][[l]][[k]] = data[[set]][[l]][[k]][indices]
      } else {
        data[[set]][[l]][[k]] = data[[set]][[l]][[k]][indices,]
      }

      if (grepl("X", l)) {
        data[[set]][[l]][[k]] = data[[set]][[l]][[k]][, 1:50]

      }

    }

  }

}

methods = c("IBMR", "IBMR_int", "IBMR_common_Gamma", "IBMR_no_Gamma",
            "IBMR_ORC_FINE", "IBMR_int_ORC_FINE", "IBMR_common_Gamma_ORC_FINE", "IBMR_no_Gamma_ORC_FINE")
debugonce(IBMR_no_Gamma)
system.time({result = evaluate_parameters(data, methods[ARRAY_ID])})
saveRDS(result, file.path(RESULT_PATH, paste0(method, ".rds")))
