library(tidyverse)
library(patchwork)

RESULT_PATH = "results/applications"
dir.create(file.path(RESULT_PATH, "figures"))

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, function(x) {
  if (file.info(x)$size > 0) {
    print(x)
    results = readRDS(x)
    data = list()
    for (result in results) {
      data = c(data, list(
        data.frame(run = result$parameters$run,
                   experiment = result$parameters$experiment,
                   value = result$parameters[[result$parameters$experiment]],
                   method = result$parameters$method,
                   replicate = result$parameters$replicate,
                   error = result$performance["error"],
                   nll = result$performance["nll"])
      ))
    }
    do.call(rbind, data)
  }
})

result = do.call(rbind, results)

methods = c("IBMR", "IBMR_int", "IBMR_common_Gamma", "IBMR_no_Gamma", "glmnet_subset", "glmnet_relabel")
methods = methods[methods %in% result$method]

dataset_names = c("hao_2020", "haniffa_2021", "tsang_2021", "blish_2020", "kotliarov_2020", "10x_sorted", "su_2020", "10x_pbmc_10k", "10x_pbmc_5k_v3", "ding_2019")
splits = expand.grid(dataset_names[-1], dataset_names[-1], stringsAsFactors = FALSE)
colnames(splits) = c("validation", "test")
splits = splits[splits[, 1] != splits[, 2], ]
rownames(splits) = 1:nrow(splits)
splits$index = as.character(1:nrow(splits))

result$value = as.character(result$value)

result = left_join(result, splits, by = c(value = "index"))

result$method = factor(result$method, levels = methods)

plasma_pal = viridis::plasma(n = length(methods) + 2)[1:length(methods)]
names(plasma_pal) = methods

pdf(file = file.path(RESULT_PATH, "figures", "application_figures_nll.pdf"), height = 3 * length(unique(result$validation)), width =  3 * length(unique(result$test)))

ggplot(
  result,
  aes(
    x = method,
    y = nll,
    group = replicate
  )
) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("Method") +
  ylab("- log likelihood") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()

pdf(file = file.path(RESULT_PATH, "figures", "application_figures_error.pdf"), height = 3 * length(unique(result$validation)), width =  3 * length(unique(result$test)))

ggplot(
  result,
  aes(
    x = method,
    y = error,
    group = replicate
  )
) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_classic() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("Method") +
  ylab("Error rate") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
