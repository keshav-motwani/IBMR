library(tidyverse)
library(patchwork)

RESULT_PATH = "results/applications_fullvaltest"
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

methods = c("IBMR", "IBMR_int", "IBMR_common_Gamma", "IBMR_no_Gamma", "subset", "relabel")
# methods = c("IBMR_no_Gamma", "subset", "relabel")
methods = methods[methods %in% result$method]

dataset_names = c("hao_2020", "haniffa_2021", "tsang_2021", "blish_2020", "kotliarov_2020", "10x_sorted", "su_2020", "10x_pbmc_10k", "10x_pbmc_5k_v3", "ding_2019")
splits = expand.grid(dataset_names[-1], dataset_names[-1], stringsAsFactors = FALSE)
colnames(splits) = c("validation", "test")
splits = splits[splits[, 1] != splits[, 2], ]
rownames(splits) = 1:nrow(splits)
splits$index = as.character(1:nrow(splits))

result$value = as.character(result$value)

result = left_join(result, splits, by = c(value = "index"))
result = result[result$method %in% methods, ]
result$method = factor(result$method, levels = methods)

plasma_pal = viridis::plasma(n = length(methods) + 2)[1:length(methods)]
names(plasma_pal) = methods

pdf(file = file.path(RESULT_PATH, "figures", "application_figures_nll.pdf"), height = 3 * length(unique(result$validation)), width =  3 * length(unique(result$test)))

for (run in gtools::mixedsort(unique(result$run))) {

p = ggplot(
  result[result$run == run, ],
  aes(
    x = method,
    y = nll,
    group = replicate
  )
) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("Method") +
  ylab("- log likelihood") +
  labs(title = run) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)

}

dev.off()

pdf(file = file.path(RESULT_PATH, "figures", "application_figures_error.pdf"), height = 3 * length(unique(result$validation)), width =  3 * length(unique(result$test)))

for (run in gtools::mixedsort(unique(result$run))) {

p = ggplot(
  result[result$run == run, ],
  aes(
    x = method,
    y = error,
    group = replicate
  )
) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("Method") +
  ylab("Error rate") +
  labs(title = run) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)

}

dev.off()

numbers = lapply(strsplit(result$run, "; "), function(x) as.numeric(sapply(strsplit(x, "= "), `[`, 2)))
result$n_sample = sapply(numbers, `[`, 1)
result$n_genes = sapply(numbers, `[`, 2)

result = result %>% group_by(run, value, method, n_sample, n_genes, validation, test) %>% summarize(mean_error = mean(error), mean_nll = mean(nll), two_se_error = 2 * sd(error)/sqrt(n()), two_se_nll = 2 * sd(nll)/sqrt(n()))
result2 = result %>% 
  ungroup() %>%
  group_by(run, method, n_sample, n_genes, test) %>% 
  summarize(mean_error = mean(mean_error), mean_nll = mean(mean_nll), two_se_error = 0, two_se_nll = 0) %>%
  mutate(validation = test)

result = rbind(result, result2)

pdf(file = file.path(RESULT_PATH, "figures", "application_figures_error_2.pdf"), height = 3 * length(unique(result$validation)), width =  3 * length(unique(result$test)))

p = ggplot(
  result[result$n_sample == 5000, ],
  aes(
    x = n_genes,
    y = mean_error,
    color = method,
    group = method,
    linetype = method,
    ymin = mean_error - two_se_error,
    ymax = mean_error + two_se_error
  )
) +
  geom_point(size = 2) +
  geom_line() + 
  geom_errorbar(width = 0.2, linetype = "solid", show.legend = FALSE) +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("# of genes") +
  ylab("Error rate") +
  labs(title = "# of genes, with 5000 cells per dataset") +
  scale_color_manual(values = plasma_pal)
print(p)

p = ggplot(
  result[result$n_genes == 1000, ],
  aes(
    x = n_sample,
    y = mean_error,
    color = method,
    group = method,
    linetype = method,
    ymin = mean_error - two_se_error,
    ymax = mean_error + two_se_error
  )
) +
  geom_point(size = 2) +
  geom_line() + 
  geom_errorbar(width = 0.2, linetype = "solid", show.legend = FALSE) +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("# of cells per dataset") +
  ylab("Error rate") +
  labs(title = "# of cells per dataset, with 1000 genes") +
  scale_color_manual(values = plasma_pal)

print(p)

dev.off()

pdf(file = file.path(RESULT_PATH, "figures", "application_figures_nll_2.pdf"), height = 3 * length(unique(result$validation)), width =  3 * length(unique(result$test)))

p = ggplot(
  result[result$n_sample == 5000, ],
  aes(
    x = n_genes,
    y = mean_nll,
    color = method,
    group = method,
    linetype = method,
    ymin = mean_nll - two_se_nll,
    ymax = mean_nll + two_se_nll
  )
) +
  geom_point(size = 2) +
  geom_line() + 
  geom_errorbar(width = 0.2, linetype = "solid", show.legend = FALSE) +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("# of genes") +
  ylab("nll rate") +
  labs(title = "# of genes, with 5000 cells per dataset") +
  scale_color_manual(values = plasma_pal)
print(p)

p = ggplot(
  result[result$n_genes == 1000, ],
  aes(
    x = n_sample,
    y = mean_nll,
    color = method,
    group = method,
    linetype = method,
    ymin = mean_nll - two_se_nll,
    ymax = mean_nll + two_se_nll
  )
) +
  geom_point(size = 2) +
  geom_line() + 
  geom_errorbar(width = 0.2, linetype = "solid", show.legend = FALSE) +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("# of cells per dataset") +
  ylab("nll rate") +
  labs(title = "# of cells per dataset, with 1000 genes") +
  scale_color_manual(values = plasma_pal)

print(p)

dev.off()

pdf(file = file.path(RESULT_PATH, "figures", "application_figures_error_3.pdf"), height = 1.3 * length(unique(result$validation)), width =  1.3 * length(unique(result$test)))

p = ggplot(
  result[result$n_sample == 5000, ],
  aes(
    x = as.factor(n_genes),
    y = mean_error,
    color = method,
    group = method,
    linetype = method,
    ymin = mean_error - two_se_error,
    ymax = mean_error + two_se_error
  )
) +
  geom_point(size = 1) +
  geom_line() + 
  geom_errorbar(width = 0.2, linetype = "solid", show.legend = FALSE) +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("# of genes") +
  ylab("Error rate") +
  labs(title = "# of genes, with 5000 cells per dataset") +
  scale_color_manual(values = plasma_pal) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)

p = ggplot(
  result[result$n_genes == 1000, ],
  aes(
    x = as.factor(n_sample),
    y = mean_error,
    color = method,
    group = method,
    linetype = method,
    ymin = mean_error - two_se_error,
    ymax = mean_error + two_se_error
  )
) +
  geom_point(size = 1) +
  geom_line() + 
  geom_errorbar(width = 0.2, linetype = "solid", show.legend = FALSE) +
  facet_grid(test ~ validation, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  scale_fill_manual(values = plasma_pal) +
  theme(legend.position = "bottom") +
  xlab("# of cells per dataset") +
  ylab("Error rate") +
  labs(title = "# of cells per dataset, with 1000 genes") +
  scale_color_manual(values = plasma_pal) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(p)

dev.off()

