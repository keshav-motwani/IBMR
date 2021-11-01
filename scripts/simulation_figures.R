library(tidyverse)
library(patchwork)

RESULT_PATH = "final_results/simulations_random_X_and_structured_Beta_final"
FIGURES_PATH = c("figures/", file.path(RESULT_PATH, "figures"))
sapply(FIGURES_PATH, function(path) dir.create(path, recursive = TRUE))

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
                   Beta_SSE = result$performance["Beta_SSE"],
                   Beta_FPR = result$performance["Beta_FPR"],
                   Beta_TPR = result$performance["Beta_TPR"],
                   KL_divergence = result$performance["KL_divergence"],
                   error = result$performance["error"])
      ))
    }
    do.call(rbind, data)
  }
})

result = do.call(rbind, results)

result$method_name = result$method
result$oracle = gsub("ORCNA", "observed", paste0("ORC", sapply(strsplit(result$method, "ORC"), `[`, 2)))
result$method = sapply(strsplit(result$method, "_ORC"), `[`, 1)

methods = c("IBMR", "IBMR_int", "IBMR_no_Gamma", "relabel", "subset")
methods = methods[methods %in% result$method]

oracles = c("observed", "ORC_fine")
oracles = oracles[oracles %in% result$oracle]

summary = result %>%
  pivot_longer(Beta_SSE:error, names_repair = "minimal", values_to = "result") %>%
  filter(method %in% methods) %>%
  mutate(value = factor(value)) %>%
  mutate(method = factor(method, levels = methods)) %>%
  mutate(oracle = factor(oracle, levels = oracles)) %>%
  mutate(group = paste0(experiment, method, oracle, name)) %>%
  mutate(experiment = factor(experiment, levels = c("N", "p", "s", "batch_effect")))

plasma_pal = rev(viridis::plasma(n = length(methods) + 2)[1:length(methods)])
names(plasma_pal) = methods

experiments = cbind(summary$run, summary$experiment)
experiments = experiments[!duplicated(experiments), , drop = FALSE]

summary = summary %>%
  group_by(run, experiment, value, method, oracle, name, method_name) %>%
  summarize(mean = mean(result), two_se = 2 * sd(result)/sqrt(n()))

experiments = cbind(summary$run, summary$name)
experiments = experiments[!duplicated(experiments), , drop = FALSE]

summary = summary %>%
  mutate(method = method_name) %>%
  filter(!(method %in% c("subset_ORC_fine", "relabel_ORC_fine")))

methods = c("IBMR", "IBMR_int", "IBMR_no_Gamma", "relabel", "subset", paste0(c("IBMR_int", "IBMR_no_Gamma"), "_ORC_fine"))
methods = methods[methods %in% summary$method]

summary = summary %>%
  mutate(method = factor(method, levels = methods))

levels(summary$method) = gsub("no_Gamma", "NG", levels(summary$method), fixed = TRUE)
levels(summary$method) = gsub("ORC_fine", "ORC", levels(summary$method), fixed = TRUE)
levels(summary$method) = gsub("_", "-", levels(summary$method), fixed = TRUE)
methods = levels(summary$method)

levels(summary$experiment) = gsub("batch_effect", "b", levels(summary$experiment), fixed = TRUE)

notorcmethods = methods[grep("ORC", methods, invert = TRUE)]
plasma_pal = rev(viridis::plasma(n = length(notorcmethods) + 2)[1:length(notorcmethods)])
names(plasma_pal) = notorcmethods
orc_pal = rep("gray", length(methods) - length(notorcmethods))
names(orc_pal) = setdiff(methods, notorcmethods)
plasma_pal = c(plasma_pal, orc_pal)

experiments = experiments[experiments[, 2] %in% c("error", "KL_divergence"), ]

names = c("error" = "Error rate", "KL_divergence" = "KL divergence")

for (glmnet in c(TRUE)) {

  for (path in FIGURES_PATH) {

    pdf(file = file.path(path, paste0("simulation_figures_line_", ifelse(glmnet, "with_glmnet", "no_glmnet"), "_1.pdf")), height = (2.5 * 2 + 1.5) * 1.1, width = length(unique(summary$experiment)) * 3 * 1.1)

    plots = list()

    for (i in 1:nrow(experiments)) {

      data = summary %>% filter(name == experiments[i, 2], run == experiments[i, 1])

      if (!glmnet) data = data %>% filter(!grepl("subset", method))

      levels = gtools::mixedsort(as.character(unique(data$value)))
      if ("int" %in% levels) levels = c("int", levels[which(levels != "int")])

      data$value = factor(data$value, levels = levels)

      for (level in levels(data$experiment)) {

        plots = c(plots, list(ggplot(
          data %>% filter(experiment == level),
          aes(
            x = value,
            y = mean,
            color = method,
            group = method,
            linetype = method,
            ymin = mean - two_se,
            ymax = mean + two_se
          )
        ) +
          geom_point(size = 1) +
          geom_line() +
          geom_errorbar(width = 0.05, linetype = "solid", show.legend = FALSE) +
          theme_bw(base_size = 16) +
          guides(color = guide_legend(nrow = 2)) +
          theme(strip.background = element_blank(), strip.placement = "outside") +
          scale_color_manual(values = plasma_pal[levels(droplevels(data$method))]) +
          theme(legend.position = "bottom", legend.key.width = grid::unit(5, "lines")) +
          xlab(level) +
          ylab(names[experiments[i, 2]]) +
          labs(subtitle = NULL, color = "Method", linetype = "Method"))) # experiments[i, 1])))

      }

    }

    for (i in setdiff(1:8, c(1, 5))) {
      plots[[i]] = plots[[i]] + ylab(NULL)
    }
    for (i in 1:4) {
      plots[[i]] = plots[[i]] + xlab(NULL)
    }

    plot = wrap_plots(plots, ncol = 4) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")

    print(plot)

    dev.off()

  }

}

write.csv(result, file.path(RESULT_PATH, "figures", "simulation_results.csv"))
write.csv(summary, file.path(RESULT_PATH, "figures", "simulation_result_summary.csv"))
