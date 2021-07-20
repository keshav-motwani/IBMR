library(tidyverse)
library(patchwork)

RESULT_PATH = "results/simulations_1_new"
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
                   Beta_SSE = result$performance["Beta_SSE"],
                   Beta_FPR = result$performance["Beta_FPR"],
                   Beta_TPR = result$performance["Beta_TPR"],
                   KL_divergence = result$performance["KL_divergence"],
                   hellinger_distance = result$performance["hellinger_distance"],
                   error = result$performance["error"],
                   best_case_Beta_SSE = result$best_case_performance["Beta_SSE"],
                   best_case_KL_divergence = result$best_case_performance["KL_divergence"],
                   best_case_hellinger_distance = result$best_case_performance["hellinger_distance"],
                   best_case_error = result$best_case_performance["error"])
      ))
    }
    do.call(rbind, data)
  }
})

result = do.call(rbind, results)

result$oracle = gsub("ORCNA", "observed", paste0("ORC", sapply(strsplit(result$method, "ORC"), `[`, 2)))
result$method = sapply(strsplit(result$method, "_ORC"), `[`, 1)

methods = c("IBMR", "IBMR_int", "IBMR_common_Gamma", "IBMR_no_Gamma", "glmnet_subset", "glmnet_relabel")
methods = methods[methods %in% result$method]

oracles = c("observed", "ORC_fine", "ORC_clean", "ORC_fine_clean")
oracles = oracles[oracles %in% result$oracle]

summary = result %>%
  pivot_longer(Beta_SSE:best_case_error, names_repair = "minimal", values_to = "result") %>%
  # filter(method != "ORACLE" | name == "error" | name == "best_case_error") %>%
  filter(method %in% methods) %>%
  mutate(value = factor(value)) %>%
  mutate(method = factor(method, levels = methods)) %>%
  mutate(oracle = factor(oracle, levels = oracles)) %>%
  mutate(group = paste0(experiment, method, oracle, name))


plasma_pal = viridis::plasma(n = length(methods) + 2)[1:length(methods)]
names(plasma_pal) = methods

experiments = cbind(summary$run, summary$experiment)
experiments = experiments[!duplicated(experiments), , drop = FALSE]

for (glmnet in c(TRUE, FALSE)) {

  pdf(file = file.path(RESULT_PATH, "figures", paste0("simulation_figures_box_", ifelse(glmnet, "with_glmnet", "no_glmnet"), ".pdf")), height = 3 * length(unique(summary$name)), width = length(unique(summary$oracle)) * 4)

  for (i in 1:nrow(experiments)) {

    data = summary %>% filter(experiment == experiments[i, 2], run == experiments[i, 1])

    if (!glmnet) data = data %>% filter(!grepl("glmnet", method))

    levels = gtools::mixedsort(as.character(unique(data$value)))
    if ("int" %in% levels) levels = c("int", levels[which(levels != "int")])

    data$value = factor(data$value, levels = levels)

    plot = ggplot(
      data,
      aes(
        x = value,
        y = result,
        fill = method
      )
    ) +
      geom_boxplot(lwd=.5, outlier.size = .4) +
      # geom_point(size = 2, position=position_dodge(0.2)) +
      # geom_line(position=position_dodge(0.2)) +
      # geom_errorbar(width = 0.2, position=position_dodge(0.2), linetype = "solid", show.legend = FALSE) +
      facet_grid(name ~ oracle, scales = "free_y") +
      theme_classic() +
      theme(strip.background = element_blank(), strip.placement = "outside") +
      scale_fill_manual(values = plasma_pal) +
      theme(legend.position = "bottom") +
      xlab(experiments[i, 2]) +
      ylab(NULL) +
      labs(subtitle = experiments[i, 1]) +
      scale_shape_manual(values=1:length(methods))

    print(plot)

  }

  dev.off()

}

summary = summary %>%
  group_by(run, experiment, value, method, oracle, name) %>%
  summarize(mean = mean(result), two_se = 2 * sd(result)/sqrt(n()))

experiments = cbind(summary$run, summary$experiment)
experiments = experiments[!duplicated(experiments), , drop = FALSE]

for (glmnet in c(TRUE, FALSE)) {

  pdf(file = file.path(RESULT_PATH, "figures", paste0("simulation_figures_line_", ifelse(glmnet, "with_glmnet", "no_glmnet"), ".pdf")), height = 3 * length(unique(summary$name)), width = length(unique(summary$oracle)) * 4)

  for (i in 1:nrow(experiments)) {

    data = summary %>% filter(experiment == experiments[i, 2], run == experiments[i, 1])

    if (!glmnet) data = data %>% filter(!grepl("glmnet", method))

    levels = gtools::mixedsort(as.character(unique(data$value)))
    if ("int" %in% levels) levels = c("int", levels[which(levels != "int")])

    data$value = factor(data$value, levels = levels)

    plot = ggplot(
      data,
      aes(
        x = value,
        y = mean,
        color = method,
        group = method,
        linetype = method,
        shape = method,
        ymin = mean - two_se,
        ymax = mean + two_se
      )
    ) +
      geom_point(size = 2, position=position_dodge(0.2)) +
      geom_line(position=position_dodge(0.2)) +
      geom_errorbar(width = 0.2, position=position_dodge(0.2), linetype = "solid", show.legend = FALSE) +
      facet_grid(name ~ oracle, scales = "free_y") +
      theme_classic() +
      theme(strip.background = element_blank(), strip.placement = "outside") +
      scale_color_manual(values = plasma_pal) +
      theme(legend.position = "bottom") +
      xlab(experiments[i, 2]) +
      ylab(NULL) +
      labs(subtitle = experiments[i, 1]) +
      scale_shape_manual(values=1:length(methods))

    print(plot)

  }

  dev.off()

}

write.csv(result, file.path(RESULT_PATH, "figures", "simulation_results.csv"))
write.csv(summary, file.path(RESULT_PATH, "figures", "simulation_result_summary.csv"))
