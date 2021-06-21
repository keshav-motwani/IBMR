library(tidyverse)
library(patchwork)

RESULT_PATH = "results/simulations_1"
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

methods = c("IBMR", "IBMR_common_Gamma", "IBMR_no_Gamma",
            "group_lasso")

summary = result %>%
  pivot_longer(Beta_SSE:best_case_error, names_repair = "minimal", values_to = "result") %>%
  # filter(method != "ORACLE" | name == "error" | name == "best_case_error") %>%
  filter(method %in% methods) %>%
  mutate(value = factor(value)) %>%
  mutate(method = factor(method, levels = methods)) %>%
  mutate(group = paste0(experiment, method, name))

plasma_pal = viridis::plasma(n = length(methods) + 2)[1:length(methods)]
names(plasma_pal) = methods

plots = list()

experiments = cbind(summary$run, summary$experiment)
experiments = experiments[!duplicated(experiments), , drop = FALSE]

for (i in 1:nrow(experiments)) {
  plots = c(
    plots,
    list(
      ggplot(
        summary %>% filter(experiment == experiments[i, 2], run == experiments[i, 1]),
        aes(
          x = as.factor(value),
          y = result,
          fill = method
        )
      ) +
        geom_boxplot(lwd=.5, outlier.size = .4) + 
        # geom_point(size = 2, position=position_dodge(0.2)) +
        # geom_line(position=position_dodge(0.2)) +
        # geom_errorbar(width = 0.2, position=position_dodge(0.2), linetype = "solid", show.legend = FALSE) +
        facet_wrap( ~ name, scales = "free", nrow = 1) +
        theme_classic() +
        theme(strip.background = element_blank(), strip.placement = "outside") +
        scale_fill_manual(values = plasma_pal) +
        theme(legend.position = "bottom") +
        xlab(experiments[i, 2]) +
        ylab(NULL) +
        labs(subtitle = experiments[i, 1]) +
        scale_shape_manual(values=1:length(methods))
    )
  )
}

wrap_plots(plots, nrow = length(plots)) + plot_layout(guides = "collect") &
  theme(legend.position='bottom')

ggsave(file.path(RESULT_PATH, "figures", "simulation_figures_box.pdf"), height = 4.2 * nrow(experiments), width = 32)

summary = summary %>%
  group_by(run, experiment, value, method, name) %>%
  summarize(mean = mean(result), two_se = 2 * sd(result)/sqrt(n()))

plots = list()

experiments = cbind(summary$run, summary$experiment)
experiments = experiments[!duplicated(experiments), , drop = FALSE]

for (i in 1:nrow(experiments)) {
  plots = c(
    plots,
    list(
      ggplot(
        summary %>% filter(experiment == experiments[i, 2], run == experiments[i, 1]),
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
        facet_wrap( ~ name, scales = "free", nrow = 1) +
        theme_classic() +
        theme(strip.background = element_blank(), strip.placement = "outside") +
        scale_color_manual(values = plasma_pal) +
        theme(legend.position = "bottom") +
        xlab(experiments[i, 2]) +
        ylab(NULL) +
        labs(subtitle = experiments[i, 1]) +
        scale_shape_manual(values=1:length(methods))
    )
  )
}

wrap_plots(plots, nrow = length(plots)) + plot_layout(guides = "collect") &
  theme(legend.position='bottom')

ggsave(file.path(RESULT_PATH, "figures", "simulation_figures_line.pdf"), height = 4.2 * nrow(experiments), width = 32)

write.csv(result, file.path(RESULT_PATH, "figures", "simulation_results.csv"))
write.csv(summary, file.path(RESULT_PATH, "figures", "simulation_result_summary.csv"))
