library(tidyverse)
library(patchwork)

RESULT_PATH = "results/application_R1"
FIGURES_PATH = c("figures/", file.path(RESULT_PATH, "figures"))
sapply(FIGURES_PATH, function(path) dir.create(path, recursive = TRUE))

files = unlist(lapply(RESULT_PATH, function(x) list.files(x, full.names = TRUE)))
files = files[grepl("rds", files)]

results = lapply(files, function(x) {
  if (file.info(x)$size > 0) {
    print(x)
    results = readRDS(x)
    data = list()
    for (result in results) {
      data = c(data, list(
        data.frame(
          run = result$parameters$run,
          experiment = result$parameters$experiment,
          n_genes = result$parameters$n_genes,
          n_sample = result$parameters$n_sample,
          value = result$parameters[[result$parameters$experiment]],
          method = result$parameters$method,
          replicate = result$parameters$replicate,
          error = result$performance["error"],
          nll = result$performance["nll"],
          time = log10(result$performance[["time"]])
        )
      ))
    }
    do.call(rbind, data)
  }
})

result = do.call(rbind, results)

methods = c("IBMR",
            "IBMR_int",
            "IBMR_common_Gamma",
            "IBMR_no_Gamma",
            "relabel",
            "subset", "Seurat", "SingleR")
methods = methods[methods %in% result$method]

dataset_names = read.csv(file.path("data/", "table_1.csv"))$dataset

splits = expand.grid(
  setdiff(dataset_names, "hao_2020"),
  setdiff(dataset_names, "hao_2020"),
  stringsAsFactors = FALSE
)
colnames(splits) = c("validation", "test")
splits = splits[splits[, 1] != splits[, 2],]
rownames(splits) = 1:nrow(splits)
splits$index = as.character(1:nrow(splits))

result$value = as.character(result$value)

result = left_join(result, splits, by = c(value = "index"))
result = result[result$method %in% methods,]
result$method = factor(result$method, levels = methods)
result$validation = factor(result$validation, levels = intersect(dataset_names, unique(result$validation)))
result$test = factor(result$test, levels = intersect(dataset_names, unique(result$test)))

levels(result$method) = gsub("no_Gamma", "NG", levels(result$method), fixed = TRUE)
levels(result$method) = gsub("ORC_fine", "ORC", levels(result$method), fixed = TRUE)
levels(result$method) = gsub("_", "-", levels(result$method), fixed = TRUE)
methods = levels(result$method)

plasma_pal = rep("black", length(methods)) # rev(viridis::plasma(n = length(methods) + 2)[1:length(methods)])
names(plasma_pal) = methods

result = result %>% group_by(run, value, method, n_sample, n_genes, validation, test) %>% summarize(
  mean_error = mean(error),
  mean_nll = mean(nll),
  mean_time = mean(time),
  se_error = sd(error) / sqrt(n()),
  se_nll = sd(nll) / sqrt(n()),
  se_time = sd(time) / sqrt(n())
)

averaged_result = result %>%
  ungroup() %>%
  group_by(run, method, n_sample, n_genes, test) %>%
  summarize(
    se_error = sd(mean_error) / sqrt(n()),
    se_nll = sd(mean_nll) / sqrt(n()),
    se_time = sd(mean_time) / sqrt(n()),
    mean_error = mean(mean_error),
    mean_nll = mean(mean_nll),
    mean_time = mean(mean_time)
  ) %>%
  mutate(validation = test)

ylabs = c("error" = "Error rate", "nll" = "Negative log-likelihood", "time" = "log10(Runtime (s))")

for (path in FIGURES_PATH) {
  for (value in c("error", "time", "nll")) {

    colors = plasma_pal
    if (value == "nll") {
      averaged_result = averaged_result %>% filter(!(method %in% c("Seurat", "SingleR"))) %>% mutate(method = factor(method, levels = setdiff(methods, c("Seurat", "SingleR"))))
      colors = colors[setdiff(methods, c("SingleR", "Seurat"))]
    }

    pdf(file = file.path(path,
                         paste0(
                           "application_figures_", value, "_1.pdf"
                         )),
        height = 9.2,
        width =  13.2 * 0.9)

    y_mean = paste0("mean_", value)
    y_se = paste0("se_", value)

    p = ggplot(
      averaged_result[averaged_result$n_sample == 10000,],
      aes(
        x = n_genes,
        y = .data[[y_mean]],
        color = method,
        group = method,
        shape = method,
        linetype = method,
        ymin = .data[[y_mean]] - .data[[y_se]],
        ymax = .data[[y_mean]] + .data[[y_se]]
      )
    ) +
      geom_point(size = 1) +
      geom_line() +
      geom_errorbar(width = 0,
                    linetype = "solid",
                    show.legend = FALSE) +
      facet_wrap(test ~ ., scales = "free", nrow = 3) +
      theme_bw(base_size = 16) +
      theme(strip.background = element_blank(),
            strip.placement = "outside") +
      theme(legend.position = "bottom",
            legend.key.width = grid::unit(5, "lines")) +
      scale_color_manual(values = colors) +
      xlab("# of genes") +
      ylab(ylabs[value]) +
      labs(color = "Method", shape = "Method", linetype = "Method")

    print(p)

    p = ggplot(
      averaged_result[averaged_result$n_genes == 1000,],
      aes(
        x = n_sample,
        y = .data[[y_mean]],
        color = method,
        group = method,
        shape = method,
        linetype = method,
        ymin = .data[[y_mean]] - .data[[y_se]],
        ymax = .data[[y_mean]] + .data[[y_se]]
      )
    ) +
      geom_point(size = 1) +
      geom_line() +
      geom_errorbar(width = 0,
                    linetype = "solid",
                    show.legend = FALSE) +
      facet_wrap(test ~ ., scales = "free", nrow = 3) +
      theme_bw(base_size = 16) +
      theme(strip.background = element_blank(),
            strip.placement = "outside") +
      theme(legend.position = "bottom",
            legend.key.width = grid::unit(5, "lines")) +
      scale_color_manual(values = colors) +
      xlab("# of cells per dataset") +
      ylab(ylabs[value]) +
      labs(color = "Method", shape = "Method", linetype = "Method")

    print(p)

    dev.off()

  }

  for (value in c("error", "time", "nll")) {
    
    colors = plasma_pal
    if (value == "nll") {
      result = result %>% filter(!(method %in% c("Seurat", "SingleR"))) %>% mutate(method = factor(method, levels = setdiff(methods, c("Seurat", "SingleR"))))
      colors = colors[setdiff(methods, c("SingleR", "Seurat"))]
    }

    pdf(
      file = file.path(path,
                       paste0(
                         "application_figures_", value, "_2.pdf"
                       )),
      height = 2 * length(unique(result$validation)),
      width =  2 * length(unique(result$test))
    )

    y_mean = paste0("mean_", value)
    y_se = paste0("se_", value)

    p = ggplot(
      result[result$n_sample == 10000,],
      aes(
        x = n_genes,
        y = .data[[y_mean]],
        color = method,
        group = method,
        shape = method,
        linetype = method,
        ymin = .data[[y_mean]] - .data[[y_se]],
        ymax = .data[[y_mean]] + .data[[y_se]]
      )
    ) +
      geom_point(size = 1) +
      geom_line() +
      geom_errorbar(width = 0,
                    linetype = "solid",
                    show.legend = FALSE) +
      facet_grid(test ~ validation, scales = "free_y") +
      theme_bw(base_size = 16) +
      theme(strip.background = element_blank(),
            strip.placement = "outside") +
      theme(legend.position = "bottom",
            legend.key.width = grid::unit(5, "lines")) +
      scale_color_manual(values = colors) +
      xlab("# of genes") +
      ylab(ylabs[value]) +
      labs(color = "Method", shape = "Method", linetype = "Method")

    print(p)

    p = ggplot(
      result[result$n_genes == 1000,],
      aes(
        x = n_sample,
        y = .data[[y_mean]],
        color = method,
        group = method,
        shape = method,
        linetype = method,
        ymin = .data[[y_mean]] - .data[[y_se]],
        ymax = .data[[y_mean]] + .data[[y_se]]
      )
    ) +
      geom_point(size = 1) +
      geom_line() +
      geom_errorbar(width = 0,
                    linetype = "solid",
                    show.legend = FALSE) +
      facet_grid(test ~ validation, scales = "free_y") +
      theme_bw(base_size = 16) +
      theme(strip.background = element_blank(),
            strip.placement = "outside") +
      theme(legend.position = "bottom",
            legend.key.width = grid::unit(5, "lines")) +
      scale_color_manual(values = colors) +
      xlab("# of cells per dataset") +
      ylab(ylabs[value]) +
      labs(color = "Method", shape = "Method", linetype = "Method")

    print(p)

    dev.off()

  }

}
