library(IBMR)
library(scran)
library(scater)
library(tidyverse)

source("scripts/application_setup.R")

CACHE_PATH = "../AnnotatedPBMC/data"
RESULT_PATH = "results/application_R1"
FIGURES_PATH = "results/application_R1/figures/"
dir.create(FIGURES_PATH, recursive = TRUE)

ding_2019 = prepare_ding_2019(CACHE_PATH, n_genes = 1000, sce = TRUE)
binning_function = ding_2019$binning_function
ding_2019 = ding_2019$sce

model = readRDS(file.path(RESULT_PATH, "n_sample_10000_n_genes_1000_split_index_65_IBMR_int_1.rds"))[[1]]
X_sd = model$fit$X_sd
model = model$fit$best_model
category_mapping = binning_function_to_category_mapping(binning_function)
Y = ding_2019$cell_type
P = predict_probabilities(model, list(as.matrix(t(logcounts(ding_2019)))))
C = predict_conditional_probabilities(P, list(Y), list(category_mapping))

predictions = list()
coarse_predictions = predict_categories(P, list(category_mapping))[[1]]
coarse_prediction = paste0("Coarse prediction\n(", round(100 * mean(Y == coarse_predictions), 1), "% accuracy)")
predictions[[coarse_prediction]] = coarse_predictions
predictions[["Fine prediction"]] = factor(predict_categories(P)[[1]], levels = colnames(model$Beta))
predictions[["Conditional prediction"]] = factor(predict_categories(C)[[1]], levels = colnames(model$Beta))

tables = list()
for (prediction in names(predictions)) {
  table = table(Y, predictions[[prediction]])
  tables[[prediction]] = table / rowSums(table)
}

data = list()
for (prediction in names(tables)) {
  result = reshape2::melt(tables[[prediction]])
  colnames(result) = c("true", "prediction", "percentage")
  result$type = prediction
  data[[prediction]] = result
}
data = do.call(rbind, data)
data$type = factor(data$type, levels = names(predictions))

library(ggplot2)

plot_data = data %>% filter(type %in% c(coarse_prediction, "Fine prediction"))
ggplot(plot_data, aes(x = prediction, y = true, fill = percentage * 100)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    breaks = c(0, 0.25, 0.5, 0.75, 1) * 100,
    limits = c(0, 1) * 100
  ) +
  facet_grid(~type, scales = "free", space = "free") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_blank()) +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(x = "Predicted category", y = "Observed label", fill = "% of observed") +
  guides(x = guide_axis(angle = 45)) +
  geom_point(data = plot_data %>% filter(percentage == 0), aes(x = prediction, y = true), size = 0.1)

ggsave(file.path(FIGURES_PATH, "application_heatmap_1.pdf"), height = 4.6, width = 13.2)

plot_data = data %>% filter(type %in% c("Conditional prediction"))
ggplot(plot_data, aes(x = prediction, y = true, fill = percentage * 100)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    breaks = c(0, 0.25, 0.5, 0.75, 1) * 100,
    limits = c(0, 1) * 100
  ) +
  facet_grid(~type, scales = "free", space = "free") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_blank()) +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  labs(x = "Predicted category", y = "Observed label", fill = "% of observed") +
  guides(x = guide_axis(angle = 45)) +
  geom_point(data = plot_data %>% filter(percentage == 0), aes(x = prediction, y = true), size = 0.1)

ggsave(file.path(FIGURES_PATH, "application_heatmap_2.pdf"), height = 4.6, width = 13.2)

Beta = model$Beta * X_sd
Beta = Beta[rowSums(Beta) != 0, ]
genes = t(apply(Beta, 2, function(x) rev(tail(names(sort((x))), 10))))
satija_genes = read.csv("data/marker_genes.csv")
satija_genes = satija_genes[!grepl("Prolif", satija_genes$Label), ]
rownames(satija_genes) = satija_genes$Label
satija_genes = apply(satija_genes, 1, function(x) unlist(strsplit(x["Markers"], ", ")))

for (category in colnames(Beta)) {
  genes[category, ] = ifelse(genes[category, ] %in% satija_genes[, category], paste0("\\textbf{", genes[category, ], "}"), genes[category, ])
}

xt = xtable::xtable(genes)
print(xt, sanitize.text = identity)
