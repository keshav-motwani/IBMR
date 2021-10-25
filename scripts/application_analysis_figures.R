library(Seurat)
library(IBMR)
library(scran)
library(scater)

source("scripts/application_setup.R")

CACHE_PATH = "../AnnotatedPBMC/data"
RESULT_PATH = "results/applications_updated_ten_not_split"

ding_2019 = prepare_ding_2019(CACHE_PATH, n_genes = 1000, sce = TRUE)
binning_function = ding_2019$binning_function
ding_2019 = as.Seurat(ding_2019$sce)

ding_2019_list = SplitObject(ding_2019, split.by = "dataset")
anchors = FindIntegrationAnchors(object.list = ding_2019_list, anchor.features = rownames(ding_2019_list[[1]]))
ding_2019 = IntegrateData(anchorset = anchors)

DefaultAssay(ding_2019) = "integrated"

ding_2019 = ScaleData(ding_2019, verbose = FALSE)
ding_2019 = RunPCA(ding_2019, npcs = 30, verbose = FALSE)
ding_2019 = RunUMAP(ding_2019, reduction = "pca", dims = 1:30)

ding_2019 = as.SingleCellExperiment(ding_2019)

model = readRDS("results/applications_updated_ten_not_split/n_sample_10000_n_genes_1000_split_index_65_IBMR_no_Gamma_1.rds")[[1]]$fit$best_model
category_mapping = binning_function_to_category_mapping(binning_function)
Y = ding_2019$cell_type
P = predict_probabilities(model, list(as.matrix(t(logcounts(ding_2019)))))
C = predict_conditional_probabilities(P, list(Y), list(category_mapping))

predictions = list()
predictions[["coarse"]] = predict_categories(P, list(category_mapping))[[1]]
predictions[["fine"]] = predict_categories(P)[[1]]
predictions[["conditional"]] = predict_categories(C)[[1]]

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

ggplot(data, aes(x = prediction, y = true, fill = percentage)) +
	geom_tile() +
	scale_fill_gradient(
        low = "white",
        high = "firebrick",
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        limits = c(0, 1)
    ) +
    facet_grid(~type, scales = "free", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.background = element_blank()) +
    theme(strip.background = element_blank(), strip.placement = "outside") +
    labs(x = "Predicted category", y = "Observed label", fill = "% of observed")

 ggsave(file.path(RESULT_PATH, "figures", "heatmap.pdf"), height = 4, width = 15)


data = as.data.frame(reducedDim(ding_2019, "UMAP"))
predictions[["observed"]] = Y
data = cbind(data, predictions)

predictions = predictions[c("observed", "coarse", "fine", "conditional")]
plots = list()
for (type in names(predictions)) {
	plot = ggplot(data = data, aes_string(x = "UMAP_1", y = "UMAP_2", color = type)) +
		geom_point(size = 0.5, alpha = 0.5) +
		theme_bw() +
		theme(legend.position = "bottom") +
		labs(subtitle = type, color = NULL) +
		scale_color_manual(values = c(ggsci::pal_d3("category10")(10), ggsci::pal_d3("category20")(20))) +
		guides(color = guide_legend(override.aes = list(size = 1, alpha = 1), ncol=3 ) )
	plots[[type]] = plot
}
patchwork::wrap_plots(plots, nrow = 2)
 ggsave(file.path(RESULT_PATH, "figures", "umap.pdf"), height = 14, width = 10)


genes = t(apply(model$Beta, 2, function(x) rev(tail(names(sort((x))), 6))))
xtable::xtable(genes)

