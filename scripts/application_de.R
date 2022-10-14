library(IBMR)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(SingleCellExperiment)

source("scripts/application_setup.R")

RESULT_PATH = "results/application_R1"
FIGURES_PATH = "results/application_R1/figures/"

library(ExperimentHub)
eh <- ExperimentHub()
query(eh, "Kang")

(sce <- eh[["EH2259"]])

# DIFFERENCE: keep all genes for now, since a fixed set of genes will be needed for our fitted model

# calculate per-cell quality control (QC) metrics
library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# DIFFERENCE: remove Megakaryocytes
sce <- sce[, sce$cell != "Megakaryocytes"]

# DIFFERENCE: replacing scran::logNormCounts with same normalization used in rest of paper
normalize_gene = function(data) {

  size_factors = Matrix::colSums(data)

  data = log1p(data %*% Matrix::Diagonal(x = 10000/size_factors))

  return(data)

}
logcounts(sce) = normalize_gene(counts(sce))


# DIFFERENCE: need to select genes required for our model
# DIFFERENCE: drop genes with duplicate Ensemble IDs (containing "ENS") and zero counts
sce <- sce[!(rowSums(counts(sce)) == 0 & grepl("ENS", rownames(sce))), ]
dim(sce)
# DIFFERENCE: rename LING02 gene to LINGO2
rownames(sce) = gsub("LING02", "LINGO2", rownames(sce))
# DIFFERENCE: remove Ensemble IDs from gene names
rownames(sce) = sapply(strsplit(rownames(sce), "_ENS"), `[`, 1)
# DIFFERENCE: remove duplicates, keeping one with more counts
sce = sce[names(sort(rowSums(counts(sce)), decreasing = TRUE)), ]
sce = sce[!duplicated(rownames(sce))]

setdiff(rownames(model$Beta), rownames(sce))
sort(rownames(sce)[grepl(paste(rownames(model$Beta), collapse = "|"), rownames(sce)) & grepl("ENS", rownames(sce))])

sce$id <- paste0(sce$stim, sce$ind)
sce$coarse_cell_type = sce$cell

binning_function = c(
  ASDC = "Dendritic cells",
  `B intermediate` = "B cells",
  `B memory` = "B cells",
  `B naive` = "B cells",
  `CD14 Mono` = "CD14+ Monocytes",
  `CD16 Mono` = "FCGR3A+ Monocytes",
  `CD4 CTL` = "CD4 T cells",
  `CD4 Naive` = "CD4 T cells",
  `CD4 TCM` = "CD4 T cells",
  `CD4 TEM` = "CD4 T cells",
  `CD8 Naive` = "CD8 T cells",
  `CD8 TCM` = "CD8 T cells",
  `CD8 TEM` = "CD8 T cells",
  cDC1 = "Dendritic cells",
  cDC2 = "Dendritic cells",
  dnT = "unobserved",
  Eryth = "unobserved",
  gdT = "unobserved",
  HSPC = "unobserved",
  ILC = "unobserved",
  MAIT = "unobserved",
  NK = "NK cells",
  NK_CD56bright = "NK cells",
  pDC = "Dendritic cells",
  Plasmablast = "B cells",
  Platelet = "unobserved",
  `Treg Memory` = "CD4 T cells",
  `Treg Naive` = "CD4 T cells"
)

model = readRDS(file.path(RESULT_PATH, "n_sample_10000_n_genes_1000_split_index_65_IBMR_int_1.rds"))[[1]]
X_sd = model$fit$X_sd
model = model$fit$best_model
category_mapping = binning_function_to_category_mapping(binning_function)
Y = as.character(sce$coarse_cell_type)
P = predict_probabilities(model, list(as.matrix(Matrix::t(logcounts(sce)[rownames(model$Beta), ]))))
C = predict_conditional_probabilities(P, list(Y), list(category_mapping))

predictions = list()
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

ggsave(file.path(FIGURES_PATH, "de_application_heatmap.pdf"), height = 4.6, width = 13.2)

sce$conditional_prediction_cell_type = predictions[["Conditional prediction"]]
table(sce$conditional_prediction_cell_type, sce$coarse_cell_type)

# saveRDS(sce, "kang_2018_with_conditional_preds.rds")
# 
# sce = readRDS("kang_2018_with_conditional_preds.rds")

# DIFFERENCE: now we remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

sce_orig = sce

# Do analysis at fine resolution

(sce <- prepSCE(sce,
                kid = "conditional_prediction_cell_type", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns

t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

t(head(assay(pb)))

res <- pbDS(pb, verbose = FALSE)
# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

# filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

de_gs_by_k <- map(tbl_fil, "gene")
results = de_gs_by_k[names(de_gs_by_k[intersect(names(de_gs_by_k), category_mapping[["CD4 T cells"]])])] # this is a list which contains the DE genes for each of the fine subcategories
names(results) = paste0(names(results), " (fine)")

sce = sce_orig

# Do analysis at coarse resolution

(sce <- prepSCE(sce,
                kid = "coarse_cell_type", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns

t(table(sce$cluster_id, sce$sample_id))

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

t(head(assay(pb)))

# (pb_mds <- pbMDS(pb))

res <- pbDS(pb, verbose = FALSE)
# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

# filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
tbl_fil <- lapply(tbl, function(u) {
  u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
  dplyr::arrange(u, p_adj.loc)
})

library(UpSetR)
de_gs_by_k <- map(tbl_fil, "gene")["CD4 T cells"] # this is a list which contains the DE genes for CD4 T cells as a whole (coarse category)
names(de_gs_by_k) = paste0(names(de_gs_by_k), " (coarse)")
results = c(de_gs_by_k, results)

pdf(file.path(FIGURES_PATH, "de_upset_plot.pdf"), height = 3, width = 5)
upset(fromList(results), nsets=length(results))
dev.off()

setdiff(results[[4]], unlist(results[c(1:3, 5)]))

