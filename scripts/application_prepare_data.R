library(AnnotatedPBMC)
library(scater)
library(ggplot2)

source("scripts/application_setup.R")

CACHE_PATH = "../AnnotatedPBMC/data"

dataset_names = c("hao_2020", "tsang_2021", "kotliarov_2020", "10x_sorted", "su_2020", "10x_pbmc_10k", "10x_pbmc_5k_v3", "ding_2019", "blish_2020", "haniffa_2021")

data = lapply(dataset_names, function(dataset) get(paste0("prepare_", dataset))(CACHE_PATH, n_genes = NA, n_sample = NA, sce = TRUE))
names(data) = dataset_names

binning_functions = lapply(data, `[[`, 2)
data = lapply(data, `[[`, 1)

select_genes = function(sce_list) {

  genes = Reduce(intersect, lapply(sce_list, rownames))

  sce_list = lapply(sce_list, function(x) x[genes, ])

  ranks = sapply(sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))

  genes = genes[order(rowMeans(ranks))]

  return(genes)

}

data_split = unlist(lapply(data, function(sce) lapply(sort(unique(sce$dataset)), function(x) sce[, sce$dataset == x])), recursive = FALSE)
genes = select_genes(data_split)
rm(data_split)
gc()
write.csv(genes, file.path(CACHE_PATH, "genes.csv"))
