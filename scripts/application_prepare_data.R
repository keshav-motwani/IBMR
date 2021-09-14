library(AnnotatedPBMC)
library(scater)

CACHE_PATH = "../AnnotatedPBMC/data"

dataset_names = c("hao_2020", "haniffa_2021", "tsang_2021", "blish_2020", "kotliarov_2020", "10x_sorted", "su_2020", "10x_pbmc_10k", "10x_pbmc_5k_v3", "ding_2019")

data = lapply(dataset_names, function(dataset) get(paste0("get_", dataset))(CACHE_PATH))
names(data) = dataset_names

data$ding_2019 = data$ding_2019[, grepl("10x", data$ding_2019$method)]

data_split = lapply(data, function(sce) lapply(sort(unique(sce$dataset)), function(x) sce[, sce$dataset == x]))

data_split = unlist(data_split, recursive = FALSE)

names(data_split) = unlist(lapply(data, function(sce) sort(unique(sce$dataset))))

rm(data)
gc()

for (i in 1:length(data_split)) {

  print(i)

  data_split[[i]] = runPCA(data_split[[i]])
  data_split[[i]] = runUMAP(data_split[[i]])

}

pdf(file.path(CACHE_PATH, "datasets.pdf"), width = 60, height = 35)
scanalysis::plot_reduced_dimensions(data_split, "UMAP", features = "cell_type", label = "cell_type", point_size = 1, facet_rows = ".sample", facet_type = "wrap")
dev.off()

select_genes = function(sce_list) {

  genes = Reduce(intersect, lapply(sce_list, rownames))

  sce_list = lapply(sce_list, function(x) x[genes, ])

  ranks = sapply(sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))

  genes = genes[order(rowMeans(ranks))]

  return(genes)

}

genes = select_genes(data_split)

write.csv(genes, file.path(CACHE_PATH, "genes.csv"))
