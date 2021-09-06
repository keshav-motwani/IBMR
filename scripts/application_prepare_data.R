library(AnnotatedPBMC)
library(scater)

CACHE_PATH = "../AnnotatedPBMC/data"

hao = get_hao_2020(CACHE_PATH)
kotliarov = get_kotliarov_2020(CACHE_PATH)
sorted = get_10x_sorted(CACHE_PATH)
pbmc_5k = get_10x_pbmc_5k_v3(CACHE_PATH)
pbmc_10k = get_10x_pbmc_10k(CACHE_PATH)
ding = get_ding_2019(CACHE_PATH)
ding = ding[, grepl("10x", ding$method)]

data = list(hao, kotliarov, sorted, pbmc_5k, pbmc_10k, ding)

data_split = lapply(data, function(sce) lapply(sort(unique(sce$dataset)), function(x) sce[, sce$dataset == x]))

data_split = unlist(data_split, recursive = FALSE)

rm(data)
gc()

for (i in 1:length(data_split)) {

  print(i)

  data_split[[i]] = runPCA(data_split[[i]])
  data_split[[i]] = runUMAP(data_split[[i]])

}

pdf(file.path(RESULT_PATH, "datasets.pdf"), width = 60, height = 35)
scanalysis::plot_reduced_dimensions(data_split, "UMAP", features = "cell_type", label = "cell_type", point_size = 1, facet_rows = ".sample", facet_type = "wrap")
dev.off()

select_genes = function(sce_list, n_genes) {

  genes = Reduce(intersect, lapply(sce_list, rownames))

  sce_list = lapply(sce_list, function(x) x[genes, ])

  ranks = sapply(sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))

  genes = genes[order(rowMeans(ranks))][1:n_genes]

  return(genes)

}

genes = select_genes(data_split, 1000)

write.csv(genes, file.path(CACHE_PATH, "genes.csv"))
