library(AnnotatedPBMC)
library(scater)
library(ggplot2)

source("scripts/application_setup.R")

CACHE_PATH = "../AnnotatedPBMC/data"

dataset_names = c("hao_2020", "tsang_2021", "kotliarov_2020", "10x_sorted", "su_2020", "10x_pbmc_10k", "10x_pbmc_5k_v3", "ding_2019")

data = lapply(dataset_names, function(dataset) get(paste0("prepare_", dataset))(CACHE_PATH, n_genes = NA, n_sample = NA, sce = TRUE))
names(data) = dataset_names

binning_functions = lapply(data, `[[`, 2)
data = lapply(data, `[[`, 1)

select_genes = function(sce_list) {

  genes = Reduce(intersect, lapply(sce_list, rownames))

  sce_list = lapply(sce_list, function(x) x[genes, ])

  ranks = sapply(sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "counts"))$vst.variance.standardized))

  genes = genes[order(rowMeans(ranks))]

  return(genes)

}

data_split = unlist(lapply(data, function(sce) lapply(sort(unique(sce$dataset)), function(x) sce[, sce$dataset == x])), recursive = FALSE)
genes = select_genes(data_split)
rm(data_split)
gc()
write.csv(genes, file.path(CACHE_PATH, "genes.csv"))

table_1 = data.frame(dataset = names(data),
                     number_of_labels = sapply(data, function(x) length(unique(x$cell_type))),
                     labels = sapply(data, function(x) paste0(sort(unique(x$cell_type)), collapse = ", ")),
                     removed_labels = sapply(data, function(x) paste0(sort(attr(x, "removed_labels")), collapse = ", ")))
table_1 = table_1[order(table_1$number_of_labels, decreasing = T), ]
write.csv(table_1, file.path(CACHE_PATH, "table_1.csv"))

set.seed(11111)

categories = c("B intermediate", "B memory", "B naive", "Plasmablast", "CD14 Mono",
               "CD16 Mono", "CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg Memory",
               "Treg Naive", "CD8 Naive", "CD8 TCM", "CD8 TEM", "dnT", "gdT",
               "MAIT", "NK", "NK_CD56bright", "ASDC", "cDC1", "cDC2", "pDC",
               "Eryth", "HSPC", "ILC", "Platelet")
dataset_names = table_1$dataset

colors = c(ggsci::pal_d3("category20")(20), ggsci::pal_d3("category20c")(20))
colors = sample(colors)
gray = "#C7C7C7FF"
colors = setdiff(colors, gray)
names(colors) = as.character(1:length(colors))
colors = c(colors, unobserved = "white")

results = list()

for (k in 1:length(binning_functions)) {

  categories_k = unique(binning_functions[[k]][categories])
  categories_k = setdiff(categories_k, "unobserved")

  common_label = c(as.character(sample(1:39, length(categories_k))), "unobserved")
  names(common_label) = c(categories_k, "unobserved")

  result = common_label[binning_functions[[k]]]
  names(result) = names(binning_functions[[k]])

  result = data.frame(category = names(result), color = result, dataset = names(binning_functions)[k])

  results = c(results, list(result))

}

data = do.call(rbind, results)
data$category = factor(data$category, levels = categories)
data$dataset = factor(data$dataset, levels = rev(dataset_names))

ggplot(data, aes(x = category, y = dataset, fill = color)) +
  geom_tile() +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none", panel.background = element_blank()) +
  geom_hline(yintercept = 0:length(binning_functions) + 0.5) +
  labs(x = "Finest resolution category", y = "Dataset")

ggsave(file.path(CACHE_PATH, "binning_functions.pdf"), height = 4, width = 8)
