library(CITEseqData)
library(SingleCellExperiment)

DATA_PATH = "data/simulation/"
dir.create(DATA_PATH, recursive = TRUE)

data = get_hao_3_prime_data(DATA_PATH)

data$split = paste0(data$donor, "_", data$time)
split = lapply(unique(data$split), function(x) data[, which(data$split == x)])
names(split) = unique(data$split)
split = split[sort(names(split))]

select_genes = function(sce_list, n_genes) {

  genes = Reduce(intersect, lapply(sce_list, rownames))
  ranks = sapply(sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))
  genes = genes[order(rowMeans(ranks))][1:n_genes]
  return(genes)

}

genes = select_genes(split, 2000)
data = data[genes, ]

data = data[, !grepl("other", data$cell_type_1)]
data = data[, !grepl("Proliferating", data$cell_type_2)]
data = data[, !grepl("ASDC|cDC1", data$cell_type_2)]
data$cell_type_2[data$cell_type_2 == "NK"] = "NK_CD56dim"

X = as.matrix(t(logcounts(data)))
Y = data$cell_type_2

require(doMC)
registerDoMC(cores = 64)
fit = glmnet::cv.glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", trace.it = 1, parallel = TRUE)
saveRDS(fit, file.path(DATA_PATH, "hao_glmnet_fit.rds"))

cell_types = colData(data)[, c("cell_type_1", "cell_type_2")]
cell_types = cell_types[!duplicated(cell_types), ]
cell_types = cell_types[order(cell_types$cell_type_1), ]

get_category_mappings = function(cell_types) {

  coarse_cell_types = sort(unique(cell_types$cell_type_1))
  fine_cell_types = sort(cell_types$cell_type_2)
  C = length(coarse_cell_types)

  inverse_category_mappings = vector("list", C)
  category_mappings = inverse_category_mappings

  replace = "cell_type_1"

  for (i in 1:C) {

    keep_fine = coarse_cell_types[i]

    inverse_category_mapping = c(cell_types[cell_types$cell_type_1 == keep_fine, "cell_type_2"], cell_types[cell_types$cell_type_1 != keep_fine, replace])
    names(inverse_category_mapping) = c(cell_types[cell_types$cell_type_1 == keep_fine, "cell_type_2"], cell_types[cell_types$cell_type_1 != keep_fine, "cell_type_2"])
    inverse_category_mapping = inverse_category_mapping[fine_cell_types]

    category_mapping = list()

    for (label in unique(inverse_category_mapping)) {

      category_mapping[[label]] = names(inverse_category_mapping)[which(inverse_category_mapping == label)]

    }

    inverse_category_mappings[[i]] = inverse_category_mapping
    category_mappings[[i]] = category_mapping

  }

  return(list(inverse_category_mappings = inverse_category_mappings, category_mappings = category_mappings, categories = fine_cell_types))

}

category_mappings = get_category_mappings(cell_types, FALSE)

saveRDS(fit, file.path(DATA_PATH, "hao_category_mappings.rds"))
