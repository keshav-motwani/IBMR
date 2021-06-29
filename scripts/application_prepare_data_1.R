library(SingleCellExperiment)
library(scran)
library(scater)
library(CITEseqData)
library(IBMR)

DATA_PATH = "data/application_1/"
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
data$time = as.factor(data$time)

data = runPCA(data)

cell_types = colData(data)[, c("cell_type_1", "cell_type_2")]
cell_types = cell_types[!duplicated(cell_types), ]
cell_types = cell_types[order(cell_types$cell_type_1), ]

get_category_mappings_train_and_validation = function(cell_types, fine = FALSE) {

  coarse_cell_types = sort(unique(cell_types$cell_type_1))
  fine_cell_types = sort(cell_types$cell_type_2)
  C = length(coarse_cell_types)

  inverse_category_mappings = vector("list", C)
  category_mappings = inverse_category_mappings

  if (!fine) {
    replace = "cell_type_1"
  } else {
    replace = "cell_type_2"
  }

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

batches = sort(unique(data$split))
train_indices = 1:6
validation_indices = 7:12
test_indices = 13:24

X = as.matrix(t(logcounts(data)))
X_list = lapply(batches, function(x) X[data$split == x, ])

category_mappings = get_category_mappings_train_and_validation(cell_types, FALSE)
category_mappings_fine = get_category_mappings_train_and_validation(cell_types, TRUE)

category_mappings_test = category_mappings_fine
category_mappings_test$inverse_category_mappings = c(category_mappings_test$inverse_category_mappings, category_mappings_test$inverse_category_mappings)
category_mappings_test$category_mappings = c(category_mappings_test$category_mappings, category_mappings_test$category_mappings)

Y = data$cell_type_2
names(Y) = Y
Y_list = lapply(batches, function(x) Y[data$split == x])
Y_list[train_indices] = mapply(Y_list[train_indices], category_mappings$inverse_category_mappings, FUN = function(y, map) map[y], SIMPLIFY = FALSE)
Y_list[validation_indices] = mapply(Y_list[validation_indices], category_mappings$inverse_category_mappings, FUN = function(y, map) map[y], SIMPLIFY = FALSE)

Z_list = compute_pca_for_Z_list(X_list[train_indices], 50)
Z_list_int = compute_pca_for_Z_list(X_list[train_indices], 0)

get_fine_categories = function(Y_list) {

  lapply(Y_list, function(Y) names(Y))

}

data = list(
  train = list(
    Y_list = Y_list[train_indices],
    Y_list_fine = get_fine_categories(Y_list[train_indices]),
    category_mappings = category_mappings,
    category_mappings_fine = category_mappings_fine,
    X_list = X_list[train_indices],
    Z_list = Z_list,
    Z_list_int = Z_list_int
  ),
  validation = list(
    Y_list = Y_list[validation_indices],
    Y_list_fine = get_fine_categories(Y_list[validation_indices]),
    category_mappings = category_mappings,
    category_mappings_fine = category_mappings_fine,
    X_list = X_list[validation_indices]
  ),
  test = list(
    Y_list_fine = get_fine_categories(Y_list[test_indices]),
    category_mappings_fine = category_mappings_test,
    X_list = X_list[test_indices]
  )
)

saveRDS(data, file.path(DATA_PATH, "data.rds"))
