prepare_real_data_application = function(split_index,
                                         n_sample,
                                         cache_path,
                                         replicate) {

  set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  dataset_names = c("hao_2020", "kotliarov_2020", "10x_sorted", "10x_pbmc_10k", "10x_pbmc_5k_v3", "ding_2019")

  datasets = lapply(dataset_names, function(dataset) get(paste0("prepare_", dataset))(cache_path, n_sample))
  names(datasets) = dataset_names

  splits = expand.grid(dataset_names[-1], dataset_names[-1], stringsAsFactors = FALSE)
  colnames(splits) = c("validation", "test")
  splits = splits[splits[, 1] != splits[, 2], ]
  rownames(splits) = 1:nrow(splits)
  split = splits[split_index, , drop = TRUE]

  train_datasets = datasets[setdiff(dataset_names, unlist(split))]
  validation_datasets = datasets[split$validation]
  test_datasets = datasets[split$test]

  categories = train_datasets[[1]]$categories

  Y_list = unlist(lapply(train_datasets, `[[`, "Y_list"), recursive = FALSE)
  X_list = unlist(lapply(train_datasets, `[[`, "X_list"), recursive = FALSE)
  category_mappings = unlist(lapply(train_datasets, `[[`, "category_mappings"), recursive = FALSE)
  inverse_category_mappings = unlist(lapply(train_datasets, `[[`, "inverse_category_mappings"), recursive = FALSE)
  category_mappings = list(categories = categories, category_mappings = category_mappings, inverse_category_mappings = inverse_category_mappings)

  Y_list_val = unlist(lapply(validation_datasets, `[[`, "Y_list"), recursive = FALSE)
  X_list_val = unlist(lapply(validation_datasets, `[[`, "X_list"), recursive = FALSE)
  category_mappings_val = unlist(lapply(validation_datasets, `[[`, "category_mappings"), recursive = FALSE)
  inverse_category_mappings_val = unlist(lapply(validation_datasets, `[[`, "inverse_category_mappings"), recursive = FALSE)
  category_mappings_val = list(categories = categories, category_mappings = category_mappings_val, inverse_category_mappings = inverse_category_mappings_val)

  Y_list_test = unlist(lapply(test_datasets, `[[`, "Y_list"), recursive = FALSE)
  X_list_test = unlist(lapply(test_datasets, `[[`, "X_list"), recursive = FALSE)
  category_mappings_test = unlist(lapply(test_datasets, `[[`, "category_mappings"), recursive = FALSE)
  inverse_category_mappings_test = unlist(lapply(test_datasets, `[[`, "inverse_category_mappings"), recursive = FALSE)
  category_mappings_test = list(categories = categories, category_mappings = category_mappings_test, inverse_category_mappings = inverse_category_mappings_test)

  output = prepare_data(Y_list = Y_list,
                        category_mappings = category_mappings,
                        category_mappings_fine = NULL,
                        X_list = X_list,
                        X_star_list = NULL,
                        Y_list_validation = Y_list_val,
                        category_mappings_validation = category_mappings_val,
                        category_mappings_fine_validation = NULL,
                        X_list_validation = X_list_val,
                        X_star_list_validation = NULL,
                        Y_list_test = Y_list_test,
                        category_mappings_test = category_mappings_test,
                        X_list_test = X_list_test,
                        alpha = NULL,
                        Beta = NULL)

}

prepare_hao_2020 = function(cache_path, n_sample = 5000) {

  data = AnnotatedPBMC::get_hao_2020(cache_path)

  SingleCellExperiment::altExp(data) = NULL
  SingleCellExperiment::counts(data) = NULL

  genes = read.csv(file.path(cache_path, "genes.csv"))[, 2]
  data = data[genes, ]

  data = data[, !grepl("Proliferating", data$cell_type_2)]

  data = data[, weighted_sample(data$cell_type, n_sample)]

  binning_function = setNames(nm = sort(unique(data$cell_type)))

  return(prepare_dataset_output(data, binning_function))

}

prepare_kotliarov_2020 = function(cache_path, n_sample = 5000) {

  data = AnnotatedPBMC::get_kotliarov_2020(cache_path)

  genes = read.csv(file.path(cache_path, "genes.csv"))[, 2]
  data = data[genes, ]

  data$cell_type = gsub("Classical monocytes|IgA\\+ monocytes", "Classical monocytes", data$cell_type, fixed = FALSE)
  data$cell_type = gsub("Transitional B|Switched B|Unswitched B", "B cells", data$cell_type, fixed = FALSE)
  data = data[, !(data$cell_type %in% c("mDC", "NKT-like", "CD8+ CD103+ T", "Unconventional CD161hi CD8+ T", "CD161+ double-negative T"))]

  data = data[, weighted_sample(data$cell_type, n_sample)]

  binning_function = c(
    ASDC = "unobserved",
    `B intermediate` = "B cells",
    `B memory` = "B cells",
    `B naive` = "B cells",
    `CD14 Mono` = "Classical monocytes",
    `CD16 Mono` = "Non-classical monocytes",
    `CD4 CTL` = "unobserved",
    `CD4 Naive` = "CD4+ naive T",
    `CD4 TCM` = "CD4+ central and transitional memory T",
    `CD4 TEM` = "CD4+ TEMRA and effector memory T",
    `CD8 Naive` = "CD8+ naive T",
    `CD8 TCM` = "CD8+ central and transitional memory T",
    `CD8 TEM` = "CD8+ TEMRA and effector memory T",
    cDC1 = "unobserved",
    cDC2 = "unobserved",
    dnT = "Double-negative T",
    Eryth = "unobserved",
    HSPC = "HSC",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "CD16++ NK",
    `NK CD56bright` = "CD56hi CD16lo NK",
    pDC = "pDC",
    Plasmablast = "B cells",
    Platelet = "unobserved",
    `Treg Memory` = "unobserved",
    `Treg Naive` = "unobserved"
  )

  return(prepare_dataset_output(data, binning_function))

}

prepare_10x_sorted = function(cache_path, n_sample = 5000) {

  data = AnnotatedPBMC::get_10x_sorted(cache_path)

  genes = read.csv(file.path(cache_path, "genes.csv"))[, 2]
  data = data[genes, ]

  data = data[, weighted_sample(data$cell_type, n_sample)]

  binning_function = c(
    ASDC = "unobserved",
    `B intermediate` = "CD19+ B cells",
    `B memory` = "CD19+ B cells",
    `B naive` = "CD19+ B cells",
    `CD14 Mono` = "CD14+ Monocytes",
    `CD16 Mono` = "unobserved",
    `CD4 CTL` = "unobserved",
    `CD4 Naive` = "CD4+/CD45RA+/CD25- Naive T cells",
    `CD4 TCM` = "CD4+/CD45RO+ Memory T Cells",
    `CD4 TEM` = "CD4+/CD45RO+ Memory T Cells",
    `CD8 Naive` = "CD8+/CD45RA+ Naive Cytotoxic T Cells",
    `CD8 TCM` = "unobserved",
    `CD8 TEM` = "unobserved",
    cDC1 = "unobserved",
    cDC2 = "unobserved",
    dnT = "unobserved",
    Eryth = "unobserved",
    HSPC = "CD34+ Cells",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "CD56+ Natural Killer Cells",
    `NK CD56bright` = "CD56+ Natural Killer Cells",
    pDC = "unobserved",
    Plasmablast = "CD19+ B cells",
    Platelet = "unobserved",
    `Treg Memory` = "CD4+/CD25+ Regulatory T Cells",
    `Treg Naive` = "CD4+/CD25+ Regulatory T Cells"
  )

  return(prepare_dataset_output(data, binning_function))

}

prepare_10x_pbmc_10k = function(cache_path, n_sample = 5000) {

  data = AnnotatedPBMC::get_10x_pbmc_10k(cache_path)

  genes = read.csv(file.path(cache_path, "genes.csv"))[, 2]
  data = data[genes, ]

  data = data[, data$cell_type != "intermediate monocyte"]

  data = data[, weighted_sample(data$cell_type, n_sample)]

  binning_function = c(
    ASDC = "unobserved",
    `B intermediate` = "B",
    `B memory` = "B",
    `B naive` = "B",
    `CD14 Mono` = "classical monocyte",
    `CD16 Mono` = "non-classical CD16+ monocyte",
    `CD4 CTL` = "unobserved",
    `CD4 Naive` = "naive CD4",
    `CD4 TCM` = "memory CD4",
    `CD4 TEM` = "memory CD4",
    `CD8 Naive` = "naive CD8",
    `CD8 TCM` = "memory CD8",
    `CD8 TEM` = "memory CD8",
    cDC1 = "unobserved",
    cDC2 = "unobserved",
    dnT = "unobserved",
    Eryth = "unobserved",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "CD16+ NK",
    `NK CD56bright` = "CD16- NK",
    pDC = "unobserved",
    Plasmablast = "B",
    Platelet = "unobserved",
    `Treg Memory` = "Treg",
    `Treg Naive` = "Treg"
  )

  return(prepare_dataset_output(data, binning_function))

}

prepare_10x_pbmc_5k_v3 = function(cache_path, n_sample = 5000) {

  data = AnnotatedPBMC::get_10x_pbmc_5k_v3(cache_path)

  genes = read.csv(file.path(cache_path, "genes.csv"))[, 2]
  data = data[genes, ]

  data = data[, data$cell_type != "intermediate monocyte"]

  data = data[, weighted_sample(data$cell_type, n_sample)]

  binning_function = c(
    ASDC = "unobserved",
    `B intermediate` = "unobserved",
    `B memory` = "memory B",
    `B naive` = "naive B",
    `CD14 Mono` = "classical monocyte",
    `CD16 Mono` = "non-classical CD16+ monocyte",
    `CD4 CTL` = "unobserved",
    `CD4 Naive` = "naive CD4",
    `CD4 TCM` = "memory CD4",
    `CD4 TEM` = "memory CD4",
    `CD8 Naive` = "naive CD8",
    `CD8 TCM` = "memory CD8",
    `CD8 TEM` = "memory CD8",
    cDC1 = "DCs",
    cDC2 = "DCs",
    dnT = "unobserved",
    Eryth = "unobserved",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "CD16+ NK",
    `NK CD56bright` = "CD16- NK",
    pDC = "DCs",
    Plasmablast = "B",
    Platelet = "unobserved",
    `Treg Memory` = "Treg",
    `Treg Naive` = "Treg"
  )

  return(prepare_dataset_output(data, binning_function))

}

prepare_ding_2019 = function(cache_path, n_sample = 5000) {

  data = AnnotatedPBMC::get_ding_2019(cache_path)

  genes = read.csv(file.path(cache_path, "genes.csv"))[, 2]
  data = data[genes, ]

  data = data[, grepl("10x", data$method) & data$cell_type != "Megakaryocyte"]

  data = data[, weighted_sample(data$cell_type, n_sample)]

  binning_function = c(
    ASDC = "Dendritic cell",
    `B intermediate` = "B cell",
    `B memory` = "B cell",
    `B naive` = "B cell",
    `CD14 Mono` = "CD14+ monocyte",
    `CD16 Mono` = "CD16+ monocyte",
    `CD4 CTL` = "CD4+ T cell",
    `CD4 Naive` = "CD4+ T cell",
    `CD4 TCM` = "CD4+ T cell",
    `CD4 TEM` = "CD4+ T cell",
    `CD8 Naive` = "Cytotoxic T cell",
    `CD8 TCM` = "Cytotoxic T cell",
    `CD8 TEM` = "Cytotoxic T cell",
    cDC1 = "Dendritic cell",
    cDC2 = "Dendritic cell",
    dnT = "unobserved",
    Eryth = "unobserved",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "Natural killer cell",
    `NK CD56bright` = "Natural killer cell",
    pDC = "Plasmacytoid dendritic cell",
    Plasmablast = "B cell",
    Platelet = "unobserved",
    `Treg Memory` = "CD4+ T cell",
    `Treg Naive` = "CD4+ T cell"
  )

  return(prepare_dataset_output(data, binning_function))

}

prepare_dataset_output = function(data, binning_function) {

  stopifnot(length(setdiff(data$cell_type, binning_function)) == 0)
  stopifnot(length(setdiff(binning_function, data$cell_type)) == 0 || setdiff(binning_function, data$cell_type) == "unobserved")

  data = split_data(data, data$dataset)

  X_list = lapply(data, function(x) t(as.matrix(SingleCellExperiment::logcounts(x))))
  Y_list = lapply(data, function(x) x$cell_type)

  category_mapping = binning_function_to_category_mapping(binning_function)

  return(list(Y_list = Y_list, X_list = X_list, categories = names(binning_function), category_mappings = replicate(length(Y_list), category_mapping, simplify = FALSE), inverse_category_mappings = replicate(length(Y_list), binning_function, simplify = FALSE)))

}

weighted_sample = function(Y, n) {

  weights = 1 / table(Y)
  indices = sample(1:length(Y), min(n, length(Y)), prob = weights[Y])

  return(indices)

}

split_data = function(sce, split_by) {

  lapply(sort(unique(split_by)), function(x) sce[, split_by == x])

}

binning_function_to_category_mapping = function(binning_function) {

  category_mapping = list()

  for (label in unique(binning_function)) {

    category_mapping[[label]] = names(binning_function)[which(binning_function == label)]

  }

  return(category_mapping)

}
