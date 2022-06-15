prepare_real_data_application = function(split_index,
                                         n_genes,
                                         n_sample,
                                         cache_path,
                                         replicate) {

  set.seed(replicate, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

  dataset_names = read.csv(file.path("data/", "table_1.csv"))$dataset
  stopifnot(length(dataset_names) == 10)

  splits = expand.grid(setdiff(dataset_names, "hao_2020"), setdiff(dataset_names, "hao_2020"), stringsAsFactors = FALSE)
  colnames(splits) = c("validation", "test")
  splits = splits[splits[, 1] != splits[, 2], ]
  rownames(splits) = 1:nrow(splits)
  stopifnot(nrow(splits) == 72)
  split = splits[split_index, , drop = TRUE]

  n_sample = rep(n_sample, length(dataset_names))
  names(n_sample) = dataset_names
  n_sample[unlist(split)] = NA

  datasets = mapply(
    dataset_names,
    n_sample,
    FUN = function(dataset, n) {
      result = get(paste0("prepare_", dataset))(cache_path, n_genes, n)
      gc()
      return(result)
    },
    SIMPLIFY = FALSE
  )
  names(datasets) = dataset_names

  train_datasets = datasets[setdiff(dataset_names, unlist(split))]
  validation_datasets = datasets[split$validation]
  test_datasets = datasets[split$test]

  categories = train_datasets[[1]]$categories

  Y_list = unlist(lapply(train_datasets, `[[`, "Y_list"), recursive = FALSE)
  X_list = unlist(lapply(train_datasets, `[[`, "X_list"), recursive = FALSE)
  category_mappings = unlist(lapply(train_datasets, `[[`, "category_mappings"), recursive = FALSE)
  inverse_category_mappings = unlist(lapply(train_datasets, `[[`, "inverse_category_mappings"), recursive = FALSE)
  category_mappings = list(categories = categories, category_mappings = category_mappings, inverse_category_mappings = inverse_category_mappings)
  sample_list = unlist(lapply(train_datasets, `[[`, "sample_list"), recursive = FALSE)

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
                        sample_list = sample_list,
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

prepare_hao_2020 = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_hao_2020(cache_path)

  SingleCellExperiment::altExp(data) = NULL
  if (!sce) SingleCellExperiment::counts(data) = NULL

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data$cell_type = ifelse(data$cell_type_2 == "Treg", data$cell_type_3, data$cell_type_2)

  removed_labels = "*Proliferating*"
  data = data[, !grepl(removed_labels, data$cell_type)]
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

  binning_function = c(
    ASDC = "ASDC",
    `B intermediate` = "B intermediate",
    `B memory` = "B memory",
    `B naive` = "B naive",
    `CD14 Mono` = "CD14 Mono",
    `CD16 Mono` = "CD16 Mono",
    `CD4 CTL` = "CD4 CTL",
    `CD4 Naive` = "CD4 Naive",
    `CD4 TCM` = "CD4 TCM",
    `CD4 TEM` = "CD4 TEM",
    `CD8 Naive` = "CD8 Naive",
    `CD8 TCM` = "CD8 TCM",
    `CD8 TEM` = "CD8 TEM",
    cDC1 = "cDC1",
    cDC2 = "cDC2",
    dnT = "dnT",
    Eryth = "Eryth",
    gdT = "gdT",
    HSPC = "HSPC",
    ILC = "ILC",
    MAIT = "MAIT",
    NK = "NK",
    NK_CD56bright = "NK_CD56bright",
    pDC = "pDC",
    Plasmablast = "Plasmablast",
    Platelet = "Platelet",
    `Treg Memory` = "Treg Memory",
    `Treg Naive` = "Treg Naive"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_kotliarov_2020 = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_kotliarov_2020(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data$cell_type = data$cell_type_1

  removed_labels = "Unconv T"
  data = data[, data$cell_type != removed_labels]
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

  binning_function = c(
    ASDC = "unobserved",
    `B intermediate` = "B",
    `B memory` = "B",
    `B naive` = "B",
    `CD14 Mono` = "Monocyte/mDC",
    `CD16 Mono` = "Non-classical monocyte",
    `CD4 CTL` = "unobserved",
    `CD4 Naive` = "CD4+ naive T/DNT",
    `CD4 TCM` = "CD4+ memory T",
    `CD4 TEM` = "CD4+ memory T",
    `CD8 Naive` = "CD8+ naive T",
    `CD8 TCM` = "CD8+ memory T",
    `CD8 TEM` = "CD8+ memory T",
    cDC1 = "unobserved",
    cDC2 = "unobserved",
    dnT = "CD4+ naive T/DNT",
    Eryth = "unobserved",
    gdT = "unobserved",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "NK",
    NK_CD56bright = "NK",
    pDC = "pDC",
    Plasmablast = "B",
    Platelet = "unobserved",
    `Treg Memory` = "CD4+ memory T",
    `Treg Naive` = "CD4+ naive T/DNT"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_haniffa_2021 = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_haniffa_2021(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data$cell_type = data$cell_type_1

  removed_labels = "*prolif*"
  data = data[, !grepl(removed_labels, data$cell_type)]
  attr(data, "removed_labels") = removed_labels

  data = data[, ]

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

  binning_function = c(
    ASDC = "DCs",
    `B intermediate` = "B_cell",
    `B memory` = "B_cell",
    `B naive` = "B_cell",
    `CD14 Mono` = "CD14",
    `CD16 Mono` = "CD16",
    `CD4 CTL` = "CD4",
    `CD4 Naive` = "CD4",
    `CD4 TCM` = "CD4",
    `CD4 TEM` = "CD4",
    `CD8 Naive` = "CD8",
    `CD8 TCM` = "CD8",
    `CD8 TEM` = "CD8",
    cDC1 = "DCs",
    cDC2 = "DCs",
    dnT = "unobserved",
    Eryth = "RBC",
    gdT = "gdT",
    HSPC = "HSC",
    ILC = "unobserved",
    MAIT = "MAIT",
    NK = "NK_16hi",
    NK_CD56bright = "NK_56hi",
    pDC = "pDC",
    Plasmablast = "Plasmablast",
    Platelet = "Platelets",
    `Treg Memory` = "Treg",
    `Treg Naive` = "Treg"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_tsang_2021 = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  # data = AnnotatedPBMC::get_tsang_2021(cache_path)
  data = AnnotatedPBMC::get_tsang_2021(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data = data[, !grepl("TCRVbeta13.1pos|TissueResMemT|double-positive T cell \\(DPT)|granulocyte|intermediate monocyte|NK_CD56loCD16lo", data$cell_type)]

  removed_labels = c("TCRVbeta13.1pos", "TissueResMemT", "double-positive T cell (DPT)", "granulocyte", "intermediate monocyte", "NK_CD56loCD16lo")
  data = data[, !(data$cell_type %in% removed_labels)]
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

  binning_function = c(
    ASDC = "conventional dendritic cell",
    `B intermediate` = "unobserved",
    `B memory` = "memory B cell",
    `B naive` = "naive B cell",
    `CD14 Mono` = "classical monocyte",
    `CD16 Mono` = "non-classical monocyte",
    `CD4 CTL` = "unobserved",
    `CD4 Naive` = "naive CD4+ T cell",
    `CD4 TCM` = "CD4-positive, alpha-beta memory T cell",
    `CD4 TEM` = "CD4-positive, alpha-beta memory T cell",
    `CD8 Naive` = "naive CD8+ T cell",
    `CD8 TCM` = "CD8-positive, alpha-beta memory T cell",
    `CD8 TEM` = "CD8-positive, alpha-beta memory T cell",
    cDC1 = "conventional dendritic cell",
    cDC2 = "conventional dendritic cell",
    dnT = "double negative T cell (DNT)",
    Eryth = "unobserved",
    gdT = "gamma-delta T cell",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "mucosal invariant T cell (MAIT)",
    NK = "NK_CD16hi",
    NK_CD56bright = "NK_CD56hiCD16lo",
    pDC = "plasmacytoid dendritic cell",
    Plasmablast = "plasmablast",
    Platelet = "platelet",
    `Treg Memory` = "regulatory T cell",
    `Treg Naive` = "regulatory T cell"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_blish_2020 = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_blish_2020(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data$cell_type = data$cell_type_1

  removed_labels = "Granulocyte"
  data = data[, data$cell_type != removed_labels]
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

  binning_function = c(
    ASDC = "DC",
    `B intermediate` = "B",
    `B memory` = "B",
    `B naive` = "B",
    `CD14 Mono` = "CD14 Monocyte",
    `CD16 Mono` = "CD16 Monocyte",
    `CD4 CTL` = "CD4 T",
    `CD4 Naive` = "CD4 T",
    `CD4 TCM` = "CD4 T",
    `CD4 TEM` = "CD4 T",
    `CD8 Naive` = "CD8 T",
    `CD8 TCM` = "CD8 T",
    `CD8 TEM` = "CD8 T",
    cDC1 = "DC",
    cDC2 = "DC",
    dnT = "unobserved",
    Eryth = "RBC",
    gdT = "gd T",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "NK",
    NK_CD56bright = "NK",
    pDC = "pDC",
    Plasmablast = "PB",
    Platelet = "Platelet",
    `Treg Memory` = "CD4 T",
    `Treg Naive` = "CD4 T"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_10x_sorted = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_10x_sorted(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  removed_labels = c()
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

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
    gdT = "unobserved",
    HSPC = "CD34+ Cells",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "CD56+ Natural Killer Cells",
    NK_CD56bright = "CD56+ Natural Killer Cells",
    pDC = "unobserved",
    Plasmablast = "CD19+ B cells",
    Platelet = "unobserved",
    `Treg Memory` = "CD4+/CD25+ Regulatory T Cells",
    `Treg Naive` = "CD4+/CD25+ Regulatory T Cells"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_10x_pbmc_10k = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_10x_pbmc_10k(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data$cell_type = data$cell_type_2

  removed_labels = "intermediate monocyte"
  data = data[, data$cell_type != removed_labels]
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

  binning_function = c(
    ASDC = "unobserved",
    `B intermediate` = "B",
    `B memory` = "B",
    `B naive` = "B",
    `CD14 Mono` = "classical monocyte",
    `CD16 Mono` = "unobserved",
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
    gdT = "unobserved",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "CD16+ NK",
    NK_CD56bright = "CD16- NK",
    pDC = "unobserved",
    Plasmablast = "B",
    Platelet = "unobserved",
    `Treg Memory` = "Treg",
    `Treg Naive` = "Treg"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_10x_pbmc_5k_v3 = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_10x_pbmc_5k_v3(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data$cell_type = data$cell_type_2

  removed_labels = "intermediate monocyte"
  data = data[, data$cell_type != removed_labels]
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

  binning_function = c(
    ASDC = "DCs",
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
    gdT = "unobserved",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "CD16+ NK",
    NK_CD56bright = "CD16- NK",
    pDC = "DCs",
    Plasmablast = "unobserved",
    Platelet = "unobserved",
    `Treg Memory` = "Treg",
    `Treg Naive` = "Treg"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_su_2020 = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_su_2020(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  data$cell_type = data$cell_type_2

  removed_labels = "intermediate monocyte"
  data = data[, data$cell_type != removed_labels]
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

  binning_function = c(
    ASDC = "myeloid DC",
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
    cDC1 = "myeloid DC",
    cDC2 = "myeloid DC",
    dnT = "unobserved",
    Eryth = "unobserved",
    gdT = "unobserved",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "CD16+ NK",
    NK_CD56bright = "CD16- NK",
    pDC = "plasmacytoid DC",
    Plasmablast = "unobserved",
    Platelet = "unobserved",
    `Treg Memory` = "Treg",
    `Treg Naive` = "Treg"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_ding_2019 = function(cache_path, n_genes = NA, n_sample = NA, sce = FALSE) {

  data = AnnotatedPBMC::get_ding_2019(cache_path)

  if (!is.na(n_genes)) {
    genes = read.csv(file.path("data/", "genes.csv"))[, 2][1:n_genes]
    data = data[genes, ]
  }

  removed_labels = "Megakaryocyte"
  data = data[, grepl("10x", data$method) & data$cell_type != removed_labels]
  attr(data, "removed_labels") = removed_labels

  data = data[, uniform_sample(data$cell_type, ifelse(is.na(n_sample), ncol(data), n_sample))]

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
    gdT = "unobserved",
    HSPC = "unobserved",
    ILC = "unobserved",
    MAIT = "unobserved",
    NK = "Natural killer cell",
    NK_CD56bright = "Natural killer cell",
    pDC = "Plasmacytoid dendritic cell",
    Plasmablast = "B cell",
    Platelet = "unobserved",
    `Treg Memory` = "CD4+ T cell",
    `Treg Naive` = "CD4+ T cell"
  )

  return(prepare_dataset_output(data, binning_function, sce))

}

prepare_dataset_output = function(data, binning_function, sce) {

  stopifnot(all(sort(names(binning_function)) == c("ASDC", "B intermediate", "B memory", "B naive", "CD14 Mono",
                                                   "CD16 Mono", "CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD8 Naive",
                                                   "CD8 TCM", "CD8 TEM", "cDC1", "cDC2", "dnT", "Eryth", "gdT",
                                                   "HSPC", "ILC", "MAIT", "NK", "NK_CD56bright", "pDC", "Plasmablast",
                                                   "Platelet", "Treg Memory", "Treg Naive")))
  stopifnot(length(setdiff(data$cell_type, binning_function)) == 0)
  if(!(length(setdiff(binning_function, data$cell_type)) == 0 || all(setdiff(binning_function, data$cell_type) == "unobserved"))) warning("Some labels from binning function not observed")

  if (sce) return(list(sce = data, binning_function = binning_function))

  # binning_function[binning_function == "unobserved"] = names(binning_function)[binning_function == "unobserved"]

  data = list(data) # split_data(data, data$dataset)

  X_list = lapply(data, function(x) t(as.matrix(SingleCellExperiment::logcounts(x))))
  Y_list = lapply(data, function(x) as.character(x$cell_type))
  sample_list = lapply(data, function(x) as.character(x$dataset))

  category_mapping = binning_function_to_category_mapping(binning_function)

  return(list(Y_list = Y_list, X_list = X_list, categories = names(binning_function), category_mappings = replicate(length(Y_list), category_mapping, simplify = FALSE), inverse_category_mappings = replicate(length(Y_list), binning_function, simplify = FALSE), sample_list = sample_list))

}

uniform_sample = function(Y, n) {

  indices = sample(1:length(Y), min(n, length(Y)))

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
