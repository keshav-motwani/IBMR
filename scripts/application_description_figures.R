library(AnnotatedPBMC)
library(scater)
library(ggplot2)
library(IBMR)

source("scripts/application_setup.R")

CACHE_PATH = "../AnnotatedPBMC/data/"
FIGURES_PATH = "figures/"
dir.create(FIGURES_PATH, recursive = TRUE)

dataset_names = c("hao_2020", "tsang_2021", "kotliarov_2020", "10x_sorted", "su_2020", "10x_pbmc_10k", "10x_pbmc_5k_v3", "ding_2019", "blish_2020", "haniffa_2021")

data = lapply(dataset_names, function(dataset) {
  x = get(paste0("prepare_", dataset))(CACHE_PATH, n_genes = 5, n_sample = NA, sce = TRUE)
  binning_function = x[[2]]
  x = x[[1]]
  number_of_labels = length(unique(x$cell_type))
  labels = paste0(sort(unique(x$cell_type)), collapse = ", ")
  removed_labels = paste0(sort(attr(x, "removed_labels")), collapse = ", ")
  return(list(number_of_labels = number_of_labels, labels = labels, removed_labels = removed_labels, binning_function = binning_function))
})
names(data) = dataset_names

references = list("hao2020integrated", "tsangliu2021time", "kotliarov2020broad", "zheng2017massively", c("su2020multi", "shasha2021superscan"), c("10x_pbmc_10k", "shasha2021superscan"), c("10x_pbmc_5k_v3", "shasha2021superscan"), "ding2019systematic", "blishwilk2020single", "haniffastephenson2021single")
references = sapply(lapply(references, function(x) paste0("\\citet{", x, "}")), function(y) paste0(y, collapse = ", "))

table_1 = data.frame(dataset = names(data),
                     number_of_labels = sapply(data, `[[`, "number_of_labels"),
                     references = references,
                     labels = sapply(data, `[[`, "labels"),
                     removed_labels = sapply(data, `[[`, "removed_labels"))
table_1 = table_1[order(table_1$number_of_labels, decreasing = T), ]
dataset_names = table_1$dataset
table_1$dataset = paste0("\\texttt{", gsub("_", "\\_", table_1$dataset, fixed = TRUE), "}")
table_1$number_of_labels = as.character(table_1$number_of_labels)
colnames(table_1) = c("Dataset", "\\# of labels", "Reference(s)", "Labels", "Removed labels")
rownames(table_1) = NULL

print(xtable::xtable(table_1[, 1:3]), sanitize.text.function = identity, sanitize.colnames.function = identity)
print(xtable::xtable(table_1[, c(1, 4, 5)]), sanitize.text.function = identity, sanitize.colnames.function = identity)

binning_functions = lapply(data[dataset_names], `[[`, "binning_function")

set.seed(11111)

categories = c("B intermediate", "B memory", "B naive", "Plasmablast", "CD14 Mono",
               "CD16 Mono", "CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg Memory",
               "Treg Naive", "CD8 Naive", "CD8 TCM", "CD8 TEM", "dnT", "gdT",
               "MAIT", "NK", "NK_CD56bright", "ASDC", "cDC1", "cDC2", "pDC",
               "Eryth", "HSPC", "ILC", "Platelet")

plot_binning_functions = function(binning_functions, categories, dataset_names = names(binning_functions)) {

  if (length(categories) > 10) colors = c(ggsci::pal_d3("category20")(20), ggsci::pal_d3("category20c")(20)) else colors = ggsci::pal_d3("category20c")(10)
  colors = sample(colors)
  gray = "#C7C7C7FF"
  colors = setdiff(colors, gray)
  names(colors) = as.character(1:length(colors))
  colors = c(colors, unobserved = "white")

  results = list()

  for (k in 1:length(binning_functions)) {

    categories_k = unique(binning_functions[[k]][categories])
    categories_k = setdiff(categories_k, "unobserved")

    common_label = c(as.character(sample(1:(length(colors) - 1), length(categories_k))), "unobserved")
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
    theme_minimal(base_size = 16) +
    theme(legend.position = "none", panel.background = element_blank()) +
    guides(x = guide_axis(angle = 45)) +
    geom_hline(yintercept = 0:length(binning_functions) + 0.5) +
    labs(x = "Finest resolution category", y = "Dataset") +
    coord_fixed()

}

plot_binning_functions(binning_functions, categories, dataset_names)
ggsave(file.path(FIGURES_PATH, "binning_functions.pdf"), height = 6, width = 13.2)

categories = c(
  "naive CD4+",
  "effector memory CD4+",
  "central memory CD4+",
  "naive CD8+",
  "effector memory CD8+",
  "central memory CD8+"
)
binning_functions = list(
  "Dataset 1" = c(rep("CD4+", 3),
                  rep("CD8+", 3)),
  "Dataset 2" = c(categories[1:3],
                  rep("CD8+", 3)),
  "Dataset 3" = c(categories[1],
                  rep("memory CD4+", 2),
                  categories[4:6])
)
binning_functions = lapply(binning_functions, function(x) {
  names(x) = categories
  return(x)
})

set.seed(34223)
plot_binning_functions(binning_functions, categories)
ggsave(
  file.path(FIGURES_PATH, "binning_functions_toy.pdf"),
  height = 4,
  width = 5
)

binning_functions = simulate_category_mappings(2, c(6, 2), c(replicate(4, c(rep(1, 5), 2), simplify = FALSE), replicate(2, rep(2, 6), simplify = FALSE)))
names(binning_functions$inv) = paste0("Dataset ", 1:6)
set.seed(11)
plot = plot_binning_functions(binning_functions$inv, binning_functions$categories) + guides(x = guide_axis(angle = 0))
plot$data$category = paste0(LETTERS[as.numeric(substr(as.character(plot$data$category), 1, 1))], substr(as.character(plot$data$category), 2, 2))
ggsave(
  file.path(FIGURES_PATH, "binning_functions_simulation.pdf"),
  plot,
  height = 4,
  width = 8
)
