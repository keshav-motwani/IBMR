DATA_PATH = "data/application_1/"
dir.create(DATA_PATH, recursive = TRUE)

data = CITEseqData::get_hao_3_prime_data(DATA_PATH)

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
registerDoMC(cores = 10)
fit = glmnet::cv.glmnet(x = X, y = Y, family = "multinomial", type.multinomial = "grouped", trace.it = 1, parallel = TRUE)
saveRDS(fit, "data/simulation/Hao_glmnet_fit.rds")
