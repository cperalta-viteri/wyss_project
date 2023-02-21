


source("../01_data_collection/lib_sce.R")
x <- loadHDF5SingleCellExperiment("../01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5")
x <- logNormCounts(x)


# predict and test performances
p <- torch_predict_cell_type(logcounts(x),"model_pruned_best_train_acc.torch")
table(x$cluster_label,colnames(p)[max.col(p)])









