


source("notebooks/01_data_collection/lib_sce.R")
x <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5")
x <- scuttle::logNormCounts(x)


# predict and test performances
source("notebooks/celltype_prediction/lib_cell_classifier.R")
p <- torch_predict_cell_type(logcounts(x),"model_pruned_best_train_acc.torch")
image(table(x$subclass_label,colnames(p)[max.col(p)]))




