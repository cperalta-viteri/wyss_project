


source("notebooks/01_data_collection/lib_sce.R")
source("notebooks/celltype_prediction/lib_cell_classifier.R")


x <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5")
x <- scuttle::logNormCounts(x)

x$p1 <- torch_predict_cell_type(counts(x),"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.torch")
x$p2 <- torch_predict_cell_type(logcounts(x),"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.torch")

x$pred1 <- colnames(x$p1)[max.col(x$p1)]
x$pred2 <- colnames(x$p2)[max.col(x$p2)]

x$y_true <- paste(x$class_label,x$neighborhood_label,x$subclass_label,sep=" --- ")
mean(x$pred1 != x$y_true)
mean(x$pred2 != x$y_true)



