#Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK=1)
#Sys.getenv("PYTORCH_ENABLE_MPS_FALLBACK")


source("notebooks/01_data_collection/lib_sce.R")


make_celltype_ds <- function() {
  # Load 10x and smart-seq cells
  allen10x <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/mouse_ctx_10x.h5");allen10x$source <- "10x"
  allenSSq <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5");allenSSq$source <- "ssv4"
  
  # Subset to common genes only
  #gplots::venn(list(SS=rownames(allenSSq),"10x"=rownames(allen10x)))
  allen10x <- allen10x[intersect(rownames(allen10x),rownames(allenSSq))]
  allenSSq <- allenSSq[rownames(allen10x)]
  
  # Merge 10x and SSq cells
  allen <- SingleCellExperiment(
    list(counts=cbind(counts(allen10x),counts(allenSSq))),
    colData = local({
      mcol <- setdiff(unique(c(names(colData(allen10x)),names(colData(allenSSq)))),"")
      m10x <- cbind(NA_character_,colData(allen10x))[match(mcol,names(colData(allen10x)),nomatch = 0) + 1L]
      mSSq <- cbind(NA_character_,colData(allenSSq))[match(mcol,names(colData(allenSSq)),nomatch = 0) + 1L]
      names(m10x) <- names(mSSq) <- mcol
      rbind(m10x,mSSq)
    })
  )
  allen$source <- factor(allen$source)
  
  # Compute target labels to predict
  allen$y_true <- paste(allen$class_label,allen$neighborhood_label,allen$subclass_label,sep=" --- ")
  allen$y_true <- factor(allen$y_true)
  stopifnot(all(!is.na(allen$y_true)))
  
  # Split cells into 3 datasets (trn/val/tst)
  set.seed(12345L)
  allen$mlset <- sample(factor(c("trn","val","tst"),c("trn","val","tst")),size = ncol(allen),replace = TRUE,prob = c(0.7,0.15,0.15))

  
  # Balance train and validation sets
  local({
    val_ds <- allen[,allen$mlset %in% "val"]
    val_ds <- val_ds[,sample(seq(ncol(val_ds)),size = 3000,replace = TRUE,prob = 1/(1+table(val_ds$y_true,val_ds$source))[cbind(val_ds$y_true,val_ds$source)])]
    table(val_ds$y_true,val_ds$source)
    val_ds <- scuttle::logNormCounts(val_ds)
    val_ds <- val_ds[,sample(ncol(val_ds))]
    list(
      x = t(logcounts(val_ds)) |> as("Matrix"),
      y = val_ds$y_true
    ) |> saveRDS("notebooks/celltype_prediction/allen_subclass_val.rds")
  })
  
  local({
    trn_ds <- allen[,allen$mlset %in% "trn"]
    trn_ds <- trn_ds[,sample(seq(ncol(trn_ds)),size = 15000,replace = TRUE,prob = 1/(1+table(trn_ds$y_true,trn_ds$source))[cbind(trn_ds$y_true,trn_ds$source)])]
    table(trn_ds$y_true,trn_ds$source)
    trn_ds <- scuttle::logNormCounts(trn_ds)
    trn_ds <- trn_ds[,sample(ncol(trn_ds))]
    list(
      x = t(logcounts(trn_ds)) |> as("Matrix"),
      y = trn_ds$y_true
    ) |> saveRDS("notebooks/celltype_prediction/allen_subclass_trn.rds")
  })
}


