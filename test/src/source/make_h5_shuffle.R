library(Seurat)
library(HDF5Array)
library(torch)


local({
  write_h5_dataset<- function(e,h5f,chunkdim=c(nrow(e),1L),level=1L) {
    e[[1]]->e_count_scale
    e[[2]]->e_meta
    
    h5write(colnames(e_count_scale),h5f,"/row_names")
    h5write(rownames(e_count_scale),h5f,"/col_names")
    #h5write(col_means_train,h5f,"/col_means_train")
    #h5write(col_stddevs_train,h5f,"/col_stddevs_train")
    h5createGroup(h5f,"/col_meta")
    h5write(e_meta$region,h5f,"/col_meta/region")
    h5write(e_meta$ML,h5f,"/col_meta/ML")
    h5write(e_meta$AP,h5f,"/col_meta/AP")
    h5write(e_meta$age,h5f,"/col_meta/age")
    h5write(e_meta$strain,h5f,"/col_meta/strain")
    h5write(e_meta$datasets,h5f,"/col_meta/datasets")
    h5write(e_meta$chemistry,h5f,"/col_meta/chemistry")
    h5write(e_meta$Phase,h5f,"/col_meta/Phase")
    h5write(as.vector(e_meta$integrated_snn_res.0.4),h5f,"/col_meta/cluster")
    h5write(rownames(e_meta),h5f,"/col_meta/rownames_meta")
    
    writeHDF5Array(DelayedArray(e_count_scale),h5f,"/x",chunkdim=chunkdim,level=level)
  }
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Make a randomized train/valid/test sets of 10000 cells
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  make_shuffled_datasets <- function() {
    # Load the data object
    e<-readRDS(file = "out/3-Purkinje_analysis/Purkinje_Telley/Telley_All_PJ_integratedl.Rds")
    
    e@project.name<-"Telley_All_PJ"
    
    # Randomize cells
    set.seed(12345)
    sample<-sample(ncol(e))
    e_count<-e@assays$RNA@counts[,sample]
    
    e_meta<-e@meta.data[sample,]
    write_h5_dataset(list(e_count,e_meta),h5f="out/3-Purkinje_analysis/Purkinje_Telley/Telley_All_PJ_integratedl.h5",level=1L,chunkdim = c(nrow(e_count),1L))
    
  }
  make_shuffled_datasets()
})





library(HDF5Array)
source("src/julien/draft/torch/lib_cell_classifier.R")
source("src/julien/draft/torch/lib_torch_dataset.R")

read_from_h5_generic <- function(h5_file) {
  meta <- as.data.frame(h5read(h5_file,"col_meta"))
  rownames(meta) <-   as.vector(h5read(h5_file,"row_names"))
  
  
  x <- t(HDF5Array(h5_file,"x",type="integer"))
  colnames(x) <- as.vector(h5read(h5_file,"col_names"))
  rownames(x) <- as.vector(h5read(h5_file,"row_names"))
  
  e <- SummarizedExperiment(
    t(x),
    colData = meta[rownames(x),,drop=FALSE],
  )
}






make_ds <- function() {
  # Load 10x and smart-seq cells
  Arlotta <- read_from_h5_generic("data/datasets/Arlotta_Bandler_model_shuffle.h5")
  #x <- readRDS("data/Natalia_AllDataQC_sumexp_2021.rds")
  
  # Subset to common genes only
  #gplots::venn(list(SS=rownames(allenSSq),"10x"=rownames(Arlotta),x=rownames(x)))
  #Arlotta <- Arlotta[intersect(rownames(Arlotta),rownames(allenSSq))]
  #allenSSq <- allenSSq[rownames(Arlotta)]
  
  # Compute target labels on 10x dataset
  Arlotta$y_true <-Arlotta$class_label
  
  #Arlotta$y_true[!Arlotta$class_label%in%c("Glutamatergic","CR")]<-NA
  Arlotta$y_true[Arlotta$y_true %in% c("","|")] <- NA
  Arlotta$y_true[!Arlotta$y_true %in% names(which(table(Arlotta$y_true)>=15))] <- NA  # remove cell-types with few cells
  Arlotta$y_true <- factor(Arlotta$y_true)
  
  
  
  
  
  # Create a training/validation set
  set.seed(12345L)
  
  splt <- sample(c("trn","val"),size=ncol(Arlotta),replace=TRUE,prob = c(0.7,0.3))
  trn_ds <- Arlotta[,splt=="trn"]
  trn_ds <- trn_ds[,balanced_bootstrap(trn_ds$y_true,2000)]
  val_ds <- Arlotta[,splt=="val"]
  val_ds <- val_ds[,balanced_bootstrap(val_ds$y_true,500)]
  
  assay(trn_ds) <- as(assay(trn_ds),"Matrix")
  assay(val_ds) <- as(assay(val_ds),"Matrix")
  
  saveRDS(trn_ds,"out/datasets/Arlotta_trn_ds.rds")
  saveRDS(val_ds,"out/datasets/Arlotta_val_ds.rds")
}

trn_ds <- readRDS("out/datasets/Arlotta_trn_ds.rds")
val_ds <- readRDS("out/datasets/Arlotta_val_ds.rds")

trn_ds <- se_dataset(trn_ds,"y_true")
val_ds <- se_dataset(val_ds,"y_true")




DEVICE <- if (cuda_is_available()) "cuda" else "cpu"
trn_dl <- dataloader(trn_ds,batch_size=64L,num_workers=0L,worker_packages = c("HDF5Array","SummarizedExperiment"))
val_dl <- dataloader(val_ds,batch_size=64L,num_workers=0L,worker_packages = c("HDF5Array","SummarizedExperiment"))



model <- cell_classifier(rownames(trn_ds$o),levels(trn_ds$y_true))
# fit1_log <- model$fit(trn_dl,val_dl,num_epoch = 30L,device=DEVICE)
# torch_save(model,"out/Arlotta_final_model_unprune.torch")
# saveRDS(fit1_log,"out/Arlotta_final_model_unprune_log.rds")

invisible(model$compute_prune_mask())
fit2_log <- model$fit(trn_dl,val_dl,num_epoch = 100L,device=DEVICE)
torch_save(model,"out/cell_type_Bandler_Arlotta_final_model_prune.torch")
saveRDS(fit2_log,"out/cell_type_Bandler_Arlotta_final_model_prune_log.rds")


#plot(fit1_log$log)
plot(fit2_log$log)

#plot(fit1_log$val_log)
plot(fit2_log$val_log)








