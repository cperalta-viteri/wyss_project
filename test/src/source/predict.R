library(torch)
library(SummarizedExperiment)
library(ggplot2)
library(Seurat)
library(cowplot)



source("lib_cell_classifier.R")
source("/lib_torch_dataset.R")
#source("/home/dnf/Desktop/final_model_julien_torch/allen20/src/julien/lib_allen.R") #removed

device <- if (cuda_is_available()) "cuda" else "cpu"

#counts<-x@assays$RNA@counts

predict_id<-function(x){
  
  dir.create(paste0(x@project.name))
  
  DefaultAssay(x) <- "RNA"
  
batch_predict <- function(m,x) {
    dl <- dataloader(x,batch_size=1024L,num_workers=0L,worker_packages = c("HDF5Array","SummarizedExperiment"))
    g <- coro::generator(function() for (o in dl) {print(".");yield(m(o$x))})
    P <- coro::collect(g())
    P <- torch_cat(P)
    P <- as_array(P)
    colnames(P) <- m$class_names
    P
  }
  
  
  # Model usage
  model2_log <- readRDS("/home/dnf/Desktop/crb_integration/crb_int/out/model/Allen_only_SS_final_model_prune_log.rds")
  plot(model2_log$log)
  plot(model2_log$val_log)
  
  model2 <- torch_load("/home/dnf/Desktop/crb_integration/crb_int/out/model/Allen_only_SS_final_model_prune.torch",device = "cpu")
  model2$eval()
  
  counts<-x@assays$RNA@counts
  e <- SummarizedExperiment(
    rbind(0,counts)[1L+match(model2$feature_names,rownames(x),nomatch = 0L),],
    colData = x@meta.data
  )
  e$P <- batch_predict(model2,se_dataset(e))
  e$pred <- colnames(e$P)[max.col(e$P)]
  e$pred_conf <- rowQ(-e$P,2) - rowQ(-e$P,1)
  
  
  predictions<-e$pred
  
  predictions<-data.frame(predictions=predictions,id_cell = colnames(x))
  
  
  predictions<-predictions[match(colnames(x),predictions$id_cell),]
  
  
  
  x$predictions<-predictions$predictions
  
  p1<-DimPlot(x,group.by = "predictions",raster=FALSE,label=TRUE,label.size = 2)+ggtitle("predictions")+NoLegend()
  p1
  
  x
  
}
