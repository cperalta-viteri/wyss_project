#############Geom flat violin
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggnewscale)
library(reticulate)
library(scales)
library(ggplot2)
library(Seurat)
#library(DropSeq.util) # not able to install/ Ludov said it's not used
library(RColorBrewer)
library(wesanderson)

#first create a .Renviron file in directory
# containing : RETICULATE_PYTHON=path_to_venv
py_config()


theme_set(theme_classic())



source('src/source/Source_new_cellHasR.R')
source('src/source/Source_allen.R')
#source('src/source/predict.R') # we need lib_allen source file


dir.create("out/1-Loading/")
dir.create("pdf/1-Loading/")




Telley_muA<-local({
  Telley_muA_mrna<-Read10X_h5("data/muA/filtered_feature_bc_matrix.h5")
  
  htos<-Read10X_h5("data/muA/muAh/raw_feature_bc_matrix.h5")
  
  
  mRNA<-Telley_muA_mrna$`Gene Expression`

  #mRNA<-mRNA[,colSums(mRNA>=1)>100]
  
  print(rowMeans(htos))
  
  # Select cell barcodes detected by both RNA and HTO In the example datasets we have already
  # filtered the cells for you, but perform this step for clarity.
  joint.bcs <- intersect(colnames(mRNA), colnames(htos))
  
  # Subset RNA and HTO counts by joint cell barcodes
  mRNA <- mRNA[, joint.bcs]
  htos <- as.matrix(htos[, joint.bcs])
  
  # Confirm that the HTO have the correct names
  rownames(htos)
  # Setup Seurat object
  
  #hashtag <- CreateSeuratObject(counts = mRNA)
  
  
  hashtag<-CreateSeuratObject(count=mRNA)
  
  hashtag[["HTO"]] <- CreateAssayObject(counts = htos)
  hashtag
  
})


Telley_muA <- Qc(Telley_muA,
                 HTO_tag=c("HTO-1","HTO-2","HTO-3","HTO-4","HTO-5","HTO-6","HTO-7","HTO-8","HTO-9","HTO-10"),
                 projet="Telley_muA")

# requires cellhashR (problematic)
# devtools::install_github(repo = 'bimberlab/cellhashR', ref = 'master', dependencies = TRUE, upgrade = 'always')
# requires python


Telley_muA<-local({

  Telley_muA$region<-NA
  Telley_muA$region[Telley_muA$consensuscall=="HTO-1"]<-"VA"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-2"]<-"VP"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-3"]<-"VC"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-4"]<-"H1A"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-5"]<-"H1P"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-6"]<-"H1C"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-7"]<-"H2A"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-8"]<-"H2P"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-9"]<-"H3"
  Telley_muA$region[Telley_muA$consensuscall=="HTO-10"]<-"DCN"
  
  
  Telley_muA$ML<-NA
  Telley_muA$ML[grep("V",Telley_muA$region)]<-"Vermis"
  Telley_muA$ML[grep("H1",Telley_muA$region)]<-"Paravermis"
  Telley_muA$ML[grep("H2",Telley_muA$region)]<-"Hemisphere_lat"
  Telley_muA$ML[grep("H3",Telley_muA$region)]<-"Hemisphere_lat"
  Telley_muA$ML[grep("DCN",Telley_muA$region)]<-"DCN"
  
  table(Telley_muA$ML)
  
  
  Telley_muA$AP<-NA
  
  Telley_muA$AP[grep("A",Telley_muA$region)]<-"Anterior"
  Telley_muA$AP[grep("P",Telley_muA$region)]<-"Posterior"
  Telley_muA$AP[grep("C",Telley_muA$region)]<-"Central"
  Telley_muA$AP[grep("H3",Telley_muA$region)]<-"Central"
  Telley_muA$AP[grep("DCN",Telley_muA$region)]<-"DCN"
  Telley_muA$sex<-NA
  Telley_muA$age<-"P0"
  Telley_muA$strain<-"CD1"
  
  
  
  plot_save_umap(Telley_muA,genes=c("Slc1a3","Calb1","Pax2","Rbfox3","Lhx9","Sox10","Pvalb","Ret","Nrgn","Calb2","Eomes"))
  
  
  saveRDS(Telley_muA, file = "out/1-Loading/Telley_muA/Telley_muA_all.Rds")
  
  Telley_muA
})


Telley_muA$htodemux

