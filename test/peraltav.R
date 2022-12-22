
library(HDF5Array)
library(Seurat)
library(SingleCellExperiment) # just for testing
library(scater)
library(cellhashR)
library(DropletUtils)

#DropletUtils::read10xCounts()

set.seed(12)

# rna-seq & atac-seq
h5ls("data/muA/raw_feature_bc_matrix.h5")
# hto tags
h5ls("data/muA/muAh/raw_feature_bc_matrix.h5")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ============= reading data ===============
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tenx <- TENxMatrix("data/muA/filtered_feature_bc_matrix.h5")
htos <- TENxMatrix("data/muA/muAh/raw_feature_bc_matrix.h5")
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## features x samples
## 'gene expresion: 32285 x 20000 sparse Matrix of class "dgCMatrix"
## 'peaks: 175091 x 20000 sparse Matrix of class "dgCMatrix"

## Expression data
gene_exp <- tenx[grep("ENS", rownames(tenx)),]
# setting gene names from feature group in h5 file. Filtering is done to get only genes.
row.names(gene_exp) <- as.character(HDF5Array("data/muA/filtered_feature_bc_matrix.h5", "/matrix/features/name")[
  grep("^G",HDF5Array("data/muA/filtered_feature_bc_matrix.h5", "/matrix/features/feature_type"))])

#'* I don't know what to do with the intervals *
ranges_mrna <- h5read("data/muA/filtered_feature_bc_matrix.h5", "/matrix/features/interval")[
  grep("^G",HDF5Array("data/muA/filtered_feature_bc_matrix.h5", "/matrix/features/feature_type"))]

## ATAC data
atac <- tenx[grep("ENS", rownames(tenx), invert = TRUE),]
# setting gene names from feature group in h5 file. Filtering is done to get only peaks.
rownames(atac) <- as.character(HDF5Array("data/muA/filtered_feature_bc_matrix.h5", "/matrix/features/name")[nrow(gene_exp)+1:nrow(atac)])


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ============= subsetting by HTO tags ===============
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

length(intersect(colnames(gene_exp), colnames(htos))) #20.000

gene_exp <- gene_exp[, intersect(colnames(gene_exp), colnames(htos))]
htos <- htos[, intersect(colnames(gene_exp), colnames(htos))]


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ============= SingleCellExperiment Object ===============
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sce <- SingleCellExperiment(assays = list(counts=gene_exp),
                            colData=list(t(htos))) # IDN why it adds a new column: "sizeFactor"
                           
sce <- scuttle::logNormCounts(sce)

#'* It's possible to add per-cell QC metrics *
#sce <- addPerCellQCMetrics(sce, subset=list(Mito=rownames(rowData(sce))[grep("mt-", rownames(rowData(sce)))]))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ============= UMAP =============== 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this worked but I had to incude logNormCounts
## umap does not accepts HDF5Array products

#sce <- runUMAP(sce) 

plotUMAP(sce)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ============= Demultiplexing ===============
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(ggplot2)
require(ggExtra)
require(dplyr)
require(cowplot)
require(Matrix)
require(gridExtra)
require(ggrepel)
require(ggpubr)
require(ggnewscale)
require(cellhashR)

reticulate::py_config()
HTO_tag=c("HTO.1","HTO.2","HTO.3","HTO.4","HTO.5","HTO.6","HTO.7","HTO.8","HTO.9","HTO.10")

# getting a list of cellbarcode
fileConn <- file(file.path("peralta_qc/out/cellbarcodeWhitelist.txt"))
writeLines(gsub("-", ".", colnames(assay(sce))), fileConn)
close(fileConn)

# writting table of hto counts filtered by HTO tag
write.table(t(colData(sce)[, colnames(colData(sce)) %in% HTO_tag]),
            file = "peralta_qc/out/hto.txt", sep = ",", col.names = TRUE, row.names = FALSE)

write.table(t(barcode_df),
            file = "peralta_qc/out/hto.txt", sep = ",", col.names = TRUE, row.names = TRUE)

pdf(file.path("peralta_qc/pdf/hto.pdf"), width = 20, height = 20)

barcodeData <- ProcessCountMatrix(rawCountData = "peralta_qc/out/hto.txt",
                                  minCountPerCell = 5,
                                  simplifyBarcodeNames = FALSE,
                                  cellbarcodeWhitelist = "peralta_qc/out/cellbarcodeWhitelist.txt",
                                  barcodeWhitelist = HTO_tag)

methods=c("multiseq","gmm_demux", "htodemux", "dropletutils")
calls <- GenerateCellHashingCalls(barcodeMatrix = barcodeData,
                                  methods = methods, 
                                  doTSNE = FALSE,
                                  majorityConsensusThreshold = 0.3)
dev.off()











