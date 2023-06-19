library(patchwork)
library(ggplot2)

library(HDF5Array)
library(SingleCellExperiment)
library(cellhashR)
library(DropletUtils)
library(scater)
library(scuttle)
library(scales)
library(ggExtra)



#' @title plot_UMIvsGene
#' @description plot_UMIvsGene displays a scatter plot of the relations between log10(UMI) on X axis and log10(genes)  on Y axis. #' a linear regression is fitted to the data and displayed.
#' @export
#' @param x A SingleCellExperiment object. It has to contain in the colData: "sum", "detected".
#' @param colour_by A character vector to set the colour of the data points. It has to be a column in colData(x). 
#' @param size_by A character vector to set the size of the data points. It has to be a column in colData(x).
#' @return A  ggplot figure.
plot_UMIvsGene <- function(x, ...) {
  p <- plotColData(x, x="sum", y="detected", ...) +
    scale_x_log10() + scale_y_log10() + 
    scale_colour_gradient(aesthetics = c("colour"),low = "gray90", high = "black") +
    scale_size_area(max_size = 3)+
    geom_smooth(method = "lm", show.legend = FALSE)+
    xlab("nCount_RNA")+
    ylab("nFeature_RNA") +
    theme_classic() + 
    theme(legend.position = "left") +
    labs(colour="percent.mito") +
    guides(size=guide_legend("percent.ribo"))
  #p <- ggMarginal(p, type = "densigram",bins=100, fill="grey", alpha=.3)
  
  p
}


#' @title plot_QC_Violin
#' @description plot_QC_Violin generates a violin + scatter + box plot for the variable desired.
#' @export
#' @param x A SingleCellExperiment object.
#' @param y Character vector from colData information of the SingleCellExperiment object to be displayed.
#' Default to `sum`which is the total count of UMIs.
#' @param ylab Character vector to be set as the label to be displayed on the Y axis.
#' @param title Character vector to be set as the title of the plot.
#' @return A graphic display of the variable specified.
plot_QC_Violin <- function(x, y="sum", ylab="nCount_RNA", title="Number of UMIs", ...){
  plotColData(x, y=y, ...) +
    geom_point(size = 0.25, alpha = 0.5, position = position_jitter(width=0.25,seed = 12345)) +
    geom_boxplot(width = 0.3, outlier.size = 0.5, outlier.alpha = 0.8) +
    xlab(NULL) +
    ylab(ylab) + 
    ggtitle(title) +
    theme(axis.text.y=element_text(size = 14,hjust = 0.5), title = element_text(size = 14,hjust = 0.5))
}



#' @title plot_QC_scatter
#' @description plot_QC_scatter produces a fully customized scatter plot for two continuous variables. Furthermore, by passing 
#' thresholds for the metrics of interest, the function will add either a vertical or horizontal line, intercepting either x-axis
#' or y-axis at the specified threshold.
#' @export
#' @param sce_object A SingleCellExperiment object.
#' @param x Continuous variable to be displayed in the x-axis.
#' @param y Continuous variable to be displayed in the y-axis.
#' @param xthreshold Value at which display a vertical line, marking the specified threshold.
#' @param ythreshold Value at which display a horizontal line, marking the specified threshold.
#' @return A customized scatter plot of class ggplot.
plot_QC_scatter <- function(sce_object, x, y, xthreshold=NULL, ythreshold=NULL, xlab=x, ylab=y, ...) {
  p <- plotColData(sce_object, x = x, y = y, ...)
  if (!is.null(ythreshold)) {p <- p + geom_hline(aes(yintercept = ythreshold), colour = "red", linetype = 2)}
  if (!is.null(xthreshold)) {p <- p + geom_vline(aes(xintercept = xthreshold), colour = "red", linetype = 2)}
  p <- p + annotate("text",x=max(sce_object[[x]]),hjust=1,y=ythreshold,vjust=0,label=sum(sce_object[[y]]>=ythreshold),size=3)
  p <- p + annotate("text",x=max(sce_object[[x]]),hjust=1,y=ythreshold,vjust=1,label=sum(sce_object[[y]]<ythreshold),size=3)
  p <- p + annotate("text",x=xthreshold,hjust=1,y=max(sce_object[[y]]),vjust=1,label=sum(sce_object[[x]]>=xthreshold),size=3,angle=90)
  p <- p + annotate("text",x=xthreshold,hjust=1,y=max(sce_object[[y]]),vjust=0,label=sum(sce_object[[x]]<xthreshold),size=3,angle=90)
  p <- p + annotate("text",x=xthreshold,y=ythreshold,label=sum((sce_object[[x]]>=xthreshold) & (sce_object[[y]]>=ythreshold)),size=3,hjust=-1,vjust=0.5,angle=45)
  p <- p + annotate("text",x=xthreshold,y=ythreshold,label=sum((sce_object[[x]]<xthreshold) & (sce_object[[y]]<ythreshold)),size=3,hjust=2,vjust=0.5,angle=45)
  p <- p + annotate("text",x=xthreshold,y=ythreshold,label=sum((sce_object[[x]]<xthreshold) & (sce_object[[y]]>=ythreshold)),size=3,hjust=2,vjust=0.5,angle=-45)
  p <- p + annotate("text",x=xthreshold,y=ythreshold,label=sum((sce_object[[x]]>=xthreshold) & (sce_object[[y]]<ythreshold)),size=3,hjust=-1,vjust=0.5,angle=-45)
  p <- p + xlab(xlab) + ylab(ylab) 
  p
}




#' @title demultiplexer
#' @description demultiplexer calls cellhashR::GenerateCellHashingCalls. Filtering of the
#' matrix of counts is performed prior call of the function.
#' @export
#' @param sce_object A SingleCellExperiment object.
#' @param methods List containing the algorithms to be used for the demultiplexing. This includes
#' "multiseq", "gmm_demux", "htodemux", "dropletutils".
#' @param project Character vector to use to name the project under the directory of pdf/1-Loading. The PDFs output from GenerateCellHashigCalls will be saved to this directory. Default to "New_project".
#' @param filter_singlet Boolean parameter to specify whether the resulting assay should be filtered to keep only cell barcodes assigned as singlets. Default to FALSE.
#' @return A SingleCellExperiment object with information of consensus call and consensus call global stored in colData(sce_object).
demultiplexer <- function(sce_object,project,methods=c("multiseq","gmm_demux", "htodemux", "dropletutils"),filter_singlet=FALSE) {
  barcodeData <- as.matrix(counts(swapAltExp(sce_object, "hto")))
  
  pdf(file.path("pdf/1-Loading",project,"hto.pdf"), width = 20, height = 20)
  calls <- cellhashR::GenerateCellHashingCalls(
    barcodeMatrix = barcodeData,
    methods = methods,
    doTSNE = FALSE,
    majorityConsensusThreshold = 0.3
  )
  # create summary plots to contrast cells based on call-status. Help inform why some cells were not called
  SummarizeCellsByClassification(calls = calls, barcodeMatrix = barcodeData)
  dev.off()
  
  calls$cellbarcode <- gsub("[.]","-",calls$cellbarcode)
  stopifnot(setequal(colnames(sce_object), calls$cellbarcode))
  colData(sce_object)$consensus_call <- calls$consensuscall
  colData(sce_object)$consensus_call_global <- calls$consensuscall.global
  # interest in filtering by calls$consensuscall.global
  if (filter_singlet) {
    sce_object <- sce_object[, sce_object$Barcode %in% calls[ calls$consensuscall.global=="Singlet", ]$cellbarcode]
  }
  sce_object
}


#' @title do_vst_transform
#' @description do_vst_transform Calls sctransform::vst on the count matrix, asking for returning the correct UMI matrix.
#' The corrected UMI matrix is then log transformed. The matrix can be obtained as it is or within a sce_object object, by passing TRUE
#' to the argument `as.AltExp`.
#' @export
#' @param sce_object A SingleCellExperiment objec.
#' @param as.AltExp Boolean parameter to specify whether add the log-transformed count matrix as alternative experiment to the 
#' sce_object object.
#' @return When `as.AltExp` is TRUE (default), the function returns a DelayedMatrix. When `as.AltExp` is TRUE, a sce_object object is returned. 
do_vst_transform <- function(sce_object, as.AltExp=TRUE){
  vst <- sctransform::vst(as(assay(sce_object),'Matrix'), return_corrected_umi = TRUE)
  vst <- DelayedArray(log1p(vst$umi_corrected)/log(2))
  if (as.AltExp) {
    altExp(sce_object, "Gene_Expression_VST") <- SummarizedExperiment(list(count_vst=vst), colData = colData(sce_object))
    sce_object
  } else {
    vst
  }
}


#' @title do_AnPearson_transform
#' @description do_AnPearson_transform carries out Analytical Pearson Residuals Normalization for UMI count matrix,
#' doi.org/10.1186/s13059-021-02451-7 .
#' @export
#' @param sce_object A SingleCellExperiment object.
#' @param theta A numerical value, default to 100. Set the measure of overdispersion to the Poisson distribution.
#' @return A SCE object with the normalized matrix added as alternative experiment under the name `Gene_Expression_PearNorm`.
do_AnPearson_transform <- function(sce_object, theta=100){
  counts_sum0 <- rowSums(assay(sce_object))
  counts_sum1 <- colSums(assay(sce_object))
  counts_sum <- sum(sce_object$sum)
  
  mu <- DelayedArray(crossprod(t(counts_sum0), counts_sum1) / counts_sum)
  z <- DelayedArray((DelayedArray(assay(sce_object)) - mu) / sqrt(mu + (mu**2/theta)))
  z <- replace(z, is.na(z), c(0))
  z <- z[, order(match(colnames(sce_object[, sce_object$qc=="pass"]), colnames(z)))]
  altExp(sce_object, "Gene_Expression_PearNorm") <- SingleCellExperiment(list(count_AP=z), colData = colData(sce_object))
  sce_object
}




normalize_sce <- function(x,method = c("logNorm","SCTransform","AnPearson")) {
  method <- match.arg(method)
  x <- x[, x$qc %in% "pass"]
  x <- switch(method,
      "logNorm" = scuttle::logNormCounts(x),
      "SCTransform" = do_vst_transform(x),
      "AnPearson" = do_AnPearson_transform(x)
  )
  x
}

make_transcriptomic_qc_pdf <- function(sce,pdf_dir,suffix="",mito_max=15,min_nfeatures=200) {
  
  # Violin of each variable independently
  local({
    p1 <- plot_QC_Violin(sce, y="subsets_Mito_percent", ylab="percent_mito", title = "Mitochondial genes")
    p2 <- plot_QC_Violin(sce, y="subsets_Ribo_percent", ylab="percent_ribo", title = "Ribosomal genes")
    p3 <- plot_QC_Violin(sce, y="sum", ylab="nCount_RNA", title = "Number of UMI")
    p4 <- plot_QC_Violin(sce, y="detected", ylab="nFeature_RNA", title = "Number of genes")
    print(p1 + p2 + p3 + p4)
    ggplot2::ggsave(
      filename = file.path(pdf_dir, sprintf("01-violin%s.pdf",suffix)),
      width = 297, height = 210, units = "mm"
    )
  })
  
  # UMI vs Detected
  local({
    p <- plot_UMIvsGene(
      sce,
      colour_by = "subsets_Mito_percent",
      size_by = "subsets_Ribo_percent"
    )
    p <- ggExtra::ggMarginal(p,type = "densigram",bins=100, fill="grey", alpha=0.3)
    print(p)
    ggplot2::ggsave(plot=p,
      filename = file.path(pdf_dir, sprintf("02-UMI-vs-detected%s.pdf",suffix)),
      width = 297, height = 210, units = "mm"
    )
  })
  
  # Mitochondria and detection criteria
  local({
    p <- plot_QC_scatter(
      sce,
      x = "detected", y = "subsets_Mito_percent",
      xthreshold = min_nfeatures,ythreshold = mito_max,
      xlab="Num. Detected Gene",ylab="% mito"
    ) +
      annotate(geom="rect",xmin=min_nfeatures,xmax=+Inf,ymin=0,ymax=mito_max,fill="green",alpha=0.3)
    p <- ggExtra::ggMarginal(p,type = "densigram",bins=100, fill="grey",  alpha=0.3)
    print(p)
    ggplot2::ggsave(plot=p,
      filename = file.path(pdf_dir,sprintf("03-MT-vs-detected%s.pdf",suffix)),
      width = 297, height = 210, units = "mm"
    )
  })
}






transcriptomic_pipeline <- function(x,
               output_dir = "qc_out",
               mito_max = 15,
               min_nfeatures = 200,
               normalization_method = "logNorm",
               ...) {
  dir.create(output_dir)
  pdf_dir <- file.path(output_dir,"pdf")
  dir.create(pdf_dir)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Compute per cell QC metrics
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  mito.genes <- grep(pattern = "^mt-", x = rownames(x), value = TRUE, ignore.case = TRUE)
  ribo.genes <- grep(pattern = "^(Rpl|Rps|Mrp)", x = rownames(x), value = TRUE, ignore.case = TRUE)
  x <- scuttle::addPerCellQCMetrics(x, subset = list(Mito = mito.genes,Ribo = ribo.genes))

  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Compute QC status of each cell and filter at each step
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  x$qc <- "pass"
  make_transcriptomic_qc_pdf(x,pdf_dir,suffix="_00-before-qc")
  
  x$qc[(x$qc %in% "pass") & (x$detected<min_nfeatures)] <- "fails_minFeat"
  make_transcriptomic_qc_pdf(x[,x$qc %in% "pass"],pdf_dir,suffix="_01-after-detect")
  
  x$qc[(x$qc %in% "pass") & x$subsets_Mito_percent > mito_max] <- "fails_maxMito"
  make_transcriptomic_qc_pdf(x[,x$qc %in% "pass"],pdf_dir,suffix="_02-after-mito")
  
  # fitting linear model on UMI vs detected only on previously QC passing cells
  lm <- lm(detected~sum,data=data.frame(detected=log10(x$detected[x$qc %in% "pass"]),sum=log10(x$sum[x$qc %in% "pass"])))
  x$qc[(x$qc %in% "pass") & (log10(x$detected) - predict(lm,data.frame(sum=log10(x$sum))) < -0.09)] <- "fails_lm"
  print(table(x$qc))
  
  x <- x[,x$qc %in% "pass"]
  make_transcriptomic_qc_pdf(x,pdf_dir,suffix="_03-after-lm")
  
  
  # Normalize
  x <- normalize_sce(x,method = normalization_method)
  
  # Compute UMAP layout
  x <- scater::runUMAP(x, BPPARAM=BiocParallel::MulticoreParam(6L))

  # display UMAP and color by numer of detected genes
  p <- scater::plotUMAP(x[, x$qc=="pass"],point_size=.05, colour_by="detected") + 
    coord_equal() + 
    guides(colour=guide_legend(NULL, override.aes = list(size=2)))
  print(p)
  ggplot2::ggsave(
    filename = file.path(pdf_dir,"XX_UMAP_detected.pdf"),
    plot = p, width = 297, height = 210, units = "mm")
  return(x)
}




# This is the parking with code to eventually integrate to the pipeline
# PARKING
function() {
  if (isTRUE(hto)) {x <- demultiplexer(x,project,...)}
  

  if (length(mito.genes) >= 1){
    p2 <- plotColData(x,
                      x="sum", y="detected",
                      colour_by="subsets_Mito_percent",
                      size_by="subsets_Ribo_percent") +
      scale_x_log10() + scale_y_log10() + 
      scale_colour_gradient(aesthetics = c("colour"),low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method = "lm", show.legend = FALSE)+
      geom_hline(aes(yintercept = min_nfeatures), colour = "darkgreen", linetype = 2)+
      geom_abline(intercept = lm$coefficients[1] - 0.09, slope = lm$coefficients[2], color="orange") +
      xlab("nCount_RNA")+
      ylab("nFeature_RNA") +
      theme_classic() + 
      theme(legend.position = "left") +
      labs(colour="percent.mito") +
      guides(size=guide_legend("percent.ribo")) +
      annotate(geom = "text",size=4, label = paste0(dim(x[, !x$qc == "fails_maxMito" & !x$qc == "fails_minFeat"])[2], " QC passed cells"),
               x = 1000, y = 5000)
    
    p2 <- ggExtra::ggMarginal(p2, type = "densigram",bins=100, fill="grey", alpha=.3)
    
    ggplot2::ggsave(filename = paste0("pdf/1-Loading/",project,"/nGenesvsnUMIs_cutoff.pdf"),
           plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ######## Cells to exclude lie below an intercept offset of -0.09
    p3 <- plotColData(x[, !x$qc == "fails_maxMito" & !x$qc == "fails_minFeat"],
                      x="sum", y="detected",
                      colour_by = "qc") +
      scale_x_log10() + scale_y_log10() +
      scale_colour_manual(values = c("darkred", "dodgerblue")) +
      geom_smooth(method = "lm", color="black",show.legend = FALSE)+
      geom_abline(intercept = lm$coefficients[1] - 0.09, slope = lm$coefficients[2], color="orange") +
      xlab("log10(nCount_RNA)")+
      ylab("log10(nFeature_RNA)") +
      theme_classic() + 
      theme(legend.position = "left") +
      labs(colour="valid cells") +
      annotate(geom = "text",size=4,
               label = paste0(
                 table(x[, !x$qc == "fails_maxMito" & !x$qc == "fails_minFeat"]$qc)[[2]], " QC passed cells\n",
                 table(x[, !x$qc == "fails_maxMito" & !x$qc == "fails_minFeat"]$qc)[[1]], " QC filtered"),
               x = 2000, y = 10000)
    
    p3 <- ggExtra::ggMarginal(p3, type = "densigram",bins=100, fill="grey", alpha=.3)
    
    p3
    
    ggplot2::ggsave(filename = paste0("pdf/1-Loading/",project,"/nGenesvsnUMIs_validCells.pdf"),
           plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  } ### close if mito genes > 1
  
  if (!length(mito.genes) >= 1){
    
    ######## Gene/UMI scatter plot before filtering
    p1 <- plot_UMIvsGene(x,
                         colour_by = "subsets_Mito_percent",
                         size_by = "subsets_Ribo_percent")
    p1 <- p1 + geom_hline(aes(yintercept = log10(min_nfeatures)), colour = "darkgreen", linetype = 2)
    p1 <- ggExtra::ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
    
    ggplot2::ggsave(filename = paste0("pdf/1-Loading/",project,"/nGenesvsnUMIs_beforefiltering.pdf"),
           plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ####### Filter cells base on min features and lm
    # line 477
    lm <- lm(
      log10(x[, !x$qc == "fails_minFeat"]$detected) ~ 
        log10(x[, !x$qc == "fails_minFeat"]$sum))
    
    p2 <- plotColData(x[, !x$qc == "fails_minFeat"],
                      x="sum", y="detected",
                      colour_by="subsets_Mito_percent",
                      size_by="subsets_Ribo_percent") +
      scale_x_log10() + scale_y_log10() + 
      scale_colour_gradient(aesthetics = c("colour"),low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method = "lm", show.legend = FALSE)+
      geom_hline(aes(yintercept = min_nfeatures), colour = "darkgreen", linetype = 2)+
      geom_abline(intercept = lm$coefficients[1] - 0.09, slope = lm$coefficients[2], color="orange") +
      xlab("log10(nCount_RNA)")+
      ylab("log10(nFeature_RNA)") +
      theme_classic() + 
      theme(legend.position = "left") +
      labs(colour="percent.mito") +
      guides(size=guide_legend("percent.ribo")) +
      annotate(geom = "text",size=4, label = paste0(dim(x[, !x$qc == "fails_minFeat"])[2], " QC passed cells"),
               x = 1000, y = 5000)
    
    p2 <- ggExtra::ggMarginal(p2, type = "densigram",bins=100, fill="grey", alpha=.3)
    
    ggplot2::ggsave(filename = paste0("pdf/1-Loading/",project,"/nGenesvsnUMIs_cutoff.pdf"),
           plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    
    ######## Cells to exclude lie below an intercept offset of -0.09
    p3 <- plotColData(x[, !x$qc == "fails_minFeat"],
                      x="sum", y="detected",
                      colour_by = "qc") +
      scale_x_log10() + scale_y_log10() +
      scale_colour_manual(values = c("darkred", "dodgerblue")) +
      geom_smooth(method = "lm", color="black",show.legend = FALSE)+
      geom_abline(intercept = lm$coefficients[1] - 0.09, slope = lm$coefficients[2], color="orange") +
      xlab("log10(nCount_RNA)")+
      ylab("log10(nFeature_RNA)") +
      theme_classic() + 
      theme(legend.position = "left") +
      labs(colour="valid cells") +
      annotate(geom = "text",size=4,
               label = paste0(
                 table(x[, !x$qc == "fails_minFeat"]$qc)[[2]], " QC passed cells\n",
                 table(x[, !x$qc == "fails_minFeat"]$qc)[[1]], " QC filtered"),
               x = 2000, y = 10000)
    
    p3 <- ggExtra::ggMarginal(p3, type = "densigram",bins=100, fill="grey", alpha=.3)
    
    p3
    
    ggplot2::ggsave(filename = paste0("pdf/1-Loading/",project,"/nGenesvsnUMIs_validCells.pdf"),
           plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    
  }  ### close if mito genes < 1
  
  ##### filtering by solely valid cells "passed"
  p1 <- plot_QC_Violin(
    x[, x$qc == "pass"],
    y="sum", ylab="nCount_RNA", title = "Number of UMI")
  
  
  p2 <- plot_QC_Violin(
    x[, x$qc == "pass"],
    y="detected", ylab="nFeature_RNA", title="Number of Genes")
  
  if (length(mito.genes)>=1) {
    p3 <- plot_QC_Violin(
      x[, x$qc == "pass"],
      y="subsets_Mito_percent", ylab="mito percent", title = "Mitochondrial Genes")
    
    p4 <- plot_QC_Violin(
      x[, x$qc == "pass"],
      y="subsets_Ribo_percent", ylab = "ribo percent", title = "Ribosomal Genes")
    
    p <- cowplot::plot_grid(plotlist = list(p3, p4, p1, p2), ncol = 2, align = "hv",
                            rel_widths = c(1, 1, 1, 1),rel_heights =  c(1, 1, 1, 1))
    
    ggsave(filename = paste0("pdf/1-Loading/",project,"/violin_plot_QC_after_filtering.pdf"),
           plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p1 <- plot_UMIvsGene(x[, x$qc == "pass"],
                        colour_by = "subsets_Mito_percent",
                        size_by = "subsets_Ribo_percent")
    p1 <- ggExtra::ggMarginal(p1, type = "densigram",bins=100, fill="grey", alpha=0.3)
    ggplot2::ggsave(filename = paste0("pdf/1-Loading/",project,"/nUMIsVSnGenes_afterfiltering.pdf"),
                    plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
  }
  
  if(!length(mito.genes)>=1){
    ggsave(filename = paste0("pdf/1-Loading/",project,"/QC_after filtering_1.pdf"),
           plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ggsave(filename = paste0("pdf/1-Loading/",project,"/QC_after filtering_2.pdf"),
           plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
  }
  
  x
}




