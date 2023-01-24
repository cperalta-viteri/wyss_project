

QC <- function(sce_object,
               organism="mouse",
               hto=TRUE,
               methods=c("multiseq","gmm_demux", "htodemux", "dropletutils"),
               project="New_projec",
               HTO_tag=NULL,
               filter_singlet=FALSE,
               mito_max=15,
               min_nfeatures=log10(200)){
  
  #' @title demultiplexer
  #' @description demultiplexer calls cellhashR::GenerateCellHashingCalls. Filtering of the
  #' matrix of counts is performed prior call of the function.
  #' @export
  #' @param sce_object A SingleCellExperiment object.
  #' @param methods List containing the algorithms to be used for the demultiplexing. This includes
  #' "multiseq", "gmm_demux", "htodemux", "dropletutils".
  #' @param project Character vector to use to name the project under the directory of pdf/1-Loading. The PDFs output from GenerateCellHashigCalls will be saved to this directory. Default to "New_project".
  #' @param HTO_tag List of character strings containing the name of the HTO tags to keep, the white-list
  #' @param filter_singlet Boolean parameter to specify whether the resulting assay should be filtered to keep only cell barcodes assigned as singlets. Default to FALSE.
  #' @return A SingleCellExperiment object with information of consensus call and consensus call global stored in colData(sce_object).
  
  demultiplexer <- function(sce_object){
    require(cellhashR)
    
    if (!file.exists(file.path("pdf/1-Loading",project))){
      dir.create(file.path("pdf/1-Loading",project))
    }
    
    
    barcode_df <- as.data.frame(
      assay(swapAltExp(sce_object, "hto"))[rownames(assay(swapAltExp(sce_object, "hto"))) %in% HTO_tag,])
    # demultiplexing
    pdf(file.path("pdf/1-Loading",project,"hto.pdf"), width = 20, height = 20)
    calls <- GenerateCellHashingCalls(barcodeMatrix = barcode_df,
                                      methods = methods, 
                                      doTSNE = FALSE,
                                      majorityConsensusThreshold = 0.3)
    # create summary plots to contrast cells based on call-status. Help inform why some cells were not called
    SummarizeCellsByClassification(calls = calls, barcodeMatrix = barcode_df)
    dev.off()
    
    stopifnot(setequal(colnames(sce_object), calls$cellbarcode))
    colData(sce_object)$consensus_call <- calls$consensuscall
    colData(sce_object)$consensus_call_global <- calls$consensuscall.global
    # interest in filtering by calls$consensuscall.global
    if(filter_singlet==TRUE){
      sce_object <- sce_object[, sce_object$Barcode %in% calls[ calls$consensuscall.global=="Singlet", ]$cellbarcode]
    }
    sce_object
  }
  
  #' @title plot_QC_Violin
  #' @description plot_QC_Violin generates a violin + scatter + box plot for the variable desired.
  #' @export
  #' @param sce_object A SingleCellExperiment object.
  #' @param y Character vector from colData information of the SingleCellExperiment object to be displayed.
  #' Default to `sum`which is the total count of UMIs.
  #' @param ylab Character vector to be set as the label to be displayed on the Y axis.
  #' @param title Character vector to be set as the title of the plot.
  #' @return A graphic display of the variable specified.
  
  plot_QC_Violin <- function(sce_object, y="sum", ylab="nCount_RNA", title="Number of UMIs", ...){
    p <- plotColData(sce_object, y=y, ...) +
      geom_point(size=0.25, alpha=0.5, position = position_jitter(width=0.25,seed = 12345))+
      geom_boxplot(width=0.025, outlier.size = 0.5, outlier.alpha = 0.8)+
      xlab(NULL)+
      ylab(ylab)+ 
      ggtitle(title)+
      theme(axis.text.y=element_text(size = 14,hjust = 0.5), title = element_text(size = 14,hjust = 0.5))
    p
  }
  
  #' @title plot_UMIvsGene
  #' @description plot_UMIvsGene displays a scatter plot of the relations between log10(UMI) on X axis and log10(genes)  on Y axis. #' a linear regression is fitted to the data and displayed.
  #' @export
  #' @param sce_object A SingleCellExperiment object. It has to contain in the colData: "sum", "detected".
  #' @param colour_by A character vector to set the colour of the data points. It has to be a column in colData(sce_object). 
  #' @param size_by A character vector to set the size of the data points. It has to be a column in colData(sce_object).
  #' @return A  ggplot figure.
  
  plot_UMIvsGene <- function(sce_object, colour_by = NULL, size_by=NULL, ...){
    
    p <- plotColData(sce_object, x="sum", y="detected",
                     colour_by=colour_by,
                     size_by=size_by, ...) +
      scale_x_log10() + scale_y_log10() + 
      scale_colour_gradient(aesthetics = c("colour"),low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method = "lm", show.legend = FALSE)+
      xlab("log10(nCount_RNA)")+
      ylab("log10(nFeature_RNA)") +
      theme_classic() + 
      theme(legend.position = "left") +
      labs(colour="percent.mito") +
      guides(size=guide_legend("percent.ribo"))
    
    #p <- ggMarginal(p, type = "densigram",bins=100, fill="grey", alpha=.3)
    
    p
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
  
  plot_QC_scatter <- function(sce_object, x, y, xthreshold=NULL, ythreshold=NULL, ...){
    p <- plotColData(sce_object, x=x, y=y, ...)
    if (!is.null(ythreshold)){p <- p + geom_hline(aes(yintercept = ythreshold), colour = "red", linetype = 2)}
    if(!is.null(xthreshold)){p <- p + geom_vline(aes(xintercept = xthreshold), colour = "red", linetype = 2)}
    if(!is.null(c(xthreshold, ythreshold))) {
      annotations <- data.frame(
        xpos = c(-Inf,-Inf,Inf,Inf),
        ypos =  c(-Inf, Inf,-Inf,Inf),
        annotateText = c(
          as.character(table(sce_object$subsets_Mito_percent <= ythreshold & sce_object$detected <= xthreshold)[2]), #BL:6
          as.character(table(sce_object$subsets_Mito_percent >= ythreshold & sce_object$detected <= xthreshold)[2]), #TL:5
          as.character(table(sce_object$subsets_Mito_percent <= ythreshold & sce_object$detected >= xthreshold)[2]), #BR: 19547
          as.character(table(sce_object$subsets_Mito_percent >= ythreshold & sce_object$detected >= xthreshold)[2])), #TR: 442
        hjustvar = c(-1,-1,5,10),
        vjustvar = c(-5,10,-5,10)
      )
      p <- p + geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
    }
    p
  }
  
  
  #' @title set_QC_status
  #' @description set_QC_status Created a column in colData indicating whether cells passed or failed three QC criterias: maximum #' mitochondrial gene percentage, minimum features (detected RNA) and laying above a linear regression line with an offset of 
  #' -0.09. For the last criteria, a linear model is fitted using as formula log10(features) ~ log10(UMIs)
  #' @export
  #' @param sce_object A SingleCellExperimet object containing QC metrics computed with scuttle::addPerCellQCMetrics. For more
  #' information ?scuttle::addPerCellQCMetrics .
  #' @param mito_percent_column Character string indicating the name of the column in colData. This column is added by
  #' scuttle::addPerCellQCMetrics.
  #' @param mito_max Numerical value indicating the threshold for the maximum percentage of mitochondrial genes. 
  #' Default to 15.
  #' @param min_nfeatures Numerical value indicating the threshold for the minimum of features. Default to log10(200).
  #' @return
  
  set_QC_status <- function(sce_object=NULL,
                            mito_max=15,
                            min_nfeatures=log10(200)){
    if(isTRUE("subsets_Mito_percent" %in% colnames(colData(sce)))){
      
      # Assessing what cells pass the threshold for mitochondrial gene percentage
      # retrieveCellInfo looks into colData for the column passed and retrieves a list. The second element in the list contains the numerical values.
      colData(sce_object)$qc <- gsub(
        "TRUE", "passed", gsub("FALSE", "fails_maxMito", retrieveCellInfo(sce_object, by="subsets_Mito_percent")[[2]] <= mito_max))
      
    }
    
    
    # filtering by passed cells from mito_percent criteria and assess whether those cells passe or no the min_feature criteria.
    colData(sce_object)$qc <-
      replace(x=colData(sce_object)$qc,
              list=which(sce_object$qc == "passed"),
              values=gsub(
                "TRUE", "passed", gsub("FALSE", "fails_minFeat",
                                       log10(sce_object[,sce_object$qc == "passed"]$detected)>=min_nfeatures)))
    
    # fitting linear model with the formula log10(features) ~ log10(UMI)
    lm <- lm(log10(sce_object[, sce_object$qc == "passed"]$detected) ~ log10(sce_object[, sce_object$qc == "passed"]$sum))
    
    # Assessing what cells passed the criteria of laying above a linear model with a offset of -0.09.
    colData(sce_object)$qc <-
      replace(x=colData(sce_object)$qc,
              list=which(sce_object$qc =="passed"),
              values=gsub("TRUE", "passed", gsub("FALSE", "fails_lm", log10(sce_object[, sce_object$qc == "passed"]$detected) > 
                                                   log10(sce_object[, sce_object$qc == "passed"]$sum) * lm$coefficients[2] + (lm$coefficients[1]-.09))))
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
    } else {vst}
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
    z <- DelayedArray((assay(sce_object) - mu) / sqrt(mu + (mu**2/theta)))
    z <- replace(z, is.na(z), c(0))
    z <- z[, order(match(colnames(sce_object[, sce_object$qc=="passed"]), colnames(z)))]
    altExp(sce_object, "Gene_Expression_PearNorm") <- 
      SingleCellExperiment(list(count_AP=z), colData = colData(sce_object))
    sce_object
  }
  
  
  #### calls ###
  
  if (isTRUE(hto)){sce_object <- demultiplexer()}
  
  switch (organism,
    mouse = sce_object <- addPerCellQCMetrics(
      sce_object, subset=list(Mito=grep("mt-", rownames(sce_object),value=TRUE),
                             Ribo=grep("(^Rpl|^Rps|Mrp)", rownames(sce_object), value = TRUE))),
    human = sce_object <- addPerCellQCMetrics(
      sce_object, subset=list(Mito=grep("MT-", rownames(sce_object),value=TRUE),
                             Ribo=grep("(^Rpl|^Rps|Mrp)", rownames(sce_object), value = TRUE)))
  )
  
  mito.genes <- grep(pattern = "^mt-", x = rownames(sce_object), value = TRUE, ignore.case = TRUE)
  
  p3 <- plot_QC_Violin(sce_object, y="sum", ylab="nCount_RNA", title = "Number of UMI")
  p4 <- plot_QC_Violin(sce_object, y="detected", ylab="nFeature_RNA", title = "Number of genes")
  
  
  if (length(mito.genes) >= 1){
    p1 <- plot_QC_Violin(sce_object, y="subsets_Mito_percent", ylab="percent_mito", title = "Mitochondial genes")
    p2 <- plot_QC_Violin(sce_object, y="subsets_Ribo_percent", ylab="percent_ribo", title = "Ribosomal genes")
    
    p <- cowplot::plot_grid(plotlist = list(p1, p2, p3, 4),
                            ncol=2,
                            align = "hv",
                            rel_widths = c(1, 1, 1, 1),
                            rel_heights =  c(1, 1, 1, 1))
    ggplot2::ggsave(filename = file.path("pdf/1-Loading",project, "violin_plot_QC_before_filtering.pdf"),
                    plot = p, device = cairo_pdf, width = 297, height = 210, units = "mm")
    
    

  
  #p <- plot_UMIvsGene(sce_object)
  #p <- ggExtra::ggMarginal(p, type = "densigram",bins=100, fill="grey", alpha=0.3)
  
    if(isTRUE(hto)){
      if (sum(sce_object$subsets_Mito_percent <= mito_max) == length(sce_object$subsets_Mito_percent <= mito_max)){
        p <- plotColData(sce_object, x="detected", y="subsets_Mito_percent",
                         color_by = "consensus_call_global")+
          geom_hline(aes(yintercept = mito_max), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          xlab("nFeature_RNA")+
          ylab("percent.mito")+
          annotate(geom = "text",size=4, 
                   label = paste0(as.numeric(table( sce_object$subsets_Mito_percent <= mito_max)[1])," cells remain"),
                   x = max(sce_object$detected)-1000, y = mito_max+5)
        p <- ggExtra::ggMarginal(p, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"),
               plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
      }
      
      if(sum(sce_object$subsets_Mito_percent > mito_max)==length(sce_object$subsets_Mito_percent <= mito_max)){
        p <- plotColData(sce_object, x="detected", y="subsets_Mito_percent",
                         color_by = "consensus_call_global")+
          geom_hline(aes(yintercept = mito_max), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          xlab("nFeature_RNA")+
          ylab("percent.mito")+
          annotate(geom = "text",size=4, 
                   label = paste0(as.numeric(table( sce_object$subsets_Mito_percent <= mito_max)[1])," cells removed"),
                   x = max(sce_object$detected)-1000, y = mito_max+5)
        p <- ggExtra::ggMarginal(p, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"),
               plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
      }
      
      if(sum(sce_object$subsets_Mito_percent <= mito_max)!=length(sce_object$subsets_Mito_percent <= mito_max)&
         sum(sce_object$subsets_Mito_percent > mito_max)!=length(sce_object$subsets_Mito_percent <= mito_max)){
        
        p <- plotColData(sce_object, x="detected", y="subsets_Mito_percent",
                         color_by = "consensus_call_global")+
          geom_hline(aes(yintercept = mito_max), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          xlab("nFeature_RNA")+
          ylab("percent.mito")+
          annotate(geom = "text",size=4, 
                   label = paste0(as.numeric(table( sce_object$subsets_Mito_percent <= mito_max)[1]),"cells removed\n",
                                  as.numeric(table( sce_object$subsets_Mito_percent <= mito_max)[2]), "cells remain"),
                   x = max(sce_object$detected)-1000, y = mito_max+5)
        p <- ggExtra::ggMarginal(p, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"),
               plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
      }
      
    } else {
      if (sum(sce_object$subsets_Mito_percent <= mito_max) == length(sce_object$subsets_Mito_percent <= mito_max)){
        p <- plotColData(sce_object, x="detected", y="subsets_Mito_percent")+
          geom_hline(aes(yintercept = mito_max), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          xlab("nFeature_RNA")+
          ylab("percent.mito")+
          annotate(geom = "text",size=4, 
                   label = paste0(as.numeric(table( sce_object$subsets_Mito_percent <= mito_max)[1])," cells remain"),
                   x = max(sce_object$detected)-1000, y = mito_max+5)
        p <- ggExtra::ggMarginal(p, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"),
               plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
      }
      
      if(sum(sce_object$subsets_Mito_percent > mito_max)==length(sce_object$subsets_Mito_percent <= mito_max)){
        p <- plotColData(sce_object, x="detected", y="subsets_Mito_percent")+
          geom_hline(aes(yintercept = mito_max), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          xlab("nFeature_RNA")+
          ylab("percent.mito")+
          annotate(geom = "text",size=4, 
                   label = paste0(as.numeric(table( sce_object$subsets_Mito_percent <= mito_max)[1])," cells removed"),
                   x = max(sce_object$detected)-1000, y = mito_max+5)
        p <- ggExtra::ggMarginal(p, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"),
               plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
      }
      
      if(sum(sce_object$subsets_Mito_percent <= mito_max)!=length(sce_object$subsets_Mito_percent <= mito_max)&
         sum(sce_object$subsets_Mito_percent > mito_max)!=length(sce_object$subsets_Mito_percent <= mito_max)){
        
        p <- plotColData(sce_object, x="detected", y="subsets_Mito_percent")+
          geom_hline(aes(yintercept = mito_max), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          xlab("nFeature_RNA")+
          ylab("percent.mito")+
          annotate(geom = "text",size=4, 
                   label = paste0(as.numeric(table( sce_object$subsets_Mito_percent <= mito_max)[1]),"cells removed\n",
                                  as.numeric(table( sce_object$subsets_Mito_percent <= mito_max)[2]), "cells remain"),
                   x = max(sce_object$detected)-1000, y = mito_max+5)
        p <- ggExtra::ggMarginal(p, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"),
               plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
      }
    }
    
    ######## filtering by mito percent
    
    ######## Gene/UMI scatter plot before filtering
    p1 <- plot_UMIvsGene(sce_object[[, !sce_object$qc == "fails_maxMito"]],
                   colour_by = "subset_Mito_percent",
                   size_by = "subset_Ribo_percent")
    p1 <- ggExtra::ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_beforefiltering.pdf"),
           plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ####### Filter cells base on mito and feature
    lm <- lm(log10(sce_object[, sce_object$qc != "fails_maxMito" & sce_object$qc != "fails_minFeat"]$detected) ~ 
               log10(sce_object[, sce_object$qc != "fails_maxMito" & sce_object$qc != "fails_minFeat"]$sum))
    
    p2 <- plotColData(sce_object[, !sce_object$qc == "fails_maxMito" & !sce_object$qc == "fails_minFeat"],
                      x="sum", y="detected",
                      colour_by="subsets_Mito_percent",
                      size_by="subsets_Ribo_percent") +
      scale_x_log10() + scale_y_log10() + 
      scale_colour_gradient(aesthetics = c("colour"),low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method = "lm", show.legend = FALSE)+
      geom_hline(aes(yintercept = 10^min_nfeatures), colour = "darkgreen", linetype = 2)+
      geom_abline(intercept = lm$coefficients[1] - 0.09, slope = lm$coefficients[2], color="orange") +
      xlab("log10(nCount_RNA)")+
      ylab("log10(nFeature_RNA)") +
      theme_classic() + 
      theme(legend.position = "left") +
      labs(colour="percent.mito") +
      guides(size=guide_legend("percent.ribo")) +
      annotate(geom = "text",size=4, label = paste0(dim(sce_object[, !sce_object$qc == "fails_maxMito" & !sce_object$qc == "fails_minFeat"])[2], " QC passed cells"),
               x = 1000, y = 5000)
    
    p2 <- ggExtra::ggMarginal(p2, type = "densigram",bins=100, fill="grey", alpha=.3)
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_cutoff.pdf"),
           plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ######## Cells to exclude lie below an intercept offset of -0.09
    p3 <- plotColData(sce_object[, !sce_object$qc == "fails_maxMito" & !sce_object$qc == "fails_minFeat"],
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
                 table(sce_object[, !sce_object$qc == "fails_maxMito" & !sce_object$qc == "fails_minFeat"]$qc)[[2]], " QC passed cells\n",
                 table(sce_object[, !sce_object$qc == "fails_maxMito" & !sce_object$qc == "fails_minFeat"]$qc)[[1]], " QC filtered"),
               x = 2000, y = 10000)
    
    p3 <- ggExtra::ggMarginal(p3, type = "densigram",bins=100, fill="grey", alpha=.3)
    
    p3
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_validCells.pdf"),
           plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  } ### close if mito genes > 1
  
  if (!length(mito.genes) >= 1){
    
    ######## Gene/UMI scatter plot before filtering
    p1 <- plot_UMIvsGene(sce_object,
                         colour_by = "subset_Mito_percent",
                         size_by = "subset_Ribo_percent")
    p1 <- p1 + geom_hline(aes(yintercept = min_nfeatures), colour = "darkgreen", linetype = 2)
    p1 <- ggExtra::ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_beforefiltering.pdf"),
           plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ####### Filter cells base on min features and lm
    # line 477
    
    
    
    
    
    
    
    
    
    
    
    
  }  ### close if mito genes < 1
  
  
  
  
########### closing bracket
}





