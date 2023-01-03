
QC_per_sample<-function(seurat,projet,hto=FALSE){
  
  
  print(projet)
  
  dir.create(paste0("out/1-Loading_Qc/",projet))
  path=paste0("projet='",projet,"'")
  
  seurat <- Qc(seurat,filter_singlet=FALSE,HTO_tag=NULL,hto=hto,vars.to.regress=c("nFeature_RNA","percent.mito"),projet=projet,py_function=py_run_string(path))
  
  saveRDS(seurat, file = paste0("out/1-Loading/",projet,"/seurat.Rds"))
  seurat
}





library(ggplot2)



count_meta<-function(meta_count,name,projet,seurat){
  
  # print(name)
  # 
  # m_meta_name<-round((table(meta_count,seurat@active.ident)/rowSums(table(meta_count,seurat@active.ident)))*100)
  # 
  # #m_meta_name<-table(meta_count,seurat@active.ident)
  # 
  # levels_meta_name<-factor(rownames(m_meta_name))
  # levels_cluster<-factor(colnames(m_meta_name),levels=sort(as.integer(colnames(m_meta_name))))
  # 
  # m_meta_name <- data.frame(cluster=as.factor(col(m_meta_name)),meta_name=as.factor(row(m_meta_name)),n=as.vector(m_meta_name))
  # 
  # 
  # m_meta_name$cluster<- levels_cluster[m_meta_name$cluster]
  # m_meta_name$meta_name<-levels_meta_name[m_meta_name$meta_name]
  # 
  # 
  # p1 <- ggplot(m_meta_name,aes(x=cluster,y=meta_name,fill=n)) + 
  #   geom_tile(aes(fill=n)) +
  #   geom_text(aes(label=n,colour="white",size=0.5)) +
  #   scale_fill_viridis_c() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_equal()+ggtitle(paste0(name))
  # 
  # 
  
  ####
  
  
  m_meta_name<-t(round(t(table(meta_count,seurat@active.ident))/colSums(table(meta_count,seurat@active.ident))*100))
  #m_meta_name<-table(meta_count,seurat@active.ident)
  
  
  levels_meta_name<-factor(rownames(m_meta_name))
  levels_cluster<-factor(colnames(m_meta_name),levels=sort(as.integer(colnames(m_meta_name))))
  
  m_meta_name <- data.frame(cluster=as.factor(col(m_meta_name)),meta_name=as.factor(row(m_meta_name)),n=as.vector(m_meta_name))
  m_meta_name$cluster<- levels_cluster[m_meta_name$cluster]
  m_meta_name$meta_name<-levels_meta_name[m_meta_name$meta_name]
  
  
  p2 <- ggplot(m_meta_name,aes(x=cluster,y=meta_name,fill=n)) + 
    geom_tile(aes(fill=n)) +
    geom_text(aes(label=n,colour="white",size=0.5)) +
    scale_fill_viridis_c() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_equal()+ggtitle(paste0(name))
  
  
  ggsave(filename = paste0("pdf/1-Loading/",projet,"/Umap_hto_.pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  
  
  m_meta_name
  
}





### QC function

Qc<-function(seurat,
             filter_singlet=FALSE,
             methods=c("multiseq","gmm_demux", "htodemux", "dropletutils"),
             projet="Telley",
             mito_max=15,
             HTO_tag=HTO_tag,
             dims=1:30,
             min_nfeatures=log10(200),
             hto=TRUE,
             py_function=py_run_string("projet='Telley'"),
             scrublet=FALSE,
             vars.to.regress=c("nFeature_RNA","percent.mito"),
             Seurat=TRUE,
             res=c(1.0),
             fill= hue_pal()(7)[1],
             color= hue_pal()(7)[1]){
  
  require(Seurat)
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
  require(gridExtra)
  
  
  dir.create(paste0("pdf/1-Loading/",projet))
  dir.create(paste0("out/1-Loading/",projet))
  
  count.data <- GetAssayData(object = seurat[["RNA"]], slot = "counts")
  meta.data<-seurat@meta.data
  Raw.data <- CreateSeuratObject(counts = count.data,meta.data=meta.data ,project=projet)
  
  if(hto==TRUE){
    # retrieving HTO tags for demultiplexing
    ## filtering based on passed argument
    htos<-seurat@assays$HTO[HTO_tag]
    
    Raw.data[["HTO"]] <- CreateAssayObject(counts = htos)
    
    #Raw.data<-NormalizeData(Raw.data,assay = "HTO")
    #Raw.data <- NormalizeData(Raw.data, assay = "HTO", normalization.method = "CLR")
    
    #Raw.data<-extract_singlets_fromHTO(Raw.data,projet=projet,init=init,HTO_tag=HTO_tag,positive.quantile=positive.quantile,subset=subset)
    
    
    
    fileConn<-file(paste0("out/1-Loading/",projet,"/cellbarcodeWhitelist.txt"))
    writeLines(gsub("-",".",colnames(htos)), fileConn)
    close(fileConn)
  
    # writing table that will be read by cellHashR
    write.table(htos,file=paste0("out/1-Loading/",projet,"/hto.txt"),sep=",")
    
    pdf(paste0("pdf/1-Loading/",projet,"/hto.pdf"),width = 20,height = 20)
    
    # parse HTO table, specifying the HTOs (whitelist)
    barcodeData <- ProcessCountMatrix(rawCountData =paste0("out/1-Loading/",projet,"/hto.txt"),
                                      minCountPerCell = 5,
                                      simplifyBarcodeNames = FALSE,
                                      cellbarcodeWhitelist=paste0("out/1-Loading/",projet,"/cellbarcodeWhitelist.txt"),
                                      barcodeWhitelist = HTO_tag)
    
    # Create QC plots of barcode normalization
    #PlotNormalizationQC(barcodeData)
    
    #methods <- c("bff_cluster", "bff_threshold", "gmm_demux", "multiseq", "htodemux", "dropletutils")
    # Generate the final cell hashing calls
    calls <- GenerateCellHashingCalls(barcodeMatrix = barcodeData,
                                      methods = methods,
                                      doTSNE=FALSE,
                                      majorityConsensusThreshold=0.3)
    
    # create summary plots to contrast cells based on call-status. Help inform why some cells were not called
    SummarizeCellsByClassification(calls = calls, barcodeMatrix = barcodeData)
    
    
    dev.off()
    
    calls$cellbarcode<-gsub("[.]","-",calls$cellbarcode)

    
    Raw.data<-subset(Raw.data,cells=calls$cellbarcode)
    
    Raw.data@meta.data<-cbind(Raw.data@meta.data,calls)
  
    
    table1<-table(Raw.data$consensuscall.global)

    if(filter_singlet==TRUE){
      Raw.data@active.ident<-factor(Raw.data$consensuscall.global)
      Raw.data <- subset(Raw.data, idents = "Singlet", invert = FALSE)
    }
  }
  
  #### when no HTO tags
  
  DefaultAssay(Raw.data) <- "RNA"
  
  
  Raw.data@meta.data$Barcodes <- rownames(Raw.data@meta.data)
  
  mytheme_classic <- theme_classic()+ 
    theme(axis.text.y = element_text(size=14, color="black"),
          axis.ticks = element_line(size=1, color="black"),
          axis.title.y = element_text(size=18, color="black"),
          axis.text.x = element_text(size=14, color="black"),
          axis.title.x=element_text(size=18),
          panel.background =  element_blank(),
          legend.text = element_text(size=14))
  
  
  
  
  ########## Percent of mitochondrial counts
  mito.genes <- grep(pattern = "^mt-", x = rownames(Raw.data), value = TRUE)
  

  if(length(mito.genes)>=1){
    percent.mito <- (Matrix::colSums(Raw.data[mito.genes, ])/Matrix::colSums(Raw.data))*100
    Raw.data <- AddMetaData(object = Raw.data, metadata = percent.mito, col.name = "percent.mito")
    
    ########### Percent of mitochondrial ribosomal
    ribo.genes <- grep(pattern = "(^Rpl|^Rps|^Mrp)", x = rownames(Raw.data), value = TRUE)
    percent.ribo <- Matrix::colSums(Raw.data[ribo.genes, ])/Matrix::colSums(Raw.data)*100
    Raw.data <- AddMetaData(object = Raw.data, metadata = percent.ribo, col.name = "percent.ribo")
    
    p1 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=nCount_RNA))+
      geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
      geom_boxplot(width=0.025, outlier.size = 0.5, outlier.alpha = 0.8)+
      geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
                 alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
      xlab(NULL)+ ggtitle("Number of UMIs")+
      mytheme_classic + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                             axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
    
    p2 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=nFeature_RNA))+
      geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
      geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
      geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
                 alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
      xlab(NULL)+ ggtitle("Number of Genes")+
      mytheme_classic + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                             axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
    
    p3 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=percent.mito))+
      geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
      geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
      geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
                 alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
      xlab(NULL)+ ggtitle("Mitochondrial genes")+
      mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                             axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
    
    p4 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=percent.ribo))+
      geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
      geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
      geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
                 alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
      xlab(NULL)+ ggtitle("Ribosomal genes")+
      mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                             axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
    
    p <- plot_grid(plotlist =  list(p1,p2, p3, p4),
                   ncol=2,
                   align = "hv",
                   rel_widths = c(1, 1, 1, 1),
                   rel_heights =  c(1, 1, 1, 1))
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/Flatviolin QC_before filtering.pdf"),
           plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ########### Relation between nUMI and nGene detected
    p1 <- ggplot(Raw.data@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + 
      geom_point(aes(color=percent.mito, size= percent.ribo)) + 
      scale_colour_gradient(low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method="lm")+ mytheme_classic + theme(legend.position = "left")
    p1 <- ggMarginal(p1, type = "densigram",bins=100, fill="grey", alpha=0.3) 
    
    
    p2 <- ggplot(Raw.data@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + 
      geom_point(aes(color=percent.mito, size= percent.ribo)) + 
      scale_colour_gradient(low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method="lm")+mytheme_classic + theme(legend.position = "left")
    p2 <- ggMarginal(p2, type = "densigram",bins=100, fill="grey", alpha=0.3)
    
    
    p <- plot_grid(plotlist = list(p1,p2), ncol=2, align = "hv", rel_widths = c(1, 1), rel_heights =  c(1, 1))
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesVSnUMIs.pdf"),
           plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    Cell.QC.Stat <- Raw.data@meta.data
    max.mito.thr <- mito_max
    
    if(hto==TRUE){
      
      
      if(sum(Cell.QC.Stat$percent.mito <= max.mito.thr)==length(Cell.QC.Stat$percent.mito <= max.mito.thr)){
        
        p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mito,color=consensuscall.global)) +
          geom_point(alpha=0.5) +
          geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          annotate(geom = "text",size=6, label = paste0(as.numeric(table( Cell.QC.Stat$percent.mito <= max.mito.thr)[1])," cells remain"), 
                   x = max(Cell.QC.Stat$nFeature_RNA)-1000, y = max.mito.thr+5)+
          mytheme_classic
        
        
        p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      if(sum(Cell.QC.Stat$percent.mito > max.mito.thr)==length(Cell.QC.Stat$percent.mito <= max.mito.thr)){
        
        p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mito,color=consensuscall.global)) +
          geom_point(alpha=0.5) +
          geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          annotate(geom = "text",size=6, label = paste0(as.numeric(table( Cell.QC.Stat$percent.mito <= max.mito.thr)[1])," cells removed"), 
                   x = max(Cell.QC.Stat$nFeature_RNA)-1000, y = max.mito.thr+5)+
          mytheme_classic
        
        
        p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      if(sum(Cell.QC.Stat$percent.mito <= max.mito.thr)!=length(Cell.QC.Stat$percent.mito <= max.mito.thr)&sum(Cell.QC.Stat$percent.mito > max.mito.thr)!=length(Cell.QC.Stat$percent.mito <= max.mito.thr)){
        
        p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mito,color=consensuscall.global)) +
          geom_point(alpha=0.5) +
          geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          annotate(geom = "text",size=6, label = paste0(as.numeric(table( Cell.QC.Stat$percent.mito <= max.mito.thr)[1])," cells removed\n",
                                                        as.numeric(table( Cell.QC.Stat$percent.mito <= max.mito.thr)[2])," cells remain"), 
                   x = max(Cell.QC.Stat$nFeature_RNA)-1000, y = max.mito.thr+5)+
          mytheme_classic
        p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      
      
    }
    
    if(hto==FALSE){
      
      if(sum(Cell.QC.Stat$percent.mito <= max.mito.thr)==length(Cell.QC.Stat$percent.mito <= max.mito.thr)){
        
        p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mito)) +
        geom_point(alpha=0.5) +
        geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
        geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          annotate(geom = "text",size=6, label = paste0(as.numeric(table( Cell.QC.Stat$percent.mito <= max.mito.thr)[1])," cells remain"), 
          x = max(Cell.QC.Stat$nFeature_RNA)-1000, y = max.mito.thr+5)+
            mytheme_classic
       
        
      p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
      
      ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      if(sum(Cell.QC.Stat$percent.mito > max.mito.thr)==length(Cell.QC.Stat$percent.mito <= max.mito.thr)){
        
        p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mito)) +
          geom_point(alpha=0.5) +
          geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          annotate(geom = "text",size=6, label = paste0(as.numeric(table( Cell.QC.Stat$percent.mito <= max.mito.thr)[1])," cells removed"), 
                   x = max(Cell.QC.Stat$nFeature_RNA)-1000, y = max.mito.thr+5)+
          mytheme_classic
        
        
        p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      if(sum(Cell.QC.Stat$percent.mito <= max.mito.thr)!=length(Cell.QC.Stat$percent.mito <= max.mito.thr)&sum(Cell.QC.Stat$percent.mito > max.mito.thr)!=length(Cell.QC.Stat$percent.mito <= max.mito.thr)){
        
        p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mito)) +
          geom_point(alpha=0.5) +
          geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
          geom_vline(aes(xintercept = 10^min_nfeatures), colour = "red", linetype = 2) +
          annotate(geom = "text",size=6, label = paste0(as.numeric(table( Cell.QC.Stat$percent.mito <= max.mito.thr)[1])," cells removed\n",
                                                        as.numeric(table( Cell.QC.Stat$percent.mito <= max.mito.thr)[2])," cells remain"), 
                   x = max(Cell.QC.Stat$nFeature_RNA)-1000, y = max.mito.thr+5)+
          mytheme_classic
        p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/percent.mitoVSnGenes.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
    }
    
    ######### Filter cells based on these thresholds
    Cell.QC.Stat <- Cell.QC.Stat  %>% filter(percent.mito <= max.mito.thr)
    
    ######## Set low and hight thresholds on the number of detected genes
    min.nFeature.thr <- min_nfeatures
    
    
    
    ######## Gene/UMI scatter plot before filtering
    p1 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point(aes(color=percent.mito, size= percent.ribo)) + 
      scale_colour_gradient(low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method="lm") +
      geom_hline(aes(yintercept = min.nFeature.thr), colour = "darkgreen", linetype = 2) +
      mytheme_classic +theme(legend.position = "left")
    p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_beforefiltering.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ####### Filter cells base on both metrics
    Cell.QC.Stat <- Cell.QC.Stat %>%  filter(log10(nFeature_RNA) >= min.nFeature.thr)
    lm.model <- lm(data = Cell.QC.Stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
    
    
    
    p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point(aes(color=percent.mito, size= percent.ribo)) + 
      scale_colour_gradient(low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method="lm") +
      geom_hline(aes(yintercept = min.nFeature.thr), colour = "darkgreen", linetype = 2) +
      geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
      annotate(geom = "text",size=6, label = paste0(dim(Cell.QC.Stat)[1], " QC passed cells"), x = 4, y = 3.8)+mytheme_classic + theme(legend.position = "left")
    p2<-ggMarginal(p2, type = "densigram",bins=100, fill="grey", alpha=0.3)
    
    
    ggsave(filename  = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_cutoff.pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    
    Cell.QC.Stat$valideCells <- (log10(Cell.QC.Stat$nFeature_RNA) > (log10(Cell.QC.Stat$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09)))
    
    ######## Cells to exclude lie below an intercept offset of -0.09 
    
    p3 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point(aes(colour = valideCells), alpha=0.5, show.legend = F) +
      scale_color_manual(values =  c("darkred", "dodgerblue"))+
      geom_smooth(method="lm", color="black") +
      geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
      theme(legend.position="none") +
      annotate(geom = "text", size=6, label = paste0(as.numeric(table(Cell.QC.Stat$valideCells)[2]), " QC passed cells\n",
                                                     as.numeric(table(Cell.QC.Stat$valideCells)[1]), " QC filtered"), x = 4, y = 3.8)+mytheme_classic
    p3<-ggMarginal(p3, type = "densigram",bins=100, fill="grey",  alpha=0.3)
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_validCells.pdf"), plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
  }
  
  
  if(!length(mito.genes)>=1){
    
    Cell.QC.Stat <- Raw.data@meta.data
    
    
    ######## Set low and hight thresholds on the number of detected genes
    min.nFeature.thr <- min_nfeatures
    
    
    
    ######## Gene/UMI scatter plot before filtering
    p1 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point() + 
      scale_colour_gradient(low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method="lm") +
      geom_hline(aes(yintercept = min.nFeature.thr), colour = "darkgreen", linetype = 2) +
      mytheme_classic +theme(legend.position = "left")
    p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey",  alpha=0.3)
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_beforefiltering.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ####### Filter cells base on both metrics
    Cell.QC.Stat <- Cell.QC.Stat %>%  filter(log10(nFeature_RNA) >= min.nFeature.thr)
    lm.model <- lm(data = Cell.QC.Stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
    
    
    
    p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point() + 
      scale_colour_gradient(low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+
      geom_smooth(method="lm") +
      geom_hline(aes(yintercept = min.nFeature.thr), colour = "darkgreen", linetype = 2) +
      geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
      annotate(geom = "text",size=6, label = paste0(dim(Cell.QC.Stat)[1], " QC passed cells"), x = 4, y = 3.8)+mytheme_classic + theme(legend.position = "left")
    p2<-ggMarginal(p2, type = "densigram",bins=100, fill="grey", alpha=0.3)
    
    
    ggsave(filename  = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_cutoff.pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    
    Cell.QC.Stat$valideCells <- (log10(Cell.QC.Stat$nFeature_RNA) > (log10(Cell.QC.Stat$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09)))
    
    ######## Cells to exclude lie below an intercept offset of -0.09 
    
    p3 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point(aes(colour = valideCells), alpha=0.5, show.legend = F) +
      scale_color_manual(values =  c("darkred", "dodgerblue"))+
      geom_smooth(method="lm", color="black") +
      geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
      theme(legend.position="none") +
      annotate(geom = "text", size=6, label = paste0(as.numeric(table(Cell.QC.Stat$valideCells)[2]), " QC passed cells\n",
                                                     as.numeric(table(Cell.QC.Stat$valideCells)[1]), " QC filtered"), x = 4, y = 3.8)+mytheme_classic
    p3<-ggMarginal(p3, type = "densigram",bins=100, fill="grey",  alpha=0.3)
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nGenesvsnUMIs_validCells.pdf"), plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
  }
  
  ######## Remove unvalid cells
  Cell.QC.Stat <- Cell.QC.Stat %>% filter(valideCells)
  Raw.data <- subset(Raw.data, cells = as.vector(rownames(Cell.QC.Stat)))
  
  ######## Plot final QC metrics
  p1 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=nCount_RNA))+
    geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
    geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
    geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
               alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
    xlab(NULL)+ ggtitle("Number of UMIs")+
    mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
  
  p2 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=nFeature_RNA))+
    geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
    geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
    geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
               alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
    xlab(NULL)+ ggtitle("Number of Genes")+
    mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
  
  if(length(mito.genes)>=1){
    p3 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=percent.mito))+
      geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
      geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
      geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
                 alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
      xlab(NULL)+ ggtitle("Mitochondrial genes")+
      mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
    
    p4 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=percent.ribo))+
      geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
      geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
      geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
                 alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
      xlab(NULL)+ ggtitle("Ribosomal genes")+
      mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
    
    p <- plot_grid(plotlist =  list(p1,p2, p3, p4), ncol=2, align = "hv", rel_widths = c(1, 1, 1, 1),rel_heights =  c(1, 1, 1, 1))
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/Flatviolin QC_after filtering.pdf"), plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p1 <- ggplot(Raw.data@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + 
      geom_point(aes(color=percent.mito, size= percent.ribo)) + 
      scale_colour_gradient(low = "gray90", high = "black") +
      scale_size_area(max_size = 3)+ 
      geom_smooth(method="lm")+mytheme_classic+theme(legend.position = "left")
    p1<-ggMarginal(p1, type = "densigram",bins=100, fill="grey", alpha=0.3)
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/nUMIsVSnGenes_afterfiltering.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    
  }
  
  if(!length(mito.genes)>=1){
    
    p1 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=nCount_RNA))+
      geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
      geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
      geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
                 alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
      xlab(NULL)+ ggtitle("Number of UMIs")+
      mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
    
    p2 <- ggplot(Raw.data@meta.data, aes(x=orig.ident, y=nFeature_RNA))+
      geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
      geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
      geom_point(aes(x=as.numeric(orig.ident)-0.3), fill="black", color="black", shape=16, size=0.25,
                 alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
      xlab(NULL)+ ggtitle("Number of Genes")+
      mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
    
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/QC_after filtering_1.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/QC_after filtering_2.pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
  }
  
  
  
  dir.create(paste0("out/1-Loading/",projet))
  
  
  if(scrublet==TRUE){
    
    
    
    ####### Export filtered matrix
    dir.create(paste0("out/1-Loading/",projet,"/QCdata/"))
    exprData <- Matrix(as.matrix(Raw.data@assays$RNA@data), sparse = TRUE)
    writeMM(exprData, paste0("out/1-Loading/",projet,"/QCdata/matrix.mtx"))
    py_run_string("import scrublet as scr")
    py_run_string("import scipy.io")
    py_run_string("import numpy as np")
    py_run_string("import os")
    py_function
    
    #py_run_string("input_dir = 'pdf'")
    py_run_string("path='out/1-Loading/'+projet")
    #Load raw counts matrix and gene list
    py_run_string("counts_matrix = scipy.io.mmread( path+'/QCdata/matrix.mtx').T.tocsc()")
    #Initialize Scrublet object
    py_run_string("scrub = scr.Scrublet(counts_matrix,
                     expected_doublet_rate=0.1,
                     sim_doublet_ratio=2,
                     n_neighbors = 8)")
    #Run the default pipeline
    py_run_string("doublet_scores, predicted_doublets = scrub.scrub_doublets()")
    
    ########## Import scrublet's doublet score
    Raw.data@meta.data$Doubletscore <- py$doublet_scores
    
    trhs_auto_doublet=mean(max(py$doublet_scores[!py$predicted_doublets]),min(py$doublet_scores[py$predicted_doublets]))
    Raw.data@meta.data$Predicted_doublets <- ifelse(py$doublet_scores > trhs_auto_doublet, "Doublet","Singlet" )
    
    print(table(Raw.data@meta.data$Predicted_doublets))
    
    ####### Plot doublet score
    p1 <- ggplot(Raw.data@meta.data, aes(x = Doubletscore, stat(ndensity))) +
      geom_histogram(bins = 200, colour ="lightgrey")+
      geom_vline(xintercept = trhs_auto_doublet, colour = "orangered", linetype = 2)+
      annotate(geom = "text",size=6, label = paste0(as.numeric(table(Raw.data@meta.data$Predicted_doublets)[2]), " QC passed cells\n",
                                                    as.numeric(table(Raw.data@meta.data$Predicted_doublets)[1]), " QC filtered"), x = 0.5, y = 1)+
      mytheme_classic+
      scale_y_continuous(labels = scales::number_format(decimal.mark = "."))
    
    
    ggsave(filename = paste0("pdf/1-Loading/",projet,"/histodoublet.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    
    
    
    
    
    if(hto==TRUE){
      hto_count<-seurat[["HTO"]]@counts[,colnames(seurat[["HTO"]]@counts)%in%c(colnames(Raw.data))]
      
      hto_count<-hto_count[,match(colnames(Raw.data),colnames(hto_count))]
      
      Raw.data[["HTO"]] <- CreateAssayObject(counts = hto_count )
      
    }
    
    if(Seurat==TRUE){
      
      
      ######### SCT transform
      if(length(mito.genes)>=1){
        
        Raw.data <- SCTransform(Raw.data, vars.to.regress = vars.to.regress)
      }
      
      if(!length(mito.genes)>=1){
        
        vars.to.regress<-vars.to.regress[-which(vars.to.regress=="percent.mito")]
        Raw.data <- SCTransform(Raw.data,vars.to.regress = vars.to.regress)
      }
      
      
      # Run the standard workflow for visualization and clustering
      Raw.data <- RunPCA(Raw.data)
      
      # Run UMAP
      Raw.data <- RunUMAP(Raw.data, dims = dims)
      
      
      # Determine the K-nearest neighbor graph
      Raw.data <- FindNeighbors(object =Raw.data,  dims = dims)
      
      # Determine the clusters for various resolutions                                
      Raw.data <- FindClusters(object =Raw.data,resolution = res)
      
      
      
      
      ######## Plot the UMAP of predicted singlet doublet
      umap.df <- as.data.frame(Raw.data@reductions$umap@cell.embeddings)
      umap.df <-umap.df %>% cbind(seurat_clusters=Raw.data@meta.data$seurat_clusters,
                                  Predicted_doublets= Raw.data@meta.data$Predicted_doublets)
      centroidsUMAP_cl <- aggregate(cbind(UMAP_1,UMAP_2)~seurat_clusters,umap.df,median)
      umap.df$Predicted_doublets <- factor(umap.df$Predicted_doublets, levels = c("Singlet", "Doublet"))
      umap.df$seurat_clusters<- factor(umap.df$seurat_clusters, levels =unique(umap.df$seurat_clusters))
      
      p1 <- ggplot(umap.df , aes(x=UMAP_1, y=UMAP_2)) + 
        geom_point(aes(color=Predicted_doublets),shape=16, size=1, show.legend = NA, alpha=0.8)+
        scale_color_manual(values=c("grey50", "orangered"), guide=guide_legend(override.aes = list(size=3),ncol=1, title = NULL ))+
        geom_text_repel(data=centroidsUMAP_cl,aes(label=seurat_clusters), color="black",size=3, show.legend = F)+
        xlab("umap_1")+ylab("umap_2")+labs(color="cluster")+
        mytheme_classic+expand_limits(x=c(min(umap.df$UMAP_1),max(umap.df$UMAP_1)),y=c(min(umap.df$UMAP_2),max(umap.df$UMAP_2)))
      
      
      ggsave(filename = paste0("pdf/1-Loading/",projet,"/UMAP_predicteddoublet_scrublet.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      
      
      if(hto==TRUE){
        
        umap.df <- as.data.frame(Raw.data@reductions$umap@cell.embeddings)
        umap.df <-umap.df %>% cbind(seurat_clusters=Raw.data@meta.data$seurat_clusters,
                                    Predicted_id= Raw.data@meta.data$consensuscall)
        centroidsUMAP_cl <- aggregate(cbind(UMAP_1,UMAP_2)~seurat_clusters,umap.df,median)
        umap.df$Predicted_id <- factor(umap.df$Predicted_id)
        umap.df$seurat_clusters<- factor(umap.df$seurat_clusters, levels =unique(umap.df$seurat_clusters))
        
        p1 <- ggplot(umap.df , aes(x=UMAP_1, y=UMAP_2)) + 
          geom_point(aes(color=Predicted_id),shape=16, size=1, show.legend = NA, alpha=0.8)+
          xlab("umap_1")+ylab("umap_2")+labs(color="cluster")+
          mytheme_classic+expand_limits(x=c(min(umap.df$UMAP_1),max(umap.df$UMAP_1)),y=c(min(umap.df$UMAP_2),max(umap.df$UMAP_2)))
        
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/UMAP_Predicted_id.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
        umap.df <- as.data.frame(Raw.data@reductions$umap@cell.embeddings)
        umap.df <-umap.df %>% cbind(seurat_clusters=Raw.data@meta.data$seurat_clusters,
                                    Predicted_id= Raw.data@meta.data$consensuscall.global)
        centroidsUMAP_cl <- aggregate(cbind(UMAP_1,UMAP_2)~seurat_clusters,umap.df,median)
        umap.df$Predicted_id <- factor(umap.df$Predicted_id)
        umap.df$seurat_clusters<- factor(umap.df$seurat_clusters, levels =unique(umap.df$seurat_clusters))
        
        p2 <- ggplot(umap.df , aes(x=UMAP_1, y=UMAP_2)) + 
          geom_point(aes(color=Predicted_id),shape=16, size=1, show.legend = NA, alpha=0.8)+
          xlab("umap_1")+ylab("umap_2")+labs(color="cluster")+
          mytheme_classic+expand_limits(x=c(min(umap.df$UMAP_1),max(umap.df$UMAP_1)),y=c(min(umap.df$UMAP_2),max(umap.df$UMAP_2)))
        
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/UMAP_Predicted_global.pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
        
      }
      
      
      
      # 
      # ########### plot log RNA UMIs & pie chart proportion
      # 
      # Raw.data@meta.data$Predicted_doublets <- factor(Raw.data@meta.data$Predicted_doublets, levels=c("Singlet", "Doublet"))
      # 
      # p1 <- ggplot(Raw.data@meta.data, aes(x=Predicted_doublets, y=log(nCount_RNA)))+
      #   geom_flat_violin(aes(fill=Predicted_doublets, color=Predicted_doublets), trim= T, show.legend = F)+
      #   scale_fill_manual(values=alpha(c("grey40", "orangered"),0.5), guide=guide_legend(title = NULL))+
      #   scale_color_manual(values=c("grey50", "orangered"), guide=guide_legend(title = NULL))+
      #   new_scale_color()+ 
      #   geom_point(aes(x=as.numeric(Predicted_doublets)-0.3 ,fill=Predicted_doublets,color=Predicted_doublets), shape=16,size=0.5, 
      #              alpha=0.5, position = position_jitter(width=0.225,seed = 12345), show.legend = F)+
      #   scale_color_manual(values=c("grey50", "orangered"))+ ylab("log RNA UMIs")+ xlab(NULL)+
      #   stat_summary(data=Raw.data@meta.data,fun.data=mean_sdl, fun.args = list(mult=1), 
      #                geom="errorbar", color=(c("grey50", "orangered")), size=1, width= 0.1)+
      #   stat_summary(data=Raw.data@meta.data,fun.y=mean, geom="point", color=(c("grey50","orangered")),size=1.5,stroke=1.5, shape=16)+
      #   mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())+
      #   stat_compare_means(data=Raw.data@meta.data,aes(label = ..p.signif..), comparisons = list(c('Singlet', 'Doublet')), method = "wilcox.test", 
      #                      size=12, tip.length= 0,label.y =max(log(Raw.data@meta.data$nCount_RNA))+0.025)
      # 
      # Singlet_doublet_all.df <- as.data.frame(table(Raw.data@meta.data$Predicted_doublets))
      # colnames(Singlet_doublet_all.df)[1] <- "Predicted_doublets"
      # 
      # Singlet_doublet_all.df <- Singlet_doublet_all.df %>%
      #   arrange(desc(Predicted_doublets)) %>%
      #   mutate(lab.ypos = cumsum(Freq) - 0.5*Freq)
      # 
      # Singlet_doublet_all.df$Predicted_doublets <- factor(Singlet_doublet_all.df$Predicted_doublets, levels = c('Singlet', 'Doublet'))
      # 
      # p2 <- ggplot(Singlet_doublet_all.df, aes(x="", y=Freq, fill=Predicted_doublets))+
      #   scale_fill_manual(values = c("grey50", "orangered"))+
      #   geom_bar(width = 1, stat = "identity", aes(color=Predicted_doublets))+
      #   scale_color_manual(values=rep("white",2), guide=F)+
      #   coord_polar("y", start=0)+labs(fill=NULL)+
      #   theme_void()+ theme(title=element_text(size=16))+
      #   geom_label_repel(aes(y=lab.ypos, label = percent(Freq/sum(Freq), accuracy = 0.01)), size=5, show.legend = F, nudge_x = 1)+
      #   theme(legend.text = element_text(size=14))
      # 
      # p <- plot_grid(plotlist = list(p1, p2), ncol=2, rel_heights = c(2,1), rel_widths = c(2,1))
      # 
      # ggsave(filename= paste0("pdf/1-Loading/",projet,"/logUMIs proportion_predicteddoublet.pdf"), plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
      # 
      # ##### plot percent in each custer
      # 
      # percent_cl_SvsD.df <- as.data.frame(table(Raw.data@meta.data$seurat_clusters, Raw.data@meta.data$Predicted_doublets))
      # colnames(percent_cl_SvsD.df) <- c("seurat_clusters", "Predicted_doublets","n_cell")
      # 
      # percent_cl_SvsD.df$Predicted_doublets  <- factor(percent_cl_SvsD.df$Predicted_doublets, levels=c("Singlet", "Doublet"))
      # 
      # 
      # percent_cl_SvsD.df <- percent_cl_SvsD.df %>% 
      #   cbind(percent_cell=percent_cl_SvsD.df$n_cell/rowSums(matrix(percent_cl_SvsD.df$n_cell,nrow=length(unique(percent_cl_SvsD.df$seurat_clusters)), ncol=2))*100)
      # 
      # p1 <- ggplot(percent_cl_SvsD.df, aes(x=seurat_clusters, y=percent_cell))+
      #   geom_col(aes(fill=Predicted_doublets))+
      #   scale_fill_manual(values=c("grey50", "orangered"), guide= guide_legend(title = NULL))+
      #   scale_y_continuous(expand=c(0,0))+
      #   xlab(NULL)+ ylab("Fraction of singlet/doublet cells (%)")+
      #   mytheme_classic+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=12), axis.line.x = element_blank(), 
      #                         axis.ticks.x = element_blank(), aspect.ratio = 0.5)
      # 
      # ggsave(filename= paste0("pdf/1-Loading/",projet,"/fractioncl_predicteddoublet.pdf"), plot =p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      # 
      # 
      
      QCFiltered.data<-Raw.data
      QCFiltered.data
    }
    
    
    if(Seurat==FALSE){QCFiltered.data <- Raw.data}
    
  }
  
  if(scrublet==FALSE){
    
    if(hto==TRUE){
      hto_count<-seurat[["HTO"]]@counts[,colnames(seurat[["HTO"]]@counts)%in%c(colnames(Raw.data))]
      
      hto_count<-hto_count[,match(colnames(Raw.data),colnames(hto_count))]
      
      Raw.data[["HTO"]] <- CreateAssayObject(counts = hto_count )
      
    }
    
    if(Seurat==TRUE){
      
      
      ######### SCT transform
      if(length(mito.genes)>=1){
        
        Raw.data <- SCTransform(Raw.data, vars.to.regress = vars.to.regress)
      }
      
      if(!length(mito.genes)>=1){
        
        vars.to.regress<-vars.to.regress[-which(vars.to.regress=="percent.mito")]
        Raw.data <- SCTransform(Raw.data,vars.to.regress = vars.to.regress)
      }
      
      
      # Run the standard workflow for visualization and clustering
      Raw.data <- RunPCA(Raw.data)
      
      # Run UMAP
      Raw.data <- RunUMAP(Raw.data, dims = dims)
      
      
      # Determine the K-nearest neighbor graph
      Raw.data <- FindNeighbors(object =Raw.data,  dims = dims)
      
      # Determine the clusters for various resolutions                                
      Raw.data <- FindClusters(object =Raw.data,resolution = res)
      
      
      # Assign identity of clusters
      #Idents(object =Raw.data) <- "SCT_snn_res.2"
      
      
      
      # ######## Plot the UMAP of predicted singlet doublet
      # umap.df <- as.data.frame(Raw.data@reductions$umap@cell.embeddings)
      # umap.df <-umap.df %>% cbind(seurat_clusters=Raw.data@meta.data$seurat_clusters)
      # centroidsUMAP_cl <- aggregate(cbind(UMAP_1,UMAP_2)~seurat_clusters,umap.df,median)
      # umap.df$seurat_clusters<- factor(umap.df$seurat_clusters, levels =unique(umap.df$seurat_clusters))
      # 
      # p1 <- ggplot(umap.df , aes(x=UMAP_1, y=UMAP_2)) + 
      #   geom_point(aes(color=seurat_clusters),shape=16, size=1, show.legend = NA, alpha=0.8)+
      #   geom_text_repel(data=centroidsUMAP_cl,aes(label=seurat_clusters), color="black",size=3, show.legend = F)+
      #   xlab("umap_1")+ylab("umap_2")+labs(color="cluster")+
      #   mytheme_classic+expand_limits(x=c(min(umap.df$UMAP_1),max(umap.df$UMAP_1)),y=c(min(umap.df$UMAP_2),max(umap.df$UMAP_2)))
      # 
      # 
      # ggsave(filename = paste0("pdf/1-Loading/",projet,"/UMAP.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      
      
      
      if(hto==TRUE){
        
        umap.df <- as.data.frame(Raw.data@reductions$umap@cell.embeddings)
        umap.df <-umap.df %>% cbind(seurat_clusters=Raw.data@meta.data$seurat_clusters,
                                    Predicted_id= Raw.data@meta.data$consensuscall)
        centroidsUMAP_cl <- aggregate(cbind(UMAP_1,UMAP_2)~seurat_clusters,umap.df,median)
        umap.df$Predicted_id <- factor(umap.df$Predicted_id)
        umap.df$seurat_clusters<- factor(umap.df$seurat_clusters, levels =unique(umap.df$seurat_clusters))
        
        p1 <- ggplot(umap.df , aes(x=UMAP_1, y=UMAP_2)) + 
          geom_point(aes(color=Predicted_id),shape=16, size=1, show.legend = NA, alpha=0.8)+
          xlab("umap_1")+ylab("umap_2")+labs(color="cluster")+
          mytheme_classic+expand_limits(x=c(min(umap.df$UMAP_1),max(umap.df$UMAP_1)),y=c(min(umap.df$UMAP_2),max(umap.df$UMAP_2)))
        
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/UMAP_Predicted_id.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
        umap.df <- as.data.frame(Raw.data@reductions$umap@cell.embeddings)
        umap.df <-umap.df %>% cbind(seurat_clusters=Raw.data@meta.data$seurat_clusters,
                                    Predicted_id= Raw.data@meta.data$consensuscall.global)
        centroidsUMAP_cl <- aggregate(cbind(UMAP_1,UMAP_2)~seurat_clusters,umap.df,median)
        umap.df$Predicted_id <- factor(umap.df$Predicted_id)
        umap.df$seurat_clusters<- factor(umap.df$seurat_clusters, levels =unique(umap.df$seurat_clusters))
        
        p2 <- ggplot(umap.df , aes(x=UMAP_1, y=UMAP_2)) + 
          geom_point(aes(color=Predicted_id),shape=16, size=1, show.legend = NA, alpha=0.8)+
          xlab("umap_1")+ylab("umap_2")+labs(color="cluster")+
          mytheme_classic+expand_limits(x=c(min(umap.df$UMAP_1),max(umap.df$UMAP_1)),y=c(min(umap.df$UMAP_2),max(umap.df$UMAP_2)))
        
        
        ggsave(filename = paste0("pdf/1-Loading/",projet,"/UMAP_Predicted_global.pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
        
        
      }
      
      
      QCFiltered.data<-Raw.data
      QCFiltered.data
    }
    
    
    if(Seurat==FALSE){QCFiltered.data <- Raw.data}
    
    
  }
  
  if(hto==TRUE){
    table2<-table(QCFiltered.data$consensuscall.global)
    
    print(table1)
    print(table2)
    
    
    capture.output(table1, file = paste0("pdf/1-Loading/",projet,"/hto_demux.txt"), append = TRUE)
    capture.output(table2, file = paste0("pdf/1-Loading/",projet,"/hto_demux_final_after_qc.txt"), append = TRUE)
    
  }
  
  QCFiltered.data
  
}

###############################################


"%||%" <- function(a, b) {
  if (!is.null(a))
    a
  else
    b
}
geom_flat_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "dodge",
           trim = TRUE,
           scale = "area",
           show.legend = NA,
           inherit.aes = TRUE,
           ...) {
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(trim = trim,
                    scale = scale,
                    ...)
    )
  }
GeomFlatViolin <-
  ggproto(
    "GeomFlatViolin",
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)
      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        dplyr::group_by(.data = ., group) %>%
        dplyr::mutate(
          .data = .,
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )
    },
    draw_group = function(data, panel_scales, coord)
    {
      # Find the points for the line to go all the way around
      data <- base::transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
      # Make sure it's sorted properly to draw the outline
      newdata <-
        base::rbind(
          dplyr::arrange(.data = base::transform(data, x = xminv), y),
          dplyr::arrange(.data = base::transform(data, x = xmaxv), -y)
        )
      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])
      ggplot2:::ggname("geom_flat_violin",
                       GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },
    draw_key = draw_key_polygon,
    default_aes = ggplot2::aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),
    required_aes = c("x", "y")
  )









seurat_marker<-function(seurat,gene_list){  
  
  seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  
  seurat.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
  
  
  top1 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
  top2 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  top5 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  top10 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  top50<- seurat.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
  
  
  print(DoHeatmap(seurat, features = top10$gene) + NoLegend())
  
  print(FeaturePlot(seurat, features = top1$gene, reduction = "umap", pt.size = 0.5))
  
  print(FeaturePlot(seurat, features = gene_list, reduction = "umap", pt.size = 0.5))
  
  return(list(seurat,seurat.markers))
}






plot_save_umap<-function(seurat,genes=c("Slc1a3","Calb1","Pax2","Rbfox3","Lhx9","Sox10","Pvalb","Ret","Nrgn","Calb2","Eomes"),hto=TRUE){
  
  

    dir.create(paste0("pdf/1-Loading/",seurat@project.name,"/all/"))
    p1<-DimPlot(seurat,group.by = "region",shuffle = TRUE)
    
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_region.pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p2<-DimPlot(seurat,group.by = "ML",shuffle = TRUE)
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_ML.pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p3<-DimPlot(seurat,group.by = "AP",shuffle = TRUE)
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_AP.pdf"), plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p4<-FeaturePlot(seurat,features = genes,max.cutoff="q99")
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_Genes_markers.pdf"), plot = p4, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p5<-FeaturePlot(seurat,features = c("nFeature_RNA","percent.mito"))
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_Nfeatures_percent_mito.pdf"), plot = p5, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
  
    if(hto==TRUE){
      
  
 
      dir.create(paste0("pdf/1-Loading/",seurat@project.name,"/split/"))
               
    p1<-DimPlot(seurat,group.by = "region",shuffle = TRUE, split.by="consensuscall.global")
    
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/split/Umap_region.pdf"),plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p2<-DimPlot(seurat,group.by = "ML",shuffle = TRUE, split.by="consensuscall.global")
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/split/Umap_ML.pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p3<-DimPlot(seurat,group.by = "AP",shuffle = TRUE, split.by="consensuscall.global")
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/split/Umap_AP.pdf"), plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p4<-FeaturePlot(seurat,features = genes,max.cutoff="q99", split.by="consensuscall.global",by.col = FALSE)
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/split/Umap_Genes_markers.pdf"), plot = p4, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
    p5<-FeaturePlot(seurat,features = c("nFeature_RNA","percent.mito"), split.by="consensuscall.global")
    ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/split/Umap_Nfeatures_percent_mito.pdf"), plot = p5, device =cairo_pdf, width = 297, height = 210, units = "mm")
    
  }
  
                  
}




plot_save_umap_database<-function(seurat,genes=c("Slc1a3","Calb1","Pax2","Rbfox3","Lhx9","Sox10","Pvalb","Ret","Nrgn","Calb2","Eomes"),meta=NULL){
  
  
  
  dir.create(paste0("pdf/1-Loading/",seurat@project.name,"/all/"))
  
  
  p1<-DimPlot(seurat,group.by = meta[1],shuffle = TRUE)
  ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_",meta[1],".pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  p2<-DimPlot(seurat,group.by = meta[2],shuffle = TRUE)
  ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_",meta[2],".pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  p3<-DimPlot(seurat,group.by = meta[3],shuffle = TRUE)
  ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_",meta[3],".pdf"), plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  p4<-DimPlot(seurat,group.by = meta[4],shuffle = TRUE)
  ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_",meta[4],".pdf"), plot = p4, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  p5<-FeaturePlot(seurat,features = genes,max.cutoff="q99")
  ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_Genes_markers.pdf"), plot = p5, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  p6<-FeaturePlot(seurat,features = c("nFeature_RNA","percent.mito"))
  ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_Nfeatures_percent_mito.pdf"), plot = p6, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  p7<-DimPlot(seurat,shuffle = TRUE)
  ggsave(filename = paste0("pdf/1-Loading/",seurat@project.name,"/all/Umap_cluster_seurat.pdf"), plot = p7, device =cairo_pdf, width = 297, height = 210, units = "mm")
  
  
}


plot_umap_singlet<-function(x,Y=NULL,ref.cl=NULL,ref.cl.df=NULL,replicat=FALSE,features=NULL,merge=FALSE,res=2,meta=NULL,fill= hue_pal()(7)[1], color= hue_pal()(7)[1],datasets,dims=1:50,vars.to.regress=c(c("nFeature_RNA","percent.mito")),genes=c("Neurod2","Gdf10","Slc1a3","Calb1","Pax2","Ptf1a","Atoh1","Pax6","Rbfox3","Lhx9","Sox10","Pvalb","Ret","Nrgn","Calb2","Eomes","Slc6a5", "Sst","Prkcd", "Sorcs3", "Ptprk","Nxph1", "Cdh22","Htr2a", "Edil3", "Slc6a5","Gdf10","Slc17a6","Slc17a7","Gad2","Gad1","Lhx9","Bcl11b","Snap25","Cx3cl1","C1qa","Slc6a5","Slc6a1","Kit")){
  x <- tryCatch(
    {
      
      genes<-unique(genes)
      x$age[x$age%in%c("P00")]<-"P0"
      x$age[x$age%in%c("P07")]<-"P7"
      x$age[x$age%in%c("P04")]<-"P4"
      
      x$age<-factor(x$age,levels = c("P0","P4","P5","P7","P10","P11","P14","P16","P23","P30"))
      
      
      if(merge==FALSE){
        
        
        if(replicat==TRUE){
          
          x$datasets<-paste(unique(x$datasets),datasets,unique(x$age),unique(x$replicat),sep="_")
          
        }
        
        if(replicat==FALSE){
          
          
          x$datasets<-paste(unique(x$datasets),datasets,unique(x$age),sep="_")
          
        }
        
      }  
      
      if(merge==TRUE){
        
        
        if(replicat==TRUE){
          x$datasets<-paste(datasets,"all",sep="_")
          
        }
        
        if(replicat==FALSE){
          
          x$datasets<-paste(datasets,"all",sep="_")
          
        }
        
      }  
      
      
      
      
      
      
      print(unique(x$datasets))
      
      dir.create(paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets)))
      dir.create(paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/"))
      
      
      
      mytheme_classic <- theme_classic()+ 
        theme(axis.text.y = element_text(size=14, color="black"), axis.ticks = element_line(size=1, color="black"), axis.title.y = element_text(size=18, color="black"),
              axis.text.x = element_text(size=14, color="black"), axis.title.x=element_text(size=18) ,panel.background =   element_blank(), legend.text = element_text(size=14))
      
      
      ######## Plot final QC metrics
      p1 <- ggplot(x@meta.data, aes(x=datasets, y=nCount_RNA))+
        geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
        geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
        geom_point(aes(x=as.numeric(factor(datasets))-0.3), fill="black", color="black", shape=16, size=0.25,
                   alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
        xlab(NULL)+ ggtitle("Number of UMIs")+
        mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
      
      p2 <- ggplot(x@meta.data, aes(x=datasets, y=nFeature_RNA))+
        geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
        geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
        geom_point(aes(x=as.numeric(factor(datasets))-0.3), fill="black", color="black", shape=16, size=0.25,
                   alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
        xlab(NULL)+ ggtitle("Number of Genes")+
        mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
      
      p3 <- ggplot(x@meta.data, aes(x=datasets, y=percent.mito))+
        geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
        geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
        geom_point(aes(x=as.numeric(factor(datasets))-0.3), fill="black", color="black", shape=16, size=0.25,
                   alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
        xlab(NULL)+ ggtitle("Mitochondrial genes")+
        mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
      
      p4 <- ggplot(x@meta.data, aes(x=datasets, y=percent.ribo))+
        geom_flat_violin(fill=alpha(fill,0.5), color=color, trim= T, show.legend = NA)+
        geom_boxplot(width=0.025, outlier.size = 0.25, outlier.alpha = 0.8)+
        geom_point(aes(x=as.numeric(factor(datasets))-0.3), fill="black", color="black", shape=16, size=0.25,
                   alpha=0.5, position = position_jitter(width=0.25,seed = 12345), show.legend = F)+
        xlab(NULL)+ ggtitle("Ribosomal genes")+
        mytheme_classic +theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), title = element_text(size = 14,hjust = 0.5))
      
      p <- plot_grid(plotlist =  list(p1,p2, p3, p4), ncol=2, align = "hv", rel_widths = c(1, 1, 1, 1),rel_heights =  c(1, 1, 1, 1))
      
      ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/QC.pdf"), plot = p, device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      
      
      
      
      DefaultAssay(x) <- "RNA"
      
      
      x <- SCTransform(x, vars.to.regress = vars.to.regress)
      
      # Run the standard workflow for visualization and clustering
      x <- RunPCA(x)
      
      # Run UMAP
      x <- RunUMAP(x, dims = dims)
      
      
      # Determine the K-nearest neighbor graph
      x <- FindNeighbors(object =x,  dims = dims)
      
      # Determine the clusters for various resolutions                                
      x <- FindClusters(object =x,resolution = res)
      
      
      
      
      
      if(sum(colnames(x@meta.data)==meta[1])==1){
        
        if(sum(is.na(x@meta.data[,meta[1]])>0)){
            x@meta.data[,meta[1]][is.na(x@meta.data[,meta[1]])]<-"unknown"
        }
        
        p1<-DimPlot(x,group.by = meta[1],shuffle = TRUE)
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Umap_",meta[1],".pdf"), plot = p1, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
      }
      
      if(sum(colnames(x@meta.data)==meta[2])==1){
        if(sum(is.na(x@meta.data[,meta[2]])>0)){
          x@meta.data[,meta[2]][is.na(x@meta.data[,meta[2]])]<-"unknown"
        }
        
        p2<-DimPlot(x,group.by = meta[2],shuffle = TRUE)
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Umap_",meta[2],".pdf"), plot = p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      if(sum(colnames(x@meta.data)==meta[3])==1){
        if(sum(is.na(x@meta.data[,meta[3]])>0)){
          x@meta.data[,meta[3]][is.na(x@meta.data[,meta[3]])]<-"unknown"
        }
        
        p3<-DimPlot(x,group.by = meta[3],shuffle = TRUE)
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Umap_",meta[3],".pdf"), plot = p3, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      if(sum(colnames(x@meta.data)==meta[4])==1){
        
        if(sum(is.na(x@meta.data[,meta[4]])>0)){
          x@meta.data[,meta[4]][is.na(x@meta.data[,meta[4]])]<-"unknown"
        }
        
        p4<-DimPlot(x,group.by = meta[4],shuffle = TRUE)
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Umap_",meta[4],".pdf"), plot = p4, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      if(sum(colnames(x@meta.data)==meta[5])==1){
        if(sum(is.na(x@meta.data[,meta[5]])>0)){
          x@meta.data[,meta[5]][is.na(x@meta.data[,meta[5]])]<-"unknown"
        }
        
        p0<-DimPlot(x,group.by = meta[5],shuffle = TRUE)
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Umap_",meta[5],".pdf"), plot = p0, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      
      if(merge==TRUE){
        

        pmerge<-DimPlot(x,group.by = "datasets_org",shuffle = TRUE)
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Umap_datasets_org.pdf"), plot = pmerge, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      
      } 
      
      
      
      
      DefaultAssay(x) <- "SCT"
      
      p5<-FeaturePlot(x,features = genes,min.cutoff="q10",max.cutoff="q99",combine=FALSE,cols = rev(brewer.pal(5,"Spectral")))
      
      pdf(paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Pie_chart_and_Umap_Genes_markers.pdf"),width = 20, height = 20)
      
      pie =data.frame(percent=table(x@active.ident)/sum(table(x@active.ident))*100)
      
      sink(paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/pie_chart_pourcentage.txt"))
      print(pie)
      
      sink()
      
      
      cols = wes_palette("GrandBudapest1", 4, type = "discrete")[4:1]
      print(pie(pie$percent.Freq, labels = pie$percent.Var1), col = cols)
      
    
      print(p5)
      dev.off()
      
      p6<-FeaturePlot(x,features = c("nFeature_RNA","percent.mito"),min.cutoff="q10",max.cutoff="q99",combine=FALSE,cols = rev(brewer.pal(5,"Spectral")))
      ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Umap__mito_nFeatures.pdf"), plot = p6[[1]]+p6[[2]], device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      p7<-DimPlot(x,shuffle = TRUE)
      ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Umap_cluster_res=",res,".pdf"), plot = p7, device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      
      # Assign cell-cycle scores based on Tirosh et al, 2015 gene list
      s.genes <- c("Mcm5", "Pcna", "Tym5", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2", "Rap2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", " Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8")
      g2m.genes <- c("Hmgb2", "Ddk1","Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", " Cenpk", "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "kif20b", "Hjurp", "Cdca3", "Hn1", "Cdc20", "Ttk", "Cdc25c", "kif2c", "Rangap1", "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")
      
      x <- CellCycleScoring(object = x,
                            s.features = s.genes,
                            g2m.features = g2m.genes)
      
      
      
      p0<-DimPlot(x,
                  reduction = "umap",
                  group.by = "Phase",
                  
                  label=T,
                  label.size = 4 )
      p0
      
      ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Cell_cycle_scores.pdf"), plot = p0, device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      
      
      
      count_meta<-function(meta_count,name,x){
        
        
        
        print(name)
        
        m_meta_name<-t(round(t(table(meta_count,x@active.ident))/colSums(table(meta_count,x@active.ident))*100))
        
        levels_meta_name<-factor(rownames(m_meta_name))
        levels_cluster<-factor(colnames(m_meta_name),levels=sort(as.integer(colnames(m_meta_name))))
        
        m_meta_name <- data.frame(cluster=as.factor(col(m_meta_name)),meta_name=as.factor(row(m_meta_name)),n=as.vector(m_meta_name))
        
        
        m_meta_name$cluster<- levels_cluster[m_meta_name$cluster]
        m_meta_name$meta_name<-levels_meta_name[m_meta_name$meta_name]
        
        
        p1 <- ggplot(m_meta_name,aes(x=cluster,y=meta_name,fill=n)) + 
          geom_tile(aes(fill=n)) +
          geom_text(aes(label=n,colour="white",size=0.2)) +
          scale_fill_viridis_c() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_equal()+ggtitle(paste0(name))
        
        
        
        ####
        
        
        m_meta_name<-round(table(meta_count,x@active.ident)/rowSums(table(meta_count,x@active.ident))*100)
        #m_meta_name<-table(meta_count,x@active.ident)
        
        
        levels_meta_name<-factor(rownames(m_meta_name))
        levels_cluster<-factor(colnames(m_meta_name),levels=sort(as.integer(colnames(m_meta_name))))
        
        m_meta_name <- data.frame(cluster=as.factor(col(m_meta_name)),meta_name=as.factor(row(m_meta_name)),n=as.vector(m_meta_name))
        m_meta_name$cluster<- levels_cluster[m_meta_name$cluster]
        m_meta_name$meta_name<-levels_meta_name[m_meta_name$meta_name]
        
        
        p2 <- ggplot(m_meta_name,aes(x=cluster,y=meta_name,fill=n)) + 
          geom_tile(aes(fill=n)) +
          geom_text(aes(label=n,colour="white",size=0.2)) +
          scale_fill_viridis_c() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_equal()+ggtitle(paste0(name))
        
        
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/heatmap_",name,".pdf"), plot = p1+p2, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
        
        p1+p2
        
      }
      
      
      
      if(sum(colnames(x@meta.data)==meta[1])==1&length(levels(factor(x@meta.data[,meta[1]])))>1){
        
        count_meta(meta_count=x@meta.data[,meta[1]],name=meta[1],x=x)
        
      }
      
      
      if(sum(colnames(x@meta.data)==meta[2])==1&length(levels(factor(x@meta.data[,meta[2]])))>1){
        
        count_meta(meta_count=x@meta.data[,meta[2]],name=meta[2],x=x)
        
      }
      
      
      
      if(sum(colnames(x@meta.data)==meta[3])==1&length(levels(factor(x@meta.data[,meta[3]])))>1){
        
        count_meta(meta_count=x@meta.data[,meta[3]],name=meta[3],x=x)
        
      }
      
      if(sum(colnames(x@meta.data)==meta[4])==1&length(levels(factor(x@meta.data[,meta[4]])))>1){
        
        count_meta(meta_count=x@meta.data[,meta[4]],name=meta[4],x=x)
        
      }
      
      
      if(sum(colnames(x@meta.data)==meta[5])==1&length(levels(factor(x@meta.data[,meta[5]])))>1){
        
        count_meta(meta_count=x@meta.data[,meta[5]],name=meta[5],x=x)
        
      }
      
      
  
      
      
      
      
      
      
      
      
      
      
       
      
      
      
      
      
      
      
      
      
      p0<-DimPlot(x,shuffle = TRUE,label=TRUE)+NoLegend()
      
      
      Vascular <-c("Mgp","Slc47a1","Dapl1","Slc6a13","Ptgds","Dcn","Fam180a","Igf2","Sema3g","Gkn3","Kcnj8","Anpep","Acta2","Mustn1","Ecscr","Gadd45a","Ly6c1","Scgb3a1","Cxcl2")
      
      
      genes.list_Vascular <- list(Vascular)
      enrich.name <- "Vascular"
      x <- AddModuleScore(x,
                          features = genes.list_Vascular,
                          pool = NULL,
                          n.bin = 5,
                          seed = 1,
                          ctrl = length(genes.list_Vascular),
                          k = FALSE,
                          name = enrich.name,
                          random.seed = 1)
      
      
      p6<-FeaturePlot(x,features = c("Vascular1"),min.cutoff = "q10", max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
      ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Vascular_scores.pdf"), plot = p0+p6, device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      
      
      
      
      
      Choroid <-c("Kcne2","Kcnj13","Cfap126","Morn5")
      
      if(sum(Choroid%in%rownames(x@assays$SCT))>=2){
        
        genes.list_Choroid <- list(Choroid)
        enrich.name <- "Choroid"
        x <- AddModuleScore(x,
                            features = genes.list_Choroid,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Choroid),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        p7<-FeaturePlot(x,features = c("Choroid1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Choroid_scores.pdf"), plot = p0+p7, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
        
      }
      
      
      
      
      endothelial <-c("Clec1b","Prl3b1","Emcn","Gimap4","Aplnr","Ecscr")
      
      if(sum(endothelial%in%rownames(x@assays$SCT))>=2){
        
        genes.list_endothelial <- list(endothelial)
        enrich.name <- "endothelial"
        x <- AddModuleScore(x,
                            features = genes.list_endothelial,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_endothelial),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p8<-FeaturePlot(x,features = c("endothelial1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/endothelial_scores.pdf"), plot = p0+p8, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      Astro <-c("Gdf10","Aqp4","Aldh1l1","Gfap","Gja1","Glul","Nfia","Mt1","S100b","Slc1a2","Slc1a3","Sox9")
      
      if(sum(Astro%in%rownames(x@assays$SCT))>=2){
        
        genes.list_Astro <- list(Astro)
        enrich.name <- "Astro"
        x <- AddModuleScore(x,
                            features = genes.list_Astro,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Astro),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p9<-FeaturePlot(x,features = c("Astro1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Astro_scores.pdf"), plot = p0+p9, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      Purkinje <-c("Calb1","Ebf1","Plcb3","Plcb4","Rorb","Sema3a","Ppp1r17")
      
      if(sum(Purkinje%in%rownames(x@assays$SCT))>=2){
        
        genes.list_Purkinje <- list(Purkinje)
        enrich.name <- "Purkinje"
        x <- AddModuleScore(x,
                            features = genes.list_Purkinje,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Purkinje),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p10<-FeaturePlot(x,features = c("Purkinje1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Purkinje_scores.pdf"), plot = p0+p10, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      Neurons_Dcn_glu <-c("Lhx9","Bcl11b","Tbr1","Slc17a6")
      if(sum(Neurons_Dcn_glu%in%rownames(x@assays$SCT))>=2){
        
        
        genes.list_Neurons_Dcn_glu <- list(Neurons_Dcn_glu)
        enrich.name <- "Neurons_Dcn_glu"
        x <- AddModuleScore(x,
                            features = genes.list_Neurons_Dcn_glu,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Neurons_Dcn_glu),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p11<-FeaturePlot(x,features = c("Neurons_Dcn_glu1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Neurons_Dcn_glu_scores.pdf"), plot = p0+p11, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      # Neurons_Dcn_gaba <-c("Slc6a1","Kit")
      # 
      # genes.list_Neurons_Dcn_gaba <- list(Neurons_Dcn_gaba)
      # enrich.name <- "Neurons_Dcn_gaba"
      # x <- AddModuleScore(x,
      #                     features = genes.list_Neurons_Dcn_gaba,
      #                     pool = NULL,
      #                     n.bin = 5,
      #                     seed = 1,
      #                     ctrl = length(genes.list_Neurons_Dcn_gaba),
      #                     k = FALSE,
      #                     name = enrich.name,
      #                     random.seed = 1)
      # 
      # 
      # 
      # p12<-FeaturePlot(x,features = c("Neurons_Dcn_gaba1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
      # ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Neurons_Dcn_gaba_scores.pdf"), plot = p0+p12, device =cairo_pdf, width = 297, height = 210, units = "mm")
      # 
      
      
      
      
      
      
      
      
      Oligo <-c("Olig2","Sox10")
      
      if(sum(Oligo%in%rownames(x@assays$SCT))>=2){
        
        genes.list_Oligo <- list(Oligo)
        enrich.name <- "Oligo"
        x <- AddModuleScore(x,
                            features = genes.list_Oligo,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Oligo),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p13<-FeaturePlot(x,features = c("Oligo1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Oligo_scores.pdf"), plot = p0+p13, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      Microglia <-c("Cx3cl1","C1qa")
      if(sum(Microglia%in%rownames(x@assays$SCT))>=2){
        
        
        genes.list_Microglia <- list(Microglia)
        enrich.name <- "Microglia"
        x <- AddModuleScore(x,
                            features = genes.list_Microglia,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Microglia),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p14<-FeaturePlot(x,features = c("Microglia1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Microglia_scores.pdf"), plot = p0+p14, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      Granule_cells_mature<-c("Rbfox3","Gabra6")
      if(sum(Granule_cells_mature%in%rownames(x@assays$SCT))>=2){
        
        
        genes.list_Granule_cells_mature<- list(Granule_cells_mature)
        enrich.name <- "Granule_cells_mature"
        x <- AddModuleScore(x,
                            features = genes.list_Granule_cells_mature,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Granule_cells_mature),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p15<-FeaturePlot(x,features = c("Granule_cells_mature1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Granule_cells_mature_scores.pdf"), plot = p0+p15, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      
      Granule_cells_immature<-c("Atoh1","Pax6")
      if(sum(Granule_cells_immature%in%rownames(x@assays$SCT))>=2){
        
        
        genes.list_Granule_cells_immature <- list(Granule_cells_immature)
        enrich.name <- "Granule_cells_immature"
        x <- AddModuleScore(x,
                            features = genes.list_Granule_cells_immature,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Granule_cells_immature),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p15bis<-FeaturePlot(x,features = c("Granule_cells_immature1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Granule_cells_immature_scores.pdf"), plot = p0+p15bis, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      
      Golgi<-c("Grm2", "Sst")
      if(sum(Golgi%in%rownames(x@assays$SCT))>=2){
        
        
        genes.list_Golgi <- list(Golgi)
        enrich.name <- "Golgi"
        x <- AddModuleScore(x,
                            features = genes.list_Golgi,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Golgi),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p16<-FeaturePlot(x,features = c("Golgi1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Golgi_scores.pdf"), plot = p0+p16, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
        
      }
      PLI<-c("Htr2a", "Edil3","Calb2","Aldh1a3")
      
      if(sum(PLI%in%rownames(x@assays$SCT))>=2){
        
        genes.list_PLI <- list(PLI)
        enrich.name <- "PLI"
        x <- AddModuleScore(x,
                            features = genes.list_PLI,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_PLI),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p17<-FeaturePlot(x,features = c("PLI1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/PLI_scores.pdf"), plot = p0+p17, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      MLI<-c("Prkcd", "Ptprk", "Pvalb","Ret")
      
      if(sum(MLI%in%rownames(x@assays$SCT))>=2){
        
        genes.list_MLI <- list(MLI)
        enrich.name <- "MLI"
        x <- AddModuleScore(x,
                            features = genes.list_MLI,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_MLI),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p18<-FeaturePlot(x,features = c("MLI1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/MLI_scores.pdf"), plot = p0+p18, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      
      Gaba_intermediate_prog<-c("Ptf1a", "Pax2")
      
      if(sum(Gaba_intermediate_prog%in%rownames(x@assays$SCT))>=2){
        
        genes.list_Gaba_intermediate_prog <- list(Gaba_intermediate_prog)
        enrich.name <- "Gaba_intermediate_prog"
        x <- AddModuleScore(x,
                            features = genes.list_Gaba_intermediate_prog,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_Gaba_intermediate_prog),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p19<-FeaturePlot(x,features = c("Gaba_intermediate_prog1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Gaba_intermediate_prog_scores.pdf"), plot = p0+p19, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      
      pan_gaba<-c("Gad1", "Gad2")
      
      if(sum(pan_gaba%in%rownames(x@assays$SCT))>=2){
        
        
        genes.list_pan_gaba <- list(pan_gaba)
        enrich.name <- "pan_gaba"
        x <- AddModuleScore(x,
                            features = genes.list_pan_gaba,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_pan_gaba),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p20<-FeaturePlot(x,features = c("pan_gaba1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/pan_gaba_scores.pdf"), plot = p0+p20, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      pan_glu<-c("Slc17a6", "Slc17a7","Rbfox3")
      
      if(sum(pan_glu%in%rownames(x@assays$SCT))>=2){
        
        
        genes.list_pan_glu <- list(pan_glu)
        enrich.name <- "pan_glu"
        x <- AddModuleScore(x,
                            features = genes.list_pan_glu,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_pan_glu),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p21<-FeaturePlot(x,features = c("pan_glu1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/pan_glu_scores.pdf"), plot = p0+p21, device =cairo_pdf, width = 297, height = 210, units = "mm")
      }
      
      pan_neuron<-c("Slc17a6","Slc6a5","Slc17a7","Rbfox3","Gad1","Gad2","Snap25")
      
      if(sum(pan_neuron%in%rownames(x@assays$SCT))>=2){
        
        
        genes.list_pan_neuron <- list(pan_neuron)
        enrich.name <- "pan_neuron"
        x <- AddModuleScore(x,
                            features = genes.list_pan_neuron,
                            pool = NULL,
                            n.bin = 5,
                            seed = 1,
                            ctrl = length(genes.list_pan_neuron),
                            k = FALSE,
                            name = enrich.name,
                            random.seed = 1)
        
        
        
        p22<-FeaturePlot(x,features = c("pan_neuron1"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
        ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/pan_neuron_scores.pdf"), plot = p0+p22, device =cairo_pdf, width = 297, height = 210, units = "mm")
        
      }
      
      
      
      
      
      
      
      
      p23<-FeaturePlot(x,features = c("Eomes"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
      ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/UBC_scores.pdf"), plot = p0+p23, device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      
      p24<-FeaturePlot(x,features = c("Slc6a5"),min.cutoff = "q10",max.cutoff="q99",cols = rev(brewer.pal(5,"Spectral")))
      ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/Pan_glycinergic.pdf"), plot = p0+p24, device =cairo_pdf, width = 297, height = 210, units = "mm")
      
      
      
      
      
      
      
      
      
      
      
      # 
      # #hicat_pipeline
      # hicat_pipeline <- function(X,Y,ref.cl,ref.cl.df,features){
      # 
      #   attr(X,"norm") <- logCPM(X@assays$RNA@counts[features,])
      #   
      #   print("it's running")
      #   
      #   attr(X,"mapping_results") <- map_sampling(train.dat    = attr(Y,"norm"),
      #                                             train.cl     = ref.cl,
      #                                             test.dat     = attr(X,"norm"),
      #                                             markers      = attr(Y,"marker_genes"),
      #                                             markers.perc = 0.8,
      #                                             iter         = 100)
      #   
      #   attr(X,"map.df") <- attr(X,"mapping_results")$map.df %>%
      #     mutate(pred.cl = as.numeric(as.character(pred.cl))) %>%
      #     left_join(ref.cl.df, by = c("pred.cl" = "class_id")) 
      #   return(X)
      # }
      # 
      # print("run_hicat")
      # #Apply hicat
      # hicat_x <- hicat_pipeline(x, Y=Y,ref.cl=ref.cl,ref.cl.df=ref.cl.df,features=features)
      # print("finish_hicat")
      # 
      # 
      # #AddMetadata to P8_wt_GlutNs
      # x <- AddMetaData(object = x,
      #                      metadata = list(attr(hicat_x,"map.df")$pred.cl,attr(hicat_x,"map.df")$prob,attr(hicat_x,"map.df")$class_label,
      #                                      attr(hicat_x,"map.df")$class_color), 
      #                      col.name = c("class_id","class_prob","class_label","class_color"))
      # 
      # 
      # 
      # p30<-DimPlot(x,group.by = "class_label",label = TRUE)+NoLegend()
      # 
      # ggsave(filename = paste0("pdf/2-integration_postnatal/",datasets,"/",unique(x$datasets),"/singlets/mapping.pdf"), plot = p30, device =cairo_pdf, width = 297, height = 210, units = "mm")
      # 
      # 
      
      
      
      
      
      
      
      
      return(x)
      
      
      
      
    },
    error=function(cond) {
      message(paste("function does not seem to exist:", unique(x$datasets)))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(x)
    }
  )    
  return(x)
  
  
}






