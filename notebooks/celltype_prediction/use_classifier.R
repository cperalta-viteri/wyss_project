


source("notebooks/01_data_collection/lib_sce.R")
source("notebooks/celltype_prediction/lib_cell_classifier.R")


x <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5")
x <- scuttle::logNormCounts(x,center_size_factors=FALSE,size.factors=x$sum/1e4)

x$p1 <- torch_predict_cell_type(counts(x),"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.torch")
x$pred1 <- colnames(x$p1)[max.col(x$p1)]

x$p2 <- torch_predict_cell_type(logcounts(x),"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.torch")
x$pred2 <- colnames(x$p2)[max.col(x$p2)]

x$y_true <- paste(x$class_label,x$neighborhood_label,x$subclass_label,sep=" --- ")
mean(x$pred1 != x$y_true)
mean(x$pred2 != x$y_true)


table(sub("(.*) --- .*","\\1",x$y_true),sub("(.*) --- .*","\\1",x$pred1))



#-#-#-#-#-#-#-#-#-#-#-#-#
# Test on human data
#-#-#-#-#-#-#-#-#-#-#-#-#
H <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/human_m1_10x.h5")
H$sum <- local({
  m <- scuttle::perCellQCMetrics(H,BPPARAM=BiocParallel::MulticoreParam(5))
  m$sum
})
H <- scuttle::logNormCounts(H,center_size_factors=FALSE,size.factors=H$sum/1e4)
H <- H[,1:1000]

P <- torch_predict_cell_type(logcounts(H),"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.torch",ignore.case = TRUE)
pred <- factor(colnames(P)[max.col(P)],colnames(P))
g <- local({
  n <- table(label=pred[drop=TRUE]) |> as.data.frame()
  pat <- "(((.*) --- (.*)) --- (.*))"
  g <- graph_from_data_frame(
    rbind(
      cbind("root",sub(pat,"\\3",n$label)),
      cbind(sub(pat,"\\3",n$label),sub(pat,"\\2",n$label)),
      cbind(sub(pat,"\\2",n$label),sub(pat,"\\1",n$label))
    ) |> unique()
  )
  V(g)$size <- n$Freq[match(V(g)$name,n$label)]
  V(g)$size[is.na(V(g)$size)] <- 0
  V(g)$label <- sub(".* --- ","",V(g)$name)
  g
})
ggraph(g,"partition",weight=V(g)$size,circular=TRUE,const.area=FALSE) +
  geom_node_arc_bar(aes(fill=depth)) +
  geom_node_text(aes(label=label),size=3) +
  coord_equal()



ggpred <- function(P) {
  library(ggplot2)
  
  pred <- factor(colnames(P)[max.col(P)],colnames(P))
  H <- data.frame(label = pred)
  pat <- "(((.*) --- (.*)) --- (.*))"
  H$lev1 <- H$lev2 <- H$lev3 <- H$lev1_lab <- H$lev2_lab <- H$lev3_lab <- H$label
  levels(H$lev1) <- sub(pat,"\\3",levels(H$lev1))
  levels(H$lev2) <- sub(pat,"\\2",levels(H$lev2))
  levels(H$lev3) <- sub(pat,"\\1",levels(H$lev3))
  levels(H$lev1_lab) <- sub(pat,"\\3",levels(H$lev1_lab))
  levels(H$lev2_lab) <- sub(pat,"\\4",levels(H$lev2_lab))
  levels(H$lev3_lab) <- sub(pat,"\\5",levels(H$lev3_lab))
  
  # Show prediction proportions
  ggplot(H,aes(fill=lev1)) + 
    geom_bar(aes(x=1,group=lev1),width=1,color="black") +
    geom_bar(aes(x=2,group=lev2),width=1,color="black") + 
    geom_bar(aes(x=3,group=lev3),width=1,color="black") +
    stat_count(geom="text",aes(x=1,label=lev1_lab,group=lev1),position = position_stack(vjust=0.5),size=3) +
    stat_count(geom="text",aes(x=2,label=lev2_lab,group=lev2),position = position_stack(vjust=0.5),size=3) +
    stat_count(geom="text",aes(x=3,label=lev3_lab,group=lev3),position = position_stack(vjust=0.5),size=3) +  
    #stat_count(geom="text_repel",aes(x=3.5,label=lev3_lab,group=lev3),position = position_stack(vjust=0.5),size=3) +
    coord_polar("y") +
    theme_void()  
}


ggpred(P)





