

library(scml)
library(ABAData)
library(SingleCellExperiment)
library(luz)
library(ggplot2)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Predict human data with mouse model
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
x <- ABAHumanSSv4("data/ABAHumanSSv4/")
x <- scuttle::logNormCounts(x,size.factors=x$sum/1e4)
fm <- luz_load("data/aba_mouse_subclass_classifier.luz")
P <- predict_delayed(t(logcounts(x)),fm,ignore.case = TRUE)
pred <- factor(colnames(P)[max.col(P)],colnames(P))

# show predictions
n <- table(x$subclass_label,pred)
(n / rowSums(n)) |>
  as.data.frame() |>
  ggplot(aes(x=pred,y=Var1,label=sprintf("%.1f",100*Freq))) + 
    geom_tile(aes(fill=Freq)) +
    geom_text(size=3) + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

# show predictions with ggplot
ggpred <- function(P) {
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



library(tidygraph)
library(igraph)
library(ggraph)
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






