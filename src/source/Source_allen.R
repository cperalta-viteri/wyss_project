library(rjson)
library(jpeg)
library(IRanges)
library(grid)


getStructures <- function(acronyms) {
  url <- "http://api.brain-map.org/api/v2/data/Structure/query.csv?criteria=[graph_id$eq17]"
  url <- paste0(url,sprintf("[acronym$in'%s']",paste(acronyms,collapse="','")))
  read.csv(url)  
}

#' @export
getSectionDatasets <- function(genes=NULL,ages=NULL,plan_of_section_id=NA) {
  url <- "http://api.brain-map.org/api/v2/data/SectionDataSet/query.csv?criteria=[failed$eqfalse],products[id$eq3]"
  url <- paste0(url,sprintf(",specimen(donor(age[name$in'%s']))",paste(ages,collapse="','")))
  url <- paste0(url,sprintf(",genes[acronym$in'%s']",paste(genes,collapse="','")))
  if (!is.na(plan_of_section_id)) url <- paste0(url,sprintf(",[plane_of_section_id$eq%d]",plan_of_section_id))
  url <- paste0(url,"&order=genes.acronym,ages.embryonic+desc,ages.name")
  url <- paste0(url,'&num_rows=all')
  #url <- paste0(url,'&tabular=plane_of_section_id,ages.embryonic,ages.days,genes.acronym')
  read.csv(url)
}

#' @import IRanges
#' @import rjson
structure2image <- function(structure.ids,dataset.ids) {
  urls <- sprintf("http://api.brain-map.org/api/v2/structure_to_image/%d.json?structure_ids=%s",
                  dataset.ids,
                  paste(structure.ids,collapse=",")
  )
  as.data.frame(unlist(SplitDataFrameList(lapply(urls,function(u) {
    json <- fromJSON(file=u)
    as.data.frame(json$msg[[1]]$image_sync)
  }))))
}

#' @import IRanges
#' @import rjson
image2image <- function(image.id,x,y,dataset.ids) {
  urls <- sprintf("http://api.brain-map.org/api/v2/image_to_image/%d.json?x=%f&y=%f&section_data_set_ids=%d",
                  image.id,x,y,dataset.ids)
  as.data.frame(unlist(SplitDataFrameList(lapply(urls,function(u) {
    json <- fromJSON(file=u)
    as.data.frame(json$msg[[1]]$image_sync)
  }))))
}


getSectionImages <- function(ids) {
  url <- sprintf("http://api.brain-map.org/api/v2/data/SectionImage/query.csv?criteria=[id$in%s]&num_rows=all",paste(ids,collapse=","))
  read.csv(url)
}



#' @import jpeg
downloadImage <- function(path,x,y,width=400,height=400,downsample=2) {
  url <- sprintf("http://api.brain-map.org/cgi-bin/imageservice?downsample=%d&path=%s&top=%f&left=%f&width=%f&height=%f",
                 downsample,
                 as.character(path),
                 y-height*0.5*2^downsample,
                 x-width*0.5*2^downsample,
                 width,height
  )
  print(url)
  img <- readBin(url,what="raw",n=10000000)
  readJPEG(img)
}


#-#-#-#-#-#-#
# Main
#-#-#-#-#-#-#

#' Retreive and plot SectionImages from allenbrain atlas for given genes at given age
#' @param genes a character vector of gene names
#' @param structure.acronym name of the structure to zoom in
#' @param plan_of_section_id identifier of the section plane: 1: sagital; 2:coronal
#' @param downsample an integer controlling image resolution
#' @param img.size an integer controlling size of the image
#' @param ref_dataset_id identifier of the reference dataset
#' @export
allenbrainplot <- function(genes,ages=c("E11.5","E13.5","E15.5","E18.5","P4","P14","P28"),
                           structure.acronym="DPall",
                           plan_of_section_id=2,downsample=3,img.size=200,
                           ref_dataset_id=getSectionDatasets("Pax6","E15.5",plan_of_section_id)$id[1]) {
  # Retreive data
  x <- getSectionDatasets(genes,ages,plan_of_section_id)
  print(x)
  
  # find the reference image for the given structure in dataset 100083728=$id[1]
  ref_img <- structure2image(getStructures(structure.acronym)$id,ref_dataset_id)
  y <- image2image(ref_img$section_image_id,ref_img$x,ref_img$y,x$id)
  
  xy <- merge(x,y,by.x="id",by.y="section_data_set_id",all=TRUE)
  z <- getSectionImages(xy$section_image_id)
  xyz <- merge(xy,z,by.x="section_image_id",by.y="id",all=TRUE,suffixes=c(".xy",".z"))
  
  xyz$img <- mapply(downloadImage,xyz$path,xyz$x.xy+xyz$x.z,xyz$y.xy+xyz$y.z,downsample=round(downsample + log2(xyz$resolution / max(xyz$resolution))),MoreArgs=list(width=img.size,height=img.size),SIMPLIFY=FALSE)
  xyz$img.exp <- mapply(downloadImage,xyz$expression_path,xyz$x.xy+xyz$x.z,xyz$y.xy+xyz$y.z,downsample=round(downsample + log2(xyz$resolution / max(xyz$resolution))),MoreArgs=list(width=img.size,height=img.size),SIMPLIFY=FALSE)
  
  xyz$acronym <- factor(xyz$acronym,genes)
  xyz$name <- factor(xyz$name,rev(ages))
  xyz <- xyz[order(xyz$acronym,xyz$name),]
  xyz <- split(xyz,xyz$acronym)
  
  for(a in names(xyz)) {
    plot(c(0,2),c(0,nlevels(xyz[[a]]$name)),type="n",asp=1,main=a,xlab="",ylab="",bty="n",axes=FALSE)
    text(0,seq_len(nlevels(xyz[[a]]$name))-0.5,levels(xyz[[a]]$name),adj=c(1,0.5),cex=0.75)
    for(i in seq_len(nrow(xyz[[a]]))) {
      rasterImage(xyz[[a]]$img[[i]],0,as.integer(xyz[[a]]$name[i])-1,1,as.integer(xyz[[a]]$name[i]))
      rasterImage(xyz[[a]]$img.exp[[i]],1,as.integer(xyz[[a]]$name[i])-1,2,as.integer(xyz[[a]]$name[i]))
    }
  }
}










# Train a multiclass classifier on x.trn, y.trn and select 2xn most powerful variables (n up + n down).

# Train a multiclass classifier on x.trn, y.trn and select 2xn most powerful variables (n up + n down).
feature.selected.multiclassSVM <- function(x.trn,y.trn,x.tst=x.trn,n=3,...) {
  y.trn <- as.factor(y.trn)
  l <- as.vector((1/table(y.trn)[y.trn])) * (1 - table(seq_along(y.trn), y.trn))
  l <- l*100/(ncol(l)-1)
  w <- nrbm(ontologyLoss(x.trn,y.trn,l=l),...)
  W <- array(w,attr(w,"model.dim"),attr(w,"model.dimname"))
  
  
  global.model <- data.frame(feature.name = rownames(W),W=W,stringsAsFactors = FALSE)
  global.model$RK <- apply(-W,2,rank,ties.method="first")
  
  global.model$Z <- t((t(W)-colMeans(W))/matrixStats::colSds(W))
  global.model$PVAL <- 2*pnorm(abs(global.model$Z),lower.tail=FALSE)
  global.model$FDR <- apply(global.model$PVAL,2,p.adjust,method="fdr")
  global.model <- cbind(global.model,trn.avg=local({
    n <- rowsum(x.trn,y.trn)
    n <- n / as.vector(table(y.trn)[rownames(n)])
    t(n) 
  }))
  
  
  # x.trn.reduced_1 <- x.trn[,rowMin(global.model$RK)==1]
  # x.trn.reduced_2 <- x.trn[,rowMin(global.model$RK)==2]$
  # x.trn.reduced<-cbind(x.trn.reduced_1,x.trn.reduced_2)
  # 
  
  x.trn.reduced <- x.trn[,rowMin(global.model$RK)<=n]
  
  
  w <- nrbm(ontologyLoss(x.trn.reduced,y.trn,l=l),...)
  W <- array(w,attr(w,"model.dim"),attr(w,"model.dimname"))
  
  
  # x.tst.reduced_1 <- x.tst[,rowMin(global.model$RK)<=n]
  # x.tst.reduced_2 <- x.tst[,rowMin(global.model$RK)==2]
  # x.tst<-cbind(x.tst.reduced_1,x.tst.reduced_2)
  # 
  
  x.tst<-x.tst[,rowMin(global.model$RK)<=n]
  
  
  y <- predict(w,x.tst)
  attr(y,"w") <- w
  attr(y,"global.model") <- global.model
  attr(y,"genes") <- colnames(x.tst)
  
  y
}