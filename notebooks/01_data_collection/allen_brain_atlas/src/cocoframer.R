

library(cocoframer)
library(rgl)

o <- get_mba_ontology()
o <- flatten_mba_ontology(o) 
o <- generate_taxons(o)

A <- get_ccf_grid_annotation()

sub <- filter_mba_ontology_children(o,"SS",include_parent = TRUE)
pos <- reshape2::melt(A) |>
  subset(value %in% sub$id)

a <- 200
clear3d()
plot_ccf_meshes(list(
  a=ccf_2017_mesh("root"),
  b=ccf_2017_mesh("SS")
),c("b"),bg_structure="a",fg_alpha = 0.4)
rgl::rgl.viewpoint(-40,phi = 50)
plot3d((pos$Var1-1)*a,(ncol(A)-pos$Var2)*a,(pos$Var3-1)*a,add = TRUE)




G <- get_exp_gene_relationships()

V_Sox10 <- get_aba_ish_data(get_gene_aba_ish_ids("Sox10"))
V_Lhx9 <- get_aba_ish_data(get_gene_aba_ish_ids("Lhx9"))
V_Pvalb <- get_aba_ish_data(get_gene_aba_ish_ids("Pvalb")[[1]])
V_Calb1 <- get_aba_ish_data(get_gene_aba_ish_ids("Calb1")[[1]])
MSK <- (V_Sox10>5) & (V_Lhx9>5) & (V_Pvalb>5) & (V_Calb1>5)


pos <- reshape2::melt(MSK) |>
  merge(reshape2::melt(A),by=1:3) |>
  subset(value.y>0)
pos$Var3_abs <- 29.5 + abs(pos$Var3 - 29.5)


a <- 200
clear3d()
plot_ccf_meshes(list(
    a=ccf_2017_mesh("root"),
    b=ccf_2017_mesh("CB"),
    c=ccf_2017_mesh("SS")
  ),c("b","c"),bg_structure="a",fg_alpha = 0.4)
plot3d((pos$Var1-1)*a,(ncol(A)-pos$Var2)*a,(pos$Var3_abs-1)*a,add = TRUE)


dim(abind::abind(V_Sox10,V_Lhx9,V_Pvalb,V_Calb1,along = 4))




