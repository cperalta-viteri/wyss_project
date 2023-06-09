library(rgl)
library(cocoframer)
library(IRanges)

# to plot a mesh manually
#shapelist3d(ccf_mesh)


# ==== Approach: biomarkers =====
# Given a list of biomarkers, find the regions where they are expressed

genes=c("Slc1a3","Calb1","Pax2","Rbfox3",
        "Lhx9","Sox10","Pvalb","Ret","Nrgn","Calb2","Eomes")

# get annotated grid
A <- get_ccf_grid_annotation()

#get root mesh, for data overlay
root_mesh <- ccf_2017_mesh("root")

# to process a list of ISH experiments for genes
# ids_l <- Map(get_gene_aba_ish_ids, genes)
# ids_l <- Map("[", ids_l[genes], 1)
# ids_l <- ids_l[!is.na(ids_l)]

#unlist(unname(Map("[", ids_l[genes], 1)))


# getting ISH experimental data for given genes
## we get, gene expression and xyz coordinates
V_Sox10 <- get_aba_ish_data(get_gene_aba_ish_ids("Sox10"))

V_Lhx9 <- get_aba_ish_data(get_gene_aba_ish_ids("Lhx9"))
V_Pvalb <- get_aba_ish_data(get_gene_aba_ish_ids("Pvalb")[[1]])
V_Calb1 <- get_aba_ish_data(get_gene_aba_ish_ids("Calb1")[[1]])

# apply different filtering according to profile
MSK <- (V_Sox10>5) & (V_Lhx9>5) & (V_Pvalb>5) & (V_Calb1>5)

#MSK <- (V_Sox10>5) & (V_Lhx9>15) & (V_Pvalb>5) & (V_Calb1>15)

pos <- reshape2::melt(MSK) |>
  subset(value)
# when merging, we get the structure ID
pos <- merge(pos,reshape2::melt(A),by=1:3)

# computing the scaling factor
sf <- range(root_mesh$it)[2]-range(root_mesh$it)[1]/67-1
sf <- sf*0.0038

# first plot mesh, as reference, and then the points for the filteres gene expression profile
plot_ccf_meshes(list(a=ccf_2017_mesh("root")),"a", fg_alpha = 0.2)
plot3d(pos$Var1*sf,pos$Var2*sf,pos$Var3*sf,add = TRUE)
#col = hcl.colors(pos$value.y)

# meshes <- Map(ccf_2017_mesh, acr) #c("root", "PFL", "CB")
# plot_ccf_meshes(meshes, fg_structure = c("PFL", "CB"), bg_structure = "root")




# ==== Can we visualize a single data point, referring to a single cell? ====
## With MSK <- sox10>=49, this structure, is missing from mesh zip file...

MSK <- (V_Sox10>=49) #10689 id
pos <- reshape2::melt(MSK) |>
  subset(value)
pos <- merge(pos,reshape2::melt(A),by=1:3)
#we retain the expression value
pos <- merge(pos,reshape2::melt(A),by=1:3)

# df <- pos[,1:3]
# df <- rbind(df, pos[,1:3])

## we can get Allen color code
o[o$id %in% pos$value.y,]$color_hex_triplet -> c

# we can get also the brain region for a gene expression profile
o[o$id %in% pos$value.y,]$acronym -> acr

plot_ccf_structure_points_v2(A, pos$value.y, scaleF = sf, col=sprintf("#%s", c))
## No, we can visualize structure, not cell



# some structures are missing from mesh zip file, so in the case we can get the origin of a single cell,
# but with ish experimental data, we can visualize the region!





# ===== Approach: structure =====
# Given a region of interest, get the ontology path to related substructures

o <- get_mba_ontology()
o <- flatten_mba_ontology(o)
o <- generate_taxons(o)

subset(o,name %in% "Cerebellum")
subset(o,name %in% "root")
cerebelum_substructures <- subset(o,any(CharacterList(strsplit(o$taxons,";")) %in% 512))$id


sub <- filter_mba_ontology_children(o,"CB",include_parent = TRUE)
mba_structure_id("CB")

plot_ccf_meshes(list(a=ccf_2017_mesh("root")),"a")

plot_ccf_meshes(list(a=ccf_2017_mesh("root"),b=ccf_2017_mesh("CB")),"b",bg_structure="a")




table(A %in% cerebelum_substructures)
table(cerebelum_substructures %in% A)


























# ==== Functions ====

plot_ccf_structure_points_v2 <- function(ccf_arr,
                                      structure_id,
                                      scaleF=NULL,
                                      ...) {

  vol_melt <- reshape2::melt(ccf_arr)
  vol_melt <- vol_melt[vol_melt$value > 0,]

  if(!structure_id %in% unique(vol_melt$value)) {
    stop("Structure not found in ccf_arr.")
  }

  str_melt <- vol_melt[vol_melt$value == structure_id,]

  if (is.null(scaleF)){
    rgl::plot3d(str_melt$Var1,str_melt$Var2,str_melt$Var3, ...)

  } else rgl::plot3d(str_melt$Var1*scaleF,str_melt$Var2*scaleF,str_melt$Var3*scaleF,
                     add=TRUE, ...)

}
