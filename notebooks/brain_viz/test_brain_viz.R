library(rgl)
library(cocoframer)
#library(httr2)
#library(xml2)

genes=c("Slc1a3","Calb1","Pax2","Rbfox3",
        "Lhx9","Sox10","Pvalb","Ret","Nrgn","Calb2","Eomes")

# allen <- request("http://api.brain-map.org/api/v2/data/query.json?criteria=model::SectionDataSet,rma::criteria,[failed$eq'false'],products[abbreviation$eq'Mouse'],plane_of_section[name$eq'coronal'],genes[acronym$eq'Adora2a']")
# out <- req_perform(allen) %>% resp_body_json(simplifyVector = TRUE)

section_id <- get_gene_aba_ish_ids("Sox10") #Adora2a
ish_3d_array <- get_aba_ish_data(section_id)
ccf_mesh <- ccf_2017_mesh(acronym ="root" ,structure_id = "CB")
#ccf_mesh <- Map(ccf_2017_mesh, ish_str$acronym, section_id)
ish_str <- get_aba_ish_structure_data(section_id)

taxon <- generate_taxons(ontology)
reshape2::melt(ish_3d_array) -> xyz

#ccf_ann <- get_ccf_annotation()
grid <- get_ccf_grid_annotation()

plot_ccf_structure_points(ccf_arr = grid, structure_id = )

# to get xyz coordinates:
xyz <- xyz[xyz$value>0, ]
points3d(xyz$Var1,xyz$Var2, xyz$Var3)

plot_ccf_meshes(ccf_mesh, fg_structure = "CB")

# to plot a mesh manually
#shapelist3d(ccf_mesh)

#================================ Julien

o <- get_mba_ontology()
o <- flatten_mba_ontology(o)
o <- generate_taxons(o)

subset(o,name %in% "Cerebellum")
subset(o,name %in% "root")
cerebelum_substructures <- subset(o,any(CharacterList(strsplit(o$taxons,";")) %in% 512))$id


sub <- filter_mba_ontology_children(o,"CB",include_parent = TRUE)
mba_structure_id("CB")



plot_ccf_meshes(list(a=ccf_2017_mesh("CB")),"a")

plot_ccf_meshes(list(a=ccf_2017_mesh("root"),b=ccf_2017_mesh("CB")),"b",bg_structure="a")
plot3d(7829.740,3703.973,5694.498,add=TRUE)

o <- ccf_2017_mesh("root")

A <- get_ccf_grid_annotation()



table(A %in% cerebelum_substructures)
table(cerebelum_substructures %in% A)





V_Sox10 <- get_aba_ish_data(get_gene_aba_ish_ids("Sox10"))
V_Lhx9 <- get_aba_ish_data(get_gene_aba_ish_ids("Lhx9"))
V_Pvalb <- get_aba_ish_data(get_gene_aba_ish_ids("Pvalb")[[1]])
V_Calb1 <- get_aba_ish_data(get_gene_aba_ish_ids("Calb1")[[1]])


MSK <- (V_Sox10>10) & (V_Lhx9>10) & (V_Pvalb>10) & (V_Calb1>10)


pos <- reshape2::melt(MSK) |>
  subset(value)
pos <- merge(pos,reshape2::melt(A),by=1:3)



plot3d(pos$Var1,pos$Var2,pos$Var3,add = TRUE)


xyz <- reshape2::melt(V) |>
  subset(!is.na(value) & value>0)
plot3d(xyz$Var1,xyz$Var2,xyz$Var3,add = TRUE)
dim(A)
