
female <- c(0, "female", "feminine", "f")
male <- c(1, "male", "masculine", "m")
transcriptomics <- c("transciptomics", "scrna", "sc-rna", "single cell rna", "single-cell rna")
epigenomics <- c("epigenomics","scatac", "sc-atac", "single cell atac", "single-cell atac")

metaD <- read.csv("geo_data/GSE213828/GSE213828_merged_sample_info.txt.gz",row.names = "qp_well_fq", sep = "\t")
