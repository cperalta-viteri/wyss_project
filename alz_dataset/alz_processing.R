

library(SingleCellExperiment)



e <- readRDS("data/sce.EC.scAlign.assigned.rds")
s <- readRDS("data/sce.SFG.scAlign.assigned.rds")

# mse loss | marginrankingloss



sce_merged <- cbind(e,s)


sce_merged




