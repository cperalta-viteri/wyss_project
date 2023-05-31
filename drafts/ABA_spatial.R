
library(EBImage)
library(reticulate)
pd <- import("pandas")
L <- pd$read_pickle("~/Downloads/405_1.stv")

plot(as.raster(L$data / max(65535)))
