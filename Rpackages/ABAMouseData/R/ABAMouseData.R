
#' @importFrom HDF5Array HDF5Array
#' @importFrom SingleCellExperiment SingleCellExperiment
.create_SCE_from_base <- function(base) {
    hub <- ExperimentHub()

    ## row and column data
    rdatapath <- paste0("TENxBrainData/", base, "_rowData.rds")
    rowData <- query(hub, rdatapath)[[1]]

    rdatapath <- paste0("TENxBrainData/", base, "_colData.rds")
    colData <- query(hub, rdatapath)[[1]]

    ## HDF5, from ExperimentHub:
    rdatapath <- paste0("TENxBrainData/", base, "_rectangular.h5")
    h5file <- query(hub, rdatapath)[[1]]

    h5array <- HDF5Array(h5file, "counts")

    SingleCellExperiment(
        list(counts = h5array), rowData = rowData, colData = colData
    )
}

#' @export
ABAMouseTENxData <- function() {
    #.create_SCE_from_base("")
}

#' @export
ABAMouseSSv4Data <- function() {
    #.create_SCE_from_base("")
}
