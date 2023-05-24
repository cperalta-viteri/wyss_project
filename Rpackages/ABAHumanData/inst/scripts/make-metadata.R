## validate with `AnnotationHub::readMetadataFromCsv("TENxBrainData")`
## (above pkg directory)

main.data <- data.frame(
    Title = paste0(
        rep(
            c(
                "Brain scRNA-seq data"
            ),
            each=4
        ),
        ", ",
        rep(
            c(
                "'CSV-based 10x Genomics' format",
                "'dense matrix' format",
                "sample (column) annotation",
                "gene (row) annotation"
            ),
            times=2
        )
    ),
    Description = paste(
        rep(
            c(
                "Single-cell RNA-seq data for 15,673 brain cells from human."
            ),
            each=4
        ),
        rep(
            c(
                "'csv-based 10x Genomics' format originally provided by Allen Brain Institute.",
                "Full rectangular, block-compressed format, 1.5GB block size.",
                "Inferred sample descriptions.",
                "Gene descriptions."
            ),
            times=2
        )
    ),
    RDataPath = file.path("ABAHumanData",
        paste0(
            rep(
                c(
                    "15k_neurons_filtered_gene_bc_matrices_csv"
                ),
                each=4
            ),
            rep(c(".h5", "_rectangular.h5", "_colData.rds", "_rowData.rds"), times=2)
        )
    ),
    BiocVersion="3.8",
    Genome="mm10",
    SourceType="CSV",
    SourceUrl=rep(
        c(
            "https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv"
        ),
        each=4
    ),
    SourceVersion=rep(c("f97b8c67e57a534816b9b2bccf31476b", "7696db0199a30a36f857e86a1ad1da8d"), each=4),
    Species="Homo Sapiens",
    TaxonomyId="9606",
    Coordinate_1_based=TRUE,
    DataProvider="Allen Brain Institute",
    Maintainer="Julien Prados <julien.prados@unige.ch>",
    RDataClass="character",
    DispatchClass=rep(c("H5File", "H5File", "Rds", "Rds"), times=2),
    stringsAsFactors = FALSE
)

write.csv(file="metadata.csv", main.data, row.names=FALSE)
