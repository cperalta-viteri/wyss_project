


# Allen Brain Atlas Mouse Cortex
.download_ABAHumanTENx <- function(dir) {
    # Download the files from GEO as well as the cell-type hierarchy
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/dend.json',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/tsne.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/medians.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv',dir)
    .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/64/6d/646d3592-aff6-4364-8c3f-9e64b902638a/human_dendrogram.rds',dir)
    .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0c/0c/0c0c882d-1c31-40a9-8039-3bf2706a77cd/sample-exp_component_mapping_human_10x_apr2020.zip',dir)
    .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/70/32/70326830-e306-4743-a02c-a8da5bf9eb56/readme-m1-10.txt',dir)
    .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/38/05/3805d61e-5c4d-446d-9331-cb61eb181198/taxonomy.txt',dir)
}

.download_ABAHumanSSv4 <- function(dir) {
    # Download the files from GEO as well as the cell-type hierarchy
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/intron.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/exon.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/trimmed_means.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/medians.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/dend.json',dir)
    .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/tsne.csv',dir)
    .download_url('https://transcriptomic-viewer-downloads.s3-us-west-2.amazonaws.com/human/taxonomy.zip',dir)
    .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/30/80/30803016-9515-45fd-a6af-8595902c176f/readme_human_ctx.txt',dir)
    .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b1/9f/b19f5e10-c6b2-4cf4-80fc-182e69a4511f/sample-exp_component_mapping_human_smart-seq_oct_2019.zip',dir)
}





#' @export
prepare_ABAHumanTENx <- function(outdir=".") {
    tmpdir <- tempdir()
    .download_ABAHumanTENx(tmpdir)
    e <- .add_qc_metrics(.load_ABAHumanTENx(tmpdir))
    .save_SCE(e,outdir,"ABAHumanTENx")
}

#' @export
ABAHumanTENx <- function(dir) {
}


#' @export
prepare_ABAHumanSSv4 <- function(outdir=".") {
    tmpdir <- tempdir()
    .download_ABAHumanSSv4(tmpdir)
    e <- .add_qc_metrics(.load_ABAHumanSSv4(tmpdir))
    .save_SCE(e,outdir,"ABAHumanSSv4")
}

#' @export
ABAHumanSSv4 <- function(dir) {
}






