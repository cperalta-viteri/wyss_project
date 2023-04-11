
#metaD <- read.csv("geo_data/GSE213828/GSE213828_merged_sample_info.txt.gz", sep = "\t")
# row.names = "qp_well_fq"

# l <- list(donor_sex="female",
#           donor_sex_male=NULL,
#           organism="human",
#           donor_sex_col="NA",
#           barcode_col="barcode")



format_df <- function(fields, dataframe) {
  fields <- Filter(Negate(is.null), fields)
  # get fields containing columns names as mask
  msk <- grep("*_col$", names(fields), value = TRUE)
  #browser()
  # get the column name removing NA values
  columns <- grep("(?m)^(?!(NA))", unlist(unname(fields[msk])), value = TRUE, perl = TRUE, ignore.case = TRUE)
  # get fields containing values
  fields[msk]=NULL

  n <- unlist(names(Map("[", fields)))
  v <- unlist(unname(Map("[", fields)))
  #browser()
  output_df <- as.data.frame(dataframe[, columns])
  colnames(output_df) <- columns

  f2 <- function(n, v){
    if(!v %in% c("NA", "na")){
      output_df[n] <- v
      output_df
    }
  }

  output_df <-Reduce(data.frame, Map(f2, n, v))
  #browser()
  output_df <- output_df[, !colnames(output_df) %in% grep("\\.[[:digit:]]", colnames(output_df), value = TRUE)]

}



















