
# metaD <- read.csv("geo_data/GSE213828/GSE213828_merged_sample_info.txt.gz", sep = "\t")
# row.names = "qp_well_fq"

# df <- data.frame(list("age"=sample(c(1,2), 5, replace = TRUE),
# "sex"=sample(c(0,1),5,replace = TRUE),
# "gender"=sample(c("m", "f"), 5, replace = TRUE)))
#
# l <- list(sex_col="gender",
#           sex_male_value="m",
#           sex_female_value="f",
#           organism="human")

# rename column
#colnames(df)[which(colnames(df)==l$sex_col)] <- "donor_sex"

# map values
# plyr::mapvalues(df$donor_sex, c(l$sex_male_value, l$sex_female_value), c("male", "female"))

# metaD <- cbind(metaD, list("gender"=sample(c("m", "f"), nrow(metaD), replace = TRUE)))

format_df <- function(fields, dataframe) {
  #browser()
  # filter out null and empty fields
  fields[fields==""] <- NA
  none <- fields[which(is.na(fields))]
  fields <- Filter(Negate(is.null), fields)
  fields <- fields[!fields == ""]
  # get fields containing columns names as mask
  msk <- grep("*_col$", names(fields), value = TRUE)

  if(!is.null(fields$sex_male_value) & !is.null(fields$sex_female_value)){
    donor_sex <- plyr::mapvalues(dataframe[[fields$sex_col]],
                               c(stringr::str_squish(fields$sex_male_value), stringr::str_squish(fields$sex_female_value)), c("male", "female"))

  }
  # browser()
  # get the column name removing NA values
  columns <- grep("(?m)^(?!(NA))", unlist(unname(fields[msk[msk!="sex_col"]])), value = TRUE, perl = TRUE, ignore.case = TRUE)
  # get fields containing values
  fields[msk]=NULL; fields$sex_male_value=NULL; fields$sex_female_value=NULL
  #
  n <- unlist(names(Map("[", fields)))
  v <- unlist(unname(Map("[", fields)))
  # #browser()
  output_df <- as.data.frame(dataframe[, columns])
  colnames(output_df) <- columns

  f2 <- function(n, v){
    output_df[n] <- v
    output_df
  }

  output_df <-Reduce(data.frame, Map(f2, n, v))
  output_df <- output_df[, !colnames(output_df) %in% grep("\\.[[:digit:]]", colnames(output_df), value = TRUE)]

  if(exists("donor_sex")){
    output_df <- cbind(output_df, donor_sex)
  } else {
    donor_sex <- NA
    output_df <- cbind(output_df, donor_sex)
  }

}



















