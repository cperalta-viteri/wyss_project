
library(IRanges)
library(igraph)



go_basic_obo_con <- function() {
  url("http://purl.obolibrary.org/obo/go/go-basic.obo")
}



# read OBO format and return an igraph with all "is_a" and "relationship"
read_obo <- function(obo_file) {
  #obo_file <- "http://purl.obolibrary.org/obo/go/go-basic.obo"
  
  pat <- "^([^:]+) *: *(.*)$"
  obo <- readLines(obo_file)
  obo <- obo[grepl("^\\[",obo) | grepl(pat,obo)]
  obo <- DataFrame(
    section = cumsum(grepl("^\\[",obo)),
    key = sub(pat,"\\1",obo),
    value = sub(pat,"\\2",obo)
  )
  obo <- subset(obo,section %in% section[key=="[Term]"])
  obo <- subset(obo,key %in% c("is_a","name","id","is_obsolete","namespace","relationship"))
  
  get_relationship <- function(obo) {
    R <- with(obo[obo$key=="relationship",],splitAsList(value,section))
    R <- sub(" *!.*$","",R) # remove comments
    pat <- "^([^ ]+) ((GO:[0-9]+)|(KW-[0-9]+))$"
    stopifnot(all(grepl(pat,unlist(R))))
    DataFrame(
      relationship_type = sub(pat,"\\1",R),
      relationship_goid = sub(pat,"\\2",R)
    )
  }
  
  obo$section <- factor(obo$section)
  obo <- DataFrame(
    goid = with(subset(obo,key=="id"),splitAsList(value,section)) |> drop(),
    term = with(subset(obo,key=="name"),splitAsList(value,section)) |> drop(),
    is_obsolete = (with(subset(obo,key=="is_obsolete"),splitAsList(value,section)) |> drop()) %in% "true",
    namespace = with(subset(obo,key=="namespace"),splitAsList(value,section)) |> drop(),
    is_a = with(subset(obo,key=="is_a"),splitAsList(value,section)),
    get_relationship(obo)
  )
  obo$is_a <- sub(" *!.*$","",obo$is_a) # remove comments
  
  is_a <- stack(setNames(obo$is_a,obo$goid))
  is_a$type <- "is_a"
  relationship <- stack(setNames(obo$relationship_goid,obo$goid))
  relationship$type <- unlist(obo$relationship_type)
  
  go <- graph_from_data_frame(rbind(relationship,is_a),vertices = as.data.frame(obo[1:4]))
  
  if (!is.dag(go)) warning("Loaded graph is not a DAG")
  return(go)
}



read_gmt <- function(con) {
  txt <- readLines(con)
  pat <- "^([^\t]*)\t([^\t]*)\t(.*)$"
  z <- DataFrame(
    STANDARD_NAME = sub(pat,"\\1",txt),
    EXTERNAL_DETAILS_URL = sub(pat,"\\2",txt),
    MEMBERS = CharacterList(strsplit(sub(pat,"\\3",txt),"\t",fixed=TRUE))
  )
}


read_gaf <- function(gaf_file) {
  txt <- readLines(gaf_file)
  stopifnot(grepl("^!gaf-version:",txt[1]))
  txt <- txt[!grepl("^!",txt)]
  gaf <- read.table(text = txt,sep="\t",quote="",comment.char = "")
  colnames(gaf) <- c("db","db_object_id","db_object_symbol","qualifier","go_id","db_ref","evidence","with_or_from","aspect","db_object_name","db_object_synonym","db_object_type","taxon","date","assigned_by","annotation_extension","gene_product_form_id")
  return(gaf)
}













