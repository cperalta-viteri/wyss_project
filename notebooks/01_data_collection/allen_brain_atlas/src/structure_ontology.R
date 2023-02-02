



read_structure_ontology <- function(json_file) {
  
  node2df_ <- function(n) {
    l <- lapply(n$children,node2df_)
    l <- do.call(rbind,l)
    n$children <- NULL
    n[sapply(n,is.null)] <- NA
    rbind(data.frame(n),l)
  }
  
  json <- rjson::fromJSON(file=json_file,simplify = FALSE)
  node2df_(json$msg[[1]])
}


human_ontology <- read_structure_ontology("http://api.brain-map.org/api/v2/structure_graph_download/10.json")
mouse_ontology <- read_structure_ontology("http://api.brain-map.org/api/v2/structure_graph_download/1.json")

# Check ontology sharing
human_ontology$id %in% mouse_ontology$id
na.exclude(human_ontology$atlas_id) %in% mouse_ontology$atlas_id
na.exclude(human_ontology$acronym) %in% mouse_ontology$acronym

