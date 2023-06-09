


#' Read structure ontology
#'
#' Consult following page for a description
#' http://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies
#' @importFrom rjson fromJSON
#' @export
#' @param json_file path to JSON file
#' @examples
#' human_ontology <- aba_read_structure_ontology(
#' "http://api.brain-map.org/api/v2/structure_graph_download/10.json")
#' mouse_ontology <- aba_read_structure_ontology(
#' "http://api.brain-map.org/api/v2/structure_graph_download/1.json")
#'
#' \donttest{
#' # Check ontology sharing
#' gplots::venn(list(human_ontology$id,mouse_ontology$id))
#' gplots::venn(list(na.exclude(human_ontology$atlas_id),na.exclude(mouse_ontology$atlas_id)))
#' gplots::venn(list(na.exclude(human_ontology$acronym),na.exclude(mouse_ontology$acronym)))
#' }
aba_read_structure_ontology <- function(json_file) {
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




