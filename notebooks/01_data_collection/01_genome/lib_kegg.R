

# Convert a json object retrieved from BRITE into a DataFrame
brite_nodes <- function(n,path="0") {
  e <- data.frame(path=path,label=n$name)
  children_paths <- file.path(path,seq_along(n$children))
  children_nodes <- mapply(brite_nodes,n$children,children_paths,SIMPLIFY = FALSE) |> 
    do.call(what = rbind)
  rbind(e,children_nodes)
}

# Get BRITE ontology
brite_ontology_from_json <- function(...,root="0") {
  V <- rjson::fromJSON(...) |>
    brite_nodes(root)
  H <- data.frame(parent=dirname(V$path),child=V$path) |>
    subset(child != root) |>
    graph_from_data_frame(vertices = V)
  V(H)$ontology_id <- V(H)$name
  V(H)$depth <- distances(H,V(H)[!.to(E(H))],mode = "out")
  V(H)$pathway_id <- sub("^([^ ]+).*","\\1",V(H)$label)
  V(H)$pathway_id[!grepl("^([^ ]+).*",V(H)$label)] <- NA
  H
}

# Read gene-pathway annotations (remove organism prefix)
read_kegg_link <- function(file) {
  Ann <- read.table(file,sep="\t",col.names = c("locus_tag","pathway_id"),header=FALSE)
  Ann$locus_tag <- sub("^[^:]*:","",Ann$locus_tag)
  Ann$locus_tag <- sub("\\.[0-9]+$","",Ann$locus_tag)
  Ann$pathway_id <- sub("^(path|md):[a-z]*_?","",Ann$pathway_id)
  Ann
}


