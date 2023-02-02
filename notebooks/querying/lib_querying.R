

library(httr2)
library(jsonlite)
library(rentrez)


# ===== ESearch =====

# base URL: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/search.fcgi

## Requieres:

# db=gds

## Parameters

# term =
# retmax, default retrieve 20 elements
# retstart, default=0, restart search

# use base URL with DB as GEO (gds) and adding retmode as json"
#esearch <- request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi")
# setting parameters
 # mouse_esearch <- esearch %>% req_url_query(
 #   "db"="gds", # set db to GEO (gds)
 #   "retmode"="json", # retrieve json format
 #   "term"="Mus musculus[Organism],brain,scRNA",
 #   "retmax" = 10000) # term Mus musculus to filter
# perform request and get the body
#uids <- req_perform(mouse_esearch) |> resp_body_json()
#  retrieving the ids as a character list
# uids <- unlist(uids$esearchresult$idlist)

 uids_entrez <- entrez_search(db="gds",
                              term="Mus musculus[Organism], brain, scRNA",
                              retmax=10000, use_history = TRUE)

# ==== ESummary ====

# base url: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi

## Required parameters

# db
# id= list of ids from the DB
#
# esummary <- request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi")
# #
#  mouse_esummary <- esummary %>% req_url_query(
#    "db"="gds",
#    "retmode"="json",
#    "term"="Mus musculus[Organism]",
#    "id"=paste(uids, collapse = ","))
# #
#  summaries <- mouse_esummary %>% req_perform() %>% resp_body_json(flatten=TRUE)
# to access info:
#summaries$result["200207249"]$`200207249`$summary

#id=uids_entrez$ids,
#sum_entrez <- entrez_summary(db="gds",  web_history = uids_entrez$web_history)
# sum_entrez$`200207249`$summary


# ==== GEO query ====
# base URL, adding GSE with the number (removing 200) brings ot the webpage

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206940

gds_quey <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="




# ==== Functions ====

investigator <- function(db="gds",
                         term1=NULL,
                         term2=NULL,
                         return_id=TRUE,
                         do_search=FALSE,
                         indx=NULL){
  environment(query_fun) <- environment()
  switch (term1,
    mouse = {
      #browser()
      organism <- "Mus musculus[Organism]"
      if (isFALSE(do_search)) return(output <- query_fun())
      if (isTRUE(do_search)) return(output <- query_fun(indx = indx))
    },
    human = {
      organism <- "Homo sapiens[Organism]"
      if (isFALSE(do_search)) return(output <- query_fun())
      if (isTRUE(do_search)) return(output <- query_fun(indx = indx))
    }
  )
}


query_fun <- function(indx=NULL){
  if(is.null(term2)){query <- sprintf(organism)} else {query <- sprintf("%s, %s", organism, term2)}
  ids <- entrez_search(db=db, term=query, retmax=10000)
  if (length(ids$ids)>100) ids <- slicer(ids$ids, jump = 100) else ids <- ids$ids
  #browser()
  if (return_id==TRUE) return(ids)
  #browser()
  if(do_search==TRUE  & is.list(ids)) return(smry <- entrez_summary(db=db, id=ids[[indx]]))
  if(do_search==TRUE) return(smry <- entrez_summary(db=db, id=ids))
}


slicer <- function(eid, jump=100){
  stop_pos <- length(eid)
  start_pos <- seq(1, stop_pos, jump)
  fun <- function(x=NULL){
    y <- if(x+jump<=stop_pos) x+(jump-1) else stop_pos
    s <- c(x:y)
    ids <- eid[s]
  }
  Map(fun, start_pos)
}



#
# library(readxl)
# astrocytes_df <- read_excel("Astrocytes_Transcriptomic_repository.xlsx")
# astrocytes_ids <- (gsub("^(?!GSE)", NA ,astrocytes_df$`Data location`, perl = TRUE))
# astrocytes_ids <- astrocytes_ids[-which(is.na(astrocytes_ids))]
#
# accession_list <- unlist(Map("[[", uids[1:length(uids)], "accession"))






























