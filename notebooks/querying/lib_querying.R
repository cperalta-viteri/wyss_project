

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
esearch <- request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi")
# setting parameters
mouse_esearch <- esearch %>% req_url_query(
  "db"="gds", # set db to GEO (gds)
  "retmode"="json", # retrieve json format
  "term"="Mus musculus[Organism], brain, scRNA") # term Mus musculus to filter
# perform request and get the body
uids <- req_perform(mouse_esearch) |> resp_body_json()
# retrieving the ids as a character list
uids <- unlist(uids$esearchresult$idlist)

uids_entrez <- entrez_search(db="gds", term="Mus musculus[Organism], brain, scRNA")

# ==== ESummary ====

# base url: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi

## Required parameters

# db
# id= list of ids from the DB

esummary <- request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi")

mouse_esummary <- esummary %>% req_url_query(
  "db"="gds",
  "retmode"="json",
  #"term"="Mus musculus[Organism]",
  "id"=paste(uids, collapse = ","))

summaries <- mouse_esummary %>% req_perform() %>% resp_body_json(flatten=TRUE)
# to access info:
summaries$result["200207249"]$`200207249`$summary


sum_entrez <- entrez_summary(db="gds", id=uids_entrez$ids)
sum_entrez$`200207249`$summary
# ==== GEO query ====
# base URL, adding GSE with the number (removing 200) brings ot the webpage

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206940

















































