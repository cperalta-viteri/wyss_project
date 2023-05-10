
library(shiny)
library(bs4Dash)
library(reactlog)
library(shinyjs)



library(httr2)
library(jsonlite)
library(rentrez)
gds_quey <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="



investigator <- function(db="gds",
                         term1=NULL,
                         term2=NULL,
                         return_id=TRUE,
                         do_search=FALSE,
                         indx=NULL){
  environment(query_fun) <- environment()
  switch (term1,
          mouse = {
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
  if(do_search==TRUE && !is.null(indx)) return(smry <- entrez_summary(db=db, id=ids[[indx]]))
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


reactlogReset()
reactlog_enable()


create_ui <- function(name, n, title){
    box(id = paste(name, sep=""),
        width=12,
        textOutput(paste(name, "summary", n, sep = "")),
        uiOutput(paste0(name, "gds_link", n)),
        uiOutput(paste0(name, "check", n)),
        collapsible = T,
        collapsed =T,
        status = "success",
        title = title)
}


ui <- dashboardPage(
    dashboardHeader(title = "GEO download"),
    dashboardSidebar(disable = TRUE),
    dashboardBody(
      fluidRow(  
        box(width = 12,
          title = "Query",
          selectInput(
            inputId = "query_org",
            label = "Select your organism",
            choices = c("mouse", "human")
          ),
          textInput(
            inputId = "query_terms",
            label = "Query Text",
            value = NULL
          ),
          actionButton(inputId = "go", label = "Submit")
        )
      ),
      fluidRow(
            box(width = 12,
                title = "Query Results",
                uiOutput("body")
            )
        )
    )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {
    esearch <- eventReactive(input$go, {
        esearch <- investigator(
          db="gds", term1 = input$query, term2 = input$term2,
          return_id = TRUE, do_search = FALSE)
    })
    sum_entrez <- eventReactive(input$go,{
        if(is.list(esearch())){
            sum_entrez <- investigator(db="gds", term1 = input$query, term2 = input$term2,
                                       return_id = FALSE, do_search = TRUE, indx = as.numeric(my_list$index))
        } else {
            sum_entrez <- investigator(db="gds", term1 = input$query, term2 = input$term2,
                                       return_id = FALSE, do_search = TRUE)
        }
        sum_entrez
    })

    #==== Output summary and GEO link ====
    observeEvent(input$go, {
        for (i in 1:length(esearch())) {
            local({
                myi <- i
                output[[paste("value","summary", i, sep = "")]] <-
                    renderText(
                        sprintf("%s\n",unlist(sum_entrez()[[myi]]["summary"]))
                        )
                output[[paste0("value", "gds_link", i)]] <-
                    renderUI(
                        tags$a(href=
                            sprintf("%s%s",gds_quey, unlist(sum_entrez()[[myi]]["accession"])),
                            sprintf("GEO: %s", unlist(sum_entrez()[[myi]]["accession"])),
                            target="_blank"
                        )
                        )
            })
        }
    }, priority = 5)

    #==== Rendering boxes for hits ====
    observeEvent(input$go, {
        output$body <- renderUI(Map(create_ui, n=1:length(esearch()),
                         name="value",
                         title=Map("[[", sum_entrez()[1:length(esearch())], "title")))
    }, priority = 10)

    #==== Setting reactive list for later queries ====
    my_list <- reactiveValues(index=1)

    observeEvent(input$go, {
        req(is.list(esearch()))
        output$count <- renderUI(
            selectInput("indices", "batches",
                        c(1:length(esearch())),
                        selected = my_list$index)
            )
    }, priority = 4)

    observeEvent(input$indices, {
        req(input$indices)
        my_list$index <- input$indices

    })

    df <- reactive(read.csv("../acronym_df.csv"))

    observeEvent(input$in_stage, {
        switch(input$query,
               mouse = {if(input$in_stage == "developing"){
                   output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri",
                                                              label = "Brain region of origin",
                                                              choices = df()$mouse_developing))
               }else{output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri",
                                                                label = "Brain region of origin",
                                                                choices = df()$mouse_adult))}
                   },
               human = {if(input$in_stage == "developing"){
                   output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri", #in
                                                                 label = "Brain region of origin",
                                                                 choices = df()$human_developing))
               }else{output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri",
                                                                   label = "Brain region of origin",
                                                                   choices = df()$human_adult))}

               }
               )
        switch (input$query,
                mouse = {if(input$in_stage == "developing"){
                    output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                                  label = "Select brain region of injection",
                                                                  choices = df()$mouse_developing))
                }else{output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                                    label = "Select brain region of injection",
                                                                    choices = df()$mouse_adult))}
                },
                human= {if(input$in_stage == "developing"){
                    output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                                  label = "Select brain region of injection",
                                                                  choices = df()$human_developing))
                }else{output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                                    label = "Select brain region of injection",
                                                                    choices = df()$human_adult))}

                }
        )
    })

    observeEvent(input$download, {

        wd <- "../geo_data"
        if(!dir.exists(wd)) dir.create(wd)
        id <- input$in_uid
        if(!dir.exists(file.path(wd, id))) dir.create(file.path(wd, id))
        link <- sum_entrez()[[sub("GSE", "200", id)]][["ftplink"]]
        # input$query, input$in_strain, input$in_acq, input$in_stage, input$in_brainOri, input$in_brainIn
        l <- list(organism=input$query, strain=input$strain,
                  acq_mode=input$in_acq, stage=input$in_stage,
                  brain_orig=input$in_brainOri, brain_inject=input$in_brainIn)
        #write(jsonlite::toJSON(l, pretty = TRUE, auto_unbox = TRUE), sprintf("%s/%s/metadata_%s.json", wd, id, id))
        jsonlite::write_json(l, pretty = TRUE, auto_unbox = TRUE, path= sprintf("%s/%s/metadata_%s.json", wd, id, id))
        system(sprintf("wget -np -nd -r -P %s %s", file.path(wd, id), file.path(link, "suppl")))

        #system("wget -np -nd -r -P test_dir/ ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206940/")
        #sprintf("wget -np -nd -r -P %s %s", file.path(wd, id), file.path(link, "suppl"))

    })


}

# Run the application
shinyApp(ui = ui, server = server)



















