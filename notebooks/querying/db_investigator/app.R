#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(bs4Dash)
library(rentrez)
library(reactlog)
library(shinyjs)
source("../lib_querying.R")

reactlogReset()
reactlog_enable()


create_ui <- function(name, n, title){
    box(id = paste(name, sep=""),
        width=12,
        #textOutput(paste(name, "summary", n, sep = "")),
        textOutput(paste(name, "summary", n, sep = "")),
        uiOutput(paste0(name, "gds_link", n)),
        uiOutput(paste0(name, "check", n)),
        collapsible = T,
        collapsed =T,
        status = "success",
        #background = "lightblue",
        title = title)
}


ui <- dashboardPage(
    dashboardHeader(title = "Welcome to GEO investigator"),
    dashboardSidebar(skin = "light",
        # selectInput(
        #     inputId = "DB",
        #     choices = "GEO",
        #     label = "Select database"
        #
        # )
        #dashboardHeader(title = "GEO Query", border = FALSE),
        selectInput(
            inputId = "query",
            label = "Select your organism",
            choices = c("mouse", "human")
        ),
        textInput(
            inputId = "term2",
            label = "Opt terms",
            value = NULL
        ),
        actionButton(inputId = "go", label = "Submit job"),
        uiOutput("count")

        ),
    dashboardBody(
        fluidRow(
            tabBox(width = 8,
                #title = "Search",
                id="tabset1",
                tabPanel(
                    id="tab1",
                    title = "Search",
                    uiOutput("body")
                )),
            tabBox(width = 4,
                tabPanel(
                    id="tab2",
                    title = "Metadata",
                    textOutput(outputId = "disp_organism"),
                    tags$head(tags$style("#disp_organism{color:blue; font-size: 20px}")),
                    textInput(inputId = "in_strain", label = "Type the strain"),
                    selectInput(inputId = "in_acq", label = "Select the data acquisition mode", choices = c("scRNA", "nuclei", "scATAC")),
                    selectInput(inputId = "in_stage", label = "Select Organism Stage", choices = c("developing", "adult"), selected = NULL),
                    uiOutput(outputId = "out_brainOri"),
                    uiOutput(outputId = "out_brainIn"),
                    textInput(inputId = "in_uid", label = "Type the accession number"),
                    actionButton(inputId = "download", label = "Click to Download")
                )
            )
        )
        #uiOutput("body")
    )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {

    esearch <- eventReactive(input$go, {
        esearch <- investigator(db="gds", term1 = input$query, term2 = input$term2,
                                return_id = TRUE, do_search = FALSE)
    })
    #browser()
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
        output$disp_organism <- renderText(stringr::str_to_title(input$query))
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
                # output[[paste0("value", "check", i)]] <-
                #     renderUI(
                #         checkboxInput(
                #             inputId = sprintf("checker_%s", unlist(sum_entrez()[[myi]]["accession"])),
                #             label="check to download",
                #             value = FALSE)
                #     )
            })
        }
    }, priority = 5)

    #==== Rendering boxes for hits ====
    observeEvent(input$go, {
        output$body <-
            renderUI(Map(create_ui, n=1:length(esearch()),
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

#options(shiny.port=8888)
























