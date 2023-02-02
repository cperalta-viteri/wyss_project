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
                    textInput(inputId = "in_organism", label = "Type the organism"),
                    textInput(inputId = "in_strain", label = "Type the strain"),
                    textInput(inputId = "in_age", label = "Type the age range"),
                    textInput(inputId = "in_organ", label = "Type the organ", value = "brain"),
                    selectInput(inputId = "in_acq", label = "Select the data acquisition mode", choices = c("scRNA", "nucleai", "scATAC")),
                    selectInput(inputId = "in_brainIn", label = "Select brain region injection", choices = c(NULL)),
                    selectInput(inputId = "in_brainOut", label = "Select brain region collection", choices = c(NULL)),
                    textInput(inputId = "in_uid", label = "Type the accession nº"),
                    textInput(inputId = "in_cmatrix", label = "Type the matrix format")
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
        # esearch <- entrez_search(db="gds",
        #                          term=sprintf("Mus musculus[Organism], %s", input$term2),
        #                          retmax=10000)
    })
    #browser()
    sum_entrez <- eventReactive(input$go,{
        if(length(esearch()) > 100){

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
        #browser()
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
        req(length(esearch()) > 100)
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
    #check_list <- reactive({grep("checker*", names(input), value = TRUE)})

    #check_values <- reactiveValues(item=reactive({grep("checker*", names(input), value = TRUE)}))

}

# Run the application
shinyApp(ui = ui, server = server)



























