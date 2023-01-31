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
        #textOutput(paste(name, n, sep = "")),
        textOutput(paste(name, "summary", n, sep = "")),
        uiOutput(paste0(name, "gds_link", n)),
        collapsible = T,
        collapsed =T,
        status = "success",
        #background = "lightblue",
        title = title)
}


ui <- dashboardPage(
    dashboardHeader(title = "Welcome to DB investigator"),
    dashboardSidebar(skin = "light",
        selectInput(
            inputId = "DB",
            choices = "GEO",
            label = "Select database"

        ),
        selectInput(
            inputId = "query",
            label = "Select your query",
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
        # fluidRow(
        #     tabBox(
        #         #title = "Search",
        #         id="tabset1",
        #         tabPanel(
        #             id="tab1",
        #             uiOutput("body")
        #         )
        #     )
        # )
        uiOutput("body")
    )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {

    session$onSessionEnded(function() {
        stopApp()
    })

    esearch <- eventReactive(input$go, {
        esearch <- investigator(db="gds", term1 = input$query, term2 = input$term2,
                                return_id = TRUE, do_search = FALSE)
    })
    #browser()
    sum_entrez <- eventReactive(input$go,{

        sum_entrez <- investigator(db="gds", term1 = input$query, term2 = input$term2,
                                   return_id = FALSE, do_search = TRUE, indx = as.numeric(my_list$index))
        sum_entrez
    })

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
            })
        }
    }, priority = 5)

    observeEvent(input$go, {
        output$body <-
            renderUI(Map(create_ui, n=1:length(esearch()),
                         name="value",
                         title=Map("[[", sum_entrez()[1:length(esearch())], "title")))
    }, priority = 10)

    my_list <- reactiveValues(index=1)

    observeEvent(input$go, {
        output$count <- renderUI(
            selectInput("indices", "batches",
                        c(1:length(esearch())),
                        selected = my_list$index)
        )
    }, priority = 4)

    observeEvent(input$indices, {
        req(input$indices)
        my_list$index <- input$indices
        #browser()
    })

}

# Run the application
shinyApp(ui = ui, server = server)



























