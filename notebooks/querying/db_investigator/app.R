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

reactlogReset()
reactlog_enable()

create_ui <- function(name, n, title){
    box(id = paste(name, "box", sep=""),
        width=12,
        textOutput(paste(name, n, sep = "")),
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
        textInput(
            inputId = "query",
            label = "Enter your query",
            value = "Mus musculus[Organism], brain, scRNA"
        ),
        actionButton(inputId = "go", label = "Submit job")
        #menuItem("Match images", tabName = "match")

        ),
    dashboardBody(
        #,
        uiOutput("body")
    )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {
    session$onSessionEnded(function() {
        stopApp()
    })

    sum_entrez <- reactive({
        uids_entrez <- entrez_search(db="gds", term=input$query)
        sum_entrez <- entrez_summary(db="gds", id=uids_entrez$ids)
        sum_entrez
    })

    observeEvent(input$go, {
        #browser()
        for (i in 1:20) {
            local({
                myi <- i
                output[[paste("value", i, sep = "")]] <-
                    renderText(unlist(sum_entrez()[[myi]]["summary"]))
            })
        }
    }, priority = 5)
    observeEvent(input$go, {
        output$body <-
            renderUI(Map(create_ui, n=1:20,
                         name="value",
                         title=Map("[[", sum_entrez()[1:20], "title")))
    }, priority = 10)
}

# Run the application
shinyApp(ui = ui, server = server)



























