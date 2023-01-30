#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rentrez)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Welcome to DB investigator"),

    # Sidebar
    sidebarLayout(
        sidebarPanel(
            selectInput(
                inputId = "DB",
                choices = "GEO",
                label = "Select database"
            ),
            textInput(
                inputId = "query",
                label = "Enter your query",
                value = ""
            ),
            submitButton()
        )
        ,

        #
        mainPanel(
            renderUI("value")
            )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    output$value <- reactive({
        uids_entrez <- rentrez::entrez_search(db="gds", term=input$query)
        sum_entrez <- rentrez::entrez_summary(db="gds", id=uids_entrez$ids)
        sum_entrez
    })
}

# Run the application
shinyApp(ui = ui, server = server)



























