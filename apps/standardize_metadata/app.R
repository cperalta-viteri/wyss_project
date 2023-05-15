
library(tidyverse)
library(shiny)
library(shinyWidgets)
library(bs4Dash)
library(DT)

inPredicateInput <- function(id) {
    span(
      "when",
      pickerInput(paste0(id,"_incol"),choices=paste0("col",1:10),inline = TRUE),
      "in",
      pickerInput(paste0(id,"_inval"),choices = c("lev1","lev2"),inline = TRUE,multiple = TRUE, options = pickerOptions(style="label label-info")),
      "then",
      pickerInput(paste0(id,"_targetval"),choices=c("male","female"),inline = TRUE)
    )
}

formulaPredicateInput <- function() {
  span(
    textInput("formula",label = "formula",value = "col1 == col2 ~ col3"),
  )
}

case_when_box <- function(target_name) {
  div(
    actionButton("btn","+ add mapping"),br(),
    inPredicateInput("id1"),br(),
    inPredicateInput("id2"),br(),
    formulaPredicateInput(),
    textInput("default","default value",value="NA"),
    footer =  "generated string"
  )
}


ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    fluidRow(
      DTOutput(outputId = "input_df")
    ),
    fluidRow(
      accordion(id="acc",
        accordionItem(title="sex",case_when_box("sex"),collapsed = FALSE),
        accordionItem(title="age",case_when_box("age"))
      )
    ),
    fluidRow(
      DTOutput(outputId = "output_df")
    )
  )
)

server <- function(input, output) {
  input_df <- reactive({iris})
  output_df <- reactive({
    input_df() |>
      filter(Sepal.Length>5)
  })
  output$input_df <- renderDT(input_df())
  output$output_df <- renderDT(output_df())
}

shinyApp(ui = ui, server = server)






# iris |> 
#   mutate(
#     sex = case_when(
#       Species %in% c("setosa","virginica") ~ "male",
#       Species %in% c("versicolor") ~ "female",
#       TRUE ~ NA
#   )) 
# 
# iris |>
# mutate(
#   sex = case_when(
#     Species %in% c("setosa","virginica") ~ Species,
#     Sepal.Length > 5 ~ "female",
#     Species %in% c("versicolor") ~ "female"
#     TRUE ~ NA
#   ))










