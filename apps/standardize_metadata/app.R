
library(tidyverse)
library(shiny)
library(shinyWidgets)
library(bs4Dash)
library(DT)



case_when_item <- function(type=c("in","text"),incol=NULL,inval=NULL,targetval=NULL,text=NULL) {
  type <- match.arg(type)
  list(
    type = type,
    incol = incol,
    inval = inval,
    targetval = targetval,
    text = text
  )
}

case_when_rule <- function(target_col,case_when_items=list()) {
  list(
    target_col = target_col,
    case_when_items = case_when_items
  )
}

standard_fields <- list(
  sex = list(suggestions=c("male","female"),validation="male|female",rules=list(case_when_rule("colA"),case_when_rule("colB"))),
  age = list(suggestions=c("E18","P0","P1","P2"),validation="[PE][0-9]+(\\.5)?"),
  strain = list(suggestions=c("C56BL6","FVB"))
)



caseWhenItemInput <- function(id,incol_choices=NULL,inval_choices=NULL,targetval_choices=NULL) {
  
  tags$table(
    tags$tr(width = "100%",
            tags$td(width = "60%", div(style = "font-size:10px;", "This is label1. I want all labels on left")),
            tags$td(width = "40%", textInput(inputId = "a", label = NULL))),
    tags$tr(width = "100%",
            tags$td(width = "60%", tags$div(style = "font-size:10pX;", "label2")),
            tags$td(width = "40%", textInput(inputId = "b", label = NULL)))
  )
  
  tags$table(
    tags$tr(
      width = "100%",
      tags$td(pickerInput(paste0(id,"_type"),choices=c("in","text"),inline=TRUE)),
      tags$td(
        conditionalPanel(
          condition = sprintf("input.%s_type == 'in'",id),
          span(
            "when",
            pickerInput(paste0(id,"_incol"),choices=incol_choices,inline = TRUE),
            "in",
            pickerInput(paste0(id,"_inval"),choices = inval_choices,inline = TRUE,multiple = TRUE, options = pickerOptions(style="label label-info")),
            "then",
            pickerInput(paste0(id,"_targetval"),choices=targetval_choices,inline = TRUE)
          )  
        ),
        conditionalPanel(
          condition = sprintf("input.%s_type == 'text'",id),
          textInput(paste0(id,"_incol"),label="",value = "col1 == col2 ~ col3"),
        )
      )
    )
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
      column(width=12,uiOutput("fields"))
    ),
    fluidRow(
      textOutput("mutate_string")
    ),    
    fluidRow(
      DTOutput(outputId = "output_df")
    )
  )
)

server <- function(input, output, session) {
  input_df <- reactive({iris})

  output_df <- reactive({
    input_df() |>
      filter(Sepal.Length>5)
  })
  
  output$mutate_string <- reactive({
    lapply(names(standard_fields),function(fn) {
        r <- lapply(seq_along(standard_fields[[fn]]$rules),function(i) {
          incol <- input[[paste0("rule_",fn,"_",i,"_incol")]]
          inval <- input[[paste0("rule_",fn,"_",i,"_inval")]]
          targetval <- input[[paste0("rule_",fn,"_",i,"_targetval")]]
          paste(incol," %in% c(",paste0('"',inval,'"',collapse=","),") ~ ",targetval)
        })
        r[length(r) + 1L] <- sprintf('TRUE ~ NA')
        sprintf("mutate(%s=case_when(%s))",fn,paste(r,collapse = ","))
    })
  })
  
  output$fields <- renderUI({
    lapply(names(standard_fields),function(fn) {
      # Create UI
      rules <- lapply(seq_along(standard_fields[[fn]]$rules),function(i) {
        caseWhenItemInput(paste0("rule_",fn,"_",i),incol_choices = colnames(input_df()))
      })

      # Add observers
      lapply(seq_along(standard_fields[[fn]]$rules),function(i) {
        incol <- paste0("rule_",fn,"_",i,"_incol")
        observeEvent(input[[incol]],{
          updatePickerInput(session, paste0("rule_",fn,"_",i,"_inval"),choices = as.character(unique(input_df()[[ input[[incol]] ]])))
        })
        updatePickerInput(session, paste0("rule_",fn,"_",i,"_targetval"),choices = as.character(standard_fields[[fn]]$suggestions))
      })
      
      # observeEvent(input[[paste0("rule_",fn,"_add")]],{
      #   standard_fields[[fn]]$rules[[length(standard_fields[[fn]]$rules) + 1L]] <- case_when_item("in")
      # })

      # Make the box
      box(
        title = fn,
        collapsible = TRUE,
        width = 12,
        rules,
        actionButton(paste0("rule_",fn,"_add"),"+"),
        textInput(paste0("rule_",fn,"_default"),"default value")
      )
    })
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










