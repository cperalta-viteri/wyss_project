
library(shiny)
library(DT)
source("lib_std.R")



dashboardPage(
  dashboardHeader(),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    fluidRow(
      DTOutput(outputId = "input_df")
    ),
    fluidRow(
      box(width = 12,collapsible = FALSE,
        column(
          width=12,
          actionButton(paste0("test_btn"),"test",status = "danger"),
          p("List of standard fields to map"),
          radioGroupButtons("current_working_field_name",choices = names(std_fields),status = "info",individual = TRUE),
          box(title = "Mapping rules",width = 12,collapsible = FALSE,
            uiOutput("cwf_rules"),
            actionButton(paste0("rule_add"),"Add rule"),
            textInput(paste0("rule_default"),"default value")
          ),
          actionButton(paste0("next_field"),"Next",status = "primary")
        )
      )
    ),
    fluidRow(
      textOutput("mutate_string")
    ),    
    fluidRow(
      #DTOutput(outputId = "output_df")
    )
  )
)



