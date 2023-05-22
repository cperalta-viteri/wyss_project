
library(DT)
source("lib_std.R")


fields <- list(
  project_id = field_item(targetval_suggestions=NULL,validation=NULL,rules=NULL),
  batch_id = field_item(targetval_suggestions=NULL,validation=NULL,rules=NULL),
  cell_name = field_item(targetval_suggestions=NULL,validation=NULL,rules=NULL),
  donor_organism = field_item(targetval_suggestions=c("mouse","human"),validation=NULL,rules=NULL),
  donor_sex = field_item(targetval_suggestions=c("M","F"),validation="^[MF]$",rules=list(case_when_item(),case_when_item())),
  donor_age = field_item(targetval_suggestions=c("E18","P0","P1","P2"),validation="[PE][0-9]+(\\.5)?"),
  donor_mouse_strain = field_item(targetval_suggestions=c("C56BL6","FVB"),rules=list(case_when_item())),
  donor_mouse_line = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  donor_mouse_genotype = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  
  donor_human_genetic_population = field_item(targetval_suggestions=c("C56BL6","FVB"),rules=list(case_when_item())),
  donor_human_phenotype = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  donor_mouse_genotype = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  
  sample_protocol_dissociation = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  sample_protocol_chemistry = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  sample_tissue_type = field_item(targetval_suggestions=c("Postnatal","Embryonic","Organoid","Neuronal Culture"),rules=list(case_when_item())),
  sample_tissue_origin = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  sample_nervous_region = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  
  default_condition = field_item(targetval_suggestions=NULL,rules=list(case_when_item())),
  disease_model = field_item(targetval_suggestions=NULL,rules=list(case_when_item()))
)



ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    fluidRow(
      DTOutput(outputId = "input_df")
    ),
    fluidRow(
      column(
        width=12,
        radioGroupButtons("current_working_field_name",choices = names(fields),status = "info",individual = TRUE),
        uiOutput("fields"),
        actionButton(paste0("rule_add"),"Add rule"),
        textInput(paste0("rule_default"),"default value"),
        actionButton(paste0("next_field"),"Next",status = "primary")
      )
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
  
  # output$mutate_string <- reactive({
  #   lapply(names(standard_fields),function(fn) {
  #   })
  # })

  
  cwf_name <- reactiveVal(1L)
  cwf <- reactiveVal(list())
  
    
  observeEvent(input$current_working_field_name,{
    cwf(fields[[input$current_working_field_name]])
  })
  observeEvent(input$rule_add,{
    v <- isolate(cwf())
    v$rules[[length(v$rules) + 1L]] <- case_when_item()
    cwf(v)
  })

  output$fields <- renderUI({
      rules <- lapply(seq_along(cwf()$rules),function(i) {
        caseWhenItemInput(paste0("rule_",i),incol_choices = colnames(input_df()))
      })

      # # Add observers
      # lapply(seq_along(standard_fields[[fn]]$rules),function(i) {
      #   incol <- paste0("rule_",fn,"_",i,"_incol")
      #   observeEvent(input[[incol]],{
      #     updatePickerInput(session, paste0("rule_",fn,"_",i,"_inval"),choices = as.character(unique(input_df()[[ input[[incol]] ]])))
      #   })
      #   updatePickerInput(session, paste0("rule_",fn,"_",i,"_targetval"),choices = as.character(standard_fields[[fn]]$suggestions))
      # })
      #
      # # observeEvent(input[[paste0("rule_",fn,"_add")]],{
      # #   standard_fields[[fn]]$rules[[length(standard_fields[[fn]]$rules) + 1L]] <- case_when_item("in")
      # # })
      #
      # Make the box
      div(
        rules
      )

  })
  
  output$input_df <- renderDT(input_df())
  output$output_df <- renderDT(output_df())
}

shinyApp(ui = ui, server = server)




