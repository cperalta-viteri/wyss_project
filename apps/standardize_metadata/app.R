
library(DT)
source("lib_std.R")


fields <- list(
  project_id = field_item(),
  batch_id = field_item(),
  cell_name = field_item(),
  donor_organism = field_item(targetval_suggestions=c("mouse","human")),
  donor_sex = field_item(targetval_suggestions=c("M","F"),validation="^[MF]$",rules=list(case_when_item(),case_when_item())),
  donor_age = field_item(targetval_suggestions=c("E18","P0","P1","P2"),validation="[PE][0-9]+(\\.5)?"),
  donor_mouse_strain = field_item(targetval_suggestions=c("C56BL6","FVB")),
  donor_mouse_line = field_item(rules=list(case_when_item())),
  donor_mouse_genotype = field_item(rules=list(case_when_item())),
  
  donor_human_genetic_population = field_item(),
  donor_human_phenotype = field_item(),
  donor_mouse_genotype = field_item(),
  
  sample_protocol_dissociation = field_item(),
  sample_protocol_chemistry = field_item(),
  sample_tissue_type = field_item(targetval_suggestions=c("Postnatal","Embryonic","Organoid","Neuronal Culture")),
  sample_tissue_origin = field_item(),
  sample_nervous_region = field_item(),
  
  default_condition = field_item(),
  disease_model = field_item()
)



ui <- dashboardPage(
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
          radioGroupButtons("current_working_field_name",choices = names(fields),status = "info",individual = TRUE),
          box(title = "Mapping rules",width = 12,collapsible = FALSE,
            uiOutput("fields"),
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

server <- function(input, output, session) {
  input_df <- reactive({iris})

  output_df <- reactive({
    input_df() |>
      filter(Sepal.Length>5)
  })
  
  cwf <- reactiveValues()
  
  observeEvent(input$current_working_field_name,{cwf$name <- input$current_working_field_name})
  
  # Copy selected field into cwf
  observeEvent(cwf$name,{
    for(nm in names(fields[[cwf$name]])) {
      cwf[[nm]] <- fields[[cwf$name]][[nm]]
    }
  })
  
  
  # Add a rule at the end of the list
  observeEvent(input$rule_add,{
    cwf$rules[[isolate(length(cwf$rules)) + 1L]] <- case_when_item()
  })
  
  observeEvent(input$test_btn,{
    print("test !")
    for(nm in names(cwf)) {
      cwf[[nm]] <- NULL
    }
    browser()
  })
  

  output$fields <- renderUI({
      rules <- lapply(seq_along(cwf$rules),function(i) {
        rule_id <- paste0("rule_",i)
        observeEvent(input[[paste0(rule_id,"_incol")]],{
          choices <- unique(as.character(input_df()[[input[[paste0(rule_id,"_incol")]]]]))
          updatePickerInput(session,paste0(rule_id,"_inval"),choices = choices)
        })
        
        caseWhenItemInput(
          rule_id,
          incol_choices = colnames(input_df()),
          targetval_choices = cwf$targetval_suggestions
        )
      })
      rules
  })

  # lapply(seq_along(cwf()$rules),function(i) {
  #   rule_id <- paste0("rule_",i)
  #   observeEvent(input[[paste0(rule_id,"_incol")]],{
  #     print(rule_id)
  #     updatePickerInput(paste0(rule_id,"_inval"),choices = unique(as.character(input_df()[[input[[paste0(rule_id,"_incol")]]]])))
  #   })
  # })
  
  
  output$input_df <- renderDT(input_df())
  output$output_df <- renderDT(output_df())
}

shinyApp(ui = ui, server = server)




