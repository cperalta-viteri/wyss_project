
library(DT)
source("lib_std.R")


std_fields <- list(
  project_id = field_item(),
  batch_id = field_item(),
  cell_name = field_item(),
  donor_organism = field_item(targetval_suggestions=c("mouse","human")),
  donor_sex = field_item(targetval_suggestions=c("M","F"),validation="^[MF]$",rules=list(case_when_item(),case_when_item())),
  donor_age = field_item(targetval_suggestions=c("E18","P0","P1","P2"),validation="[PE][0-9]+(\\.5)?"),
  donor_mouse_strain = field_item(targetval_suggestions=readLines("data/jax_inbred_strains.txt")),
  donor_mouse_line = field_item(targetval_suggestions=readLines("data/mouse_line.txt"),rules=list(case_when_item())),
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

server <- function(input, output, session) {
  input_df <- reactive({iris})
  
  output_df <- reactive({
    input_df() |>
      filter(Sepal.Length>5)
  })
  
  # This reactive variable store the current working field (cwf)
  cwf <- reactiveVal()
  
  # Copy selected field into cwf
  observeEvent(input$current_working_field_name,{
      cwf(std_fields[[input$current_working_field_name]])
  })
  
  # Add a rule at the end of the list
  observeEvent(input$rule_add,{
    x <- isolate(cwf())
    x$rules[[length(x$rules) + 1L]] <- case_when_item()
    print("add")
    cwf(x)
  })

  output$cwf_rules <- renderUI({
      rules <- lapply(seq_along(cwf()$rules),function(i) {
        rule_id <- paste0("rule_",i)
        
        # When "_incol" change => update _inval
        observeEvent(input[[paste0(rule_id,"_incol")]],{
          choices <- unique(as.character(input_df()[[input[[paste0(rule_id,"_incol")]]]]))
          updatePickerInput(session,paste0(rule_id,"_inval"),choices = choices)
        })
        
        # When the "_rm" button is pressed
        # observeEvent(input[[paste0(rule_id,"_rm")]],{
        #   print("rm")
        #   print(i)
        #   x <- isolate(cwf())
        #   x$rules <- x$rules[-i]
        #   cwf(x)
        # })
        
        caseWhenItemInput(
          rule_id,
          incol_choices = colnames(input_df()),
          targetval_choices = cwf()$targetval_suggestions
        )
      })
      rules
  })
  
  output$input_df <- renderDT(input_df())
  output$output_df <- renderDT(output_df())
}

shinyApp(ui = ui, server = server)




