
library(shiny)
library(DT)
source("lib_std.R")



std_fields <- list(
  field_item(
    field_name = "batch_id",
    description = "The batch group together cells from the same experiment."
  ),
  field_item(
    field_name = "cell_name",
    description = "Identify a cell within an experiment (e.g. barcode)."
  ),
  field_item(
    field_name = "donor_id",
    description = "Identifier of the donor where the cell is coming from."
  ),
  field_item(
    field_name = "donor_sex",
    description = "Sex of the donor where the cell is coming from.",
    targetval_suggestions = c("M","F"),
    validation = "^[MF]$",rules=list(case_when_item(),case_when_item())
  ),
  field_item(
    field_name = "donor_age",
    description = "Age of the donor where the cell is coming from.",
    targetval_suggestions = c("E18","P0","P1","P2"),
    validation = "[PE][0-9]+(\\.5)?"
  ),
  field_item(
    field_name = "donor_mouse_strain",
    description = "Strain of the donor where the cell is coming from.",
    targetval_suggestions = readLines("data/jax_inbred_strains.txt")
  )
  # donor_mouse_line = field_item(targetval_suggestions=readLines("data/mouse_line.txt"),rules=list(case_when_item())),
  # donor_mouse_genotype = field_item(rules=list(case_when_item())),
  # 
  # donor_human_genetic_population = field_item(),
  # donor_human_phenotype = field_item(),
  # donor_mouse_genotype = field_item(),
  # 
  # sample_protocol_dissociation = field_item(),
  # sample_protocol_chemistry = field_item(),
  # sample_tissue_type = field_item(targetval_suggestions=c("Postnatal","Embryonic","Organoid","Neuronal Culture")),
  # sample_tissue_origin = field_item(),
  # sample_nervous_region = field_item(),
  # 
  # default_condition = field_item(),
  # disease_model = field_item()
)



fieldInput <- function(f) {
  tags$div(
    id = f$field_name,
    tags$p(tags$b("Description:"),helpText(f$description)),
    lapply(seq_along(f$rules),\(i) caseWhenItemInput(sprintf("%s_%d",f$field_name,i),targetval_choices = f$targetval_suggestions)),
    actionBttn(sprintf("%s_add",f$field_name),icon=icon("add"),size="xs")
  )
}


ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    fluidRow(
      box(width=12,title = "Input data.frame",
        #style = "font-size: 75%; height: 50%",
        DTOutput(outputId = "input_df")
      )
    ),
    fluidRow(
      box(width = 12,title = "Metadata columns mapping rules",collapsible = FALSE,
        fluidRow(
            column(width=4,
                   textInput("project_id","project_id","my_project"),
                   selectInput("organism","organism",c("Mouse","Human"),selectize = TRUE)
            )
        ),
        fluidRow(
            do.call(what=\(...) tabBox(...,width=12,collapsible = FALSE),lapply(std_fields,\(f) tabPanel(f$field_name,fieldInput(f))))
        )
      )
    ),
    fluidRow(
      textOutput("mutate_string")
    ),    
    fluidRow(
      tabBox(height = "400px",
        width = 12,
        tabPanel(title="Standardized Output",
          DTOutput(outputId = "output_df")
        ),
        tabPanel(title="Transformation",
          verbatimTextOutput("transformation")
        )
      )
    )
  )
)



server <- function(input, output, session) {
  input_df <- reactive({
    iris
  })
  
  output_df <- reactive({
    input_df() |> filter(Sepal.Length>5)
  })
  
  output$input_df <- renderDT(
    input_df() |> 
      datatable(options = list(dom="t",pageLength=5L,scrollY="200px",paging=FALSE,scrollCollapse=TRUE),selection = "none",filter="top",rownames = FALSE)
  )
  output$output_df <- renderDT({
    output_df() |>
      datatable(options = list(dom="t",pageLength=5L,scrollY="200px",paging=FALSE,scrollCollapse=TRUE),selection = "none",filter="top",rownames = FALSE)
  })
  
  output$transformation <- renderText({
    paste0(
      "input_df |>\n",
      paste0("    mutate(",sapply(std_fields,"[[","field_name")," = case_when(\n","    ))",collapse=" |>\n")
    )
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
      
      caseWhenItemInput(
        rule_id,
        incol_choices = colnames(input_df()),
        targetval_choices = cwf()$targetval_suggestions
      )
    })
    rules
  })
  
}

shinyApp(ui, server)
