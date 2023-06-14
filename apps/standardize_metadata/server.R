
library(DT)
source("lib_std.R")


function(input, output, session) {
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





