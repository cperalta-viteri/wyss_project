

library(shiny)
library(bs4Dash)
library(DT)


ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    fluidRow(
      tabBox(width = 12,
             id="tabset1",
             tabPanel(
               id="tab_dt",
               title = "Data frame viz",
               fileInput(inputId = "in_file",
                         label = "Select your metadata file",
                         buttonLabel = "Browse",
                         placeholder = "No file selected"),
               actionButton(inputId = "in_go",
                            label = "Submit file"),
               DTOutput(outputId = "out_table")
             )
      ),
      box(width = 4,
          title = "Sex",
          id="b_sex",
          uiOutput("out_sex")),
      box(width = 4,
          title = "Modality",
          id="b_modality",
          uiOutput("out_modality"))
    )
))


server <- function(input, output) {
   observeEvent(input$in_go, {
     #browser()
     f <- input$in_file
     f <- read.csv(f$datapath, sep = "\t")
     output$out_table <- renderDT(head(f))
     # ==== sex ====
     output$out_sex <- renderUI(
       box(width = 7,
           id="box_sex",
        selectInput(inputId = "in_sex_col", label = "What column refers to sex?", choices=colnames(f)),
        uiOutput("out_sex_male"),
        uiOutput("out_sex_female")
       ))

     observeEvent(input$in_sex_col, {
       #browser()
       output$out_sex_male <- renderUI(
         selectizeInput(inputId = "in_sex_male", label = "What value refers to male?",
                     choices=unlist(strsplit(gsub('["c(|)"]',"", levels(as.factor(f[input$in_sex_col]))), ",")))
         )
       output$out_sex_female <- renderUI(
         selectizeInput(inputId = "in_sex_male", label = "What value refers to female?",
                        choices=unlist(strsplit(gsub('["c(|)"]',"", levels(as.factor(f[input$in_sex_col]))), ",")))
       )
       })
     # ==== modality ====
     output$out_modality <- renderUI(
       box(width = 7, id="box_modality",
           selectInput(inputId = "in_mod_col", label = "What column refers to modality?", choices = colnames(f)),
           uiOutput("out_mod_type")

     ))

     observeEvent(input$in_mod_col, {
       output$out_mod_type <- renderUI(
         selectInput(inputId = "in_mod_type", label="Choose the modality:",
                        choices = unlist(strsplit(gsub('["c(|)"]', "", levels(as.factor(f[input$in_mod_col]))), ",")))
       )
     })



    })



}

# Run the application
shinyApp(ui = ui, server = server)
