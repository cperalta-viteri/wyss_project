

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
               # fileInput(inputId = "in_file",
               #           label = "Select your metadata file",
               #           buttonLabel = "Browse",
               #           placeholder = "No file selected"),
               selectInput(inputId = "in_file",
                              label="Select your metadata",
                           choices = ls(".GlobalEnv")[sapply(ls(".GlobalEnv"), function(i) class(get(i))) =="data.frame"]),
               textInput(inputId = "in_project_name", label = "Introduce the name of this projec", value = "New"),
               actionButton(inputId = "in_go",
                            label = "Submit file"),
               DTOutput(outputId = "out_table")
             )
      ),
      box(width = 4,
          id="b_sex",
          title = "Sex",
          uiOutput("out_sex")
          ),
      box(width = 4,
          id="b_modality",
          title = "Modality",
          uiOutput("out_modality")
          ),
      box(width = 4,
          id="b_brain",
          title = "Brain info",
          #textOutput(outputId = "disp_organism"),
          #tags$head(tags$style("#disp_organism{color:blue; font-size: 20px}")),
          #textInput(inputId = "in_strain", label = "Type the genotype"),
          #selectInput(inputId = "in_acq", label = "Select the data acquisition mode", choices = c("scRNA", "nuclei", "scATAC")),
          selectInput(inputId = "in_stage", label = "Select Organism Stage", choices = c("developing", "adult"), selected = NULL),
          uiOutput(outputId = "out_brainOri"),
          uiOutput(outputId = "out_brainIn"),
          actionButton(inputId = "act_save",label = "Save metadata")
          #textInput(inputId = "in_uid", label = "Type the accession number"),
          #actionButton(inputId = "download", label = "Click to Download")
          )
    )
))


server <- function(input, output, session) {

  observeEvent(input$in_go, {
     #browser()
     f <- get(input$in_file)
     #f <- read.csv(f$datapath, sep = "\t")
     output$out_table <- renderDT(head(f))
     # ==== sex ====
     output$out_sex <- renderUI(
       box(width = 7,
           id="box_sex",
        selectizeInput(inputId = "in_sex_col", label = "What column refers to sex?", choices=colnames(f),
                       options=list(create=TRUE)),
        uiOutput("out_sex_male"),
        uiOutput("out_sex_female")
       ))

     observeEvent(input$in_sex_col, {
       #browser()
       output$out_sex_male <- renderUI(
         selectizeInput(inputId = "in_sex_male", label = "What value refers to male?",
                     choices=unlist(strsplit(gsub('["c(|)"]',"", levels(as.factor(f[input$in_sex_col]))), ",")),
                     options=list(create=TRUE))
         )
       output$out_sex_female <- renderUI(
         selectizeInput(inputId = "in_sex_female", label = "What value refers to female?",
                        choices=unlist(strsplit(gsub('["c(|)"]',"", levels(as.factor(f[input$in_sex_col]))), ",")),
                        options=list(create=TRUE))
       )
       })
     # ==== modality ====
     output$out_modality <- renderUI(
       box(width = 7, id="box_modality",
           selectizeInput(inputId = "in_mod_col", label = "What column refers to modality?", choices = colnames(f),
                          options=list(create=TRUE)),
           uiOutput("out_mod_type")

     ))

     observeEvent(input$in_mod_col, {
       output$out_mod_type <- renderUI(
         selectizeInput(inputId = "in_mod_type", label="Choose the modality:",
                        choices = unlist(strsplit(gsub('["c(|)"]', "", levels(as.factor(f[input$in_mod_col]))), ",")),
                        options=list(create=TRUE))
       )
     })

     # === brain info ====
     df <- reactive(read.csv("../acronym_df.csv"))
     observeEvent(input$in_stage, {

       if(input$in_stage == "developing"){
         output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri",
                                                        label = "Brain region of origin",
                                                        choices = df()$human_developing))
       }else{output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri",
                                                            label = "Brain region of origin",
                                                            choices = df()$human_adult))}
       if(input$in_stage == "developing"){
         output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                       label = "Select brain region of injection",
                                                       choices = df()$human_developing))
       }else{output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                           label = "Select brain region of injection",
                                                           choices = df()$human_adult))}
     })


     observeEvent(input$act_save, {
       l <- list(sex_col=input$in_sex_col,
                 sex_male_value=input$in_sex_male,
                 sex_female_value=input$in_sex_female,
                 modality_col=input$in_mod_col,
                 modality=input$in_mod_type)
       wd <- "../metadata"
       if(!dir.exists(wd)) dir.create(wd)
       jsonlite::write_json(l, pretty=TRUE, auto_unbox=TRUE, path = file.path(wd, input$in_project_name))
     })


    })

}

# Run the application
shinyApp(ui = ui, server = server)




















