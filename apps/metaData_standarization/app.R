

library(shiny)
library(bs4Dash)
library(DT)
source("lib_standards.R")

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    fluidRow(
      tabBox(width = 12,
             id="tabset1",
             tabPanel(
               id="tab_dt",
               title = "Raw data frame",
               selectInput(inputId = "in_file",
                              label="Select your metadata",
                           choices = ls(".GlobalEnv")[sapply(ls(".GlobalEnv"), function(i) class(get(i))) =="data.frame"]),
               textInput(inputId = "in_project_name", label = "Introduce the name of this project", value = "New"),
               actionButton(inputId = "in_go",
                            label = "Submit file"),
               DTOutput(outputId = "out_table")
             ),
             tabPanel(title = "Standardized data frame", DTOutput(outputId = "out_table_2"))
      ),
      box(id="b_donorInfo", title  = "Donor Information", width = 4,
             tabsetPanel(id="set_donorInfo",
                         tabPanel(
                                  title = "Organism & ID",
                                  uiOutput("out_D_org_id")
                                  ),
                         tabPanel(
                                 id="tab_2",
                                 title = "Donor sex",
                                 uiOutput("out_D_sex")
                                 ),
                         tabPanel(
                                 id="tab_3",
                                 title = "Donor age",
                                 uiOutput("out_D_age")
                                 )
                         )
      ),
      box(width = 4,
          id="b_experiment",
          title = "Experiment Information",
          tabsetPanel(id="exper_information",
                      tabPanel(
                        title = "Data acquisition mode",
                        selectizeInput(inputId = "in_acq", label = "Select the data acquisition mode",
                                    choices = c("scRNA", "nuclei", "scATAC"),
                                    options=list(create=TRUE))
                        ),
                      tabPanel(
                        title = "Experimental conditions",
                        uiOutput("out_condition_col")
                      )
                      )
          ),
      box(width = 4,
          title = "Cell Information",
          tabsetPanel(id="set_cellInfo",
                      tabPanel(
                        title = "Brain Region",
                        selectInput(inputId = "in_stage", label = "Select Organism Stage",
                                    choices = c("developing", "adult")),
                        uiOutput(outputId = "out_brainOri"),
                        uiOutput(outputId = "out_brainIn"),
                      ),
                      tabPanel(
                        title = "Cell type",
                        uiOutput("out_cellType_col")
                      ),
                      tabPanel(
                        title = "Barcode",
                        uiOutput("out_barcode_col")
                      )

          )),
      box(id = "bbb", width = 3,
          textInput(inputId = "in_varName",
                    label = "Introduce the name of the variable storing the generated data frame", value = "New"),
          actionButton(status = "success",inputId = "act_save",label = "Standardize data")
          )

    )
))


server <- function(input, output) {

  #df <- eventReactive(input$in_go, {f <- get(input$in_file)})

  observeEvent(input$in_go, {
     f <- get(input$in_file)
     #browser()

     output$out_table <- renderDT(head(f))
     #-----------------------
     # Donor information
     #-----------------------

     # ==== Organism & age ====
     output$out_D_org_id <- renderUI(
       tabPanel(title = "title",
                box(width = 7,
                    id="box_organism",
                    title = "Organism",
                    selectInput(inputId = "in_organism", label = "Select donor organism", choices = c("human", "mouse")),
                    collapsible = FALSE
                    ),
                box(width = 7,
                    id="box_id",
                    title = "Dono ID",
                    selectizeInput(inputId = "in_id", label = "Type donor ID",  choices=colnames(f),
                                   options=list(create=TRUE)),
                    collapsible = FALSE
                    )
       ))

     # ==== Sex ====
     output$out_D_sex <- renderUI(
       box(width = 7,
           id="box_Dsex",
           title = "Donor sex",
           collapsible = FALSE,
        selectizeInput(inputId = "in_sex_col", label = "Select the column that refers to sex (if none, type 'NA')", choices=colnames(f),
                       options=list(create=TRUE), multiple=TRUE),
        uiOutput("out_sex_male"),
        uiOutput("out_sex_female")
        ))
     observeEvent(input$in_sex_col, {
       #browser()
       if (input$in_sex_col %in% colnames(f)){
         output$out_sex_male <- renderUI(
           selectizeInput(inputId = "in_sex_male", label = "What value refers to male?",
                          choices=unlist(strsplit(gsub('["c(|)"]',"", levels(as.factor(f[input$in_sex_col]))), ",")),
                          options=list(create=TRUE)))

         output$out_sex_female <- renderUI(
           selectizeInput(inputId = "in_sex_female", label = "What value refers to female?",
                          choices=unlist(strsplit(gsub('["c(|)"]',"", levels(as.factor(f[input$in_sex_col]))), ",")),
                          options=list(create=TRUE)))

       } else {output$out_sex_male <- renderUI(selectizeInput(inputId = "in_sex",
                                                              label = "Select default donor sex:",
                                                              choices=c("female", "male")))}
       }, priority = 0)

     # ==== Age ====
     output$out_D_age <- renderUI(

       box(width = 7,
           id="box_age",
           title = "Donor age",
           collapsible = FALSE,
           selectizeInput(inputId = "in_age_col", label = "Select the column (if none, type 'NA')", choices=colnames(f),
                          options=list(create=TRUE)
                          ),
           uiOutput(outputId = "out_age_val")
      ))

     observeEvent(input$in_age_col,{

       if (!input$in_age_col %in% colnames(f)){
         output$out_age_val <- renderUI(textInput(inputId = "in_age_value", label = "Type donor age"))}
     })

     #-----------------------
     # Cell information
     #-----------------------

     # === Brain region ====
     df <- reactive(read.csv("../acronym_df.csv"))


     observeEvent(input$in_organism,{

       observeEvent(input$in_stage, {
         switch(input$in_organism,
                mouse = {if(input$in_stage == "developing"){
                  output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri",
                                                                 label = "Brain region of origin",
                                                                 choices = df()$mouse_developing))
                }else{output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri",
                                                                     label = "Brain region of origin",
                                                                     choices = df()$mouse_adult))}
                },
                human = {if(input$in_stage == "developing"){
                  output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri", #in
                                                                 label = "Brain region of origin",
                                                                 choices = df()$human_developing))
                }else{output$out_brainOri <- renderUI(selectizeInput(inputId = "in_brainOri",
                                                                     label = "Brain region of origin",
                                                                     choices = df()$human_adult))}

                }
         )
         switch (input$in_organism,
                 mouse = {if(input$in_stage == "developing"){
                   output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                                 label = "Select brain region of injection",
                                                                 choices = df()$mouse_developing))
                 }else{output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                                     label = "Select brain region of injection",
                                                                     choices = df()$mouse_adult))}
                 },
                 human= {if(input$in_stage == "developing"){
                   output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                                 label = "Select brain region of injection",
                                                                 choices = df()$human_developing))
                 }else{output$out_brainIn <- renderUI(selectizeInput(inputId = "in_brainIn",
                                                                     label = "Select brain region of injection",
                                                                     choices = df()$human_adult))}

                 }
         )

       })
     })

     # === Cell type ====
     output$out_cellType_col <- renderUI(
       box(id = "box_cell_type",
           width = 12,
           selectizeInput(inputId = "in_cell_type_col", label="Select the column containing the cell type (if none, type 'NA')",
                      choices = colnames(f),
                      options=list(create=TRUE)),
                      uiOutput(outputId = "out_cell_type")
       )
     )

     observeEvent(input$in_cell_type_col, {
       if (!input$in_cell_type_col %in% colnames(f)) {
         output$out_cell_type <- renderUI(textInput(inputId = "in_cell_type", label = "Type default cell type, if any"))}
     })


     # === Barcode ====
     output$out_barcode_col <- renderUI(
       selectizeInput(inputId = "in_barcode_col", label="Select the column containing barcodes (if none, type 'NA')",
                      choices=colnames(f), options=list(create=TRUE))
     )

     # === Condition ====
     output$out_condition_col <- renderUI(
       selectizeInput(inputId = "in_condition_col", label="Select the column containing the experimental conditions",
                      choices=colnames(f), options=list(create=TRUE))
     )
    })
  #in_organism, in_id, in_sex_col, in_sex_male, in_sex_female, in_sex, in_age_col,in_acq,
  #in_brainOri, in_brainIn, in_cell_type, in_barcode_col, in_condition_col

  observeEvent(input$act_save, {
    #browser()
    l <- list(organism=input$in_organism,
              donor_id=input$in_id,
              sex_col=input$in_sex_col,
              sex_male_value=input$in_sex_male,
              sex_female_value=input$in_sex_female,
              sex_default=input$in_sex,
              donor_age_col=input$in_age_col,
              donor_age=input$in_age_value,
              data_mod=input$in_acq,
              brain_injec=input$in_brainIn,
              brain_origin=input$in_brainOri,
              cell_type_col=input$in_cell_type_col,
              cell_type=input$in_cell_type,
              barcode_col=input$in_barcode_col,
              condition_col=input$in_condition_col
    )
    #browser()

    out_df <- format_df(l, get(input$in_file))

    output$out_table_2 <- renderDT(head(out_df))

    wd <- "../metadata_processed"
    if(!dir.exists(wd)) dir.create(wd)
    jsonlite::write_json(l, pretty=TRUE, auto_unbox=TRUE, path = file.path(wd, paste(input$in_project_name, ".json", sep = "")))

    write.csv(out_df, file.path(wd, paste0(input$in_project_name, ".csv")), row.names = FALSE)
    assign(input$in_varName, out_df, envir = .GlobalEnv)
    #test <<- out_df
     })



}

# Run the application
shinyApp(ui = ui, server = server)




















