
library(shiny)
library(shinyWidgets)
library(bs4Dash)
library(rentrez)


# Handle pressing return
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'

download_hit <- function(h) {
  outdir <- file.path("downloads",h$accession)
  dir.create(outdir)
  system(sprintf("wget -np -nd -r -P %s %s", outdir, file.path(h$ftplink, "suppl")))
}

hit_box <- function(h) {
  btn_id <- paste0("download_",h$accession)
  box(
    width = 12,status = "success",collapsible = FALSE,
    title = h$title,
    fluidRow(tags$b("Abstract:"),h$summary),
    footer = 
      fluidRow(
        column(width=4,align="center",h$pdat),
        column(width=4,align="center",tags$b(h$n_samples,"samples")),
        column(width=4,align="center",actionButton(inputId = btn_id, label = tagList(icon("download"),h$accession), status = "primary")),
        column(width=12,align="center",tags$a("See record on NCBI",href=sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s",h$accession),target="_blank"))
      )
  )
}


ui <- dashboardPage(
    dashboardHeader(
      title = "GEO download",
      tags$head(tags$script(HTML(jscode)))
    ),
    dashboardSidebar(disable = TRUE),
    dashboardBody(
      fluidRow(
        box(width = 12,
          title = "Query",
          pickerInput("query_org",label="Organism",choices=c("Mus musculus", "Homo sapiens"),multiple = FALSE),
          tagAppendAttributes(
            textInput(
              inputId = "query_text",
              label = "Query Text",
              value = ""
            ),
            `data-proxy-click` = "submit_query"
          ),
          infoBoxOutput("queryIB",width = 12),
          actionButton(inputId = "submit_query", label = "Submit", status = "primary")
        )
      ),
      fluidRow(
            box(width = 12,
                title = "Results",
                valueBoxOutput("resultCountIB"),
                uiOutput("results"),
                footer = column(offset=2,width=8,align="center",
                  textOutput("current_result_page",inline=TRUE),"/",textOutput("total_result_page",inline=TRUE),br(),
                  actionButton("prev_page","<"),actionButton("next_page",">")
                )
            )
        )
    )
)



# Define server logic
server <- function(input, output, session) {
    query <- reactive({
      c(
        paste0(input$query_org,"[Organism]"),
        input$query_text
      ) |>
      paste0(collapse = " AND ")
    })
    
    output$queryIB <- renderInfoBox({infoBox("Full query string",query(),color = "info",fill = TRUE)})
    
    page_size <- 20L
    current_page <- reactiveVal(1L)
    observeEvent(input$next_page,{current_page(isolate(current_page()) + 1L)})
    observeEvent(input$prev_page,{current_page(max(isolate(current_page()) - 1L,1L))})
    output$current_result_page <- reactive({current_page()})

    observeEvent(input$submit_query,current_page(1L))
    
    search <- eventReactive(input$submit_query,{
      entrez_search("gds",query(), use_history=TRUE,retmax=0,retmode = "xml",rettype="xml")
    })
    hits <- reactive({
      entrez_summary("gds",web_history=search()$web_history,retmax=page_size,retstart=(current_page()-1L)*page_size,retmode = "json")
    })
    
    output$total_result_page <- reactive({((search()$count-1L) %/% page_size) + 1L})
    output$resultCountIB <- renderValueBox(valueBox(search()$count,"query hits in total"))
    output$selectedResultCountIB <- renderValueBox(valueBox(0,"selected"))

    output$results <- renderUI({
      lapply(hits(),function(h) {
        observeEvent(input[[paste0("download_",h$accession)]],download_hit(h))
        hit_box(h)
      })
    })
}

# Run the application
shinyApp(ui = ui, server = server)



















