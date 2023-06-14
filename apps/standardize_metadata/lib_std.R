
library(tidyverse)
library(shiny)
library(shinyWidgets)
library(bs4Dash)



case_when_item <- function(type=c("when","text"),incol=NULL,inval=NULL,targetval=NULL,text=NULL) {
  type <- match.arg(type)
  list(
    type = type,
    incol = incol,
    inval = inval,
    targetval = targetval,
    text = text
  )
}

field_item <- function(targetval_suggestions=NULL,validation=NULL,rules=list()) {
  list(
    targetval_suggestions = targetval_suggestions,
    validation = validation,
    rules = rules
  )
}

caseWhenItemInput <- function(id,incol_choices=NULL,inval_choices=NULL,targetval_choices=NULL) {
  tags$table(
    tags$tr(
      width = "100%",
      tags$td(pickerInput(paste0(id,"_type"),choices=c("when","text"),inline=TRUE)),
      tags$td(
        conditionalPanel(
          condition = sprintf("input.%s_type == 'when'",id),
          span(
            pickerInput(paste0(id,"_incol"),choices=incol_choices,inline = TRUE),
            "in",
            pickerInput(paste0(id,"_inval"),choices = inval_choices,inline = TRUE,multiple = TRUE, options = pickerOptions(style="label label-info")),
            "then",
            selectizeInput(paste0(id,"_targetval"),label=NULL,choices=targetval_choices,options = list(create=TRUE))
          )  
        ),
        conditionalPanel(
          condition = sprintf("input.%s_type == 'text'",id),
          textInput(paste0(id,"_incol"),label="",value = "col1 == col2 ~ col3"),
        )
      ),
      tags$td(
        actionButton(paste0(id,"_up"),icon("circle-arrow-up")),
        actionButton(paste0(id,"_dn"),icon("circle-arrow-down")),
        actionButton(paste0(id,"_rm"),icon("trash"))
      )
    )
  )
  
}


mutate_string <- function(f) {
  r <- lapply(seq_along(standard_fields[[fn]]$rules),function(i) {
    incol <- input[[paste0("rule_",fn,"_",i,"_incol")]]
    inval <- input[[paste0("rule_",fn,"_",i,"_inval")]]
    targetval <- input[[paste0("rule_",fn,"_",i,"_targetval")]]
    paste(incol," %in% c(",paste0('"',inval,'"',collapse=","),") ~ ",targetval)
  })
  r[length(r) + 1L] <- sprintf('TRUE ~ NA')
  sprintf("mutate(%s=case_when(%s))",fn,paste(r,collapse = ","))
}




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





