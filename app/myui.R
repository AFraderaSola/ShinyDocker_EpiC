############################
######## Libraries #########
############################

library(shinythemes)
library(mailtoR)
library(DT)
library(shinyBS)
library(plotly)
library(tidyverse)
library(rMQanalysis)
set.seed(666)

############################
####### Functions ##########
############################

Enriched2 <- function (x, y = NULL, c = 1, s0 = 1, pvalue = 0.05, ...) 
{
  if (class(x) == "data.frame") {
    if (length(x) == 2) {
      y <- x[[2]]
      x <- x[[1]]
    }
    else {
      stop("x hast to be a data frame with length 2!")
    }
  }
  results <- list()
  results$curve <- createEnrichmentCurve2(c, s0, pvalue, ...)
  results$positive <- (c/(x - s0) - y + -log10(pvalue) <= 0 & 
                         x >= s0)
  results$enriched <- results$positive 
  results$background <- !results$enriched
  class(results) <- append(class(results), "Enriched")
  return(results)
}

createEnrichmentCurve2 <-  function (c, s0, pvalue = 1, step = 0.01, xlim = 50, linetype = 2, 
                                     size = 0.1) 
{
  x <- seq(s0, xlim, step)
  y <- c/(x - s0) + -log10(pvalue)
  return(annotate("path", x = c(NA, x), y = c(NA, y), 
                  linetype = linetype, size = size))
}


############################
####### Variables ##########
############################

load("./InputFiles/00_DNA_5mCG.RData")

############################
######### Script ###########
############################

ui <-  navbarPage(title = HTML("<b>EpiC</b>: Identifying epigenetic regulators in the plantae kingdom"),
                  theme = shinytheme("lumen"),
                  windowTitle = "EpiC database",
                  tabPanel("Welcome",
                           fluidRow(column(12,
                                           h2("Abstract"),
                                           p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore 
                                             magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                             consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                             Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum. 
                                             Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore 
                                             magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo 
                                             consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. 
                                             Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
                                             Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore 
                                             magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                             consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                             Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."))),
                           fluidRow(column(12,
                                           h3("Experimental design"),
                                           br(),
                                           img(src = "./Welcome/Image_01.png",  height = "75%", width = "75%"),
                                           br(),
                                           br(),
                                           br(),
                                           p(strong("Epigenetic modification affinity-purification screen design."),
                                             "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore 
                                             magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                             consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                             Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."))),
                           fluidRow(column(12,
                                           tags$hr(style="border-color: darkgrey;"))),
                           fluidRow(column(12,
                                           strong("Identifying epigenetic regulators in the plantae kingdom"),
                                           br(),
                                           br(),
                                           p("Lars Teschke",tags$sup("1"),",",
                                             "Albert Fradera-Sola",tags$sup("1"),",",
                                             "Falk Butter",tags$sup("1,#"), style = "font-size:9pt;"),
                                           em("# Indicates correspondance,
                                           1 Quantitative Proteomics, Institute of Molecular Biology, Mainz, Germany",
                                              style = "font-size:8pt;"),
                                           br(),
                                           br(),
                                           p("App created by Albert Fradera-Sola in November 2023", style = "font-size:9pt;"),
                                           p("Comments and bug reports to the following e-mail: ", mailtoR(email = "A.FraderaSola@imb-mainz.de",
                                                                                                           subject = "Comments and bugs: RBP Interactome shiny app",
                                                                                                           text = "A.FraderaSola@imb-mainz.de"), style = "font-size:9pt;"),
                                           em("Last update: November 2023", style = "font-size:8pt;")))),
                  
                  tabPanel("DNA modifications",
                             tags$style(HTML(".tabbable > .nav > li > a                  {background-color: #f8f8f8}
                                              .tabbable > .nav > li[class=active]    > a {background-color: #f8f8f8}")),
                           tabsetPanel(
                             tabPanel(HTML("<p style=color:#56B4E9 ><sup style=color:#56B4E9 >5</sup>mCG</p>"),
                                      navbarPage(title = HTML("<p>AP-MS analysis</p>"),
                                                 tabPanel(title = "Overview",
                                                          br(),
                                                          fluidRow(
                                                            column(2,
                                                                   img(src = "./5mCG/Image_01.png",  height = "100%", width = "100%"),
                                                                   br(),
                                                                   br(),
                                                                   p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore 
                                                                     magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                                                     consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                                                     Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum. 
                                                                     Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore 
                                                                     magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo 
                                                                     consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. 
                                                                     Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
                                                                     Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore 
                                                                     magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                                                     consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                                                     Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
                                                                   ),
                                                            column(10,
                                                                   br(),
                                                                   plotlyOutput("DNA_5mCG_Overview_Barplot", height = "1000px")
                                                                   )
                                                          )),
                                                 tabPanel(title = "Identified Interactors",
                                                          br(),
                                                          fluidRow(
                                                            column(2,
                                                                   h3("Species"),
                                                                   br(),
                                                                   radioButtons(inputId = "DNA_5mCG_species",
                                                                                label = "Which species do want to display?",
                                                                                choices =sort(unique(gsub(pattern = "_",
                                                                                                          replacement = " ",
                                                                                                          x = names(DNA_5mCG_quantified_proteins_df)))),
                                                                                selected = sort(unique(gsub(pattern = "_",
                                                                                                            replacement = " ",
                                                                                                            x = names(DNA_5mCG_quantified_proteins_df))))[1]
                                                                   )
                                                            ),
                                                            column(3,
                                                                   h3("Quantified proteins"),
                                                                   br(),
                                                                   dataTableOutput("DNA_5mCG_QP_DT"),
                                                                   br(),
                                                                   h4("Selected proteins:"),
                                                                   selectizeInput('DNA_5mCG_QP_DT_ui_selectedIDs',
                                                                                  NULL,
                                                                                  choices=NULL,
                                                                                  multiple=TRUE),
                                                                   actionButton("update_DNA_5mCG_QP_DT",
                                                                                "Update selection"),
                                                                   br(),
                                                                   br()),
                                                            column(1,
                                                                   br()),
                                                            column(6,
                                                                   column(6,
                                                                          bsCollapsePanel("Volcano plot settings",
                                                                                          column(4,
                                                                                                 selectInput(inputId = "DNA_5mCG_Pvalue", 
                                                                                                             label = h4("Select a p-value:"), 
                                                                                                             choices = c(0.05,0.01))),
                                                                                          column(8,
                                                                                                 sliderInput(inputId = "DNA_5mCG_LFC",
                                                                                                             label = "Select a log2(Fold change):",
                                                                                                             min = 0, max = 5,
                                                                                                             value = 1, step = 0.1)))),
                                                                   column(12,
                                                                          plotlyOutput("DNA_5mCG_Volcanoplot", height = "1000px")))),
                                                          br(),
                                                          fluidRow(
                                                            column(6,
                                                                   br()),
                                                            column(6,
                                                                   p(strong("Volcano plot"),
                                                                     "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore
                                                        magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                                        consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                                 Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")))))),
                             tabPanel(HTML("<sup>5</sup>mCHG"),
                                      br(),
                                      fluidRow(
                                        column(2,
                                               h3("Species"),
                                               br(),
                                               # radioButtons(inputId = "DNA_5mCHG_species",
                                               #               label = "Which species do want to display?",
                                               #               choices =sort(unique(names(DNA_5mCHG_quantified_proteins_df))),
                                               #               selected = sort(unique(names(DNA_5mCHG_quantified_proteins_df)))[1]
                                               #               )
                                               ),
                                               column(10,
                                                      img(src = "Figure1.png",  height = "100%", width = "100%"))),
                                        br(),
                                        fluidRow(
                                          column(2,
                                                 br()),
                                          column(10,
                                                 p(strong("Volcano plot"),
                                                   "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore
                                                        magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                                        consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                                 Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")))))),
                  tabPanel("RNA modifications",
                           tabsetPanel(
                             tabPanel(HTML("m<sup>6</sup>A"),
                                      tabsetPanel(
                                        tabPanel("Nuclear extract",
                                                 br(),
                                                 fluidRow(
                                                   column(2,
                                                          h3("Species"),
                                                          br(),
                                                          # radioButtons(inputId = "RNA_m6A_nuclear_species",
                                                          #               label = "Which species do want to display?",
                                                          #               choices =sort(unique(names(RNA_m6A_nuclear_species))),
                                                          #               selected = sort(unique(names(RNA_m6A_nuclear_species)))[1]
                                                          #               )
                                                          ),
                                                          column(10,
                                                                 img(src = "Figure1.png",  height = "100%", width = "100%"))),
                                                   br(),
                                                   fluidRow(
                                                     column(2,
                                                            br()),
                                                     column(10,
                                                            p(strong("Volcano plot"),
                                                              "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore
                                                        magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                                        consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                                 Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")))),
                                        tabPanel("Total extract",
                                                 br(),
                                                 fluidRow(
                                                   column(2,
                                                          h3("Species"),
                                                          br(),
                                                          # radioButtons(inputId = "RNA_m6A_total_species",
                                                          #               label = "Which species do want to display?",
                                                          #               choices =sort(unique(names(RNA_m6A_total_species))),
                                                          #               selected = sort(unique(names(RNA_m6A_total_species)))[1])
                                                          ),
                                                          column(10,
                                                                 img(src = "Figure1.png",  height = "100%", width = "100%"))),
                                                   br(),
                                                   fluidRow(
                                                     column(2,
                                                            br()),
                                                     column(10,
                                                            p(strong("Volcano plot"),
                                                              "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore
                                                        magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
                                                        consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
                                                 Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."))))))))
                  
)

