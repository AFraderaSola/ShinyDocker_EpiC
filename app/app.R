############################
######## Libraries #########
############################

library(shinythemes)
library(mailtoR)
library(DT)
library(shinyBS)
library(plotly)
# library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(rMQanalysis)
set.seed(666)

############################
######### Script ###########
############################

source('myui.R', local = TRUE)
source('myserver.R')

shinyApp(
  ui = ui,
  server = server
)

