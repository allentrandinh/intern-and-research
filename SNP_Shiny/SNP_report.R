library(stringr)
library(tidyverse)
library(shiny)
library(dplyr)

#Loading Data
gemini_data <- read.delim("./Galaxy_GEMINI_query.tabular")
unique_data <- read.table("./cleaned_data.tabular")

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Human-diseased-related SNPs Report"),
  
  # Sidebar with a slider input for number of bins 
  tabsetPanel(
    tabPanel("Raw_Data",
             sidebarLayout(
               sidebarPanel(textInput("disease_keyword", label = h4("Diseases Key Word")),submitButton(text = "Submit"),
                            checkboxGroupInput("columns", "Select columns to display",choices = names(gemini_data))),
               mainPanel(tableOutput("raw_data_key_word")
               ))),
    tabPanel("Unique_Data", dataTableOutput("unique_data"))
    ))

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$raw_data_key_word <- renderTable({
    c <- tolower(input$disease_keyword)
    table_filtered <-filter(gemini_data,grepl(c,gemini_data$clinvar_disease_name,ignore.case=TRUE))
    columns = names(gemini_data)
    if (!is.null(input$columns)){
      columns = input$columns
    }
    table_filtered[,columns,drop=FALSE]})
  output$unique_data <- renderDataTable(unique_data)
  }

# Run the application 
shinyApp(ui = ui, server = server)
