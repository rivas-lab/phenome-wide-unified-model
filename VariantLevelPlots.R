library(shiny)
library(plotly)
library(readr)
library(dplyr)
library(DT)
library(vroom)

# File path for the local data
file_path <- "variant_level_results.tsv.gz"  # Replace with your actual path

# UI
ui <- fluidPage(
  titlePanel(
    tags$div(
      tags$h1("Unified Meta Regression Model Results 400k Exomes")
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      
      # Add dropdown menus for gene and phenotype
      selectInput("geneSelect", "Select Gene:", choices = NULL),  # Populated dynamically
      selectInput("phenotypeSelect", "Select Phenotype:", choices = NULL),  # Populated dynamically
      
      # Add output for displaying selected gene and phenotype
      textOutput("geneLabel"),
      textOutput("phenotypeLabel")
    ),
    
    mainPanel(
      plotlyOutput("constraintPlot"),     # Coefficient plot
      tags$br(), tags$br(),
      plotlyOutput("pathogenicityPlot"),  # Constraint plot
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Load the data once
  data <- vroom(file_path, delim = "\t", show_col_types = FALSE)

  # Populate dropdown menu choices for genes
  updateSelectInput(session, "geneSelect", choices = unique(data$gene))

  # Observe changes in the selected gene and update phenotype dropdown
  observeEvent(input$geneSelect, {
    # Filter data to get available phenotypes for the selected gene
    available_phenotypes <- data %>%
      filter(gene == input$geneSelect) %>%
      pull(description) %>%
      unique()
    
    # Update the phenotype dropdown choices based on the selected gene
    updateSelectInput(session, "phenotypeSelect", choices = available_phenotypes)
  })
  
  # Reactive filtering based on selected gene and phenotype
  filtered_data <- reactive({
    req(input$geneSelect, input$phenotypeSelect)  # Ensure selections are made
    data %>%
      filter(gene == input$geneSelect, description == input$phenotypeSelect)  # Ensure correct column name
  })

  # Constraint Plot
  output$constraintPlot <- renderPlotly({
    plot_ly(filtered_data(), x = ~log_constraint, y = ~BETA, type = 'scatter', mode = 'markers', 
            text = ~paste("Variant ID:", markerID, "<br>Effect size:", BETA),
            hoverinfo = "text"
    ) %>%
        layout(
            title = paste("Constraint Analysis"),
            xaxis = list(title = "-log(1 - Constraint Probability)"), 
            yaxis = list(title = "Effect Size on Phenotype")
        )
  })

   # Pathogenicity Plot
   output$pathogenicityPlot <- renderPlotly({
        plot_ly(filtered_data(), x = ~log_pathogenicity, y = ~BETA, type = 'scatter', mode = 'markers',
                text = ~paste("Variant ID:", markerID, "<br> Effect size:", BETA),
                hoverinfo = "text"
        ) %>%
        layout(
            title = paste("Pathogenicity Analysis"),
            xaxis = list(title = "-log(1 - Pathogenicity Probability)"), 
            yaxis = list(title = "Effect Size on Phenotype")
        )
    })
}

# Run the app
shinyApp(ui = ui, server = server)
