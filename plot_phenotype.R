library(shiny)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(RColorBrewer)

# Get the list of files from the 'phenotypes' folder
file_list <- list.files("phenotypes", pattern = "*.gz", full.names = FALSE)

# Helper function to preprocess data
preprocess_data <- function(file_path) {
  # Load the data and rename columns for clarity
  results_data <- read.delim(gzfile(file_path), header = TRUE, sep = "\t")
  results_data <- results_data %>%
    rename(
      `p_value unified model` = p_model,
      `p_value log constraint` = p_log_constraint,
      `p_value log pathogenicity` = p_log_pathogenicity,
      `p_value pLoF indicator` = p_pLoF_indicator,
      `coef constant`=coef_constant,
      `coef log constraint`=coef_log_constraint,
      `coef log pathogenicity` =coef_log_pathogenicity,
      `coef pLoF indicator` = coef_pLoF_indicator,
      `p_value missense indicator` = p_missense_indicator,
      `coef missense indicator` = coef_missense_indicator
    )

  list(
    data = results_data,
    coef_cols = c("coef log constraint", "coef log pathogenicity", "coef pLoF indicator","coef missense indicator"),
    p_value_cols = c("p_value unified model", "p_value log constraint", "p_value log pathogenicity", "p_value pLoF indicator","p_value missense indicator")
  )
}

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
      
      # Initial dropdown menu for file selection
      selectInput("fileSelect", "Choose File:", choices = file_list, selected = NULL),
      
      # UI elements generated dynamically after file is selected
      uiOutput("dynamicFilters")
    ),

    mainPanel(
      plotlyOutput("coefPlot"),
      tags$br(), tags$br(),
      plotlyOutput("genePlot"),
      tags$br(),
      DTOutput("filteredResults")
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Load the data when a file is selected
  data_info <- reactive({
    req(input$fileSelect)
    preprocess_data(paste0("phenotypes/", input$fileSelect))
  })

  # for side panel:
  # Generate dynamic UI for p-value filters based on selected file's columns
  output$dynamicFilters <- renderUI({
    req(data_info())  # guards against null data input
    p_value_cols <- data_info()$p_value_cols
    
  tagList(
    textInput("geneSearch", "Search Gene:", ""),
    actionButton("searchButton", "Search"),
    lapply(seq_along(p_value_cols), function(i) {
      col <- p_value_cols[i]
      default_value <- if (i == 1) "0.00001" else ""  # Set a specific default for the first input
    
    textInput(
      inputId = paste0("filter_", col), 
      label = paste("Set α for", col), 
      value = default_value, 
      placeholder = "Enter a number"
    )
  })
)

  })

  # Reactive value to store the filtered data
  filteredData <- reactive({
    data <- data_info()$data
    p_value_cols <- data_info()$p_value_cols
    
    # Filter by gene search
    if (input$geneSearch != "") {
      data <- data %>% filter(grepl(input$geneSearch, gene, ignore.case = TRUE))
    }
    
    # Filter by text input fields
    for (col in p_value_cols) {
      filter_value <- as.numeric(input[[paste0("filter_", col)]])
      if (!is.na(filter_value)) {
        data <- data %>% filter(data[[col]] <= filter_value)
      }
    }
    
    data
  })

  selectedGene <- reactiveVal(NULL)
  
  output$filteredResults <- renderDT({
    data <- filteredData()
    numeric_cols <- sapply(data, is.numeric)
    numeric_colnames <- names(data)[numeric_cols]
    
    datatable(
      data,  
      selection = 'single',
      options = list(
        pageLength = 20,  
        lengthChange = TRUE,
        searching = FALSE,
        scrollX = TRUE
      )
    ) %>% 
      formatSignif(columns = numeric_colnames, digits = 5)
  })

  # Plot the p-values
  output$genePlot <- renderPlotly({
    data <- filteredData()
    selected_gene <- selectedGene()
    p_value_cols <- data_info()$p_value_cols
    
    if (nrow(data) == 0) return(NULL)

    gene_data <- data %>% filter(grepl(input$geneSearch, gene, ignore.case = TRUE))
    if (nrow(gene_data) == 0) return(NULL)
    
    plot_data <- gene_data %>%
      select(gene, all_of(p_value_cols)) %>%
      pivot_longer(cols = all_of(p_value_cols), names_to = "p_value_type", values_to = "p_value") %>%
      mutate(log_p_value = -log10(p_value)) %>%
      mutate(p_value_type = factor(p_value_type, levels = rev(p_value_cols)))
    
    plot_ly(
      plot_data, 
      x = ~log_p_value, 
      y = ~p_value_type, 
      type = 'scatter', 
      mode = 'markers',
      text = ~paste('Gene:', gene, '<br>p-value:', p_value),
      hoverinfo = 'text'
    ) %>%
      layout(
        title = paste("Plot of P-Values", input$geneSearch),
        xaxis = list(title = "-log(p-value)", type = 'linear', autorange = TRUE),
        yaxis = list(title = "P-Values", categoryorder = 'array', categoryarray = rev(p_value_cols)),
        plot_bgcolor = 'white',
        paper_bgcolor = 'white',
        font = list(color = 'black'),
        margin = list(l = 100)
      )
  })
  
  # Plot the coefficients
  output$coefPlot <- renderPlotly({
    data <- filteredData()
    coef_cols <- data_info()$coef_cols
    
    if (nrow(data) == 0) return(NULL)

    gene_data <- data %>% filter(grepl(input$geneSearch, gene, ignore.case = TRUE))
    if (nrow(gene_data) == 0) return(NULL)

    plot_data <- gene_data %>%
      select(gene, all_of(coef_cols)) %>%
      pivot_longer(cols = all_of(coef_cols), names_to = "coef_type", values_to = "coefficient") %>%
      mutate(coefficient = ifelse(coefficient > 1, 1, coefficient)) %>%
      mutate(coefficient = ifelse(coefficient < -1, -1, coefficient)) %>%
      mutate(coef_type = factor(coef_type, levels = rev(coef_cols)))
    
    plot_ly(
      plot_data, 
      x = ~coefficient, 
      y = ~coef_type, 
      type = 'scatter', 
      mode = 'markers',
      text = ~paste('Gene:', gene, '<br>coefficient:', coefficient),
      hoverinfo = 'text'
    ) %>%
      layout(
        title = paste("Plot of Coefficients", input$geneSearch),
        xaxis = list(title = "Coefficient Value", type = 'linear', autorange = TRUE),
        yaxis = list(title = "Coefficients", categoryorder = 'array', categoryarray = rev(coef_cols)),
        plot_bgcolor = 'white',
        paper_bgcolor = 'white',
        font = list(color = 'black'),
        margin = list(l = 100)
      )
  })
}

# Run the app
shinyApp(ui, server)
