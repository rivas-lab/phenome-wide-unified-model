# variantLevelResults.R
library(shiny)
library(plotly)
library(readr)
library(dplyr)
library(DT)
library(vroom)
library(tidyr)
library(stringr)
library(webshot2)
library(htmlwidgets)


# File path for the local data
file_path <- "variant_level_results_filtered.tsv.gz"  # path to the variant level results
gene_level_results <- "filtered_results.csv" # path to the gene level results

# UI
ui <- fluidPage(
  titlePanel(
    tags$div(
      tags$h2("Unified Meta Regression Model Results 400k Exomes"),
      tags$h4("Larissa Lauer, Manuel Rivas")
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      
      # Add dropdown menus for gene and phenotype
      selectInput("geneSelect", "Select Gene:", choices = NULL),  # Populated dynamically
      selectInput("phenotypeSelect", "Select Phenotype:", choices = NULL)  # Populated dynamically
    ),
    
    mainPanel(
      # Display the selected gene, phenotype, and p_model below the search
      div(style = "text-align: center;",
      tags$h3(textOutput("selectedGenePhenotype"))),

      div(style = "text-align: center;",
      tags$h4(textOutput("pModelValue"))),
      tags$br(), tags$br(),

      plotlyOutput("constraintPlot"),     # Coefficient plot
      tags$br(), tags$br(),
      plotlyOutput("pathogenicityPlot"),  # Constraint plot
      tags$br(), tags$br(),
      plotlyOutput("lofPlot"),     # Coefficient plot
      tags$br(), tags$br(),
      plotlyOutput("missensePlot"),
      tags$br(), tags$br()
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Load the data once
  data <- vroom(file_path, delim = "\t", show_col_types = FALSE)
  gene_results <- read.csv(gene_level_results)

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

  # pull the coefficients and p_values out of the gene level file
  gene_results_filtered <- reactive({
    req(input$geneSelect, input$phenotypeSelect)
    gene_results %>%
      filter(gene == input$geneSelect, description == input$phenotypeSelect) %>%
      select(p_model, 
            coef_log_constraint, p_log_constraint, 
            coef_log_pathogenicity,p_log_pathogenicity,
            coef_pLoF_indicator,p_pLoF_indicator,
            coef_missense_indicator,p_missense_indicator,
            coef_constant)
  })

  output$selectedGenePhenotype <- renderText({
  req(input$geneSelect, input$phenotypeSelect)  # Ensure selections are made
  paste(input$geneSelect, "-", input$phenotypeSelect)
})

  output$pModelValue <- renderText({
  req(gene_results_filtered())
  paste("Unified Model: p =", formatC(gene_results_filtered()$p_model, format = "e", digits = 3))
})

  
# Constraint Plot
output$constraintPlot <- renderPlotly({
    req(gene_results_filtered())  # Ensure gene_results_filtered is ready
    req(filtered_data())  # Ensure filtered data is ready
    
    # Pull coefficient and constant terms
    coef_constraint <- gene_results_filtered()$coef_log_constraint[1]
    constant <- gene_results_filtered()$coef_constant[1]
    p_constraint <- paste("(p =", formatC(gene_results_filtered()$p_log_constraint, format = "e", digits = 3), ")")
    
    # Calculate the y-values for the best fit line
    line_y <- ~coef_constraint * log_constraint + constant  # Element-wise multiplication
    
    # crop graphs
    plotdata <- filtered_data() %>%
        mutate(log_constraint = ifelse(log_constraint > 1.5, 1.5, log_constraint))

    # Plot data
    plot_ly(plotdata, x = ~log_constraint, y = ~BETA, type = 'scatter', mode = 'markers',
            error_y = list(
                array = ~SE,
                color = "gray",
                thickness = 0.5
            ), 
            text = ~paste("Variant ID:", markerID, "<br>Effect size:", BETA),
            hoverinfo = "text"
    ) %>%
       add_trace(
            x = ~log_constraint,
            y = line_y,  # Use precomputed y-values
            type = 'scatter',
            mode = 'lines',
            line = list(color = 'black'),
            name = 'Best Fit Line',
            error_y = NULL,
            showlegend = FALSE
       ) %>%
       layout(
            title = list(
              text = paste("Constraint Analysis<br><sub>", p_constraint, "</sub>"),
              x = 0.5  # Center the title and subtitle
            ),
            xaxis = list(title = "-log(1 - Constraint Probability)"), 
            yaxis = list(title = "Effect Size on Phenotype")
       )
})

# Pathogenicity Plot
output$pathogenicityPlot <- renderPlotly({
    req(gene_results_filtered())
    req(filtered_data())

    # Pull coefficient, constant, and p-value for pathogenicity
    coef_pathogenicity <- gene_results_filtered()$coef_log_pathogenicity[1]
    constant <- gene_results_filtered()$coef_constant[1]
    p_pathogenicity <- paste("p =", formatC(gene_results_filtered()$p_log_pathogenicity, format = "e", digits = 3), ")")
    
    # Calculate the y-values for the best fit line
    line_y <- ~coef_pathogenicity * log_pathogenicity + constant

    # Plot data with best fit line
    plot_ly(filtered_data(), x = ~log_pathogenicity, y = ~BETA, type = 'scatter', mode = 'markers',
            error_y = list(
                array = ~SE,
                color = "gray",
                thickness = 0.5
            ),
            text = ~paste("Variant ID:", markerID, "<br> Effect size:", BETA),
            hoverinfo = "text"
    ) %>%
       add_trace(
            x = ~log_pathogenicity,
            y = line_y,
            type = 'scatter',
            mode = 'lines',
            line = list(color = 'black'),
            name = 'Best Fit Line',
            showlegend = FALSE
       ) %>%
       layout(
            title = list(
              text = paste("Pathogenicity Analysis<br><sub>", p_pathogenicity, "</sub>"),
              x = 0.5
            ),
            xaxis = list(title = "-log(1 - Pathogenicity Probability)"), 
            yaxis = list(title = "Effect Size on Phenotype")
       )
})


   
    # LoF Analysis Plot
    output$lofPlot <- renderPlotly({
    req(gene_results_filtered())

    # Pull p-value for LoF
    p_lof <- paste("(p =", formatC(gene_results_filtered()$p_pLoF_indicator, format = "e", digits = 3), ")")

    plot_ly(filtered_data(), x = ~factor(pLoF_indicator), y = ~BETA, type = 'box') %>%
        layout(
            title = list(
              text = paste("Loss-of-Function Analysis<br><sub>", p_lof, "</sub>"),
              x = 0.5
            ),
            xaxis = list(title = "Loss-of-Function Indicator"), 
            yaxis = list(title = "Effect Size on Phenotype")
        )
})

# Missense Indicator Plot
output$missensePlot <- renderPlotly({
    req(gene_results_filtered())

    # Pull p-value for missense indicator
    p_missense <- paste("(p =", formatC(gene_results_filtered()$p_missense_indicator, format = "e", digits = 3), ")")

    plot_ly(filtered_data(), x = ~factor(missense_indicator), y = ~BETA, type = 'box') %>%
        layout(
            title = list(
              text = paste("Missense Analysis<br><sub>", p_missense, "</sub>"),
              x = 0.5
            ),
            xaxis = list(title = "Missense Indicator"), 
            yaxis = list(title = "Effect Size on Phenotype")
        )
})
}

# Run the app
shinyApp(ui = ui, server = server)
