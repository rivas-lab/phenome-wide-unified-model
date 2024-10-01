# Load necessary libraries
library(shiny)
library(DT)
library(ggplot2)
library(reshape2)
library(plotly)  # for interactive plots

# Define UI for application
ui <- fluidPage(
  titlePanel("Unified Meta Regression Model Results 400k Exomes"),
  
  sidebarLayout(
    sidebarPanel(
      # Search inputs
      textInput("search_gene", "Search Gene:", value = ""),
      textInput("search_phenocode", "Search Phenocode:", value = ""),
      textInput("search_description", "Search Description:", value = ""),
      
      # Alpha filters for p-values
      numericInput("alpha_p_model", "Set α for p_value unified model:", value = 1e-6, min = 0, max = 1),
      numericInput("alpha_p_log_constraint", "Set α for p_value log constraint:", value = 1, min = 0, max = 1),
      numericInput("alpha_p_log_pathogenicity", "Set α for p_value log pathogenicity:", value = 1, min = 0, max = 1),
      numericInput("alpha_p_pLoF_indicator", "Set α for p_value pLoF indicator:", value =1, min = 0, max = 1),
      numericInput("alpha_p_missense_indicator", "Set α for p_value missense indicator:", value = 1, min = 0, max = 1),
      width = 3
          ),
    mainPanel(
      fluidRow(
        column(12,
               plotlyOutput("coefficientsPlot", height = "300px")  # plotly for interactivity
        ),
        column(12,
               plotlyOutput("pvaluesPlot", height = "300px")  # plotly for interactivity
        )
      ),
      fluidRow(
        DTOutput("dataTable", width = "100%")
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  # Read and process the data
  data <- reactive({
    df <- read.csv("filtered_results.csv", stringsAsFactors = FALSE)
    
    # Apply search filters
    if (input$search_gene != "") {
      df <- df[grep(input$search_gene, df$gene, ignore.case = TRUE), ]
    }
    if (input$search_phenocode != "") {
      df <- df[grep(input$search_phenocode, df$phenocode, ignore.case = TRUE), ]
    }
    if (input$search_description != "") {
      df <- df[grep(input$search_description, df$description, ignore.case = TRUE), ]
    }
    
    # Apply p-value filters
    df <- df[
      df$p_model < input$alpha_p_model &
        df$p_log_constraint < input$alpha_p_log_constraint &
        df$p_log_pathogenicity < input$alpha_p_log_pathogenicity &
        df$p_pLoF_indicator < input$alpha_p_pLoF_indicator &
        df$p_missense_indicator < input$alpha_p_missense_indicator, 
    ]
    
    return(df)
  })
  
  # Coefficients Plot (Interactive with Tooltips)
  output$coefficientsPlot <- renderPlotly({
    df <- data()
    if (nrow(df) == 0) {
      return(NULL)
    }
    # Melt the dataframe for coefficients
    coef_cols <- c("coef_log_constraint", "coef_log_pathogenicity", "coef_pLoF_indicator", "coef_missense_indicator", "coef_constant")
    df_coef <- melt(df[, c("gene", coef_cols)], id.vars = "gene", variable.name = "Coefficient", value.name = "Value")
    
    # Create a scatter plot for coefficients with plotly for interactivity
    p <- ggplot(df_coef, aes(x = Value, y = Coefficient, text = paste("Gene:", gene, "<br>Value:", round(Value, 4)))) +
      geom_point(color = "#8C1515") +
      theme_minimal() +
      labs(title = "Plot of Coefficients", x = "Coefficients", y = "Coefficient Value")
    
    # Convert ggplot to plotly for interactive tooltips
    ggplotly(p, tooltip = "text")
  })
  
  # P-values Plot (Interactive with Tooltips)
  output$pvaluesPlot <- renderPlotly({
    df <- data()
    if (nrow(df) == 0) {
      return(NULL)
    }
    # Melt the dataframe for p-values
    pvalue_cols <- c("p_model", "p_log_constraint", "p_log_pathogenicity", "p_pLoF_indicator", "p_missense_indicator", "p_constant")
    df_pvalues <- melt(df[, c("gene", pvalue_cols)], id.vars = "gene", variable.name = "PValue", value.name = "Value")
    
    # Transform p-values to -log10(p-value)
    df_pvalues$NegLogPValue <- -log10(df_pvalues$Value)
    
    # Create a scatter plot for p-values with plotly for interactivity
    p <- ggplot(df_pvalues, aes(x = NegLogPValue, y = PValue, text = paste("Gene:", gene, "<br>-log10(P):", round(NegLogPValue, 4)))) +
      geom_point(color = "#4D4F53") +
      theme_minimal() +
      labs(title = "Plot of P-Values", x = "P-Values", y = "-log10(P-Value)")
    
    # Convert ggplot to plotly for interactive tooltips
    ggplotly(p, tooltip = "text")
  })
  
  # Data Table (with formatting)
  output$dataTable <- renderDT({
    df <- data()
    datatable(df,
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                searchHighlight = TRUE,
                columnDefs = list(list(className = 'dt-center', targets = "_all")) # Center align all columns
              )
    ) %>%
      # Format numeric columns in scientific notation with 4 significant figures
      formatSignif(columns = c("p_model", "p_log_constraint", "p_log_pathogenicity", "p_pLoF_indicator", "p_missense_indicator","p_constant"), digits = 4) %>%
      formatSignif(columns = c("coef_log_constraint", "coef_log_pathogenicity", "coef_pLoF_indicator", "coef_missense_indicator", "coef_constant"), digits = 4) %>%
      formatStyle(
        columns = c("p_model", "p_log_constraint", "p_log_pathogenicity", "p_pLoF_indicator", "p_missense_indicator"),
        color = styleInterval(c(1e-6, 0.05), c('green', 'blue', 'black'))  # Apply color formatting based on value
      ) %>%
      formatStyle(
        'gene',
        fontWeight = 'bold'  # Make the 'gene' column bold
      )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
