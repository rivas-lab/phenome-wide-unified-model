library(shiny)
library(plotly)
library(readr)
library(tidyr)
library(stringr)
library(dplyr)
library(webshot2)
library(htmlwidgets)

# Define the file paths
file_path <- "gene_phenotype/CASZ1_102.tsv"
pheno_file <- "phenotypes/102_unifiedmodel_MAF.05_updated.tsv.gz"

# Extract the base name without the extension
base_name <- tools::file_path_sans_ext(basename(file_path))
gene_name <- str_split(base_name, "_")[[1]][1]
pheno_name <- str_split(base_name, "_")[[1]][2]

# Load the reference data
reference_file <- "precomputing_files_for_plots/pheno_results.txt.bgz"
reference_df <- read.table(gzfile(reference_file), sep = "\t", header = TRUE, fill = TRUE)


# Map the input gene to the phenotype description
description <- reference_df$description[reference_df$phenocode == pheno_name]

# Check if the description is empty
if (length(description) == 0) {
    pheno_description <- pheno_name
} else {
    pheno_description <- description[1]
}

# Load data
data <- read_tsv(file_path, show_col_types = FALSE)

# Get the P-values for this gene from the model
model_data <- read.delim(gzfile(pheno_file), header = TRUE, sep = "\t")
model_data <- model_data %>%
  filter(gene == gene_name)
# Extract values and format in scientific notation with four significant figures
p_constraint <- formatC(model_data$p_log_constraint[1], format = "e", digits = 3)
p_pathogenicity <- formatC(model_data$p_log_pathogenicity[1], format = "e", digits = 3)
p_pLoF <- formatC(model_data$p_pLoF_indicator[1], format = "e", digits = 3)
p_missense <- formatC(model_data$p_missense_indicator[1], format = "e", digits = 3)

coef_constraint <- formatC(model_data$coef_log_constraint, format = "e", digits = 3)
coef_pathogenicity <- formatC(model_data$coef_log_pathogenicity, format = "e", digits = 3)
coef_pLoF_indicator <- formatC(model_data$coef_pLoF_indicator, format = "e", digits = 3)
coef_missense_indicator <- formatC(model_data$coef_missense_indicator, format = "e", digits = 3)

constant <- model_data$coef_constant[1]

# Define UI
ui <- fluidPage(
    titlePanel(
        HTML(
            paste(
                "<h1>Gene-Phenotype Analysis</h1>",
                "<h3>Gene: ", gene_name, 
                " | Phenotype: ", '"', pheno_description, '"', 
                "</h3>",
                "<br><br>"
            )
        )
    ),
    
    mainPanel(
        plotlyOutput("constraintPlot"),
        tags$br(), tags$br(),
        plotlyOutput("pathogenicityPlot"),
        tags$br(), tags$br(),
        plotlyOutput("lofPlot"),
        tags$br(), tags$br(),
        plotlyOutput("missensePlot"),
        tags$br(), tags$br(), tags$br(), tags$br()
    )
)

# Define server logic
server <- function(input, output) {

# Constraint Plot
    # crop outliers out of plot
    plotdata <- data %>%
        mutate(log_constraint = ifelse(log_constraint > 1.5, 1.5, log_constraint)) 

    output$constraintPlot <- renderPlotly({
        plot_ly(plotdata, x = ~log_constraint, y = ~BETA, type = 'scatter', mode = 'markers',
                error_y = list(
                    array = ~SE, 
                    color = "gray",  
                    thickness = 0.5
                ), text = ~paste("Variant ID:", markerID, "<br> Effect size:", BETA, "±", SE),
                hoverinfo = "text"
        ) %>%
        add_trace(
            x = ~log_constraint,
            y = ~ model_data$coef_log_constraint[1] * log_constraint + constant,
            type = 'scatter',
            mode = 'lines',
            line = list(color = 'black'),
            name = 'Best Fit Line',
            error_y = NULL,
            showlegend = FALSE
        ) %>%
        layout(title = list(
            text = paste(gene_name,"Constraint Analysis<br><span style='font-size: 14px;'>(p-value:", p_constraint, "| Coefficient:", coef_constraint, ")</span>")), 
            xaxis = list(title = "log(1 - Constraint Probability)"), 
            yaxis = list(title = "Effect Size")
        )
    })

    # Alpha Missense Analysis Plot
    output$pathogenicityPlot <- renderPlotly({
        plot_ly(data, x = ~log_pathogenicity, y = ~BETA, type = 'scatter', mode = 'markers',
                error_y = list(
                    array = ~SE, 
                    color = "gray", 
                    thickness = 0.5
                ), text = ~paste("Variant ID:", markerID, "<br> Effect size:", BETA, "±", SE),
                hoverinfo = "text"
        ) %>%
        # Overlay the best fit line
        add_trace(
            x = ~log_pathogenicity,
            y = ~model_data$coef_log_pathogenicity[1] * log_pathogenicity + constant,
            type = 'scatter',
            mode = 'lines',
            line = list(color = 'black'),
            name = 'Best Fit Line',
            error_y = NULL,
            showlegend = FALSE
        ) %>%
        layout(title = list(
            text = paste(gene_name, "Pathogenicity Analysis<br><span style='font-size: 14px;'>(p-value:", p_pathogenicity, "| Coefficient:", coef_pathogenicity, ")</span>")), 
            xaxis = list(title = "log(1 - Pathogenicity Probability)"), 
            yaxis = list(title = "Effect Size")
        )
    })

    # LoF Analysis Plot
    output$lofPlot <- renderPlotly({
        plot_ly(data, x = ~factor(pLoF_indicator), y = ~BETA, type = 'box') %>%
        layout(title = paste(gene_name, "Loss-of-Function Analysis<br><span style='font-size: 14px;'>(p-value:", p_pLoF, "| Coefficient:", coef_pLoF_indicator, ")</span>"), 
               xaxis = list(title = "LoF Indicator"), 
               yaxis = list(title = "Effect Size")
        )
    })

    # Missense Indicator Plot
    output$missensePlot <- renderPlotly({
        plot_ly(data, x = ~factor(missense_indicator), y = ~BETA, type = 'box') %>%
        layout(title = paste(gene_name, "Missense Analysis<br><span style='font-size: 14px;'>(p-value:", p_missense, "| Coefficient:", coef_missense_indicator, ")</span>"), 
               xaxis = list(title = "Missense Indicator"), 
               yaxis = list(title = "Effect Size")
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
