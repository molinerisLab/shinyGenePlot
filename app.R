# app.R
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(DT)
library(sortable)  # For sortable list

ui <- fluidPage(
  titlePanel("Gene Expression Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      # File inputs
      fileInput("expressionFile", "Upload Gene Expression Matrix (CSV/TSV):",
                accept = c("text/csv", "text/comma-separated-values", "text/tab-separated-values", 
                           ".csv", ".tsv", ".gz", "application/gzip")),
      fileInput("metadataFile", "Upload Sample Metadata (CSV/TSV):",
                accept = c("text/csv", "text/comma-separated-values", "text/tab-separated-values", 
                           ".csv", ".tsv", ".txt", "text/plain")),
      
      # Only show these UI elements after files are loaded
      conditionalPanel(
        condition = "output.filesLoaded",
        
        # Select genes
        selectizeInput("selectedGenes", "Select Genes:", choices = NULL, multiple = TRUE),
        
        # Select grouping variable from metadata
        selectInput("groupVar", "Group Samples by:", choices = NULL),
        
        # Select color variable from metadata
        selectInput("colorVar", "Color by:", choices = NULL),
        
        # Reorder factor levels
        uiOutput("factorLevelsUI"),
        
        # Update button
        actionButton("updatePlot", "Update Plot", class = "btn-primary")
      )
    ),
    
    mainPanel(
      # Show content in tabs
      tabsetPanel(
        # Gene expression plot tab
        tabPanel("Plot", 
                 plotOutput("expressionPlot", height = "700px")
        ),
        # Data preview tabs
        tabPanel("Expression Data", DTOutput("expressionPreview")),
        tabPanel("Metadata", DTOutput("metadataPreview"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive values to store data
  data <- reactiveValues(
    expression = NULL,
    metadata = NULL,
    plotData = NULL
  )
  
  # Function to detect file delimiter
  detectDelimiter <- function(file) {
    first_line <- readLines(file, n = 1)
    if (grepl("\t", first_line)) {
      return("\t")
    } else {
      return(",")
    }
  }
  
  # Load expression data
  observeEvent(input$expressionFile, {
    req(input$expressionFile)
    
    # Check if file is gzipped
    is_gzipped <- grepl("\\.gz$", input$expressionFile$name, ignore.case = TRUE)
    
    if (is_gzipped) {
      # Read gzipped file
      con <- gzfile(input$expressionFile$datapath, "rt")
      first_line <- readLines(con, n = 1)
      close(con)
      
      # Determine delimiter
      if (grepl("\t", first_line)) {
        delimiter <- "\t"
      } else {
        delimiter <- ","
      }
      
      # Re-open and read full file
      con <- gzfile(input$expressionFile$datapath, "rt")
      df <- read.delim(con, sep = delimiter, check.names = FALSE)
      close(con)
    } else {
      # Detect delimiter for non-gzipped file
      delimiter <- detectDelimiter(input$expressionFile$datapath)
      
      # Read expression file (assuming first column contains gene names)
      df <- read.delim(input$expressionFile$datapath, sep = delimiter, check.names = FALSE)
    }
    
    # Store gene expression data
    data$expression <- df
    
    # Update gene selection choices
    updateSelectizeInput(session, "selectedGenes", choices = df[[1]], server = TRUE)
  })
  
  # Load metadata
  observeEvent(input$metadataFile, {
    req(input$metadataFile)
    
    # Detect delimiter
    delimiter <- detectDelimiter(input$metadataFile$datapath)
    
    # Read metadata file
    metadata <- read.delim(input$metadataFile$datapath, sep = delimiter, check.names = FALSE)
    
    # Store metadata
    data$metadata <- metadata
    
    # Update grouping/coloring variable choices (excluding first column which should be sample ID)
    choices <- colnames(metadata)[-1]
    updateSelectInput(session, "groupVar", choices = choices)
    updateSelectInput(session, "colorVar", choices = choices)
  })
  
  # Signal that files are loaded
  output$filesLoaded <- reactive({
    return(!is.null(data$expression) && !is.null(data$metadata))
  })
  outputOptions(output, "filesLoaded", suspendWhenHidden = FALSE)
  
  # Generate UI for reordering factor levels
  output$factorLevelsUI <- renderUI({
    req(input$groupVar, data$metadata)
    
    # Get unique levels for the selected grouping variable
    group_levels <- unique(data$metadata[[input$groupVar]])
    
    # Create a true drag-and-drop sortable list
    tagList(
      tags$label("Drag to Reorder Group Levels:"),
      tags$div(
        style = "margin-bottom: 15px;",
        "Drag items to reorder them. The order here will be used in the plot."
      ),
      rank_list(
        input_id = "factorOrder",
        labels = group_levels,
        text = "",
        css_id = "factor-rank-list",
        class = "default-sortable",
        orientation = "vertical"
      )
    )
  })
  
  # Prepare data for plotting
  observeEvent(input$updatePlot, {
    req(data$expression, data$metadata, input$selectedGenes, input$groupVar, input$colorVar, input$factorOrder)
    
    # Get gene names column (first column)
    gene_col <- colnames(data$expression)[1]
    
    # Get sample IDs from metadata (first column)
    sample_col <- colnames(data$metadata)[1]
    
    sC<-colnames(data$expression)[2:length(data$expression)]
    sM<-data$metadata[,1]

    #test if there are missing samples in the metadata
    missing_samples <- setdiff(sC,sM)
    
    
    # Filter expression data to selected genes
    gene_data <- data$expression %>%
      filter(!!sym(gene_col) %in% input$selectedGenes)
    
    # Reshape from wide to long format
    long_data <- gene_data %>%
      pivot_longer(cols = -all_of(gene_col), 
                   names_to = "sample_id", 
                   values_to = "expression")
    
    
    
    # Trim whitespace from sample IDs to avoid join issues
    long_data$sample_id <- trimws(long_data$sample_id)
    data$metadata[[sample_col]] <- trimws(data$metadata[[sample_col]])
    
    # Join with metadata (fixed join syntax)
    plot_data <- long_data %>%
      left_join(data$metadata, by = c("sample_id" = sample_col))
    
    if(length(missing_samples) > 0) {
      warning_msg <- paste("Warning: Some samples in expression data are not in metadata and are removed from the plot:",
                           paste(missing_samples, collapse=", "))
      showNotification(warning_msg, type = "warning", duration = 10)
      
      # Remove rows with NA grouping variable to avoid NA in plot
      plot_data <- plot_data %>% filter(!is.na(!!sym(input$groupVar)))
    }
    
    # Important: Ensure ALL levels are properly ordered according to the sortable list
    # First convert to character to avoid factor level issues
    plot_data[[input$groupVar]] <- as.character(plot_data[[input$groupVar]])
    
    # Then create a new factor with explicit levels from the sortable list
    # This ensures the exact order from the UI is respected
    factor_levels <- input$factorOrder
    
    # Check if there are any levels in the data that aren't in the sortable
    all_levels_in_data <- unique(plot_data[[input$groupVar]])
    missing_levels <- setdiff(all_levels_in_data, factor_levels)
    
    # If there are missing levels, add them to the end of the factor levels
    if(length(missing_levels) > 0) {
      factor_levels <- c(factor_levels, missing_levels)
    }
    
    # Now set the factor with the complete ordered levels
    plot_data[[input$groupVar]] <- factor(plot_data[[input$groupVar]], 
                                          levels = factor_levels)
    
    # Store plot data
    data$plotData <- plot_data
  })
  
  # Generate expression plot
  output$expressionPlot <- renderPlot({
    req(data$plotData, input$groupVar, input$colorVar)
    
    # Get gene column name
    gene_col <- colnames(data$expression)[1]
    
    # Create the plot
    ggplot(data$plotData, aes_string(x = input$groupVar, 
                                     y = "expression",
                                     group = input$groupVar,
                                     color = input$colorVar)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(height = 0) +
      facet_wrap(~ .data[[gene_col]], scales = "free_y") +
      theme_bw() +
      labs(
        title = "Gene Expression by Group",
        x = input$groupVar,
        y = "Normalized Expression"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(face = "bold")
      )
  })
  
  # Data previews
  output$expressionPreview <- renderDT({
    req(data$expression)
    datatable(head(data$expression, 10), options = list(scrollX = TRUE))
  })
  
  output$metadataPreview <- renderDT({
    req(data$metadata)
    datatable(data$metadata, options = list(scrollX = TRUE))
  })
}

shinyApp(ui = ui, server = server)