library(shiny)
library(shinyFiles)
library(shinythemes)
library(glue)
library(tidyverse)
library(seqinr)

# Get the directory where this GUI script is located
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmd_args)
  if (length(match) > 0) {
    return(normalizePath(sub(needle, "", cmd_args[match])))
  } else {
    # If not run via Rscript, try to find the path another way
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

script_dir <- dirname(get_script_dir())
setwd(script_dir)

# Source the required R scripts
source("process_unaligned.R")
source("test_normality.R")
source("test_quantResponse_categPredictor.R")
source("process_aligned.R")
source("benchmarking.R")

ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Sequence Analysis Pipeline GUI"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("File Inputs"),
      fileInput("fasta_file", "FASTA File", accept = c(".fasta", ".fa", ".phy")),
      fileInput("metadata_file", "Metadata File", accept = c(".tsv", ".txt", ".csv")),
      textInput("outdir", "Output Directory", value = script_dir),
      radioButtons("sequence_type", "Sequence Type:", 
                   choices = c("Unaligned" = "unaligned", "Aligned" = "aligned"),
                   selected = "unaligned"),
      actionButton("process_btn", "Process Data", class = "btn-primary"),
      br(), br(),
      conditionalPanel(
        condition = "output.processing_done == true",
        checkboxInput("subset_fasta", "Subset input multifasta to categorized sequences"),
        conditionalPanel(
          condition = "input.subset_fasta == true",
          actionButton("subset_btn", "Run Subsetting", class = "btn-success")
        )
      )
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        
        # Tab 1: Data Import
        tabPanel("Data Import",
                 h3("Data Processing Results"),
                 verbatimTextOutput("process_output"),
                 
                 conditionalPanel(
                   condition = "output.processing_done == true",
                   h3("Benchmarking Results"),
                   verbatimTextOutput("benchmarking_output"),
                   h3("Normality Test Results"),
                   verbatimTextOutput("normality_output"),
                   h3("Category Summary"),
                   tableOutput("category_table"),
                   h3("QQ Plot"),
                   imageOutput("qqplot")
                 )
        ),
        
        # Tab 2: Figures
        tabPanel("Figures",
                 fluidRow(
                   column(4,
                          h4("Figure Generation"),
                          selectInput("figure_script", "Select Figure Script:",
                                      choices = c("make_histograms.R" = "make_histograms.R",
                                                  "make_scatterplots.R" = "make_scatterplots.R",
                                                  "make_heatmaps.R" = "make_heatmaps.R")),
                          conditionalPanel(
                            condition = "input.figure_script == 'make_histograms.R'",
                            selectInput("hist_group_var", "Grouping Variable:",
                                        choices = c("category", "genus"),
                                        selected = "category")
                          ),
                          actionButton("generate_fig_btn", "Generate Figure", class = "btn-primary")
                   ),
                   column(8,
                          h4("Generated Figures"),
                          uiOutput("figure_browser")
                   )
                 )
        ),
        
        # Tab 3: Statistical Tests
        tabPanel("Statistical Tests",
                 fluidRow(
                   column(4,
                          h4("Statistical Analysis"),
                          selectInput("stats_script", "Select Test Script:",
                                      choices = c("test_quantResponse_categPredictor.R" = "test_quantResponse_categPredictor.R",
                                                  "test_categResponse_categPredictor.R" = "test_categResponse_categPredictor.R")),
                          
                          conditionalPanel(
                            condition = "input.stats_script == 'test_quantResponse_categPredictor.R'",
                            selectInput("response_var", "Response Variable:",
                                        choices = c("sequence_length", "pident", "evalue"),
                                        selected = "sequence_length"),
                            selectInput("predictor_var", "Predictor Variable:",
                                        choices = c("category", "genus"),
                                        selected = "category"),
                            selectInput("test_type", "Test Type:",
                                        choices = c("parametric", "non-parametric"),
                                        selected = "non-parametric"),
                            selectInput("adjust_method", "P-value Adjustment:",
                                        choices = c("holm", "hochberg", "hommel", "bonferroni", 
                                                   "BH", "BY", "fdr", "none"),
                                        selected = "fdr")
                          ),
                          actionButton("run_stats_btn", "Run Statistical Test", class = "btn-primary")
                   ),
                   column(8,
                          h4("Test Results"),
                          verbatimTextOutput("stats_output")
                   )
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive values to store processed data and state
  values <- reactiveValues(
    df = NULL,
    processing_done = FALSE,
    figures_dir = NULL,
    benchmarking_result = NULL,
    process_output_text = NULL
  )
  
  # Process data when button is clicked
  observeEvent(input$process_btn, {
    req(input$fasta_file, input$metadata_file, input$outdir)
    
    tryCatch({
      showNotification("Processing data...", type = "message")
      
      # Run benchmarking first
      if (file.exists(input$metadata_file$datapath)) {
        values$benchmarking_result <- run_benchmarking(input$metadata_file$datapath)
      }
      
      # Create output directory
      dir.create(input$outdir, showWarnings = FALSE, recursive = TRUE)
      
      # Run the appropriate processing script based on sequence type
      if (input$sequence_type == "unaligned") {
        # Use the new process_unaligned_shiny function
        result <- process_unaligned_shiny(
          input$fasta_file$datapath,
          input$metadata_file$datapath,
          input$outdir
        )
        
        values$df <- result
        values$process_output_text <- capture.output({
          cat("Unaligned sequence processing completed successfully.\n")
          cat(glue("Total sequences processed: {nrow(result)}\n"))
          if ("category" %in% colnames(result)) {
            cat("Category distribution:\n")
            print(table(result$category))
          }
        })
        
      } else {
        # Run process_aligned.R (placeholder for now)
        result <- process_aligned_pipeline(
          input$fasta_file$datapath,
          input$metadata_file$datapath,
          input$outdir
        )
        values$df <- result$df
        values$process_output_text <- result$output
      }
      
      values$processing_done <- TRUE
      values$figures_dir <- file.path(input$outdir, "figures")
      
      # Run normality test if we have data
      if (!is.null(values$df) && nrow(values$df) > 0 && "sequence_length" %in% colnames(values$df)) {
        fig_dir <- file.path(input$outdir, "figures")
        dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
        values$normality_test <- test_normality(values$df, "sequence_length", fig_dir)
      }
      
      showNotification("Data processing completed!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Render process output
  output$process_output <- renderPrint({
    req(values$process_output_text)
    cat(values$process_output_text)
  })
  
  # Render benchmarking output
  output$benchmarking_output <- renderPrint({
    req(values$benchmarking_result)
    values$benchmarking_result
  })
  
  # Render normality test results
  output$normality_output <- renderPrint({
    req(values$normality_test)
    values$normality_test
  })
  
  # Render category table
  output$category_table <- renderTable({
    req(values$df, values$processing_done)
    if ("category" %in% colnames(values$df)) {
      values$df %>%
        count(category) %>%
        rename(Category = category, Count = n) %>%
        arrange(desc(Count))
    }
  })
  
  # Display QQ plot if it exists
  output$qqplot <- renderImage({
    req(values$processing_done, input$outdir)
    
    qqplot_file <- file.path(input$outdir, "figures", "qqplot_sequence_length.pdf")
    if (file.exists(qqplot_file)) {
      list(src = qqplot_file,
           contentType = 'application/pdf',
           width = 400,
           height = 300,
           alt = "QQ Plot")
    } else {
      # Return a blank image if file doesn't exist
      list(src = "", width = 0, height = 0)
    }
  }, deleteFile = FALSE)
  
  # Subset FASTA when button is clicked
  observeEvent(input$subset_btn, {
    req(values$processing_done, input$fasta_file, input$outdir)
    
    tryCatch({
      multifasta <- input$fasta_file$datapath
      ids_file <- file.path(input$outdir, "categorized_ids.txt")
      outfile <- file.path(input$outdir, "multifasta_categorized.fasta")
      
      if (file.exists(ids_file)) {
        system2("seqtk", args = c("subseq", multifasta, ids_file), stdout = outfile)
        showNotification(glue("FASTA subsetting completed! Output: {outfile}"), type = "message")
      } else {
        showNotification("categorized_ids.txt not found. Please process data first.", type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error in subsetting:", e$message), type = "error")
    })
  })
  
  # Generate figures
  observeEvent(input$generate_fig_btn, {
    req(values$df, input$outdir)
    
    tryCatch({
      figures_dir <- file.path(input$outdir, "figures")
      dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Source the selected figure script
      figure_script_path <- file.path("figures", input$figure_script)
      if (file.exists(figure_script_path)) {
        source(figure_script_path)
        
        # Call the appropriate function based on the script
        if (input$figure_script == "make_histograms.R") {
          if (!is.null(input$hist_group_var)) {
            make_histograms(values$df, figures_dir, group_var = input$hist_group_var)
          } else {
            make_histograms(values$df, figures_dir)
          }
        } else if (input$figure_script == "make_scatterplots.R") {
          make_scatterplots(values$df, figures_dir)
        } else if (input$figure_script == "make_heatmaps.R") {
          make_heatmaps(values$df, figures_dir)
        }
        
        showNotification("Figure generation completed!", type = "message")
        
        # Refresh the figure browser
        output$figure_browser <- renderUI({
          figure_files <- list.files(figures_dir, pattern = "\\.(png|jpg|jpeg|pdf)$", 
                                    full.names = TRUE, recursive = TRUE)
          
          if (length(figure_files) > 0) {
            tagList(
              selectInput("selected_figure", "Select Figure to Display:",
                         choices = setNames(figure_files, basename(figure_files))),
              imageOutput("displayed_figure", height = "400px")
            )
          } else {
            p("No figures found in the output directory.")
          }
        })
      } else {
        showNotification("Figure script not found!", type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error in figure generation:", e$message), type = "error")
    })
  })
  
  # Display selected figure
  output$displayed_figure <- renderImage({
    req(input$selected_figure)
    
    figure_path <- input$selected_figure
    if (file.exists(figure_path)) {
      # Determine content type based on file extension
      ext <- tolower(tools::file_ext(figure_path))
      content_type <- switch(ext,
                            "png" = "image/png",
                            "jpg" = "image/jpeg",
                            "jpeg" = "image/jpeg",
                            "pdf" = "application/pdf",
                            "image/png") # default
      
      list(src = figure_path,
           contentType = content_type,
           width = "100%",
           height = "auto",
           alt = basename(figure_path))
    }
  }, deleteFile = FALSE)
  
  # Run statistical tests
  observeEvent(input$run_stats_btn, {
    req(values$df)
    
    tryCatch({
      if (input$stats_script == "test_quantResponse_categPredictor.R") {
        # Run the quantitative response test
        result <- quantResponse_categPredictor(
          df = values$df,
          response = input$response_var,
          predictor = input$predictor_var,
          test_type = input$test_type,
          adjust = input$adjust_method
        )
        
        output$stats_output <- renderPrint({
          result
        })
      }
      
      showNotification("Statistical test completed!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error in statistical test:", e$message), type = "error")
    })
  })
  
  # Output for conditional panels
  output$processing_done <- reactive({
    values$processing_done
  })
  outputOptions(output, "processing_done", suspendWhenHidden = FALSE)
  
  # Wrapper function for process_aligned.R (unchanged)
  process_aligned_pipeline <- function(alignment_file, metadata_file, outdir) {
    # Create a temporary environment to run the script in
    temp_env <- new.env()
    
    # Capture output
    output_text <- capture.output({
      # Source the main processing in the temporary environment
      source("process_aligned.R", local = temp_env)
      
      # Call the function (you might need to adjust parameters based on your process_aligned.R)
      if (exists("process_alignment_and_metadata", envir = temp_env)) {
        alignment_type <- ifelse(grepl("\\.phy$", alignment_file), "phylip", "fasta")
        temp_env$df <- temp_env$process_alignment_and_metadata(
          alignment_file = alignment_file,
          alignment_type = alignment_type,
          metadata_file = metadata_file,
          sequence_name_col = "locus_tags"
        )
      }
    }, type = "message")
    
    return(list(
      df = if(exists("df", envir = temp_env)) temp_env$df else data.frame(),
      output = paste(output_text, collapse = "\n"),
      normality_test = "non-parametric"  # Placeholder - you might want to run normality test here too
    ))
  }
  
  # Wrapper function for benchmarking.R
  run_benchmarking <- function(metadata_file) {
    if (!file.exists("benchmarking.R")) {
      return("benchmarking.R not found in script directory")
    }
    
    temp_env <- new.env()
    
    output_text <- capture.output({
      source("benchmarking.R", local = temp_env)
      
      # Assuming benchmarking.R has a main function or can be run directly
      # You might need to adjust this based on your benchmarking.R structure
      if (exists("run_benchmarking", envir = temp_env)) {
        temp_env$run_benchmarking(metadata_file)
      }
    }, type = "message")
    
    return(paste(output_text, collapse = "\n"))
  }
}

shinyApp(ui = ui, server = server)