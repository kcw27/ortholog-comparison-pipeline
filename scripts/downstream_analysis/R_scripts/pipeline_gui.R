library(shiny)

# Get script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grep(file_arg, args)])
  if (length(script_path) == 0) {
    stop("Cannot determine script path. Are you running via Rscript?")
  }
  normalizePath(dirname(script_path))
}

script_dir <- get_script_dir()

# Source processing scripts
source(file.path(script_dir, "process_unaligned.R"))
source(file.path(script_dir, "process_aligned.R"))
source(file.path(script_dir, "benchmarking.R"))

# Strip quotes from input
strip_quotes <- function(x) gsub('^["\']|["\']$', '', x)

# Count lines in a file
count_lines <- function(path) {
  if (!file.exists(path)) return(NA)
  length(readLines(path, warn = FALSE))
}

ui <- navbarPage("Pipeline GUI",
  tabPanel("Data import",
    fluidPage(
      radioButtons("fasta_type", "FASTA type:",
                   choices = c("Unaligned", "Aligned"),
                   selected = "Unaligned",
                   inline = TRUE),
      textInput("fasta_file", "FASTA file path:", ""),
      textInput("metadata_file", "Metadata file path:", ""),
      textInput("outdir", "Output directory:", ""),
      actionButton("submit_btn", "Submit"),
#      conditionalPanel("output.loading",
#        tags$div(class = "spinner-border text-primary", role = "status",
#                 tags$span(class = "visually-hidden", "Loading..."))
#      ),
      verbatimTextOutput("output_text"),
      
      # QQ plot path and image
      textOutput("qqplot_path"),
      imageOutput("qqplot_image")
    )
  ),
  tabPanel("Figures",
    conditionalPanel("output.tabs_ready",
      fluidPage(
        selectInput("fig_script", "Select figure script:",
                    choices = c("make_histograms.R")),
        conditionalPanel("input.fig_script == 'make_histograms.R'",
          selectInput("group_var", "Group by:",
                      choices = c("category", "genus"),
                      selected = "category")
        ),
        actionButton("run_figures", "Run figure script"),
        uiOutput("image_selector"),
        imageOutput("selected_image")
      )
    )
  ),
  tabPanel("Statistical tests",
    conditionalPanel("output.tabs_ready",
      fluidPage(
        selectInput("test_script", "Select statistical test script:",
                    choices = c("test_quantResponse_categPredictor.R",
                                   "test_alleleFreqs_within_categories.R",
                                   "test_categFreqs_within_quantQuadrants.R",
                                   "test_resampled_distribution_properties.R",
                                   "test_sequence_diversity.R")),
        
        # Inputs for test_quantResponse_categPredictor.R
        conditionalPanel("input.test_script == 'test_quantResponse_categPredictor.R'",
          selectInput("response_col", "Response (numeric):", choices = NULL),
          selectInput("predictor_col", "Predictor:", choices = c("category", "genus")),
          selectInput("adjust_method", "Adjustment method:",
                      choices = c("", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                      selected = "fdr"),
          selectInput("test_type", "Test type:", choices = c("parametric", "non-parametric"))
        ),
        
        # Inputs for other scripts can be added similarly
        # Example: test_alleleFreqs_within_categories.R
        conditionalPanel("input.test_script == 'test_alleleFreqs_within_categories.R'",
          selectInput("allele_predictor", "Predictor:", choices = c("category", "genus"))
        ),
        
        actionButton("run_test", "Run statistical test"),
        verbatimTextOutput("test_output"),
      )
    )
  )

)


server <- function(input, output, session) {
  log_store <- reactiveValues(lines = character())
  df_store <- reactiveVal(NULL)
  tabs_ready <- reactiveVal(FALSE)
  test_type_store <- reactiveVal(NULL)
  show_output <- reactiveVal(FALSE)
  
  #loading <- reactiveVal(FALSE)
  #output$loading <- reactive({ loading() })
  #outputOptions(output, "loading", suspendWhenHidden = FALSE)


  append_log <- function(msg) {
    log_store$lines <- c(log_store$lines, msg)
  }

  output$tabs_ready <- reactive({ tabs_ready() })
  outputOptions(output, "tabs_ready", suspendWhenHidden = FALSE)
  
  output$show_output <- reactive({ show_output() })
  outputOptions(output, "show_output", suspendWhenHidden = FALSE)

  observeEvent(input$submit_btn, {
    #loading(TRUE)
    show_output(TRUE)
    
    # Sanitize inputs
    fasta_file <- normalizePath(as.character(strip_quotes(input$fasta_file)), mustWork = FALSE)
    metadata_file <- normalizePath(as.character(strip_quotes(input$metadata_file)), mustWork = FALSE)
    outdir <- normalizePath(as.character(strip_quotes(input$outdir)), mustWork = FALSE)

    #append_log(glue("FASTA file has {count_lines(fasta_file)} lines"))
    #append_log(glue("Metadata file has {count_lines(metadata_file)} lines"))
    benchmark_output <- capture.output(benchmark(metadata_file))
    lapply(benchmark_output, append_log)

    append_log("\nProcessing input data...")
    result <- if (input$fasta_type == "Unaligned") {
      process_unaligned_shiny(fasta_file, metadata_file, outdir, script_dir, log_fn = append_log)
    } else {
      process_aligned_shiny(fasta_file, metadata_file, outdir, script_dir, log_fn = append_log)
    }


    df_store(result$df)
    test_type_store(result$test_type)
    
    numeric_cols <- names(result$df)[sapply(result$df, is.numeric)]
    updateSelectInput(session, "response_col",
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])

    tabs_ready(TRUE)

    categorized_path <- file.path(outdir, "categorized_ids.txt")
    #append_log(glue("\ncategorized_ids.txt has {count_lines(categorized_path)} lines"))
    #loading(FALSE)
    

    # Run normality test and capture output
    test_output <- capture.output({
      invisible(test_type <- test_normality(result$df, "sequence_length", file.path(outdir, "figures")))
    })
    lapply(test_output, append_log)
    test_type_store(test_type)
    
    # Show QQ plot image and path
    qqplot_path <- file.path(outdir, "figures", "qqplot_sequence_length.png")
    
    output$qqplot_path <- renderText({
      if (file.exists(qqplot_path)) {
        glue("QQ plot saved to: {qqplot_path}")
      } else {
        "QQ plot not found."
      }
    })
    
    output$qqplot_image <- renderImage({
      if (file.exists(qqplot_path)) {
        list(src = qqplot_path, contentType = "image/png", alt = "QQ Plot")
      } else {
        NULL
      }
    }, deleteFile = FALSE)

    
      output$output_text <- renderText({
        invalidateLater(500, session)
        paste(log_store$lines, collapse = "\n")
      })
  })

#  output$output_text <- renderText({
#    invalidateLater(500, session)
#    paste(log_store$lines, collapse = "\n")
#  })
  
  observeEvent(input$run_test, {
    df <- df_store()
    if (is.null(df)) return()
    
    script_name <- input$test_script
    if (script_name == "") return()
    
    test_script <- file.path(script_dir, "hypothesis_testing", script_name)
    if (!file.exists(test_script)) {
      output$test_output <- renderText({ glue("Script not found: {script_name}") })
      return()
    }
  
    source(test_script, local = TRUE)
  
    result <- NULL
    if (script_name == "test_quantResponse_categPredictor.R") {
      args <- list(df, input$response_col, input$predictor_col, input$test_type)
      if (input$adjust_method != "") {
        args <- c(args, input$adjust_method)
      }
  
      # Run and capture output
      test_output <- capture.output({
        invisible(result <- do.call(quantResponse_categPredictor, args))
      })
      output$test_output <- renderText({ paste(test_output, collapse = "\n") })
  
      # Render QQ plot image
      qqplot_path <- file.path(input$outdir, "figures", glue("qqplot_{input$response_col}.png"))
      output$qqplot_image <- renderImage({
        list(src = qqplot_path, contentType = "image/png", alt = "QQ Plot")
      }, deleteFile = FALSE)
    }
  })


  # Run figure script
  observeEvent(input$run_figures, {
    df <- df_store()
    if (is.null(df)) return()
    fig_script <- input$fig_script
    group_var <- input$group_var
    script_path <- file.path(script_dir, "figures", fig_script)

    if (file.exists(script_path)) {
      append_log(glue("Running {fig_script}..."))
      if (fig_script == "make_histograms.R") {
        source(script_path, local = TRUE)
        if (group_var != "") {
          make_histograms(df, outdir, group_var)
        } else {
          make_histograms(df, outdir)
        }
        append_log("Histograms generated.")
      }
    }
  })

  # Browse images
  output$image_selector <- renderUI({
    req(tabs_ready())
    img_dir <- file.path(input$outdir, "figures")
    imgs <- list.files(img_dir, pattern = "\\.(png|jpg|jpeg)$", full.names = TRUE)
    div(style = "max-width: 275px; overflow: hidden;",
      selectInput("image_file", "Select image:", choices = imgs)
    )

  })

  output$selected_image <- renderImage({
    req(input$image_file)
    list(src = input$image_file, contentType = "image/png", alt = "Figure")
  }, deleteFile = FALSE)

  observe({
    updateSelectInput(session, "test_type",
                      choices = c("parametric", "non-parametric"),
                      selected = test_type_store())
  })


  observeEvent(input$run_test, {
    df <- df_store()
    if (is.null(df)) return()
    
    script_name <- input$test_script
    if (script_name == "") return()
    
    test_script <- file.path(script_dir, "hypothesis_testing", script_name)
    if (!file.exists(test_script)) {
      output$test_output <- renderText({ glue("Script not found: {script_name}") })
      return()
    }
  
    source(test_script, local = TRUE)
  
    # Dispatch based on selected script
    result <- NULL
    if (script_name == "test_quantResponse_categPredictor.R") {
      args <- list(df, input$response_col, input$predictor_col, input$test_type)
      if (input$adjust_method != "") {
        args <- c(args, input$adjust_method)
      }
      result <- do.call(quantResponse_categPredictor, args)
  
    } else if (script_name == "test_alleleFreqs_within_categories.R") {
      result <- test_alleleFreqs_within_categories(df, input$allele_predictor)
  
    } else if (script_name == "test_categFreqs_within_quantQuadrants.R") {
      result <- test_categFreqs_within_quantQuadrants(df)
  
    } else if (script_name == "test_resampled_distribution_properties.R") {
      result <- test_resampled_distribution_properties(df)
  
    } else if (script_name == "test_sequence_diversity.R") {
      result <- test_sequence_diversity(df)
    }
  
    output$test_output <- renderText({ paste(result, collapse = "\n") })
  })

}





shinyApp(ui = ui, server = server, options = list(host = "0.0.0.0", port = 3838))
