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
      # Conditional name map input for aligned sequences
      conditionalPanel(
        "input.fasta_type == 'Aligned'",
        textInput("name_map_file", "Name map file path (optional):", "")
      ),
      actionButton("submit_btn", "Submit Data"),
      br(), br(),
      # Reference sequence ID section
      textInput("ref_seqID", "Reference sequence ID (optional):", ""),
      actionButton("submit_ref_btn", "Submit Reference Sequence ID"),
      # Add subset FASTA button - only show when data is loaded
      conditionalPanel(
        "output.tabs_ready",
        actionButton("subset_fasta_btn", "Subset FASTA to Categorized Sequences")
      ),
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
          selectInput("hist_column", "Column for histogram:",
                      choices = NULL),
          selectInput("group_var", "Group by:",
                      choices = c("category", "genus"),
                      selected = "category"),
          selectInput("hist_width", "Bin width:",
                      choices = c(1, 5, 10, 25, 50, 100),
                      selected = 10),
          checkboxInput("consistent_y_axis", "Consistent Y-axis across plots", 
                       value = TRUE)
        ),
        actionButton("run_figures", "Run figure script"),
        actionButton("refresh_images", "Refresh Image List", icon = icon("refresh")),
        uiOutput("image_selector"),
        uiOutput("selected_image")
      )
    )
  ),
  tabPanel("Statistical tests",
    conditionalPanel("output.tabs_ready",
      fluidPage(
        selectInput("test_script", "Select statistical test script:",
                    choices = c("test_normality.R",
                                "test_quantResponse_categPredictor.R",
                                "test_alleleFreqs_within_categories.R",
                                "test_categFreqs_within_quantQuadrants.R",
                                "test_resampled_distribution_properties.R",
                                "test_sequence_diversity.R")),
        
        # Inputs for test_normality.R
        conditionalPanel("input.test_script == 'test_normality.R'",
          selectInput("normality_col", "Column to test:", choices = NULL)
        ),
        
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
  ref_seqID_store <- reactiveVal(NULL)
  ref_length_store <- reactiveVal(NULL)
  ref_sequence_store <- reactiveVal(NULL)
  
  # Store the resource path for serving files
  resource_path_added <- reactiveVal(FALSE)
  # Add reactive value to trigger image list refresh
  image_refresh_trigger <- reactiveVal(0)
  
  append_log <- function(msg) {
    log_store$lines <- c(log_store$lines, msg)
  }

  output$tabs_ready <- reactive({ tabs_ready() })
  outputOptions(output, "tabs_ready", suspendWhenHidden = FALSE)
  
  output$show_output <- reactive({ show_output() })
  outputOptions(output, "show_output", suspendWhenHidden = FALSE)

  observeEvent(input$submit_btn, {
    show_output(TRUE)
    
    # Sanitize inputs
    fasta_file <- normalizePath(as.character(strip_quotes(input$fasta_file)), mustWork = FALSE)
    metadata_file <- normalizePath(as.character(strip_quotes(input$metadata_file)), mustWork = FALSE)
    outdir <- normalizePath(as.character(strip_quotes(input$outdir)), mustWork = FALSE)
    
    # Handle name map file for aligned sequences
    name_map_file <- NULL
    if (input$fasta_type == "Aligned" && input$name_map_file != "") {
      name_map_file <- normalizePath(as.character(strip_quotes(input$name_map_file)), mustWork = FALSE)
    }

    benchmark_output <- capture.output(benchmark(metadata_file))
    lapply(benchmark_output, append_log)

    append_log("\nProcessing input data...")
    result <- if (input$fasta_type == "Unaligned") {
      process_unaligned_shiny(fasta_file, metadata_file, outdir, script_dir, log_fn = append_log)
    } else {
      process_aligned_shiny(fasta_file, metadata_file, outdir, script_dir, name_map_file = name_map_file, log_fn = append_log)
    }

    df_store(result$df)
    test_type_store(result$test_type)
                      
    # Update column selectors for statistical tests AND figures
    numeric_cols <- names(result$df)[sapply(result$df, is.numeric)]
    updateSelectInput(session, "response_col",
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])
    updateSelectInput(session, "normality_col",
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])
    updateSelectInput(session, "hist_column",
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])
  
      tabs_ready(TRUE)

#    # Run normality test and capture output
#    test_output <- capture.output({
#      invisible(test_type <- test_normality(result$df, "sequence_length", file.path(outdir, "figures")))
#    })
#    lapply(test_output, append_log)
#    test_type_store(test_type)
#    
#    # Show QQ plot image and path
#    qqplot_path <- file.path(outdir, "figures", "qqplot_sequence_length.png")
#    
#    output$qqplot_path <- renderText({
#      if (file.exists(qqplot_path)) {
#        glue("QQ plot saved to: {qqplot_path}")
#      } else {
#        "QQ plot not found."
#      }
#    })
#    
#    output$qqplot_image <- renderImage({
#      if (file.exists(qqplot_path)) {
#        list(src = qqplot_path, contentType = "image/png", alt = "QQ Plot")
#      } else {
#        NULL
#      }
#    }, deleteFile = FALSE)
    
    output$output_text <- renderText({
      invalidateLater(500, session)
      paste(log_store$lines, collapse = "\n")
    })
  })

  # Handle reference sequence ID submission - look in FASTA file
  observeEvent(input$submit_ref_btn, {
    req(input$fasta_file)
    
    ref_seqID <- strip_quotes(input$ref_seqID)
    
    if (ref_seqID == "") {
      append_log("No reference sequence ID provided.")
      ref_seqID_store(NULL)
      ref_length_store(NULL)
      ref_sequence_store(NULL)
      return()
    }
    
    fasta_file <- normalizePath(as.character(strip_quotes(input$fasta_file)), mustWork = FALSE)
    
    if (!file.exists(fasta_file)) {
      append_log(glue("FASTA file not found: {fasta_file}"))
      return()
    }
    
    # Read the FASTA file to find the reference sequence
    tryCatch({
      library(seqinr)
      alignment <- read.alignment(fasta_file, format = "fasta")
      
      seqs <- alignment[["seq"]] |>
        unlist() |> # convert from list to vector
        toupper() # originally in lowercase; convert to uppercase
      
      alignment_names = (alignment["nam"][[1]]) 
      allseqs_df <- data.frame(sequence_id = alignment_names, sequence = seqs)
      
      # Check if reference sequence ID exists
      matching_rows <- which(allseqs_df$sequence_id == ref_seqID)
      
      if (length(matching_rows) == 0) {
        append_log(glue("WARNING: Reference sequence ID '{ref_seqID}' not found in FASTA file."))
        append_log("Please submit a sequence ID that exists in the FASTA file.")
        ref_seqID_store(NULL)
        ref_length_store(NULL)
        ref_sequence_store(NULL)
      } else if (length(matching_rows) > 1) {
        append_log(glue("WARNING: Reference sequence ID '{ref_seqID}' appears {length(matching_rows)} times in FASTA file."))
        append_log("The first matching sequence will be used.")
        ref_seqID_store(ref_seqID)
        # Use first occurrence
        first_match <- matching_rows[1]
        ref_sequence_store(allseqs_df$sequence[first_match])
        ref_length_store(nchar(allseqs_df$sequence[first_match]))
      } else {
        append_log(glue("Reference sequence ID '{ref_seqID}' found in FASTA file."))
        ref_seqID_store(ref_seqID)
        ref_sequence_store(allseqs_df$sequence[matching_rows])
        ref_length_store(nchar(allseqs_df$sequence[matching_rows]))
      }
      
      if (!is.null(ref_length_store())) {
        append_log(glue("Reference sequence length: {ref_length_store()}"))
      }
      if (!is.null(ref_sequence_store())) {
        append_log(glue("Reference sequence stored (first 50 chars): {substr(ref_sequence_store(), 1, 50)}..."))
      }
    }, error = function(e) {
      append_log(glue("Error reading FASTA file: {e$message}"))
    })
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
  
    result <- NULL
    if (script_name == "test_normality.R") {
      # Run normality test on selected column
      colname <- input$normality_col
      if (is.null(colname) || colname == "") {
        output$test_output <- renderText({ "Please select a column to test." })
        return()
      }
      
      # For test_normality: capture print statements and show them, then show return value separately
      test_output <- capture.output({
        result <- test_normality(df, colname, file.path(input$outdir, "figures"))
      })
      
      # Combine captured output and return value clearly separated
      full_output <- c(test_output, "", "Final test type:", result)
      output$test_output <- renderText({ paste(full_output, collapse = "\n") })
      
    } else if (script_name == "test_quantResponse_categPredictor.R") {
      args <- list(df, input$response_col, input$predictor_col, input$test_type)
      if (input$adjust_method != "") {
        args <- c(args, input$adjust_method)
      }
  
      # For quantResponse_categPredictor: the function returns formatted output, so just capture and display it
      test_output <- capture.output({
        result <- do.call(quantResponse_categPredictor, args)
      })
      
      # Since quantResponse_categPredictor returns formatted strings, just display the result
      output$test_output <- renderText({ result })
  
      # Render QQ plot image
      qqplot_path <- file.path(input$outdir, "figures", glue("qqplot_{input$response_col}.png"))
      output$qqplot_image <- renderImage({
        list(src = qqplot_path, contentType = "image/png", alt = "QQ Plot")
      }, deleteFile = FALSE)
    } else if (script_name == "test_alleleFreqs_within_categories.R") {
      result <- test_alleleFreqs_within_categories(df, input$allele_predictor)
      output$test_output <- renderText({ paste(result, collapse = "\n") })
  
    } else if (script_name == "test_categFreqs_within_quantQuadrants.R") {
      result <- test_categFreqs_within_quantQuadrants(df)
      output$test_output <- renderText({ paste(result, collapse = "\n") })
  
    } else if (script_name == "test_resampled_distribution_properties.R") {
      result <- test_resampled_distribution_properties(df)
      output$test_output <- renderText({ paste(result, collapse = "\n") })
  
    } else if (script_name == "test_sequence_diversity.R") {
      result <- test_sequence_diversity(df)
      output$test_output <- renderText({ paste(result, collapse = "\n") })
    }
  })

  observeEvent(input$run_figures, {
    df <- df_store()
    if (is.null(df)) return()
    fig_script <- input$fig_script
    group_var <- input$group_var
    column_name <- input$hist_column
    width <- as.numeric(input$hist_width)
    consistent_y_axis <- input$consistent_y_axis  # Get the toggle value
    script_path <- file.path(script_dir, "figures", fig_script)
  
    if (file.exists(script_path)) {
      append_log(glue("Running {fig_script}..."))
      if (fig_script == "make_histograms.R") {
        source(script_path, local = TRUE)
        
        # Prepare arguments for histograms_by_source
        args <- list(
          df = df,
          outdir = file.path(input$outdir, "figures"),
          column_name = column_name,
          group_var = group_var,
          width = width,
          consistent_y_axis = consistent_y_axis  # Pass the toggle value
        )
        
        # Add reference_label if available (but NOT reference_length)
        if (!is.null(ref_seqID_store())) {
          args$reference_label <- ref_seqID_store()
          append_log(glue("Using reference sequence: {ref_seqID_store()}"))
        }
        
        # Call histograms_by_source with the appropriate arguments
        do.call(histograms_by_source, args)
        append_log(glue("Histograms generated for column: {column_name}"))
        append_log(glue("Consistent Y-axis: {consistent_y_axis}"))
        
        # Trigger image list refresh after generating figures
        image_refresh_trigger(image_refresh_trigger() + 1)
        append_log("Image list refreshed automatically.")
      }
    }
  })
  
  # Add the subset FASTA functionality
  observeEvent(input$subset_fasta_btn, {
    req(input$fasta_file, input$outdir)
    
    fasta_file <- normalizePath(as.character(strip_quotes(input$fasta_file)), mustWork = FALSE)
    outdir <- normalizePath(as.character(strip_quotes(input$outdir)), mustWork = FALSE)
    categorized_ids_file <- file.path(outdir, "categorized_ids.txt")
    
    if (!file.exists(categorized_ids_file)) {
      append_log("Error: categorized_ids.txt not found. Please submit data first.")
      return()
    }
    
    # Read the categorized IDs
    categorized_ids <- readLines(categorized_ids_file)
    append_log(glue("Found {length(categorized_ids)} categorized sequence IDs"))
    
    # Read the original FASTA file
    tryCatch({
      library(seqinr)
      alignment <- read.alignment(fasta_file, format = "fasta")
      
      seqs <- alignment[["seq"]] |>
        unlist() |> # convert from list to vector
        toupper() # originally in lowercase; convert to uppercase
      
      alignment_names = (alignment["nam"][[1]]) 
      allseqs_df <- data.frame(sequence_id = alignment_names, sequence = seqs)
      
      # Filter to only categorized sequences
      subset_df <- allseqs_df[allseqs_df$sequence_id %in% categorized_ids, ]
      
      if (nrow(subset_df) == 0) {
        append_log("Warning: No sequences from the FASTA file match the categorized IDs.")
        return()
      }
      
      append_log(glue("Subsetting FASTA: {nrow(subset_df)} sequences match categorized IDs"))
      
      # Write the subset to a new FASTA file
      subset_fasta_file <- file.path(outdir, "categorized_sequences.fasta")
      
      # Create FASTA format
      fasta_lines <- character()
      for (i in 1:nrow(subset_df)) {
        fasta_lines <- c(fasta_lines, 
                         paste0(">", subset_df$sequence_id[i]),
                         subset_df$sequence[i])
      }
      
      writeLines(fasta_lines, subset_fasta_file)
      append_log(glue("Subset FASTA saved to: {subset_fasta_file}"))
      
    }, error = function(e) {
      append_log(glue("Error subsetting FASTA file: {e$message}"))
    })
  })
  
  # Refresh button handler
  observeEvent(input$refresh_images, {
    image_refresh_trigger(image_refresh_trigger() + 1)
    append_log("Image list manually refreshed.")
  })
  
  # Image selector with refresh trigger dependency
  output$image_selector <- renderUI({
    req(tabs_ready())
    # Depend on refresh trigger to update when triggered
    image_refresh_trigger()
    
    img_dir <- file.path(input$outdir, "figures")
    if (!dir.exists(img_dir)) return(div("No figures directory"))
    
    imgs <- list.files(img_dir, pattern = "\\.(png|jpg|jpeg|pdf)$", 
                       full.names = TRUE, recursive = TRUE)
    
    if (length(imgs) == 0) return(div("No images found"))
    
    selectInput("image_file", "Select image:", 
                choices = setNames(imgs, basename(imgs)))
  })
  
  # Fixed image display using renderUI with proper resource path - show default message
  output$selected_image <- renderUI({
    # Always show the default message when no image is selected
    # This will reset when images are refreshed since input$image_file becomes NULL
    if (is.null(input$image_file)) {
      return(div("Please select an image to view.", 
                 style = "border: 1px solid #ccc; margin-top: 10px; padding: 20px; text-align: center; color: #666;"))
    }
    
    # Validate file exists
    if (!file.exists(input$image_file)) {
      return(div("Selected file not found"))
    }
    
    # Set up resource path for serving files
    img_dir <- file.path(input$outdir, "figures")
    if (!resource_path_added()) {
      addResourcePath("figures", img_dir)
      resource_path_added(TRUE)
    }
    
    file_ext <- tolower(tools::file_ext(input$image_file))
    
    # Get relative path from the figures directory
    relative_path <- sub(paste0("^", img_dir, "/?"), "", input$image_file)
    
    if (file_ext == "pdf") {
      # For PDFs, use tags$embed to display
      tags$div(
        tags$embed(
          src = file.path("figures", relative_path),
          type = "application/pdf",
          width = "100%",
          height = "600px"
        ),
        style = "border: 1px solid #ccc; margin-top: 10px;"
      )
    } else {
      # For images, use img tag with proper styling and resource path
      tags$div(
        tags$img(
          src = file.path("figures", relative_path),
          style = "max-width: 100%; height: auto; display: block;",
          alt = "Figure"
        ),
        style = "border: 1px solid #ccc; margin-top: 10px; text-align: center;"
      )
    }
  })

  observe({
    updateSelectInput(session, "test_type",
                      choices = c("parametric", "non-parametric"),
                      selected = test_type_store())
  })

}

shinyApp(ui = ui, server = server, options = list(host = "0.0.0.0", port = 3838))