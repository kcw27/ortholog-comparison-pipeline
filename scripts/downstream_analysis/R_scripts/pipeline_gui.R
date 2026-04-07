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
      
      conditionalPanel(
        "output.tabs_ready",
        wellPanel(
          h4("Subset data by group"),
          selectInput("subset_group_var", "Group variable:", 
                     choices = c("category", "subcategory", "genus", "species"), 
                     selected = "category"),
          uiOutput("subset_levels_ui"),
          fluidRow(
            column(6, actionButton("apply_subset", "Apply Subset", style = "width: 100%;")),
            column(6, actionButton("clear_subset", "Clear Subset", style = "width: 100%;"))
          ),
          br(),
          fluidRow(
            column(6, actionButton("reverse_selection", "Reverse Selection", style = "width: 100%;")),
            column(6, actionButton("deselect_all", "Deselect All", style = "width: 100%;"))
          )
        )
      ),
      
      # UPDATED: Subset FASTA and Save Data buttons
      conditionalPanel(
        "output.tabs_ready",
        fluidRow(
          column(6, 
                 actionButton("subset_fasta_btn", "Subset FASTA", style = "width: 100%;")
          ),
          column(6,
                 actionButton("save_data_btn", "Save current data to TSV", style = "width: 100%;")
          )
        )
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
                    choices = c("make_histograms.R",
                                "make_violin_plots.R")),
        # inputs for make_histograms.R
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
        
        # inputs for make_violin_plots.R
        conditionalPanel("input.fig_script == 'make_violin_plots.R'",
          selectInput("violin_column", "Column for violin plot:",
                      choices = NULL),
          selectInput("group_var", "Group by:",
                      choices = c("category", "genus"),
                      selected = "category"),
          selectInput("additional", "Overlay:",
                      choices = c("boxplot", "mean", "none"),
                      selected = "boxplot"),
        ),
        
        # NEW: inputs for make_sequence_logos.R
        conditionalPanel("input.fig_script == 'make_sequence_logos.R' && input.fasta_type == 'Aligned'",
          # Reference sequence input - only show if not already submitted
          conditionalPanel("output.ref_seqID_available == false",
            textInput("seq_logo_ref_seqID", "Reference sequence ID:", ""),
            actionButton("submit_seq_logo_ref_btn", "Submit Reference Sequence ID")
          ),
          conditionalPanel("output.ref_seqID_available == true",
            textOutput("current_ref_seqID")
          ),
          textInput("sites_str", "Enter sites to annotate in the sequence logo (comma-delimited):", "1, 2, 3"),
          selectInput("seq_colname", "Sequence column (must all be the same length):",
                      choices = NULL),
          selectInput("group_var", "Group by:",
                      choices = c("category", "genus"),
                      selected = "category"),
          numericInput("logo_width", "Width (inches):", value = 8, min = 1, max = 20, step = 0.5),
          numericInput("logo_height", "Height (inches):", value = 6, min = 1, max = 20, step = 0.5)
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
          selectInput("normality_col", "Column to test:", choices = NULL),
          selectInput("normality_group", "Group to test for variance (optional):", choices = NULL)
        ),
        
        # Inputs for test_quantResponse_categPredictor.R
        conditionalPanel("input.test_script == 'test_quantResponse_categPredictor.R'",
          selectInput("response_col", "Response (numeric):", choices = NULL),
          selectInput("predictor_col", "Predictor:", choices = c("category", "subcategory", "genus", "species")),
          selectInput("adjust_method", "Adjustment method:",
                      choices = c("", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                      selected = "fdr"),
          selectInput("test_type", "Test type:", choices = c("parametric", "non-parametric"))
        ),
        
        # NEW: Inputs for test_alleleFreqs_within_categories.R
        conditionalPanel("input.test_script == 'test_alleleFreqs_within_categories.R' && input.fasta_type == 'Aligned'",
          # Reference sequence input - only show if not already submitted
          conditionalPanel("output.ref_seqID_available == false",
            textInput("allele_ref_seqID", "Reference sequence ID:", ""),
            actionButton("submit_allele_ref_btn", "Submit Reference Sequence ID")
          ),
          conditionalPanel("output.ref_seqID_available == true",
            textOutput("current_allele_ref_seqID")
          ),
          numericInput("allele_site", "Site (reference coordinate):", value = 1, min = 1, step = 1),
          selectInput("allele_seq_colname", "Sequence column:",
                      choices = NULL),
          selectInput("allele_group_var", "Group by:",
                      choices = c("category", "genus"),
                      selected = "category"),
          uiOutput("allele_residue_ui"),  # Dynamic residue selection
          selectInput("adjust_method", "Adjustment method:",
                      choices = c("", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                      selected = "fdr")
        ),
        
        # Inputs for other scripts can be added similarly
        # Example: test_alleleFreqs_within_categories.R (old version - keep for backward compatibility)
        conditionalPanel("input.test_script == 'test_alleleFreqs_within_categories.R' && input.fasta_type != 'Aligned'",
          selectInput("allele_predictor", "Predictor:", choices = c("category", "genus"))
        ),
        
        # Inputs for other scripts
        conditionalPanel("input.test_script == 'test_categFreqs_within_quantQuadrants.R'",
          # Add inputs if needed
        ),
        # Inputs for test_resampled_distribution_properties.R
        conditionalPanel("input.test_script == 'test_resampled_distribution_properties.R'",
          selectInput("resample_group_var", "Group by:", 
                      choices = c("category", "genus"), 
                      selected = "category"),
          selectInput("resample_data_col", "Data column:", choices = NULL),
          selectInput("resample_adjust_method", "Adjustment method:",
                      choices = c("", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                      selected = "fdr"),
          selectInput("resample_test_type", "Test type:", 
                      choices = c("non-parametric", "parametric"),
                      selected = "non-parametric"),
          numericInput("n_resamples", "Number of resamples:", value = 1000, min = 10, max = 10000, step = 10),
          numericInput("resample_seed", "Random seed:", value = 42, min = 1, max = 10000, step = 1),
          checkboxInput("downsample", "Downsample to smallest group", value = TRUE)
        ),
        
        # Inputs for test_sequence_diversity.R
        conditionalPanel("input.test_script == 'test_sequence_diversity.R'",
          # Reference sequence input - only show if not already submitted
          conditionalPanel("output.ref_seqID_available == false",
            textInput("diversity_ref_seqID", "Reference sequence ID:", ""),
            actionButton("submit_diversity_ref_btn", "Submit Reference Sequence ID")
          ),
          conditionalPanel("output.ref_seqID_available == true",
            textOutput("current_diversity_ref_seqID")
          ),
          selectInput("diversity_group_var", "Group by:", 
                      choices = c("category", "genus"), 
                      selected = "category"),
          selectInput("diversity_seq_colname", "Sequence column:", choices = NULL),
          selectInput("diversity_adjust_method", "Adjustment method:",
                      choices = c("", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                      selected = "fdr"),
          selectInput("diversity_test_type", "Test type:", 
                      choices = c("non-parametric", "parametric"),
                      selected = "non-parametric")
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
  selected_df_store <- reactiveVal(NULL)  # NEW: Store for subsetted data
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
  
  # UPDATED: Get current dataframe (subsetted or full) - MOVED UP
  get_current_df <- reactive({
    if (!is.null(selected_df_store())) {
      return(selected_df_store())
    } else {
      return(df_store())
    }
  })
  
  append_log <- function(msg) {
    log_store$lines <- c(log_store$lines, msg)
  }

  output$tabs_ready <- reactive({ tabs_ready() })
  outputOptions(output, "tabs_ready", suspendWhenHidden = FALSE)
  
  output$show_output <- reactive({ show_output() })
  outputOptions(output, "show_output", suspendWhenHidden = FALSE)

  output$subset_levels_ui <- renderUI({
    req(df_store())
    df <- get_current_df()  # Use current df (subsetted or original)
    group_var <- input$subset_group_var
    
    if (!group_var %in% names(df)) {
      return(div("Selected group variable not found in data"))
    }
    
    # Calculate counts for each level
    level_counts <- table(df[[group_var]], useNA = "no")
    levels <- sort(names(level_counts))
    
    # Create labels with counts
    level_labels <- glue("{levels} ({level_counts[levels]})")
    
    # Get currently selected levels (preserve selection)
    current_selection <- input$subset_levels %||% levels
    
    checkboxGroupInput("subset_levels", "Select levels to include:",
                       choices = setNames(levels, level_labels),
                       selected = current_selection[current_selection %in% levels])
  })
  
  # Update figure script choices based on FASTA type
  observe({
    if (input$fasta_type == "Aligned") {
      updateSelectInput(session, "fig_script",
                        choices = c("make_histograms.R",
                                    "make_violin_plots.R",
                                    "make_sequence_logos.R"))
    } else {
      updateSelectInput(session, "fig_script",
                        choices = c("make_histograms.R",
                                    "make_violin_plots.R"))
    }
  })

  # NEW: Apply subset - MODIFIED for sequential subsetting
  observeEvent(input$apply_subset, {
    req(df_store(), input$subset_group_var, input$subset_levels)
    
    df <- get_current_df()  # CHANGED: Use current df instead of always starting from original
    group_var <- input$subset_group_var
    
    if (!group_var %in% names(df)) {
      append_log("Error: Selected group variable not found in data")
      return()
    }
    
    # Filter to selected levels
    selected_df <- df[df[[group_var]] %in% input$subset_levels, ]
    
    if (nrow(selected_df) == 0) {
      append_log("Warning: No data matches the selected levels")
      return()
    }
    
    selected_df_store(selected_df)
    append_log(glue("Data subset applied: {nrow(selected_df)} rows selected from {group_var} levels: {paste(input$subset_levels, collapse=', ')}"))
  })
  
  # NEW: Clear subset
  observeEvent(input$clear_subset, {
    selected_df_store(NULL)
    append_log("Data subset cleared - using full dataset")
  })
  
  # NEW: Reverse selection
  observeEvent(input$reverse_selection, {
    req(df_store(), input$subset_group_var)
  
    df <- df_store()
    group_var <- input$subset_group_var
  
    if (!group_var %in% names(df)) return()
  
    all_levels <- sort(unique(df[[group_var]]))
    all_levels <- all_levels[!is.na(all_levels)]
  
    current <- input$subset_levels %||% character(0)
  
    reversed <- setdiff(all_levels, current)
  
    updateCheckboxGroupInput(
      session, "subset_levels",
      selected = reversed
    )
  })
  
  # NEW: Deselect all levels
  observeEvent(input$deselect_all, {
    updateCheckboxGroupInput(session, "subset_levels",
                            selected = character(0))
  })
  
  # NEW: Output to check if reference sequence ID is available
  output$ref_seqID_available <- reactive({
    !is.null(ref_seqID_store())
  })
  outputOptions(output, "ref_seqID_available", suspendWhenHidden = FALSE)
  
  # NEW: Display current reference sequence ID
  output$current_ref_seqID <- renderText({
    if (!is.null(ref_seqID_store())) {
      glue("Using reference sequence: {ref_seqID_store()}")
    } else {
      "No reference sequence set"
    }
  })

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
    
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive=TRUE)
      print(glue("{outdir} created"))
    }
    
    result <- if (input$fasta_type == "Unaligned") {
      process_unaligned_shiny(fasta_file, metadata_file, outdir, script_dir, log_fn = append_log)
    } else {
      process_aligned_shiny(fasta_file, metadata_file, outdir, script_dir, name_map_file = name_map_file, log_fn = append_log)
    }
  
    df_store(result$df)
    selected_df_store(NULL)  # Reset subset when new data is loaded
    test_type_store(result$test_type)
                      
    # Update column selectors for statistical tests AND figures
    numeric_cols <- names(result$df)[sapply(result$df, is.numeric)]
    
    # Find columns with "sequence" in the name that are non-numeric
    all_cols <- names(result$df)
    sequence_cols <- all_cols[grepl("sequence", all_cols, ignore.case = TRUE) & 
                              !sapply(result$df, is.numeric) &
                              !grepl("sequence_id", all_cols, ignore.case = TRUE)]
    
    if (length(sequence_cols) == 0) {
      # If no sequence columns found, show all non-numeric columns
      sequence_cols <- all_cols[!sapply(result$df, is.numeric)]
    }
    
    
    updateSelectInput(session, "response_col",
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])
    updateSelectInput(session, "normality_col",
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])
    updateSelectInput(session, "normality_group",  
                      choices = c("none", "category", "genus"),
                      selected = "none")
    updateSelectInput(session, "hist_column",
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])
    updateSelectInput(session, "violin_column",  
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])
    updateSelectInput(session, "resample_data_col",
                      choices = numeric_cols,
                      selected = if ("sequence_length" %in% numeric_cols) "sequence_length" else numeric_cols[1])
    updateSelectInput(session, "diversity_seq_colname",
                      choices = sequence_cols,
                      selected = if ("sequence" %in% sequence_cols) "sequence" else sequence_cols[1])
    
    # NEW: Update sequence column selectors for both sequence logos AND allele frequencies
    if (input$fasta_type == "Aligned") {
      # Update sequence logos column selector
      updateSelectInput(session, "seq_colname",
                        choices = sequence_cols,
                        selected = if ("sequence" %in% sequence_cols) "sequence" else sequence_cols[1])
      
      # NEW: Update allele frequencies column selector
      updateSelectInput(session, "allele_seq_colname",
                        choices = sequence_cols,
                        selected = if ("sequence" %in% sequence_cols) "sequence" else sequence_cols[1])
    }
  
    tabs_ready(TRUE)
    
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
  
  # NEW: Handle reference sequence ID submission for sequence logos
  observeEvent(input$submit_seq_logo_ref_btn, {
    req(input$fasta_file)
    
    ref_seqID <- strip_quotes(input$seq_logo_ref_seqID)
    
    if (ref_seqID == "") {
      append_log("No reference sequence ID provided for sequence logos.")
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
  
  # NEW: Display current reference sequence ID for allele frequencies
  output$current_allele_ref_seqID <- renderText({
    if (!is.null(ref_seqID_store())) {
      glue("Using reference sequence: {ref_seqID_store()}")
    } else {
      "No reference sequence set"
    }
  })
  
  # NEW: Handle reference sequence ID submission for allele frequencies
  observeEvent(input$submit_allele_ref_btn, {
    req(input$fasta_file)
    
    ref_seqID <- strip_quotes(input$allele_ref_seqID)
    
    if (ref_seqID == "") {
      append_log("No reference sequence ID provided for allele frequency test.")
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
  
  # NEW: Dynamic UI for residue selection based on site and sequence column
  output$allele_residue_ui <- renderUI({
    req(df_store(), input$allele_site, input$allele_seq_colname, ref_seqID_store())
    
    # Add selected_df_store() as a dependency to force update when data is subsetted
    selected_df_store()
    
    # Use FULL dataframe for coordinate calibration (to ensure reference sequence is available)
    full_df <- df_store()
    site <- input$allele_site
    seq_colname <- input$allele_seq_colname
    
    if (!seq_colname %in% names(full_df)) {
      return(div("Selected sequence column not found in data"))
    }
    
    # Calculate adjusted site using reference sequence from FULL dataframe
    matching_rows <- which(full_df$sequence_id == ref_seqID_store())
    if (length(matching_rows) == 0) {
      return(div(glue("Reference sequence '{ref_seqID_store()}' not found in dataset")))
    }
    
    reference_seq <- full_df[[seq_colname]][matching_rows[1]]
    sites_adj <- calibrate_coords(reference_seq, site)
    
    if (length(sites_adj) == 0) {
      return(div(glue("Site {site} not found in reference sequence (may be in gapped region)")))
    }
    
    adj_site <- sites_adj[1]
    
    # Use CURRENT (possibly subsetted) dataframe for residue counts
    current_df <- get_current_df()
    
    # Extract characters at the adjusted site from CURRENT dataframe
    characters <- substr(current_df[[seq_colname]], adj_site, adj_site)
    characters <- characters[!is.na(characters) & characters != ""]
    
    if (length(characters) == 0) {
      return(div("No characters found at this site in current data subset"))
    }
    
    # Get unique residues and their frequencies from CURRENT dataframe
    residue_counts <- table(characters)
    residue_options <- names(residue_counts)
    
    # Create labels with frequencies from CURRENT dataframe
    residue_labels <- glue("{residue_options} ({residue_counts[residue_options]} sequences)")
    
    selectInput("allele_residue", "Residue to test:",
                choices = setNames(residue_options, residue_labels))
  })
  
  # NEW: Display current reference sequence ID for sequence diversity
  output$current_diversity_ref_seqID <- renderText({
    if (!is.null(ref_seqID_store())) {
      glue("Using reference sequence: {ref_seqID_store()}")
    } else {
      "No reference sequence set"
    }
  })
  
  # NEW: Handle reference sequence ID submission for sequence diversity
  observeEvent(input$submit_diversity_ref_btn, {
    req(input$fasta_file)
    
    ref_seqID <- strip_quotes(input$diversity_ref_seqID)
    
    if (ref_seqID == "") {
      append_log("No reference sequence ID provided for sequence diversity test.")
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
  
  # NEW: Update test script choices based on FASTA type
  observe({
    if (input$fasta_type == "Aligned") {
      # Show all tests including allele frequencies
      updateSelectInput(session, "test_script",
                        choices = c("test_normality.R",
                                    "test_quantResponse_categPredictor.R",
                                    "test_alleleFreqs_within_categories.R",
                                    "test_categFreqs_within_quantQuadrants.R",
                                    "test_resampled_distribution_properties.R",
                                    "test_sequence_diversity.R"))
    } else {
      # Hide allele frequency test for unaligned data
      updateSelectInput(session, "test_script",
                        choices = c("test_normality.R",
                                    "test_quantResponse_categPredictor.R",
                                    "test_categFreqs_within_quantQuadrants.R",
                                    "test_resampled_distribution_properties.R",
                                    "test_sequence_diversity.R"))
    }
  })

  observeEvent(input$run_test, {
    df <- get_current_df()  # UPDATED: Use current df (subsetted or full)
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
      groupname <- input$normality_group
      if (is.null(colname) || colname == "") {
        output$test_output <- renderText({ "Please select a column to test." })
        return()
      }
      
      # For test_normality: capture print statements and show them, then show return value separately
      test_output <- capture.output({
        if (groupname == "none") {
          result <- test_normality(df, colname, file.path(input$outdir, "figures"))
        } else {
          result <- test_normality(df, colname, file.path(input$outdir, "figures"), groupname)
        }
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
      # NEW: Handle allele frequency test
      if (input$fasta_type == "Aligned") {
        # Use the new UI inputs for aligned data
        if (is.null(ref_seqID_store())) {
          output$test_output <- renderText({ "Error: Reference sequence ID required for allele frequency test" })
          return()
        }
        if (is.null(input$allele_residue)) {
          output$test_output <- renderText({ "Error: Please select a residue to test" })
          return()
        }
        
        # For coordinate calibration, use the FULL dataframe to ensure reference sequence is available
        full_df <- df_store()
        
        # Calculate adjusted site using full dataframe
        matching_rows <- which(full_df$sequence_id == ref_seqID_store())
        if (length(matching_rows) == 0) {
          output$test_output <- renderText({ 
            glue("Error: Reference sequence '{ref_seqID_store()}' not found in dataset")
          })
          return()
        }
        
        reference_seq <- full_df[[input$allele_seq_colname]][matching_rows[1]]
        sites_adj <- calibrate_coords(reference_seq, input$allele_site)
        
        if (length(sites_adj) == 0) {
          output$test_output <- renderText({ 
            glue("Error: Site {input$allele_site} not found in reference sequence (may be in gapped region)")
          })
          return()
        }
        
        adj_site <- sites_adj[1]
        
        # Pass the adjusted site to the function (not the original site)
        args <- list(
          df = df,  # Use current (possibly subsetted) df for the actual statistical test
          site = adj_site,  # Pass the pre-adjusted site
          residue = input$allele_residue,
          group_var = input$allele_group_var,
          seq_colname = input$allele_seq_colname
        )
        if (input$adjust_method != "") {
          args <- c(args, input$adjust_method)
        }
        
        test_output <- capture.output({
          result <- do.call(test_alleleFreqs_within_categories, args)
        })
        
        output$test_output <- renderText({ paste(test_output, collapse = "\n") })
      }
    } else if (script_name == "test_categFreqs_within_quantQuadrants.R") {
      result <- test_categFreqs_within_quantQuadrants(df)
      output$test_output <- renderText({ paste(result, collapse = "\n") })
    } else if (script_name == "test_resampled_distribution_properties.R") {
      # Run resampled distribution properties test
      if (is.null(input$resample_data_col) || input$resample_data_col == "") {
        output$test_output <- renderText({ "Please select a data column." })
        return()
      }
      
      test_output <- capture.output({
        # Run bootstrap_stats
        bs_df <- bootstrap_stats(
          df = df,
          group_var = input$resample_group_var,
          data_col = input$resample_data_col,
          n = input$n_resamples,
          seed = input$resample_seed,
          downsample = input$downsample
        )
        
        # Run compare_bootstrapped_stats
        compare_bootstrapped_stats(
          bs_df = bs_df,
          group_var = input$resample_group_var,
          adjust = if (input$resample_adjust_method != "") input$resample_adjust_method else "fdr",
          test_type = input$resample_test_type
        )
        
        # Run plot_bootstrapped_stats
        plot_bootstrapped_stats(
          original_df = df,
          bootstrapped_df = bs_df,
          outdir = file.path(input$outdir, "figures"),
          group_var = input$resample_group_var,
          data_col = input$resample_data_col,
          downsample = input$downsample
        )
      })
      
      output$test_output <- renderText({ paste(test_output, collapse = "\n") })
    } else if (script_name == "test_sequence_diversity.R") {
      # Run sequence diversity test with new functions
      if (is.null(input$diversity_seq_colname) || input$diversity_seq_colname == "") {
        output$test_output <- renderText({ "Please select a sequence column." })
        return()
      }
      
      if (is.null(ref_seqID_store())) {
        output$test_output <- renderText({ "Error: Reference sequence ID required for sequence diversity test" })
        return()
      }
      
      test_output <- capture.output({
        # Run get_mutational_richness with reference sequence
        dist_df <- get_mutational_richness(
          df = df,
          reference_sequence = ref_sequence_store(),
          group_var = input$diversity_group_var,
          seq_colname = input$diversity_seq_colname
        )
        
        # Run compare_mutational_richness
        compare_mutational_richness(
          dist_df = dist_df,
          test_type = input$diversity_test_type,
          adjust = if (input$diversity_adjust_method != "") input$diversity_adjust_method else "fdr"
        )
        
        # Run plot_mutational_richness with reference label
        plot_mutational_richness(
          dist_df = dist_df,
          outdir = file.path(input$outdir, "figures"),
          reference_label = ref_seqID_store(),
          group_var = input$diversity_group_var,
          test_type = input$diversity_test_type
        )
      })
      
      output$test_output <- renderText({ paste(test_output, collapse = "\n") })
    }
  })  # <<< ADDED THIS CLOSING BRACE for observeEvent(input$run_test, ...)

  observeEvent(input$run_figures, {
    df <- get_current_df()
    if (is.null(df)) return()
    fig_script <- input$fig_script
    group_var <- input$group_var
    width <- as.numeric(input$hist_width)
    consistent_y_axis <- input$consistent_y_axis
    additional <- input$additional
    script_path <- file.path(script_dir, "figures", fig_script)
  
    # Get the appropriate column name based on which script is selected
    if (fig_script == "make_histograms.R") {
      column_name <- input$hist_column  # Use histogram column selector
    } else if (fig_script == "make_violin_plots.R") {
      column_name <- input$violin_column  # Use violin column selector
    } else if (fig_script == "make_sequence_logos.R") {
      column_name <- NULL  # Not used for sequence logos
    } else {
      column_name <- NULL
    }
    
    if (fig_script != "make_sequence_logos.R" && is.null(column_name)) {
      append_log("Error: No column selected")
      return()
    }
  
    if (file.exists(script_path)) {
      append_log(glue("Running {fig_script}..."))
      if (fig_script == "make_histograms.R") { # histograms
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
        
      } else if (fig_script == "make_violin_plots.R") { # violin plots
        source(script_path, local = TRUE)
        
        # Prepare arguments for violin_plots_by_source
        args <- list(
          df = df,
          outdir = file.path(input$outdir, "figures"),
          column_name = column_name,
          group_var = group_var,
          additional = additional
        )
        
        # Add reference_label if available (but NOT reference_length)
        if (!is.null(ref_seqID_store())) {
          args$reference_label <- ref_seqID_store()
          append_log(glue("Using reference sequence: {ref_seqID_store()}"))
        }
        
        # Call violin_plots_by_source with the appropriate arguments
        do.call(violin_plots_by_source, args)
        append_log(glue("Violin plots generated for column: {column_name}"))
        
        # Trigger image list refresh after generating figures
        image_refresh_trigger(image_refresh_trigger() + 1)
        append_log("Image list refreshed automatically.")
        
      } else if (fig_script == "make_sequence_logos.R") { # sequence logos
        source(script_path, local = TRUE)
        
        # Check if reference sequence is available
        if (is.null(ref_seqID_store())) {
          append_log("Error: No reference sequence ID available. Please submit one in the Data Import tab or on this page.")
          return()
        }
        
        # Check if sites string is provided
        if (input$sites_str == "") {
          append_log("Error: Please enter sites to annotate.")
          return()
        }
        
        # Prepare arguments for plot_local_logo
        args <- list(
          df = df,
          outdir = file.path(input$outdir, "figures"),
          reference_seq = ref_sequence_store(),
          sites_str = input$sites_str,
          seq_colname = input$seq_colname,
          width = input$logo_width,
          height = input$logo_height,
          group_var = group_var
        )
        
        # Call plot_local_logo with the appropriate arguments
        do.call(plot_local_logo, args)
        append_log(glue("Sequence logo generated for reference: {ref_seqID_store()}"))
        append_log(glue("Annotated sites: {input$sites_str}"))
        append_log(glue("Sequence column: {input$seq_colname}"))
        
        # Trigger image list refresh after generating figures
        image_refresh_trigger(image_refresh_trigger() + 1)
        append_log("Image list refreshed automatically.")
      }
    }
  })
  
  # UPDATED: Subset FASTA functionality
  observeEvent(input$subset_fasta_btn, {
    req(input$fasta_file, input$outdir)
    
    df <- get_current_df()  # UPDATED: Use current df (subsetted or full)
    if (is.null(df)) {
      append_log("Error: No data available for subsetting")
      return()
    }
    
    fasta_file <- normalizePath(as.character(strip_quotes(input$fasta_file)), mustWork = FALSE)
    outdir <- normalizePath(as.character(strip_quotes(input$outdir)), mustWork = FALSE)
    
    # Get sequence IDs from current dataframe
    if ("sequence_id" %in% names(df)) {
      sequence_ids <- df$sequence_id
    } else if ("sequence_id.x" %in% names(df)) {
      sequence_ids <- df$sequence_id.x
    } else {
      # Try to find any column that might contain sequence IDs
      id_cols <- names(df)[grepl("sequence|id", names(df), ignore.case = TRUE)]
      if (length(id_cols) > 0) {
        sequence_ids <- df[[id_cols[1]]]
        append_log(glue("Using column '{id_cols[1]}' for sequence IDs"))
      } else {
        append_log("Error: No sequence ID column found in data")
        return()
      }
    }
    
    sequence_ids <- unique(sequence_ids[!is.na(sequence_ids)])
    append_log(glue("Found {length(sequence_ids)} unique sequence IDs in current data"))
    
    # Read the original FASTA file
    tryCatch({
      library(seqinr)
      alignment <- read.alignment(fasta_file, format = "fasta")
      
      seqs <- alignment[["seq"]] |>
        unlist() |> # convert from list to vector
        toupper() # originally in lowercase; convert to uppercase
      
      alignment_names = (alignment["nam"][[1]]) 
      allseqs_df <- data.frame(sequence_id = alignment_names, sequence = seqs)
      
      # Filter to only sequences in current data
      subset_df <- allseqs_df[allseqs_df$sequence_id %in% sequence_ids, ]
      
      if (nrow(subset_df) == 0) {
        append_log("Warning: No sequences from the FASTA file match the IDs in current data.")
        return()
      }
      
      append_log(glue("Subsetting FASTA: {nrow(subset_df)} sequences match current data IDs"))
      
      # Write the subset to a new FASTA file
      subset_fasta_file <- file.path(outdir, "subset_sequences.fasta")
      
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
  
  # NEW: Save current data to TSV
  observeEvent(input$save_data_btn, {
    req(input$outdir)
    
    df <- get_current_df()  # Get current df (subsetted or full)
    if (is.null(df)) {
      append_log("Error: No data available to save")
      return()
    }
    
    outdir <- normalizePath(as.character(strip_quotes(input$outdir)), mustWork = FALSE)
    
    # Create filename with timestamp
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- glue("current_data_{timestamp}.tsv")
    filepath <- file.path(outdir, filename)
    
    tryCatch({
      write_tsv(df, filepath)
      append_log(glue("Current data saved to: {filepath}"))
      append_log(glue("Data dimensions: {nrow(df)} rows, {ncol(df)} columns"))
    }, error = function(e) {
      append_log(glue("Error saving data: {e$message}"))
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
    
    # Include both image files and TSV files
    files <- list.files(img_dir, pattern = "\\.(png|jpg|jpeg|pdf|tsv|fasta)$", 
                       full.names = TRUE, recursive = TRUE)
    
    if (length(files) == 0) return(div("No images or data files found"))
    
    selectInput("image_file", "Select image or data file:", 
                choices = setNames(files, basename(files)))
  })
  
  # Fixed image display using renderUI with proper resource path - show default message
  output$selected_image <- renderUI({
    # Always show the default message when no image is selected
    # This will reset when images are refreshed since input$image_file becomes NULL
    if (is.null(input$image_file)) {
      return(div("Please select an image or data file to view.", 
                 style = "border: 1px solid #ccc; margin-top: 10px; padding: 20px; text-align: center; color: #666;"))
    }
    
    # Validate file exists
    if (!file.exists(input$image_file)) {
      return(div("Selected file not found"))
    }
    
    file_ext <- tolower(tools::file_ext(input$image_file))
    
    if (file_ext == "tsv") {
      # Handle TSV files - display as a table
      tryCatch({
        # Read the TSV file
        tsv_data <- read.delim(input$image_file, stringsAsFactors = FALSE)
        
        # Create a styled table
        tags$div(
          h4(basename(input$image_file)),
          tags$p(glue("Rows: {nrow(tsv_data)}, Columns: {ncol(tsv_data)}")),
          tags$div(
            style = "max-height: 600px; overflow: auto; border: 1px solid #ccc;",
            renderTable({
              tsv_data
            }, striped = TRUE, hover = TRUE, bordered = TRUE, 
            align = 'c', digits = 3)
          ),
          style = "border: 1px solid #ccc; margin-top: 10px; padding: 10px;"
        )
      }, error = function(e) {
        div(glue("Error reading TSV file: {e$message}"),
            style = "border: 1px solid #ccc; margin-top: 10px; padding: 20px; color: red;")
      })
    } else if (file_ext == "pdf") {
      # For PDFs, use tags$embed to display
      # Set up resource path for serving files
      img_dir <- file.path(input$outdir, "figures")
      if (!resource_path_added()) {
        addResourcePath("figures", img_dir)
        resource_path_added(TRUE)
      }
      
      # Get relative path from the figures directory
      relative_path <- sub(paste0("^", img_dir, "/?"), "", input$image_file)
      
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
      # Set up resource path for serving files
      img_dir <- file.path(input$outdir, "figures")
      if (!resource_path_added()) {
        addResourcePath("figures", img_dir)
        resource_path_added(TRUE)
      }
      
      # Get relative path from the figures directory
      relative_path <- sub(paste0("^", img_dir, "/?"), "", input$image_file)
      
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