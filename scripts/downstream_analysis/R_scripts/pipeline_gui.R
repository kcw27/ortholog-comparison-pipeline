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

ui <- fluidPage(
  titlePanel("Pipeline Input GUI"),
  
  radioButtons("fasta_type", "FASTA type:",
               choices = c("Unaligned", "Aligned"),
               selected = "Unaligned",
               inline = TRUE),
  
  textInput("fasta_file", "FASTA file path:", ""),
  textInput("metadata_file", "Metadata file path:", ""),
  textInput("outdir", "Output directory:", ""),
  
  actionButton("submit_btn", "Submit"),
  verbatimTextOutput("output_text")
)

server <- function(input, output, session) {
  log_store <- reactiveValues(lines = character())

  append_log <- function(msg) {
    log_store$lines <- c(log_store$lines, msg)
  }

  observeEvent(input$submit_btn, {
    # Step 1: sanitize inputs
    fasta_file <- normalizePath(as.character(strip_quotes(input$fasta_file)), mustWork = FALSE)
    metadata_file <- normalizePath(as.character(strip_quotes(input$metadata_file)), mustWork = FALSE)
    outdir <- normalizePath(as.character(strip_quotes(input$outdir)), mustWork = FALSE)

    # Step 2: print line counts and benchmarking
    append_log(glue("FASTA file has {count_lines(fasta_file)} lines"))
    append_log(glue("Metadata file has {count_lines(metadata_file)} lines"))
    append_log("\n")
    benchmark_output <- capture.output(benchmark(metadata_file))
    lapply(benchmark_output, append_log)

    # Step 3: run processing
    append_log("Processing input data...")
    if (input$fasta_type == "Unaligned") {
      capture.output(process_unaligned_shiny(fasta_file, metadata_file, outdir, log_fn = append_log))
    } else {
      capture.output(process_aligned_shiny(fasta_file, metadata_file, outdir, log_fn = append_log))
    }

    # Step 4: count output
    categorized_path <- file.path(outdir, "categorized_ids.txt")
    append_log(glue("categorized_ids.txt has {count_lines(categorized_path)} lines"))
  })

  output$output_text <- renderText({
    invalidateLater(500, session)  # refresh every 500ms
    paste(log_store$lines, collapse = "\n")
  })
}




shinyApp(ui = ui, server = server, options = list(host = "0.0.0.0", port = 3838))
