###################################################################################
################# HELPER FUNCTIONS FOR SAVING FILES AUTOMATICALLY ###############
###################################################################################

# Initialize analysis environment with folder structure and plotting defaults
init_analysis <- function(analysis_folder, base_results_path = "data/results") {
  
  # Set Arial font for all base R plots
  par(family = "Arial")
  
  # Store analysis folder name globally for use in helper functions
  assign("ANALYSIS_FOLDER", analysis_folder, envir = .GlobalEnv)
  
  # Define full paths for results and figures
  results_path <- here::here(base_results_path, analysis_folder)
  figures_path <- here::here(base_results_path, analysis_folder, "00_figures")
  
  # Create directory structure if it doesn't exist
  if (!dir.exists(results_path)) {
    dir.create(results_path, recursive = TRUE)
  }
  if (!dir.exists(figures_path)) {
    dir.create(figures_path, recursive = TRUE)
  }
  
  # Set knitr chunk options for consistent output
  knitr::opts_chunk$set(
    echo = TRUE,                    # Show code in output
    message = FALSE,                # Hide package loading messages
    warning = FALSE,                # Hide warnings
    fig.width = 10,                 # Default figure width
    fig.height = 8,                 # Default figure height
    fig.align = 'center',           # Center figures
    dev = c("svg", "pdf"),          # Save figures in SVG and PDF formats
    dev.args = list(svg = list(family = "Arial")), # Arial font for SVG figures
    fig.path = here::here(base_results_path, analysis_folder, "00_figures", "") # Use here() for absolute path
  )
}

# Save data files in various formats to the analysis results folder
save_data <- function(data, filename, format = "csv") {
  
  # Construct full file path using global analysis folder
  full_path <- here::here("data/results", ANALYSIS_FOLDER, paste0(filename, ".", format))
  
  # Save file based on specified format
  if (format == "csv") {
    readr::write_csv(data, full_path)
  } else if (format == "tsv") {
    readr::write_tsv(data, full_path)
  } else if (format == "rds") {
    saveRDS(data, full_path)
  } else if (format == "txt") {
    readr::write_delim(data, full_path, delim = " ")
  } else if (format == "xlsx") {
    openxlsx::write.xlsx(data, full_path)
  } else {
    stop("Unsupported file format")
  }
  
  # Return data invisibly for piping
  return(invisible(data))
}

# Save ggplot objects as SVG and PDF to the figures folder
save_plot <- function(plot, filename, width = 10, height = 8, dpi = 300) {
  
  # Define figures directory path
  figures_path <- here::here("data/results", ANALYSIS_FOLDER, "00_figures")
  
  # Save as scalable SVG vector format
  svg_path <- file.path(figures_path, paste0(filename, ".svg"))
  ggplot2::ggsave(svg_path, plot, width = width, height = height, 
                  device = "svg")
  
  # Save as vector PDF for publications
  pdf_path <- file.path(figures_path, paste0(filename, ".pdf"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = "pdf")
  
  # Return plot invisibly for piping
  return(invisible(plot))
}