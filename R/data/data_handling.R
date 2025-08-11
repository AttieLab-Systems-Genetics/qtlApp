# =============================================================================
# DATA HANDLING AND VALIDATION MODULE
# =============================================================================
# Comprehensive utility module providing robust data validation, caching,
# formatting, and file operations for the QTL mapping application.
# 
# This module ensures data integrity and provides defensive programming
# patterns throughout the application by validating inputs, managing caches,
# and handling file operations safely.

# =============================================================================
# CACHING INFRASTRUCTURE
# =============================================================================
# Functions for creating and managing in-memory caches to improve app performance

#' Create a new environment for caching
#' 
#' Creates an isolated environment for storing cached data. Environments
#' provide efficient key-value storage and automatic garbage collection.
#' 
#' @param parent The parent environment (default: emptyenv())
#'   - emptyenv(): Creates isolated cache with no parent scope
#'   - globalenv(): Cache can access global variables (use with caution)
#' @return A new environment object for caching operations
#' @examples
#' # Create isolated cache
#' peak_cache <- create_cache_env()
#' 
#' # Store and retrieve cached data
#' peak_cache$my_key <- expensive_computation()
#' result <- peak_cache$my_key
create_cache_env <- function(parent = emptyenv()) {
  new.env(parent = parent)
}

#' Clear all caches
#' 
#' Safely removes all objects from multiple cache environments.
#' Used for memory management and ensuring fresh data after updates.
#' 
#' @param ... Variable number of environment objects to clear
#' @examples
#' # Clear multiple caches at once
#' clear_all_caches(peak_cache, trait_cache, plot_cache)
clear_all_caches <- function(...) {
  envs <- list(...)
  for (env in envs) {
    if (is.environment(env)) {
      # Remove all objects from this environment
      rm(list = ls(envir = env), envir = env)
    }
  }
}

#' Create cache key from multiple components
#' 
#' Generates a consistent string key for cache storage by combining
#' multiple parameters. Handles NULL values and ensures reproducible keys.
#' 
#' @param ... Components to combine into a cache key (trait, dataset, chromosome, etc.)
#' @return Cache key string suitable for environment storage
#' @examples
#' # Create key for specific analysis
#' key <- create_cache_key("insulin", "clinical_study", "chr2", "LOD_6")
#' # Result: "insulin_clinical_study_chr2_LOD_6"
create_cache_key <- function(...) {
  components <- list(...)
  # Convert all components to character and remove NULL values
  components <- lapply(components[!sapply(components, is.null)], as.character)
  # Combine with underscore separator
  paste(unlist(components), collapse = "_")
}

# =============================================================================
# INPUT VALIDATION FUNCTIONS
# =============================================================================
# Comprehensive validation functions ensuring data integrity and user input safety

#' Validate trait name against available options
#' 
#' Ensures the selected trait exists in the dataset and provides meaningful
#' error messages for debugging. Critical for preventing downstream errors.
#' 
#' @param trait The trait name to validate (user input from UI)
#' @param gene_symbols Vector of valid gene symbols from annotation data
#' @return Validated trait name (unchanged if valid)
#' @throws Error if trait is empty; Warning if trait not found in symbols
#' @examples
#' # Validate gene selection
#' valid_trait <- validate_trait("Actb", available_genes)
validate_trait <- function(trait, gene_symbols) {
  # Check for empty or NULL input
  if (is.null(trait) || nchar(trait) == 0) {
    stop("Trait name cannot be empty")
  }
  
  # Warn if trait not found (but don't stop - might be valid in some contexts)
  if (!trait %in% gene_symbols) {
    warning(sprintf("Trait '%s' not found in gene symbols", trait))
  }
  
  return(trait)
}

#' Validate dataset selection against file directory
#' 
#' Ensures the selected dataset exists in the imported file directory.
#' Prevents attempts to analyze non-existent datasets.
#' 
#' @param dataset The selected dataset identifier (from UI dropdown)
#' @param file_directory Data frame containing valid datasets with 'group' column
#' @return Validated dataset identifier
#' @throws Error if dataset is empty or not found in directory
#' @examples


#' Validate LOD score threshold within reasonable bounds
#' 
#' Ensures LOD threshold is numeric and within scientifically reasonable
#' ranges for QTL analysis. Prevents extreme values that could cause issues.
#' 
#' @param threshold The LOD threshold value (from UI slider/input)
#' @param min Minimum allowed value (default: 4 - typical significance threshold)
#' @param max Maximum allowed value (default: 120 - computational limit)
#' @return Validated threshold value
#' @throws Error if threshold is non-numeric or outside bounds
#' @examples


#' Validate chromosome selection
#' 
#' Ensures chromosome selection is valid for mouse genome analysis.
#' Handles both "All" option and specific chromosome identifiers.
#' 
#' @param chr The selected chromosome (from UI dropdown)
#' @param valid_chr Vector of valid chromosome values (autosomes + sex + mito)
#' @return Validated chromosome selection
#' @throws Error if chromosome is invalid; Returns "All" if empty
#' @examples


# =============================================================================
# DATA FORMATTING AND DISPLAY UTILITIES
# =============================================================================
# Functions for preparing data for display and creating user-friendly output

## Removed: format_peak_info()
## Rationale: superseded by `peakInfoModule.R` which renders peak details and
## founder effects. No call sites remained in the codebase.

#' Safe number formatting with error handling
#' 
#' Robustly formats numeric values for display, handling edge cases
#' like NULL, NA, or non-numeric inputs without throwing errors.
#' 
#' @param x Number to format (can be NULL, NA, or non-numeric)
#' @param digits Number of digits after decimal point (default: 2)
#' @return Formatted number string, or NA if input invalid
#' @examples
#' safe_number_format(3.14159, 2)    # "3.14"
#' safe_number_format(NULL)          # NA
#' safe_number_format("not_number")  # NA
safe_number_format <- function(x, digits = 2) {
  # Handle edge cases gracefully
  if (is.null(x) || !is.numeric(x) || is.na(x)) {
    return(NA)
  }
  # Format with specified decimal places
  format(round(x, digits), nsmall = digits)
}

#' Create HTML formatted message for UI display
#' 
#' Generates styled HTML messages for different types of user feedback.
#' Provides consistent styling across the application.
#' 
#' @param message Error message text to display
#' @param type Type of message affecting color and styling:
#'   - "error": Red, bold text for critical issues
#'   - "warning": Orange text for cautions  
#'   - "info": Blue text for informational messages
#' @return HTML object suitable for Shiny UI display
#' @examples
#' # Create error message for UI
#' error_msg <- create_message("Invalid trait selection", "error")
#' # Create info message
#' info_msg <- create_message("Analysis complete", "info")
#' @importFrom htmltools HTML
create_message <- function(message, type = "error") {
  # Define color scheme for different message types
  color <- switch(type,
    error = "#e74c3c",    # Red for errors
    warning = "#f39c12",  # Orange for warnings
    info = "#3498db",     # Blue for info
    "#2c3e50"            # Dark gray default
  )
  
  # Generate styled HTML div
  HTML(sprintf(
    '<div style="color: %s; padding: 10px; margin: 10px 0; font-weight: %s;">%s</div>',
    color,
    if(type == "error") "bold" else "normal",
    message
  ))
}

# =============================================================================
# PLOT AND EXPORT UTILITIES
# =============================================================================
# Functions for managing plot dimensions, downloads, and file operations

#' Validate plot dimensions within reasonable bounds
#' 
#' Ensures plot dimensions are numeric and within practical limits for
#' both display and file export. Provides sensible defaults for invalid inputs.
#' 
#' @param width Plot width in pixels
#' @param height Plot height in pixels  
#' @param min_width Minimum width (default: 400px - readable minimum)
#' @param max_width Maximum width (default: 2000px - practical limit)
#' @param min_height Minimum height (default: 300px - readable minimum)
#' @param max_height Maximum height (default: 1200px - practical limit)
#' @return Named list with validated width and height values
#' @examples
#' # Validate user-specified dimensions
#' dims <- validate_plot_dimensions(1200, 800)  # Valid - returns as-is
#' dims <- validate_plot_dimensions(50, 2500)   # Invalid - returns defaults
validate_plot_dimensions <- function(width, height,
                                   min_width = 400, max_width = 2000,
                                   min_height = 300, max_height = 1200) {
  # Validate width with fallback to default
  if (!is.numeric(width) || width < min_width || width > max_width) {
    width <- 1000  # Default width for most plots
  }
  
  # Validate height with fallback to default
  if (!is.numeric(height) || height < min_height || height > max_height) {
    height <- 600  # Default height for most plots
  }
  
  list(width = width, height = height)
}

#' Get preset plot dimensions based on common aspect ratios
#' 
#' Provides standardized plot dimensions for common use cases.
#' Ensures consistent appearance across different plot types.
#' 
#' @param preset Preset aspect ratio name:
#'   - "1:1": Square plots (ideal for correlation matrices)
#'   - "3:2": Classic photo ratio (good for most scientific plots)
#'   - "16:9": Widescreen ratio (good for timeseries, genomic plots)
#' @param base_size Base size for calculation (width for non-square ratios)
#' @return Named list with width and height values
#' @examples
#' # Get dimensions for different plot types
#' square_dims <- get_preset_dimensions("1:1", 600)      # 600x600
#' wide_dims <- get_preset_dimensions("16:9", 800)       # 800x450
#' classic_dims <- get_preset_dimensions("3:2", 900)     # 900x600
get_preset_dimensions <- function(preset, base_size = 800) {
  switch(preset,
    "1:1" = list(width = base_size, height = base_size),
    "3:2" = list(width = base_size, height = round(base_size * 2/3)),
    "16:9" = list(width = base_size, height = round(base_size * 9/16)),
    list(width = base_size, height = round(base_size * 0.6))  # Default to 3:2-ish
  )
}

#' Create standardized download filename
#' 
#' Generates consistent, descriptive filenames for exported files.
#' Includes relevant analysis parameters and timestamps for organization.
#' 
#' @param prefix Filename prefix (usually plot type: "manhattan", "effect", etc.)
#' @param trait Trait name being analyzed
#' @param chr Chromosome identifier (optional, excluded if "All")
#' @param ext File extension (default: "png")
#' @return Formatted filename string
#' @examples
#' # Create filename for Manhattan plot
#' filename <- create_download_filename("manhattan", "insulin", "chr2", "png")
#' # Result: "manhattan_insulin_chr2_20241205.png"
#' 
#' # Genome-wide plot (no chromosome specified)
#' filename <- create_download_filename("manhattan", "insulin", "All", "pdf")
#' # Result: "manhattan_insulin_20241205.pdf"
create_download_filename <- function(prefix, trait, chr = NULL, ext = "png") {
  components <- c(
    prefix,
    trait,
    # Only include chromosome if it's specific (not "All")
    if(!is.null(chr) && chr != "All") paste0("chr", chr),
    format(Sys.time(), "%Y%m%d")  # Add date stamp
  )
  paste0(paste(components, collapse = "_"), ".", ext)
}

# =============================================================================
# FILE SYSTEM UTILITIES
# =============================================================================
# Safe file operations with comprehensive error handling

#' Safe file path joining with validation
#' 
#' Combines path components while handling NULL, NA, and empty values.
#' Prevents path construction errors and ensures cross-platform compatibility.
#' 
#' @param ... Path components to join (directories, filenames, extensions)
#' @return Properly joined file path string
#' @examples
#' # Safe path construction
#' path <- safe_file_path("/data", "study1", NULL, "peaks.fst")
#' # Result: "/data/study1/peaks.fst" (NULL component ignored)
safe_file_path <- function(...) {
  components <- list(...)
  # Remove NULL, NA, and empty strings to prevent path issues
  components <- components[!sapply(components, function(x) is.null(x) || is.na(x) || x == "")]
  # Join with appropriate separator (handles Windows vs Unix)
  do.call(file.path, components)
}

#' Check if file exists and is readable
#' 
#' Safely verifies file accessibility without throwing errors.
#' Essential for defensive programming when dealing with user-specified files.
#' 
#' @param path File path to check
#' @return TRUE if file exists and is readable, FALSE otherwise
#' @examples
#' # Check before attempting to read
#' if (check_file_accessible("/path/to/data.fst")) {
#'   data <- fst::read_fst("/path/to/data.fst")
#' }
check_file_accessible <- function(path) {
  tryCatch({
    # Check both existence and read permissions (mode = 4)
    file.exists(path) && file.access(path, mode = 4) == 0
  }, error = function(e) {
    FALSE  # Return FALSE for any file system errors
  })
}

#' Safe read of CSV/FST file with comprehensive error handling
#' 
#' Robustly reads data files with automatic format detection and error recovery.
#' Supports both CSV and FST formats commonly used in genomics applications.
#' 
#' @param path File path to read
#' @param as_dt Convert to data.table format (default: TRUE for performance)
#' @return Data frame/data.table object, or NULL on error
#' @examples
#' # Safe file reading with error handling
#' peak_data <- safe_read_file("/data/peaks.fst")
#' if (is.null(peak_data)) {
#'   stop("Failed to load peak data")
#' }
#' 
#' # Read CSV as regular data frame
#' annotations <- safe_read_file("/data/genes.csv", as_dt = FALSE)
safe_read_file <- function(path, as_dt = TRUE) {
  tryCatch({
    # Pre-flight check for file accessibility
    if (!check_file_accessible(path)) {
      warning(sprintf("File not accessible: %s", path))
      return(NULL)
    }
    
    # Format-specific reading with appropriate packages
    if (grepl("\\.fst$", path)) {
      # FST format - fast binary format for R data frames
      fst::read_fst(path, as.data.table = as_dt)
    } else if (grepl("\\.csv$", path)) {
      # CSV format - use fast reader or base R depending on preference
      if (as_dt) {
        data.table::fread(path)  # Fast CSV reader
      } else {
        read.csv(path)           # Base R CSV reader
      }
    } else {
      warning(sprintf("Unsupported file type: %s", path))
      NULL
    }
  }, error = function(e) {
    # Comprehensive error logging for debugging
    warning(sprintf("Error reading file %s: %s", path, e$message))
    NULL
  })
}

# =============================================================================
# MODULE SUMMARY
# =============================================================================
# This module provides the foundational utilities for:
# 
# 1. **Caching**: Performance optimization through environment-based storage
# 2. **Validation**: Comprehensive input checking for all user interactions  
# 3. **Formatting**: Consistent data presentation and number formatting
# 4. **File Operations**: Safe, robust file system interactions
# 5. **Plot Management**: Standardized plotting dimensions and export filenames
# 6. **Error Handling**: Graceful degradation and informative error messages
#
# Key Design Principles:
# - **Defensive Programming**: Every function handles edge cases gracefully
# - **Consistent Interface**: Similar parameter patterns across functions
# - **User-Friendly**: Clear error messages and sensible defaults
# - **Performance-Aware**: Efficient caching and file operations
# - **Cross-Platform**: Compatible path handling and file operations
#
# =============================================================================