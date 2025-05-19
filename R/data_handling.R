# Data handling and validation module


#' Create a new environment for caching
#' @param parent The parent environment (default: emptyenv())
#' @return A new environment for caching
create_cache_env <- function(parent = emptyenv()) {
  new.env(parent = parent)
}

#' Clear all caches
#' @param ... List of environments to clear
clear_all_caches <- function(...) {
  envs <- list(...)
  for (env in envs) {
    if (is.environment(env)) {
      rm(list = ls(envir = env), envir = env)
    }
  }
}

#' Validate trait name
#' @param trait The trait name to validate
#' @param gene_symbols Vector of valid gene symbols
#' @return Validated trait name or error
validate_trait <- function(trait, gene_symbols) {
  if (is.null(trait) || nchar(trait) == 0) {
    stop("Trait name cannot be empty")
  }
  
  if (!trait %in% gene_symbols) {
    warning(sprintf("Trait '%s' not found in gene symbols", trait))
  }
  
  return(trait)
}

#' Validate dataset selection
#' @param dataset The selected dataset
#' @param file_directory Data frame containing valid datasets
#' @return Validated dataset or error
validate_dataset <- function(dataset, file_directory) {
  if (is.null(dataset) || nchar(dataset) == 0) {
    stop("Dataset cannot be empty")
  }
  
  if (!dataset %in% unique(file_directory$group)) {
    stop(sprintf("Dataset '%s' not found", dataset))
  }
  
  return(dataset)
}

#' Validate LOD threshold
#' @param threshold The LOD threshold value
#' @param min Minimum allowed value (default: 4)
#' @param max Maximum allowed value (default: 120)
#' @return Validated threshold or error
validate_lod_threshold <- function(threshold, min = 4, max = 120) {
  if (is.null(threshold) || !is.numeric(threshold)) {
    stop("LOD threshold must be a numeric value")
  }
  
  if (threshold < min || threshold > max) {
    stop(sprintf("LOD threshold must be between %d and %d", min, max))
  }
  
  return(threshold)
}

#' Validate chromosome selection
#' @param chr The selected chromosome
#' @param valid_chr Vector of valid chromosome values
#' @return Validated chromosome or error
validate_chromosome <- function(chr, valid_chr = c("All", as.character(1:19), "X", "Y", "M")) {
  if (is.null(chr) || nchar(chr) == 0) {
    return("All")
  }
  
  if (!chr %in% valid_chr) {
    stop(sprintf("Invalid chromosome: %s", chr))
  }
  
  return(chr)
}

#' Convert chromosome label to numeric
#' @param chr The chromosome label
#' @return Numeric chromosome value
chr_to_numeric <- function(chr) {
  if (chr == "X") return(20)
  if (chr == "Y") return(21)
  if (chr == "M") return(22)
  return(as.numeric(chr))
}

#' Convert numeric chromosome to label
#' @param chr_num The numeric chromosome value
#' @return Chromosome label
numeric_to_chr <- function(chr_num) {
  if (chr_num == 20) return("X")
  if (chr_num == 21) return("Y")
  if (chr_num == 22) return("M")
  return(as.character(chr_num))
}

#' Format peak information
#' @param peak Peak data frame
#' @return Formatted peak information
format_peak_info <- function(peak) {
  if (is.null(peak) || nrow(peak) == 0) {
    return(NULL)
  }
  
  # Format chromosome
  chr_label <- if(peak$chr %in% c(20,21,22)) {
    c("X","Y","M")[peak$chr-19]
  } else {
    peak$chr
  }
  
  # Create formatted info
  info <- list(
    marker = peak$marker,
    chromosome = chr_label,
    position = round(peak$pos, 2),
    lod = round(peak$lod, 2)
  )
  
  # Add strain effects if available
  strain_cols <- c("A", "B", "C", "D", "E", "F", "G", "H")
  if (all(strain_cols %in% colnames(peak))) {
    info$strain_effects <- list(
      AJ = round(peak$A, 3),
      B6 = round(peak$B, 3),
      `129` = round(peak$C, 3),
      NOD = round(peak$D, 3),
      NZO = round(peak$E, 3),
      CAST = round(peak$F, 3),
      PWK = round(peak$G, 3),
      WSB = round(peak$H, 3)
    )
  }
  
  return(info)
}

#' Create cache key
#' @param ... Components to combine into a cache key
#' @return Cache key string
create_cache_key <- function(...) {
  components <- list(...)
  # Convert all components to character and remove NULL values
  components <- lapply(components[!sapply(components, is.null)], as.character)
  # Combine with underscore
  paste(unlist(components), collapse = "_")
}

#' Safe number formatting
#' @param x Number to format
#' @param digits Number of digits after decimal point
#' @return Formatted number string
safe_number_format <- function(x, digits = 2) {
  if (is.null(x) || !is.numeric(x) || is.na(x)) {
    return(NA)
  }
  format(round(x, digits), nsmall = digits)
}

#' Create error message
#' @param message Error message text
#' @param type Type of message ("error", "warning", or "info")
#' @return HTML formatted message
create_message <- function(message, type = "error") {
  color <- switch(type,
    error = "#e74c3c",
    warning = "#f39c12",
    info = "#3498db",
    "#2c3e50"
  )
  
  HTML(sprintf(
    '<div style="color: %s; padding: 10px; margin: 10px 0; font-weight: %s;">%s</div>',
    color,
    if(type == "error") "bold" else "normal",
    message
  ))
}

#' Validate plot dimensions
#' @param width Plot width
#' @param height Plot height
#' @param min_width Minimum width (default: 400)
#' @param max_width Maximum width (default: 2000)
#' @param min_height Minimum height (default: 300)
#' @param max_height Maximum height (default: 1200)
#' @return List of validated width and height
validate_plot_dimensions <- function(width, height,
                                   min_width = 400, max_width = 2000,
                                   min_height = 300, max_height = 1200) {
  # Validate width
  if (!is.numeric(width) || width < min_width || width > max_width) {
    width <- 1000  # Default width
  }
  
  # Validate height
  if (!is.numeric(height) || height < min_height || height > max_height) {
    height <- 600  # Default height
  }
  
  list(width = width, height = height)
}

#' Get preset dimensions
#' @param preset Preset name ("1:1", "3:2", or "16:9")
#' @param base_size Base size for calculation
#' @return List of width and height
get_preset_dimensions <- function(preset, base_size = 800) {
  switch(preset,
    "1:1" = list(width = base_size, height = base_size),
    "3:2" = list(width = base_size, height = round(base_size * 2/3)),
    "16:9" = list(width = base_size, height = round(base_size * 9/16)),
    list(width = base_size, height = round(base_size * 0.6))  # Default to 3:2
  )
}

#' Create download filename
#' @param prefix Filename prefix
#' @param trait Trait name
#' @param chr Chromosome (optional)
#' @param ext File extension
#' @return Formatted filename
create_download_filename <- function(prefix, trait, chr = NULL, ext = "png") {
  components <- c(
    prefix,
    trait,
    if(!is.null(chr) && chr != "All") paste0("chr", chr),
    format(Sys.time(), "%Y%m%d")
  )
  paste0(paste(components, collapse = "_"), ".", ext)
}

#' Safe file path joining
#' @param ... Path components to join
#' @return Joined file path
safe_file_path <- function(...) {
  components <- list(...)
  # Remove NULL, NA, and empty strings
  components <- components[!sapply(components, function(x) is.null(x) || is.na(x) || x == "")]
  # Join with appropriate separator
  do.call(file.path, components)
}


#' @param path File path to check
#' @return TRUE if file exists and is readable, FALSE otherwise
check_file_accessible <- function(path) {
  tryCatch({
    file.exists(path) && file.access(path, mode = 4) == 0
  }, error = function(e) {
    FALSE
  })
}

#' Safe read of CSV/FST file
#' @param path File path
#' @param as_dt Convert to data.table (default: TRUE)
#' @return Data frame or NULL on error
safe_read_file <- function(path, as_dt = TRUE) {
  tryCatch({
    if (!check_file_accessible(path)) {
      warning(sprintf("File not accessible: %s", path))
      return(NULL)
    }
    
    if (grepl("\\.fst$", path)) {
      fst::read_fst(path, as.data.table = as_dt)
    } else if (grepl("\\.csv$", path)) {
      if (as_dt) {
        data.table::fread(path)
      } else {
        read.csv(path)
      }
    } else {
      warning(sprintf("Unsupported file type: %s", path))
      NULL
    }
  }, error = function(e) {
    warning(sprintf("Error reading file %s: %s", path, e$message))
    NULL
  })
} 