#' @export
create_cache <- function(caches = c("peaks", "trait")) {
  
}
#' @export
get_trait_choices <- function(import, selected_dataset = NULL) {
  trait_type <- get_trait_type(import, selected_dataset)
  if(is.null(trait_type)) {
    return(NULL)
  }
  trait_id <- get_trait_id(trait_type)
  trait_list = get_trait_list(import, trait_type)
  
  annotation_list <- import$annotation_list
  if(!(trait_type %in% c("genes", "isoforms"))) {
    choices <- annotation_list[[trait_type]][[trait_id]]
  } else { # "genes", "isoforms"
    choices <- join_symbol_id(annotation_list, trait_type, trait_id)
  }
  choices
}
#' @importFrom stringr str_remove
#' @export
join_symbol_id <- function(annotation_list,
                           trait_type = c("genes","isoforms"),
                           trait_id = get_trait_id(trait_type)) {
  
  trait_type <- match.arg(trait_type)
  
  id <- stringr::str_remove(annotation_list[[trait_type]][[trait_id]],
                            pattern = "ENSMUS[GT]0+")
  paste(annotation_list[[trait_type]]$symbol, id, sep = "_")
}
#' @importFrom dplyr arrange desc filter
#' @importFrom rlang .data
#' @export
highest_peaks <- function(peak_table, LOD_thr) {
  filtered_peaks <- dplyr::filter(peak_table, .data$lod >= LOD_thr) |>
    dplyr::arrange(dplyr::desc(.data$lod))
  if (nrow(filtered_peaks) == 0) return(NULL)
  filtered_peaks
}
#' @importFrom stringr str_split
#' @export
get_selected_trait <- function(import, which_trait, selected_dataset) {
 
  trait_type <- get_trait_type(import, selected_dataset)
  if(trait_type %in% c("genes", "isoforms")) {
    which_trait <- stringr::str_split(which_trait, pattern="_")[[1]][1]
  }
  which_trait
}
#' @importFrom reshape2 melt
#' @export
pivot_peaks <- function(peaks, which_peak) {
  
  if (peaks$intcovar[1] == "none") {
    # set peaks
    peak <- subset(peaks, marker == which_peak)  
    # Check if we have data after filtering
    if (nrow(peak) > 0) {
      # Select and rename columns
      peak <- peak[c("marker","A","B","C","D","E","F","G","H")]
      colnames(peak)[2:9] <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
      # Reshape data
      peak <- reshape2::melt(peak, id.vars = "marker")
    }
  } else {
    
    peak <- NULL
  }
  return(peak)
}
# internal helper functions
get_trait_type <- function(import, selected_dataset = NULL) {
  if(is.null(selected_dataset)) {
    return(NULL)
  }
  file_directory <- import$file_directory
  # Filter for the selected dataset
  file_directory_subset <- subset(file_directory, group == selected_dataset)
  
  if(nrow(file_directory_subset) == 0) {
    warning(paste("No entry found for dataset:", selected_dataset, "in file directory. Cannot determine trait_type."))
    return(NULL)
  }
  
  # Get trait_type from the first entry of the subset and convert to lowercase
  trait_type_from_file <- tolower(file_directory_subset$trait_type[1])
  
  # Standardize "clinical traits" to "clinical"
  if (trait_type_from_file == "clinical traits") {
    return("clinical")
  }
  
  return(trait_type_from_file)
}
get_trait_list <- function(import, trait_type) {
  annotation_list <- import$annotation_list
  annotation_list[[trait_type]]
}
get_trait_id <- function(trait_type) {
  switch(trait_type,
    genes    = "gene.id",
    isoforms = "transcript.id",
    "data_name")
}
chr_XYM <- function(chr) {
  if(is.null(chr)) return(NA)
  
  if (chr %in% c(20, 21, 22)) {
    c("X", "Y", "M")[chr - 19]
  } else {
    chr
  }
}
