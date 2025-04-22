#' @export
create_cache <- function(caches = c("peaks", "trait")) {
  for(cache_name in caches) {
    assign(paste0(cache_name, "_cache"), new.env(parent = emptyenv()),
      pos = globalenv())
  }
}
#' @export
get_trait_choices <- function(import, selected_dataset = NULL) {
  trait_type <- get_trait_type(import, selected_dataset)
  if(is.null(trait_type)) {
    return(NULL)
  }
  trait_id <- get_trait_id(trait_type)
  trait_list = get_trait_list(import, trait_type)
  # Choices come from elements of `import$annotation_list`
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
  # `trait_type` should be "genes" or "isoforms"
  # trait_id should be "gene.id" or "transcript.id"
  trait_type <- match.arg(trait_type)
  # Remove "ENSMUSG" (genes) or "ENSMUST" (transcripts)
  # and leading zeros from the trait ID
  # and create choices by combining the symbol and simplified ID
  # e.g. "A1bg_1".
  id <- stringr::str_remove(annotation_list[[trait_type]][[trait_id]],
                            pattern = "ENSMUS[GT]0+")
  paste(annotation_list[[trait_type]]$symbol, id, sep = "_")
}
#' @importFrom stringr str_split
#' @export
get_selected_trait <- function(import, which_trait, selected_dataset) {
  # Split the trait string by "_" to separate the symbol and ID.
  # Useful for gene or transcript names and IDs.
  # Caution: symbols may not be unique.
  # Caution: some clinical or other traits may have "_" in name.
  trait_type <- get_trait_type(import, selected_dataset)
  if(trait_type %in% c("genes", "isoforms")) {
    which_trait <- stringr::str_split(which_trait, pattern="_")[[1]][1]
  }
  which_trait
}
#' @importFrom reshape2 melt
#' @export
pivot_peaks <- function(peaks, which_peak) {
  # Check if scan data exists and is additive
  if (peaks$intcovar[1] == "none") {
    # set peaks
    peak <- subset(peaks, marker == which_peak)  # Changed from marker.id to marker
    # Check if we have data after filtering
    if (nrow(peak) > 0) {
      # Select and rename columns
      peak <- peak[c("marker","A","B","C","D","E","F","G","H")]
      colnames(peak)[2:9] <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
      # Reshape data
      peak <- reshape2::melt(peak, id.vars = "marker")
    }
  } else {
    # If not additive, return NULL
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
  file_directory <- subset(file_directory, group == selected_dataset)
  trait_type <- tolower(file_directory$trait_type[1])
  trait_type
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
