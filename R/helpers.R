#' @export
create_cache <- function(caches = c("peaks", "trait")) {
  for(cache_name in caches) {
    assign(paste0(cache_name, "_cache"), new.env(parent = emptyenv()),
      pos = globalenv())
  }
}
#' @export
get_trait_choices <- function(import, selected_group) {
  trait_type <- get_trait_type(import, selected_group)
  trait_id <- get_trait_id(trait_type)
  trait_list = get_trait_list(import, trait_type)
  # Choices come from elements of `import$annotation_list`
  annotation_list <- import$annotation_list
  if(!(trait_type %in% c("genes", "isoforms"))) {
    choices <- annotation_list[[trait_type]][[trait_id]]
  } else { # "genes", "isoforms"
    # Remove "ENSMUSG" (genes) or "ENSMUST" (transcripts)
    # and leading zeros from the trait ID
    # and create choices by combining the symbol and simplified ID
    # e.g. "A1bg_1".
    id <- stringr::str_remove(annotation_list[[trait_type]][[trait_id]],
      pattern = "ENSMUS[GT]0+")
    choices <- paste(annotation_list[[trait_type]]$symbol, id, sep = "_")
  }
  choices
}
#' @export
get_selected_trait <- function(which_trait, selected_group) {
  # Split the trait string by "_" to separate the symbol and ID.
  # Useful for gene or transcript names and IDs.
  # Caution: symbols may not be unique.
  # Caution: some clinical or other traits may have "_" in name.
  trait_type <- get_trait_type(import, selected_group)
  if(trait_type %in% c("genes", "isoforms")) {
    which_trait <- stringr::str_split(which_trait, pattern="_")[[1]][1]
  }
  which_trait
}
# internal helper functions
get_trait_type <- function(import, selected_group) {
  file_directory <- import$file_directory
  file_directory <- subset(file_directory, group == selected_group)
  trait_type <- tolower(file_directory$trait_type[1])
  trait_type
}
get_trait_list <- function(import, trait_type) {
  annotation_list <- import$annotation_list
  annotation_list[[trait_type]]
}
get_trait_id <- function(trait_type) {
  switch(trait_type, # ** Expand for new data types ** #
    genes    = "gene.id",
    isoforms = "transcript.id",
    clinical = "data_name")
}
