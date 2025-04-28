#' Find trait scan
#' 
#' @param file_dir data frame with file directory information
#' @param selected_dataset character string
#' @param selected_trait character string
#' 
#' @importFrom fst read_fst
#' @importFrom stringr str_detect str_replace
#' @importFrom data.table rbindlist setnames data.table
#' @export
trait_scan <- function(file_dir, selected_dataset, selected_trait) {
   # Create a more unique cache key that includes the full dataset name
   cache_key <- paste(selected_dataset, tolower(selected_trait), sep = "_")
   
   # Check if trait_cache exists, create if not
   if (!exists("trait_cache", envir = .GlobalEnv)) {
     assign("trait_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
     message("Created trait_cache environment.")
   }
   
   if (!is.null(trait_cache[[cache_key]])) {
       message("Using cached data for trait: ", selected_trait, " in dataset: ", selected_dataset)
       return(trait_cache[[cache_key]])
   }
  
   message("Searching for trait: ", selected_trait, " in dataset: ", selected_dataset)
  
   # Filter the file directory for the selected dataset and scan type
   file_dir <- subset(file_dir, group == selected_dataset & file_type == "scans")
  
   # Check if we found any matching files
   if (nrow(file_dir) == 0) {
       stop("No matching files found for the selected dataset: ", selected_dataset)
   }
  
   # Initialize list to store data from each file
   all_data <- list()
  
   # Process each FST file (one per chromosome)
   for (i in 1:nrow(file_dir)) {
       # Extract chromosome number from ID_code
       chr_num <- file_dir$ID_code[i]
       fst_path <- file_dir$File_path[i]
      
       if (!file.exists(fst_path)) {
           warning("File does not exist: ", fst_path)
           next
       }
      
       # Ensure we are working with an FST file
       if (!stringr::str_detect(fst_path, "fst$")) {
          warning("Path does not seem to be an FST file: ", fst_path, ". Skipping.")
          next # Skip if not an fst file path
       }
       
       message("Checking chromosome ", chr_num, " for trait: ", selected_trait, " in dataset: ", selected_dataset)
       
       # Create row index if it doesn't exist or get the path if it does
       row_index_path <- fst_rows(fst_path)
       
       # Check if row index creation was successful
       if (is.null(row_index_path)) {
         warning("Skipping chromosome ", chr_num, " due to missing or failed row index creation.")
         next # Skip to the next chromosome
       }
       
       tryCatch({
           # Read the row index to find the trait
           trait_index <- fst::read_fst(row_index_path, as.data.table = TRUE)
           trait_index[, Phenotype := tolower(Phenotype)]
           trait_rows <- trait_index[Phenotype == tolower(selected_trait), ]
          
           if (nrow(trait_rows) > 0) {
               message("Found trait in index for chromosome ", chr_num, " at rows ", trait_rows$from, "-", trait_rows$to)
               data <- fst::read_fst(fst_path,
                                    from = trait_rows$from,
                                    to = trait_rows$to,
                                    as.data.table = TRUE)
               
               # Standardize columns (ensure LOD, marker exist)
               if (!"LOD" %in% colnames(data)) {
                  possible_lod_cols <- grep("lod|LOD|score", colnames(data), ignore.case = TRUE, value = TRUE)
                  if (length(possible_lod_cols) > 0) { data.table::setnames(data, possible_lod_cols[1], "LOD") }
               }
               if (!"marker" %in% colnames(data)) {
                  possible_marker_cols <- grep("marker|id|snp", colnames(data), ignore.case = TRUE, value = TRUE)
                  if (length(possible_marker_cols) > 0) { data.table::setnames(data, possible_marker_cols[1], "marker") }
               }
               
               # Verify trait match and column presence
               if ("Phenotype" %in% colnames(data) && "LOD" %in% colnames(data) && "marker" %in% colnames(data)) {
                   data <- data[tolower(Phenotype) == tolower(selected_trait)]
                   if (nrow(data) > 0) {
                      message("Adding ", nrow(data), " rows from chromosome ", chr_num, " (indexed read)")
                      all_data[[length(all_data) + 1]] <- data
                   }
               } else {
                 warning("Required columns (Phenotype, LOD, marker) missing after indexed read for ", fst_path)
               }
           } else {
               message("Trait not found in index for chromosome ", chr_num)
           }
       }, error = function(e) {
          warning("Error processing chromosome ", chr_num, " using index: ", e$message)
          # Potentially fallback to full read here if index read fails unexpectedly
          row_index_path <- NULL # Mark index as unusable for this file
       })
       
       # If data was found via index, skip to next iteration
       if (!is.null(row_index_path)) next 
       
       # --- Fallback to FULL read if index failed or wasn't available --- 
       if (is.null(row_index_path)) { 
          message("Index reading failed or unavailable for chr ", chr_num, ". Attempting full file read.")
          tryCatch({
             data <- fst::read_fst(fst_path, as.data.table = TRUE)
             
             # Standardize columns
              if (!"LOD" %in% colnames(data)) {
                 possible_lod_cols <- grep("lod|LOD|score", colnames(data), ignore.case = TRUE, value = TRUE)
                 if (length(possible_lod_cols) > 0) { data.table::setnames(data, possible_lod_cols[1], "LOD") }
              }
              if (!"marker" %in% colnames(data)) {
                 possible_marker_cols <- grep("marker|id|snp", colnames(data), ignore.case = TRUE, value = TRUE)
                 if (length(possible_marker_cols) > 0) { data.table::setnames(data, possible_marker_cols[1], "marker") }
              }
              
             # Filter and add data
             if ("Phenotype" %in% colnames(data) && "LOD" %in% colnames(data) && "marker" %in% colnames(data)) {
                 data <- data[tolower(Phenotype) == tolower(selected_trait)]
                 if (nrow(data) > 0) {
                     message("Adding ", nrow(data), " rows from chromosome ", chr_num, " (full read)")
                     all_data[[length(all_data) + 1]] <- data
                 } else {
                     message("Trait not found in chromosome ", chr_num, " (full read)")
                 }
             } else {
                warning("Required columns missing after full read for ", fst_path)
             }
          }, error = function(e) {
             warning("Error reading full FST file for chromosome ", chr_num, ": ", e$message)
          })
       } # End fallback full read
       
   } # End chromosome loop
  
   # Check if we found any data
   if (length(all_data) == 0) {
       stop("Trait '", selected_trait, "' not found in any chromosome for dataset: ", selected_dataset)
   }
  
   # Combine all data
   combined_data <- data.table::rbindlist(all_data, fill = TRUE)
   message("Total rows in combined data: ", nrow(combined_data))
  
   # Cache the result
   trait_cache[[cache_key]] <- combined_data
   return(combined_data)
} 