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
           fst_path <- stringr::str_replace(fst_path, "csv$", "fst")
           if (!file.exists(fst_path)) {
               # Check if fst_rows function exists before calling
               if (exists("fst_rows", mode = "function")) {
                  # Try creating the FST file and index if missing
                  warning("FST file not found: ", fst_path, ". Attempting to create.")
                  # Need csv2fst function defined or loaded
                  if(exists("csv2fst", mode="function")) {
                    original_csv <- stringr::str_replace(fst_path, "fst$", "csv")
                    if(file.exists(original_csv)) {
                      fst_path <- tryCatch(csv2fst(original_csv), error = function(e) NULL)
                      if(is.null(fst_path) || !file.exists(fst_path)) {
                         warning("Failed to create FST file: ", fst_path)
                         next
                      }
                    } else {
                      warning("Original CSV not found to create FST: ", original_csv)
                      next
                    }
                  } else {
                    warning("csv2fst function not available to create missing FST file.")
                    next
                  }
               } else {
                 warning("FST file not found: ", fst_path, " and fst_rows function not available.")
                 next
               }
            }
       
      
       message("Checking chromosome ", chr_num, " for trait: ", selected_trait, " in dataset: ", selected_dataset)
      
       # Check if fst_rows function exists before calling 
         # Create row index if it doesn't exist
         row_index_path <- fst_rows(fst_path) 
         tryCatch({
            data <- fst::read_fst(fst_path, as.data.table = TRUE)
            data <- data[tolower(Phenotype) == tolower(selected_trait)]
            if (nrow(data) > 0) {
                message("Adding ", nrow(data), " rows from chromosome ", chr_num, " (full read)")
                all_data[[length(all_data) + 1]] <- data
            } else {
                message("Trait not found in chromosome ", chr_num, " (full read)")
            }
         }, error = function(e) {
            warning("Error reading full FST file for chromosome ", chr_num, ": ", e$message)
         })
         next # Skip indexed reading part
       }
      
       tryCatch({
           # Read the row index to find the trait
           trait_index <- fst::read_fst(row_index_path, as.data.table = TRUE)
          
           # Convert Phenotype column to lowercase for case-insensitive matching
           trait_index[, Phenotype := tolower(Phenotype)]
          
           # Check if the trait is present in this chromosome (case-insensitive)
           trait_rows <- trait_index[Phenotype == tolower(selected_trait), ]
          
           if (nrow(trait_rows) > 0) {
               message("Found trait in chromosome ", chr_num, " at rows ", trait_rows$from, "-", trait_rows$to)
              
               # Read only the rows for this trait
               data <- fst::read_fst(fst_path,
                                    from = trait_rows$from,
                                    to = trait_rows$to,
                                    as.data.table = TRUE)
              
               # Ensure required columns are present
               if (!"LOD" %in% colnames(data)) {
                   possible_lod_cols <- grep("lod|LOD|score", colnames(data), ignore.case = TRUE, value = TRUE)
                   if (length(possible_lod_cols) > 0) {
                       data.table::setnames(data, possible_lod_cols[1], "LOD")
                   } else {
                       warning("LOD column not found in file: ", fst_path)
                       next
                   }
               }
              
               if (!"marker" %in% colnames(data)) {
                   possible_marker_cols <- grep("marker|id|snp", colnames(data), ignore.case = TRUE, value = TRUE)
                   if (length(possible_marker_cols) > 0) {
                       data.table::setnames(data, possible_marker_cols[1], "marker")
                   } else {
                       warning("marker column not found in file: ", fst_path)
                       next
                   }
               }
              
               # Verify that we have the correct trait data (case-insensitive)
               if ("Phenotype" %in% colnames(data)) {
                   # Double-check that all rows are for the requested trait
                   data <- data[tolower(Phenotype) == tolower(selected_trait)]
                   message("Verified ", nrow(data), " rows for trait: ", selected_trait)
               }
              
               if (nrow(data) > 0) {
                   message("Adding ", nrow(data), " rows from chromosome ", chr_num)
                   all_data[[length(all_data) + 1]] <- data
               }
           } else {
               message("Trait not found in chromosome ", chr_num)
           }
       }, error = function(e) {
           warning("Error processing chromosome ", chr_num, " using index: ", e$message)
       })
   }
  
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