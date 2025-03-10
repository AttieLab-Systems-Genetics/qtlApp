# load libraries=============================================================
require("rstudioapi")
require("dplyr")
require("stringr")
require("tidyverse")
require("BiocManager")
require("ggplot2")
require("qtl2")
require("grid")
require("ggrepel")
require("gridGraphics")
require("ggpubr")
require("shiny")
require("shinyFiles")
require("bslib")
require("spsComps")
require("DT")
require("shinyjs")
require("shinycssloaders")
require("data.table")
require("reshape2")


# set shiny options===========================================================
options(shiny.maxRequestSize = 20000*1024^2)  # Increase to 20GB needed for Genoprobs, etc


# load data===================================================================
file_directory <- read.csv("/data/dev/miniViewer_3.0/file_index.csv")
chr_breaks <- read.csv("/data/dev/miniViewer_3.0/chromosomal_sep_mm11.csv")
annotation_list <- readRDS("/data/dev/miniViewer_3.0/annotation_list.rds")
markers <- readRDS(file.path("/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"))


file_directory$group <- paste0(file_directory$diet, " ", file_directory$trait_compartment, " ",
                              file_directory$trait_type, ", ", file_directory$scan_type)


# Add this at the start of your script with other global variables
file_path_cache <- new.env(parent = emptyenv())


trait_scan <- function(file_dir, selected_dataset, selected_trait) {
 # Create cache key
 cache_key <- paste(selected_dataset, "path")
  # Get file path from cache or compute it
 if (is.null(file_path_cache[[cache_key]])) {
   file_dir <- subset(file_dir, group == selected_dataset)
   file_dir <- subset(file_dir, file_type == "scans")
  
   if (nrow(file_dir) == 0) {
     stop("No matching files found for the selected dataset")
   }
  
   if (!file.exists(file_dir$File_path[1])) {
     stop("File does not exist: ", file_dir$File_path[1])
   }
  
   # Cache the FST path
   fst_path <- csv2fst(file_dir$File_path[1])
   row_path <- fst_rows(fst_path)
   file_path_cache[[cache_key]] <- list(
     fst_path = fst_path,
     row_path = row_path
   )
 }
  # Use cached paths
 paths <- file_path_cache[[cache_key]]
 scan_data <- tryCatch({
   all_rows <- fst::read_fst(paths$row_path)
   rows <- dplyr::filter(all_rows, Phenotype == selected_trait)
  
   if (nrow(rows) == 0) {
     similar_traits <- all_rows$Phenotype[grep(selected_trait, all_rows$Phenotype, ignore.case = TRUE)]
     if (length(similar_traits) > 0) {
       stop("Trait '", selected_trait, "' not found exactly as written. Did you mean one of these? ",
            paste(head(unique(similar_traits), 5), collapse = ", "))
     } else {
       stop("Trait '", selected_trait, "' not found in the dataset")
     }
   }
  
   data <- fst::read_fst(paths$fst_path, from = rows$from, to = rows$to)
   if (nrow(data) == 0) {
     stop("No data found for phenotype: ", selected_trait)
   }
   data
 }, error = function(e) {
   stop("Error reading file: ", e$message)
 })
  return(scan_data)
}


# Also add caching for peak_finder
peaks_cache <- new.env(parent = emptyenv())


peak_finder <- function(file_dir, selected_dataset) {
 cache_key <- selected_dataset
  if (is.null(peaks_cache[[cache_key]])) {
   file_dir <- subset(file_dir, group == selected_dataset)
   file_dir <- subset(file_dir, file_type == "peaks")
   peaks <- read.csv(file_dir$File_path)
   peaks <- peaks %>% relocate("marker", "lodcolumn", "chr", "pos", "lod")
   colnames(peaks)[which(colnames(peaks) == "lodcolumn")] <- "trait"
   peaks_cache[[cache_key]] <- peaks
 }
 return(peaks_cache[[cache_key]])
}
# set microfunctions==========================================================
QTL_plot_visualizer <- function(qtl.temp, phenotype, LOD_thr, mrkrs) {
 # The LOD values are already in the "LOD" column
 qtl.temp <- qtl.temp[, c("marker", "LOD")]
 colnames(qtl.temp) <- c("markers", "LOD")
  mrkrs2 <- mrkrs[c("marker", "chr", "bp_grcm39")]
 colnames(mrkrs2) <- c("markers", "chr", "position")
 mrkrs2$position <- as.numeric(mrkrs2$position) / (10^6)
  # Join the data
 qtl.temp <- merge(qtl.temp, mrkrs2, by = "markers", all.x = FALSE)
 qtl.temp <- na.omit(qtl.temp)
  # Convert chr to numeric
 qtl.temp$chr <- as.character(qtl.temp$chr)
 qtl.temp$order <- sapply(qtl.temp$chr, function(x) {
   if (x == "X") return(20)
   else if (x == "Y") return(21)
   else if (x == "M") return(22)
   else return(as.numeric(x))
 })
  qtl.temp <- qtl.temp[!is.na(qtl.temp$order), ]
 qtl.temp$chr <- qtl.temp$order
  # Create plot object
 qtl_plot_obj <- qtl.temp %>%
   group_by(chr) %>%
   summarise(chr_len = max(position)) %>%
   mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
   select(-chr_len) %>%
   left_join(qtl.temp, ., by = c("chr" = "chr")) %>%
   arrange(order, position) %>%
   mutate(BPcum = position + tot)
  # Create axis labels
 axisdf <- qtl_plot_obj %>%
   group_by(order) %>%
   summarize(center = (max(BPcum) + min(BPcum)) / 2)
  axisdf$order[axisdf$order == 20] <- "X"
 axisdf$order[axisdf$order == 21] <- "Y"
 axisdf$order[axisdf$order == 22] <- "M"
 axisdf$order <- factor(axisdf$order, levels = c(as.character(1:19), "X", "Y", "M"))
  # Create plot with x-axis
 plot_QTL <- ggplot(qtl_plot_obj, aes(x = BPcum, y = LOD)) +
   geom_line(aes(color = as.factor(chr)), alpha = 0.8, size = 0.5) +
   scale_color_manual(values = rep(c("black", "darkgrey"), 22)) +
   scale_x_continuous(
     label = axisdf$order,
     breaks = axisdf$center,
     expand = expansion(mult = 0.15)
   ) +
   ylim(0, max(qtl_plot_obj$LOD, na.rm = TRUE) * 1.25) +
   theme_bw() +
   theme(
     legend.position = "none",
     panel.border = element_blank(),
     panel.grid.major.x = element_blank(),
     panel.grid.minor.x = element_blank(),
     axis.line = element_line(colour = "black"),
     axis.text = element_text(size = 18),
     axis.title = element_text(size = 20),
     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 10)),
     plot.margin = margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
   ) +
   xlab("Chromosome") +
   ylab("LOD") +
   geom_hline(aes(yintercept = LOD_thr), color = "black", linetype = "dashed")
  plots <- list(plot_QTL, qtl_plot_obj)
 return(plots)
}


csv2fst <- function(csv_path, chunk_size = 50000) {
 if (stringr::str_detect(csv_path, "csv$")) {
   fst_path <- stringr::str_replace(csv_path, "csv$", "fst")
   if (!file.exists(fst_path)) {
     warning("Writing FST file in chunks: ", fst_path)
     all_chunks <- list()
     skip <- 0
     repeat {
       chunk <- data.table::fread(csv_path, skip = skip, nrows = chunk_size, showProgress = TRUE)
       if (nrow(chunk) == 0) break
       all_chunks[[length(all_chunks) + 1]] <- chunk
       skip <- skip + chunk_size
     }
     full_data <- data.table::rbindlist(all_chunks)
     if ("Phenotype" %in% colnames(full_data)) {
       data.table::setorder(full_data, Phenotype)
     } else {
       warning("Column 'Phenotype' not found. Data will not be sorted.")
     }
     fst::write_fst(full_data, path = fst_path, compress = 50)
   }
 } else {
   if (!stringr::str_detect(csv_path, "fst$"))
     stop("No CSV or FST name provided: ", csv_path)
   fst_path <- csv_path
 }
 return(fst_path)
}


fst_rows <- function(fst_path) {
 row_path <- stringr::str_replace(fst_path, ".fst$", "_row.fst")
 if (!file.exists(row_path)) {
   rows <- fst::read_fst(fst_path) |>
     dplyr::select(Phenotype) |>
     dplyr::mutate(rown = dplyr::row_number()) |>
     dplyr::group_by(Phenotype) |>
     dplyr::slice(c(1, dplyr::n())) |>
     dplyr::mutate(set = c("from", "to")) |>
     tidyr::pivot_wider(names_from = "set", values_from = "rown")
   fst::write_fst(rows, row_path)
 }
 return(row_path)
}




# set UI=======================================================================
ui <- fluidPage(
 useShinyjs(),
 titlePanel("Pre-scanned QTL Visualizer for Diet DO Study"),  # Fixed quotation marks
 sidebarLayout(


 sidebarPanel(
   selectizeInput("selected_dataset", "Choose a dataset",
                  choices = unique(file_directory$group), multiple = FALSE),
   sliderInput("LOD_thr", "LOD Threshold", min = 4, max = 20, value = 7.5, round = TRUE),
   selectizeInput("which_trait", "Choose a trait",
                 choices = NULL,
                 multiple = FALSE,
                 options = list(
                   placeholder = 'Type to search traits...',
                   persist = TRUE,
                   createOnBlur = FALSE
                 )),
     helpText("Select a peak to see strain effects (Only for additive scans)."),
     selectizeInput("which_peak", "Choose peak", choices = NULL, multiple = FALSE)
   ),
   mainPanel(
     # 1. Datasets available table
     DTOutput("available_data"),
     br(),
     # 2. LOD profile plot with click output
     plotOutput("scan_plot", click = "plot_click"),
     verbatimTextOutput("scan_points"),
     br(),
     # 3. Peaks datatable
     DTOutput("peaks"),
     br(),
     # 4. Strain effects graph
     plotOutput("allele_effects"),
     br(),
     # Error message display
     verbatimTextOutput("error_message")
   )
 )
)


# Server=======================================================================
server <- function(input, output, session) {
  # Show available datasets (file_directory) for reference
 output$available_data <- renderDT(
   DT::datatable(file_directory, options = list(pageLength = 5, scrollX = TRUE))
 )
  # Reactive expression to load scan data based on dataset and trait
 scan_data_reactive <- reactive({
   req(input$selected_dataset, input$which_trait)
   tryCatch({
     trait_scan(file_directory, input$selected_dataset, input$which_trait)
   }, error = function(e) {
     # Check if the error is related to trait not found
     if (grepl("Trait.*not found", e$message)) {
       # Get available traits for suggestions
       dataset_files <- file_directory %>%
         filter(group == input$selected_dataset, file_type == "scans")
       fst_path <- csv2fst(dataset_files$File_path[1])
       scan_data_sample <- fst::read_fst(fst_path, columns = "Phenotype")
       available_traits <- unique(na.omit(scan_data_sample$Phenotype))
      
       # Find similar traits
       similar_traits <- available_traits[grep(input$which_trait, available_traits, ignore.case = TRUE)]
       if (length(similar_traits) > 0) {
         output$error_message <- renderText({
           paste("Did you mean one of these traits? ",
                 paste(head(similar_traits, 5), collapse = ", "))
         })
       } else {
         output$error_message <- renderText({
           paste("No matching traits found for:", input$which_trait)
         })
       }
     } else {
       output$error_message <- renderText({ paste("Error: ", e$message) })
     }
     NULL
   })
 })


   plot_obj <- reactive({
  req(scan_data_reactive())
  QTL_plot_visualizer(scan_data_reactive(), input$which_trait, input$LOD_thr, markers)
})


  # Reactive expression for QTL plot and underlying data
 plot_base <- reactive({
   req(plot_obj())
   plot_data <- plot_obj()[[2]]
  
   ggplot(plot_data, aes(x = BPcum, y = LOD)) +
     geom_line(aes(color = as.factor(chr)), alpha = 0.8, size = 0.5) +
     scale_color_manual(values = rep(c("black", "darkgrey"), 22)) +
     scale_x_continuous(
       label = plot_data %>%
         group_by(chr) %>%
         summarise(center = mean(BPcum)) %>%
         pull(chr),
       breaks = plot_data %>%
         group_by(chr) %>%
         summarise(center = mean(BPcum)) %>%
         pull(center),
       expand = expansion(mult = 0.15)
     ) +
     ylim(0, max(plot_data$LOD, na.rm = TRUE) * 1.25) +
     theme_bw() +
     theme(
       legend.position = "none",
       panel.border = element_blank(),
       panel.grid.major.x = element_blank(),
       panel.grid.minor.x = element_blank(),
       axis.line = element_line(colour = "black"),
       axis.text = element_text(size = 18),
       axis.title = element_text(size = 20),
       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, margin = margin(t = 10)),
       plot.margin = margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
     ) +
     xlab("Chromosome") +
     ylab("LOD")
 })


 # Render the plot with the threshold line


 # Render the plot with metadata title and threshold line
 output$scan_plot <- renderPlot({
   req(plot_base(), input$LOD_thr, input$which_trait)
  
   # Get peak information for the current trait
   peaks_info <- peak_finder(file_directory, input$selected_dataset) %>%
     filter(trait == input$which_trait) %>%
     arrange(desc(lod)) %>%
     slice(1)  # Get the highest peak for this trait
  
   # Create metadata text
   trait_text <- paste0("Trait: ", input$which_trait)
   peak_text <- if(nrow(peaks_info) > 0) {
     paste0("Peak Marker: ", peaks_info$marker,
           " (Chr", peaks_info$chr, ":", round(peaks_info$pos, 2), ")")
   } else {
     "No significant peaks"
   }
   lod_text <- if(nrow(peaks_info) > 0) {
     paste0("Max LOD: ", peaks_info$lod)
   } else {
     ""
   }
  
   # Add metadata to plot
   plot_base() +
     geom_hline(yintercept = input$LOD_thr, color = "black", linetype = "dashed") +
     ggtitle(label = trait_text,
             subtitle = paste0(peak_text, "\n", lod_text)) +
     theme(
       plot.title = element_text(color = "#2C3E50", size = 18, face = "bold"),
       plot.subtitle = element_text(color = "#E74C3C", size = 16),
       plot.title.position = "plot",
       plot.margin = margin(t = 40, r = 20, b = 20, l = 20)
     )
 })


  # Show near-point metadata on click
 output$scan_points <- renderPrint({
   req(plot_obj())
   near_points <- nearPoints(plot_obj()[[2]], input$plot_click,
                             xvar = "BPcum", yvar = "LOD", threshold = 10, maxpoints = 1, addDist = TRUE)
   if (nrow(near_points) > 0) {
     print(near_points)
   } else {
     "Click on a point to see metadata"
   }
 })




 # Update trait list based on selected dataset with improved search functionality
 observeEvent(input$selected_dataset, {
   req(input$selected_dataset)
   tryCatch({
     dataset_files <- file_directory %>%
       filter(group == input$selected_dataset, file_type == "scans")
    
     if (nrow(dataset_files) > 0) {
       fst_path <- csv2fst(dataset_files$File_path[1])
       scan_data_sample <- fst::read_fst(fst_path, columns = "Phenotype")
       available_traits <- unique(na.omit(scan_data_sample$Phenotype))
      
       updateSelectizeInput(session, "which_trait",
                          choices = available_traits,
                          server = TRUE,
                          options = list(
                            maxOptions = 10000,
                            create = FALSE,  # Don't allow creating new options
                            persist = TRUE,  # Keep the typed text when backspace
                            createOnBlur = FALSE,
                            selectOnTab = TRUE,
                            placeholder = 'Type to search traits...',
                            searchField = 'value'
                          ))
     }
   }, error = function(e) {
     updateSelectizeInput(session, "which_trait", choices = character(0))
   })
 })




 # Update Peaks Table when dataset is chosen
 observeEvent(input$selected_dataset, {
   req(input$selected_dataset)
  
   tryCatch({
     # Get peaks data
     peaks <- peak_finder(file_directory, input$selected_dataset)
    
     # Get available traits from scan data
     dataset_files <- file_directory %>%
       filter(group == input$selected_dataset, file_type == "scans")
    
     if (nrow(dataset_files) > 0) {
       fst_path <- csv2fst(dataset_files$File_path[1])
       scan_data_sample <- fst::read_fst(fst_path, columns = "Phenotype")
       available_scan_traits <- unique(na.omit(scan_data_sample$Phenotype))
      
       # Filter peaks to only show traits that exist in scan data
       peaks_to_show <- peaks %>%
         filter(trait %in% available_scan_traits) %>%
         select(marker, trait, chr, pos, lod) %>% # Keep only relevant columns
         arrange(trait, desc(lod)) # Sort by trait and LOD score
      
       output$peaks <- renderDT({
         datatable(peaks_to_show,
                  options = list(
                    pageLength = 5,
                    scrollX = TRUE,
                    dom = 'lftip', # Added 'f' for search
                    searchCols = list(
                      list(searchable = TRUE), # marker
                      list(searchable = TRUE), # trait
                      list(searchable = TRUE), # chr
                      list(searchable = TRUE), # pos
                      list(searchable = TRUE)  # lod
                    ),
                    columnDefs = list(
                      list(
                        targets = c(3, 4), # pos and lod columns
                        searchable = TRUE,
                        type = 'numeric'
                      )
                    )
                  ),
                  filter = 'top', # Adds filter row at top of each column
                  rownames = FALSE)
       })
     }
   }, error = function(e) {
     output$error_message <- renderText({
       paste("Error loading peaks data:", e$message)
     })
   })
 })


 # Update Peaks Dropdown based on selected trait




# Add this observer to update peak options when trait is selected
observeEvent(input$which_trait, {
 req(input$selected_dataset, input$which_trait)
 peaks <- peak_finder(file_directory, input$selected_dataset)
 if (!is.null(peaks) && "trait" %in% colnames(peaks)) {
   filtered_peaks <- peaks %>% filter(trait == input$which_trait)
   peak_options <- unique(filtered_peaks$marker)
   updateSelectizeInput(session, "which_peak", choices = peak_options, server = TRUE)
 } else {
   updateSelectizeInput(session, "which_peak", choices = character(0), server = TRUE)
 }
})


observeEvent(scan_data_reactive(), {
   if (!is.null(scan_data_reactive())) {
     output$error_message <- renderText({ "" })
   }
 })


# Plot Strain Effects when a peak is selected
observeEvent(input$which_peak, {  # Simplified trigger
 req(input$selected_dataset)
 req(input$which_trait)
 req(input$which_peak)
  peaks <- peak_finder(file_directory, input$selected_dataset)
  # Check if scan data exists and is additive
 if (!is.null(peaks) && "intcovar" %in% colnames(peaks) && peaks$intcovar[1] == "none") {
   # Filter peaks for the selected trait
   peaks <- peaks %>% filter(trait == input$which_trait, marker == input$which_peak)
  
   # Check if we have data after filtering
   if (nrow(peaks) > 0) {
     # Select and rename columns
     peaks <- peaks[c("marker","A","B","C","D","E","F","G","H")]
     colnames(peaks)[2:9] <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
    
     # Reshape data
     peaks <- reshape2::melt(peaks, id.vars = "marker")
    
     # Define colors
     newClrs <- c(
       "AJ" = "#000000",
       "B6" = "#96989A",
       "129" = "#E69F00",
       "NOD" = "#0072B2",
       "NZO" = "#619BFF",
       "CAST" = "#009E73",
       "PWK" = "#D55E00",
       "WSB" = "#CC79A7"
     )
    
     # Create plot
     plot_alleles <- ggplot(data = peaks, aes(x = marker, y = value, color = variable)) +
       geom_point(size = 10) +
       scale_color_manual(values = newClrs) +
       theme_bw() +
       theme(
         legend.text = element_text(size = 18),
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.line = element_line(colour = "black"),
         axis.text = element_text(size = 18),
         axis.title = element_text(size = 20)
       ) +
       labs(x = "Marker ID", y = "Founder allele effect", color = "Strain") +
       geom_hline(yintercept = 0, color = "black")
    
     output$allele_effects <- renderPlot({
       plot_alleles
     })
   } else {
     output$allele_effects <- renderPlot({ NULL })
   }
 } else {
   output$allele_effects <- renderPlot({ NULL })
 }
})
}


shinyApp(ui = ui, server = server)
