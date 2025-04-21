# Server logic for the Shiny application

server <- function(input, output, session) {
  
  # --- Reactive Values and States --- 
  scan_data <- reactiveVal(NULL)
  clicked_data <- reactiveVal(NULL)
  current_trait <- reactiveVal("") 
  official_gene_symbol <- reactiveVal("")
  use_alternating_colors <- reactiveVal(TRUE)
  strain_plot <- reactiveVal(NULL)
  
  # Debounced reactive for the trait text input
  trait_debounced <- reactive({
    input$which_trait
  }) %>% debounce::debounce(800) # Use explicit prefix
  
  # --- Observers and Event Handlers --- 
  
  # Automatic search based on debounced trait input
  observeEvent(trait_debounced(), {
    req(input$selected_dataset)
    trait_input <- trimws(trait_debounced())
    
    if (nchar(trait_input) > 0 && trait_input != current_trait()) {
      message("Debounced search triggered for trait: ", trait_input)
      current_trait(trait_input)
      
      # Clear previous state
      output$error_message <- renderText("") 
      shinyjs::hide("error_message_container")
      clicked_data(NULL)
      scan_data(NULL)
      
      # Perform trait scan (function sourced from R/data_access.R)
      tryCatch({
        # Pass cache environments explicitly
        data <- trait_scan(file_directory, input$selected_dataset, trait_input, trait_cache)
        scan_data(data)
        message("Trait scan successful.")
      }, error = function(e) {
        error_msg <- paste("Error searching for trait ", sQuote(trait_input), ": ", e$message)
        output$error_message <- renderText(error_msg)
        shinyjs::show("error_message_container")
        scan_data(NULL)
        message(error_msg)
      })
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # Manual search button observer (optional, retained)
  observeEvent(input$search_trait, {
      req(input$selected_dataset, input$which_trait)
      # This directly triggers the same logic as the debounced observer would
      updateTextInput(session, "which_trait", value = input$which_trait) # Ensure debounce triggers
  })
  
  # Actions when selected dataset changes
  observeEvent(input$selected_dataset, { 
    req(input$selected_dataset) 
    message("Dataset changed to: ", input$selected_dataset)
    
    # Clear caches
    rm(list = ls(envir = trait_cache), envir = trait_cache)
    rm(list = ls(envir = peaks_cache), envir = peaks_cache)
    message("Cleared trait and peaks caches.")
    
    # Reset UI elements and reactive values
    scan_data(NULL)
    clicked_data(NULL)
    current_trait("") 
    output$error_message <- renderText("")
    shinyjs::hide("error_message_container")
    updateSelectizeInput(session, "which_peak", choices = character(0), selected = character(0))
    updateTextInput(session, "which_trait", value = "")
    
    # Update dataset choices (in case they are dynamic, though loaded in global)
    # This might be redundant if choices are static from global.R
    # updateSelectizeInput(session, "selected_dataset", 
    #                      choices = unique(file_directory$group), 
    #                      selected = input$selected_dataset)
                         
  }, ignoreInit = TRUE)

  # Update available peaks dropdown based on scan data and LOD threshold
  observe({ 
      req(scan_data(), input$LOD_thr, input$selected_dataset)
      trait_val <- current_trait()
      req(trait_val)
      message("Updating peak list for trait: ", trait_val)
      
      # Get peaks data, filtering by LOD threshold
      peaks_data <- peak_finder(file_directory, input$selected_dataset, trait_val, peaks_cache)
      
      # Ensure necessary columns exist and filter
      # Use data.table directly for efficiency
      if (!is.null(peaks_data) && nrow(peaks_data) > 0 && "lod" %in% colnames(peaks_data)) {
          # Filter using data.table syntax
          filtered_peaks <- peaks_data[lod >= input$LOD_thr]
          data.table::setorder(filtered_peaks, -lod) # Sort by LOD descending
          
          if (nrow(filtered_peaks) > 0) {
              # Check for columns needed for labels
              if(all(c("marker", "chr", "pos", "lod") %in% colnames(filtered_peaks))){
                  # Create named vector using data.table syntax for efficiency
                  marker_choices <- filtered_peaks[, marker]
                  names(marker_choices) <- filtered_peaks[, paste0(marker,
                                               " (Chr", chr, 
                                               ": ", round(pos, 2),
                                               " Mb, LOD: ", round(lod, 2), ")")]
                  
                  updateSelectizeInput(session, "which_peak",
                                       choices = marker_choices,
                                       selected = marker_choices[1]) 
                  message("Peak list updated. Selected highest peak: ", names(marker_choices)[1])
              } else {
                   warning("Required columns missing in filtered peaks data for dropdown.")
                   updateSelectizeInput(session, "which_peak", choices = character(0), selected = character(0))
              }
          } else {
              message("No peaks found above LOD threshold: ", input$LOD_thr)
              updateSelectizeInput(session, "which_peak", choices = character(0), selected = character(0))
          }
      } else {
          message("No peaks data available or 'lod' column missing.")
          updateSelectizeInput(session, "which_peak", choices = character(0), selected = character(0))
      }
  })

  # Color toggle observer
  observeEvent(input$toggle_colors, {
      req(input$toggle_colors)
      use_alternating_colors(input$toggle_colors)
      message("Alternating colors set to: ", input$toggle_colors)
  })

  # Plot dimension preset observers
  observeEvent(input$preset_1to1, { updateNumericInput(session, "plot_width", value = 800); updateNumericInput(session, "plot_height", value = 800) })
  observeEvent(input$preset_3to2, { updateNumericInput(session, "plot_width", value = 900); updateNumericInput(session, "plot_height", value = 600) })
  observeEvent(input$preset_16to9, { updateNumericInput(session, "plot_width", value = 1600); updateNumericInput(session, "plot_height", value = 900) })

  # Plotly double-click observer (for zoom out)
  observeEvent(input$plotly_scan_doubleclick, {
      message("Plotly double-click detected, resetting chromosome view to 'All'.")
      updateSelectInput(session, "selected_chr", selected = "All")
  })
  
  # --- Data Processing Reactives --- 

  # Prepare data for plotting (using function from R/plotting.R)
  plot_obj_data <- reactive({ 
      req(scan_data()) 
      trait_val <- current_trait()
      req(trait_val)
      message("Preparing plot data for trait: ", trait_val)
      QTL_plot_visualizer(scan_data(), markers) # Expects data.table, returns data.table
  })

  # Get filtered peaks (using function from R/data_access.R)
  filtered_peaks_data <- reactive({ 
      req(input$selected_dataset, current_trait(), input$LOD_thr)
      trait_val <- current_trait()
      req(trait_val)
      
      peaks_dt <- peak_finder(file_directory, input$selected_dataset, trait_val, peaks_cache)
      
      if (!is.null(peaks_dt) && nrow(peaks_dt) > 0 && "lod" %in% colnames(peaks_dt)) {
           peaks_dt[lod >= input$LOD_thr] # Filter using data.table
      } else {
          # Return an empty data.table with expected columns
          data.table::data.table(marker=character(), chr=character(), pos=numeric(), lod=numeric())
      }
  })
  
  # Reactive for the highest peak (avoids recalculating in multiple places)
  highest_peak_info <- reactive({ 
      peaks <- filtered_peaks_data()
      if (nrow(peaks) > 0) {
          data.table::setorder(peaks, -lod) # Ensure sorted by LOD
          return(peaks[1, ]) # Return data.table row for the highest peak
      } else {
          return(NULL)
      }
  })

  # --- Reactive for Base ggplot Object --- 
  # This generates the core ggplot structure, shared by plotly output and downloads
  plot_base_gg <- reactive({ 
      req(plot_obj_data()) 
      plot_data <- plot_obj_data() # This is a data.table
      trait_val <- current_trait()
      req(trait_val)
      
      message("Generating base ggplot structure.")
      
      # Filter plot data for selected chromosome
      current_chr_view <- input$selected_chr
      if (current_chr_view != "All") {
          chr_num_selected <- data.table::fcase(
              current_chr_view == "X", 20L, current_chr_view == "Y", 21L, 
              current_chr_view == "M", 22L, default = suppressWarnings(as.integer(current_chr_view))
          )
          if (!is.na(chr_num_selected)) {
              plot_data <- plot_data[chr_numeric == chr_num_selected]
              message("Filtered plot data for chromosome: ", current_chr_view)
          } else {
              warning("Invalid chromosome selected: ", current_chr_view)
              plot_data <- plot_data[0,] # Return empty if invalid chr
          }
      }
      
      # Prepare axis info only if needed (All view)
      axis_df <- NULL
      x_aes_name <- "Position (Mb)"
      if (current_chr_view == "All") {
          x_aes_name <- "Chromosome"
          # Calculate centers using data.table for efficiency
          axis_df <- plot_data[, .(center = mean(BPcum)), keyby = chr_numeric]
          axis_df[, chr_label := data.table::fcase( chr_numeric == 20L, "X", chr_numeric == 21L, "Y", 
                                                  chr_numeric == 22L, "M", default = as.character(chr_numeric))]
      }
      
      # Create base ggplot object
      gg <- ggplot(plot_data, aes(y = LOD)) # Define y aesthetic upfront
      
      # Add line geom based on view
      if (current_chr_view == "All") {
          gg <- gg + aes(x = BPcum)
          if (use_alternating_colors()) {
              gg <- gg + geom_line(aes(color = factor(chr_numeric)), linewidth = 0.75)
          } else {
              gg <- gg + geom_line(color = "#3498db", linewidth = 0.75)
          }
      } else {
          gg <- gg + aes(x = position)
          gg <- gg + geom_line(color = "#3498db", linewidth = 0.75)
      }
      
      # Add common elements
      gg <- gg + 
          geom_hline(yintercept = input$LOD_thr, color = "#e74c3c", linetype = "dashed", linewidth = 0.8) +
          scale_y_continuous(name = "LOD Score", expand = expansion(mult = c(0.02, 0.1))) +
          theme_minimal(base_family = "Lato") +
          theme(
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "#ecf0f1", linewidth = 0.2),
              axis.line = element_line(color = "#2c3e50", linewidth = 0.5),
              axis.text = element_text(size = 11, color = "#2c3e50"),
              axis.title = element_text(size = 12, face = "bold", color = "#2c3e50", family="Montserrat"),
              axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
              legend.position = "none",
              plot.title = element_text(size = 16, face = "bold", color = "#2c3e50", family="Montserrat"),
              plot.subtitle = element_text(size = 12, color = "#7f8c8d"),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA),
              plot.margin = margin(t = 10, r = 20, b = 20, l = 20) 
          )
          
      # Apply x-axis scales based on view
      if (current_chr_view == "All" && !is.null(axis_df)) {
          gg <- gg + scale_x_continuous(name = x_aes_name,
                                         breaks = axis_df$center,
                                         labels = axis_df$chr_label,
                                         expand = expansion(mult = c(0.01, 0.01)))
          # Apply alternating colors if needed
          if (use_alternating_colors()) {
              gg <- gg + scale_color_manual(values = rep(c("#3498db", "#2c3e50"), length.out = nrow(axis_df)))
          }
      } else {
          gg <- gg + scale_x_continuous(name = x_aes_name, expand = expansion(mult = c(0.02, 0.02)))
      }
          
      # Return list containing plot object and the data used (filtered or full)
      return(list(gg = gg, data_in_view = plot_data))
  })

  # --- Plot Rendering --- 

  # Render the main LOD scan plot (interactive Plotly)
  output$scan_plot <- renderPlotly({ 
      req(plot_base_gg(), current_trait()) 
      plot_info <- plot_base_gg()
      gg <- plot_info$gg
      trait_val <- current_trait()
      req(trait_val)
      
      message("Rendering Plotly plot for trait: ", trait_val)
      
      # Get highest peak (already filtered by LOD and possibly chromosome)
      peak <- highest_peak_info() 
      
      # Add peak marker point to ggplot object *before* conversion
      if (!is.null(peak)) {
          # Find the point corresponding to the peak in the data currently plotted
          peak_point_plot <- plot_info$data_in_view[marker == peak$marker]
          
          if (nrow(peak_point_plot) > 0) {
              message("Adding peak marker to plot: ", peak$marker)
              x_aes <- if (input$selected_chr == "All") sym("BPcum") else sym("position")
              gg <- gg + geom_point(data = peak_point_plot, 
                                   aes(x = !!x_aes, y = LOD, text = paste("Peak:", marker)),
                                   color = "#e74c3c", size = 3, shape = 18)
          }
      }
      
      # Create title and subtitle
      official_symbol <- official_gene_symbol()
      plot_title_text <- if (nchar(official_symbol) > 0) official_symbol else trait_val
      plot_subtitle_text <- if (!is.null(peak)) {
          peak_chr_label <- data.table::fcase(peak$chr == "X", "X", peak$chr == "Y", "Y", peak$chr == "M", "M", default = as.character(peak$chr))
          paste0("Highest Peak: ", peak$marker, " (Chr", peak_chr_label, ": ", round(peak$pos, 2), " Mb) | LOD: ", round(peak$lod, 2))
      } else {
          if(input$selected_chr != "All") "No significant peaks in this chromosome" else "No significant peaks > threshold"
      }
      
      # Convert to plotly
      plt <- plotly::ggplotly(gg, 
                      source = "scan_plot", 
                      width = input$plot_width, 
                      height = input$plot_height, 
                      tooltip = c("x", "y")) %>% 
          plotly::layout(
              title = list(
                  text = paste0("<b style='font-size: 1.2em;'>", plot_title_text, "</b><br><span style='font-size: 0.9em;'>", plot_subtitle_text, "</span>"),
                  font = list(family = "Montserrat, sans-serif"),
                  x = 0.05, xanchor = 'left', y = 0.98, yanchor = 'top'
              ),
              margin = list(t = 80, l = 60, b = 60, r = 20), 
              hoverlabel = list(
                  bgcolor = "white",
                  font = list(family = "Lato, sans-serif", size = 12, color = "#2c3e50"),
                  bordercolor = "#bdc3c7"
              ),
              xaxis = list(titlefont = list(family="Montserrat, sans-serif")), 
              yaxis = list(titlefont = list(family="Montserrat, sans-serif"))
          ) %>% 
          plotly::config(
              displaylogo = FALSE,
              modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines")
          )
      
      # Add double-click event listener for zoom out
      if (input$selected_chr != "All") {
          plt <- plotly::event_register(plt, 'plotly_doubleclick')
      }
          
      return(plt)
  })
  
  # Render the allele effects plot (static ggplot)
  output$allele_effects <- renderPlot({ 
      req(input$which_peak) 
      peak_marker <- input$which_peak
      trait_val <- current_trait()
      req(trait_val)
      message("Rendering allele effects plot for peak: ", peak_marker)
      
      peaks <- peak_finder(file_directory, input$selected_dataset, trait_val, peaks_cache)
      selected_peak_data <- peaks[marker == peak_marker] # data.table subset
      
      allele_cols <- LETTERS[1:8]
      if (nrow(selected_peak_data) > 0 && all(allele_cols %in% colnames(selected_peak_data))) {
          # Use data.table for reshaping (potentially faster)
          effects_data <- selected_peak_data[1, c("marker", ..allele_cols)]
          effects_plot_data <- data.table::melt(effects_data, 
                                                id.vars = "marker", 
                                                measure.vars = allele_cols, 
                                                variable.name = "Strain", 
                                                value.name = "Effect")
          # Set strain names
          strain_map <- setNames(c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB"), LETTERS[1:8])
          effects_plot_data[, Strain := strain_map[as.character(Strain)]]
          
          strain_colors <- c("AJ"="#F0E442", "B6"="#000000", "129"="#E69F00", "NOD"="#56B4E9", 
                             "NZO"="#009E73", "CAST"="#0072B2", "PWK"="#D55E00", "WSB"="#CC79A7")
          
          official_symbol <- official_gene_symbol()
          plot_title_trait <- if (nchar(official_symbol) > 0) official_symbol else trait_val
          
          p_alleles <- ggplot(effects_plot_data, aes(x = Strain, y = Effect, color = Strain)) +
              geom_point(size = 4, alpha = 0.8) +
              scale_color_manual(values = strain_colors) +
              theme_minimal(base_family = "Lato") +
              theme(
                  legend.position = "none", 
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=10),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size=11, family="Montserrat"),
                  plot.title = element_text(size=14, face="bold", family="Montserrat"),
                  plot.subtitle = element_text(size=11)
              ) +
              labs(y = "Founder Allele Effect Estimate",
                   title = paste("Strain Effects at Peak:", peak_marker),
                   subtitle = paste("Trait:", plot_title_trait))
          
          strain_plot(p_alleles) 
          return(p_alleles)
          
      } else {
          message("Allele effect columns (A-H) not found or peak data missing.")
          p_empty <- ggplot() + 
                       annotate("text", x=0.5, y=0.5, label="Allele effects data not available\nfor this peak or dataset.", hjust=0.5, vjust=0.5, size=4) + 
                       theme_void()
          strain_plot(NULL)
          return(p_empty)
      }
  })
  
  # --- Table Rendering --- 
  
  # Render table with info about the clicked/selected point/peak
  output$clicked_point_info <- renderDT({
      event_data <- event_data("plotly_click", source = "scan_plot") # Use event_data() for clicks
      selected_marker <- input$which_peak
      trait_val <- current_trait()
      req(trait_val)
      
      point_info_dt <- NULL # Use data.table
      source_info <- "none"
      
      # Priority: Use selected peak
      if (!is.null(selected_marker) && nchar(selected_marker) > 0) {
          message("Getting info for selected peak: ", selected_marker)
          peaks <- peak_finder(file_directory, input$selected_dataset, trait_val, peaks_cache)
          peak_match <- peaks[marker == selected_marker] # data.table subset
          
          if (nrow(peak_match) > 0) {
              point_info_dt <- peak_match[1, ] 
              source_info <- "selected_peak"
          }
      }
      
      # Fallback: Use clicked point if no peak selected or selected peak not found
      if (is.null(point_info_dt) && !is.null(event_data) && length(event_data) > 0) {
          message("Getting info for clicked point near x=", round(event_data$x, 2), ", y=", round(event_data$y, 2))
          plot_data <- plot_obj_data() 
          req(plot_data)
          
          x_coord_col <- if(input$selected_chr == "All") "BPcum" else "position"
          if(!x_coord_col %in% colnames(plot_data)) return(NULL)
          
          # Find nearest point using data.table (potentially faster for large data)
          plot_data[, distance := sqrt((get(x_coord_col) - event_data$x)^2 + (LOD - event_data$y)^2)]
          nearest_point <- plot_data[which.min(distance)]
          plot_data[, distance := NULL] # Remove temporary column
          
          if (nrow(nearest_point) > 0) {
              nearest_marker <- nearest_point[, marker]
              message("Nearest marker from click: ", nearest_marker)
              
              peaks <- peak_finder(file_directory, input$selected_dataset, trait_val, peaks_cache)
              peak_match <- peaks[marker == nearest_marker]
              
              if (nrow(peak_match) > 0) {
                  point_info_dt <- peak_match[1, ]
                  source_info <- "clicked_peak"
              } else {
                  point_info_dt <- nearest_point
                  source_info <- "clicked_non_peak"
              }
          }
      }
      
      # Format data for table
      if (!is.null(point_info_dt)) {
          display_cols_map <- list(
              Marker = "marker", Chromosome = "chr", Position_Mb = "pos", LOD = "lod",
              Trait = "trait", CI_Low_Mb = "ci_lo", CI_High_Mb = "ci_hi", Cis = "cis",
              AJ = "A", B6 = "B", `129` = "C", NOD = "D", NZO = "E", 
              CAST = "F", PWK = "G", WSB = "H"
          )
          
          table_data <- list()
          # Handle cases where point_info_dt might be from scan data (LOD, position)
          if (!"lod" %in% colnames(point_info_dt) && "LOD" %in% colnames(point_info_dt)) {
              data.table::setnames(point_info_dt, "LOD", "lod")
          }
           if (!"pos" %in% colnames(point_info_dt) && "position" %in% colnames(point_info_dt)) {
              data.table::setnames(point_info_dt, "position", "pos")
          }
          
          for (display_name in names(display_cols_map)) {
              source_col <- display_cols_map[[display_name]]
              if (source_col %in% colnames(point_info_dt)) {
                  value <- point_info_dt[[source_col]]
                  # Formatting
                  if (source_col == "chr") {
                      value <- data.table::fcase(value == "20", "X", value == "21", "Y", value == "22", "M", default = as.character(value))
                  } else if (is.numeric(value) && source_col != "cis") { 
                      value <- round(value, 3)
                  }
                  table_data[[display_name]] <- value
              }
          }
          
          final_dt <- data.table::as.data.table(table_data)
          caption_text <- switch(source_info,
              selected_peak = "Information for Selected Peak",
              clicked_peak = "Information for Peak Near Clicked Point",
              clicked_non_peak = "Information for Point Near Click (Not a Peak)",
              "Point Information"
          )
          
          DT::datatable(final_dt,
              options = list(
                  dom = 't',        
                  ordering = FALSE,
                  pageLength = 1,   
                  scrollX = TRUE    
              ),
              rownames = FALSE,
              class = 'compact hover nowrap', 
              caption = htmltools::tags$caption(
                  style = 'caption-side: top; text-align: left; color: #2c3e50; font-weight: 600; font-size: 1em; margin-bottom: 5px;',
                  caption_text
              )
          )
      } else {
          return(NULL) 
      }
  })

  # --- Download Handlers --- 
  # (generate_filename and download handlers remain similar, 
  #  but ensure they use the reactive plot_base_gg() for plot generation)
  
  # Helper function remains the same
  generate_filename <- function(...) { ... }
  
  # Download handler for strain effects plot (PNG)
  output$download_effects_plot_png <- downloadHandler(
      filename = function() { generate_filename(...) },
      content = function(file) { ... }
  )
  
  # Download handler for strain effects plot (PDF)
  output$download_effects_plot_pdf <- downloadHandler(
      filename = function() { generate_filename(...) },
      content = function(file) { ... }
  )
  
  # Download handler for main QTL plot (PNG/PDF)
  # Consolidate plot generation logic for downloads
  generate_static_plot <- reactive({ 
      plot_info <- plot_base_gg()
      gg <- plot_info$gg
      trait_val <- current_trait()
      req(gg, trait_val)
      
      peak <- highest_peak_info()
      if (!is.null(peak)) {
          peak_point_plot <- plot_info$data_in_view[marker == peak$marker]
          if (nrow(peak_point_plot) > 0) {
              x_aes <- if (input$selected_chr == "All") sym("BPcum") else sym("position")
              gg <- gg + geom_point(data = peak_point_plot, aes(x = !!x_aes, y = LOD), color = "#e74c3c", size = 3, shape = 18)
          }
      }
      
      official_symbol <- official_gene_symbol()
      plot_title_text <- if (nchar(official_symbol) > 0) official_symbol else trait_val
      chr_info_text <- if(input$selected_chr != "All") paste0(" (Chr ", input$selected_chr, ")") else ""
      plot_subtitle_text <- if (!is.null(peak)) { ... } else { ... } # Same subtitle logic
      
      gg <- gg + ggtitle(label = paste0(plot_title_text, chr_info_text), subtitle = plot_subtitle_text) +
               theme(plot.title = element_text(size=16), plot.subtitle = element_text(size=12))
      return(gg)
  })
  
  output$download_qtl_plot_png <- downloadHandler(
      filename = function() { generate_filename(...) },
      content = function(file) {
          gg <- generate_static_plot()
          req(gg)
          ggsave(file, plot = gg, width = input$plot_width/96, height = input$plot_height/96, dpi = 300, device = "png")
      }
  )
  
  output$download_qtl_plot_pdf <- downloadHandler(
      filename = function() { generate_filename(...) },
      content = function(file) {
          gg <- generate_static_plot()
          req(gg)
          ggsave(file, plot = gg, width = input$plot_width/72, height = input$plot_height/72, device = cairo_pdf)
      }
  )
  
  # --- Initial Setup --- 
  observe({ 
      req(file_directory)
      dataset_choices <- unique(file_directory$group)
      if (is.null(input$selected_dataset) || input$selected_dataset == "") {
          updateSelectizeInput(session, "selected_dataset", choices = dataset_choices, selected = dataset_choices[1])
      }
  }, priority = 10) 

}
# --- End of server.R ---
