server <- function(input, output, session) {
  # Reactive value to store scan data
  scan_data <- reactiveVal(NULL)
  
  # Add reactive value to track click events
  clicked_data <- reactiveVal(NULL)
  
  # Add reactive value to track the current trait being searched
  current_trait <- reactiveVal("")
  
  # Create a debounced version of the trait input to prevent too many searches while typing
  # This will trigger the search 800ms after the user stops typing
  trait_debounced <- reactive({
    input$which_trait
  }) %>% debounce(800)
  
  # Observe the debounced trait input and trigger search automatically
  observeEvent(trait_debounced(), {
    req(input$selected_dataset)
    trait <- trait_debounced()
    
    # Only search if the trait is not empty and different from the current one
    if (nchar(trait) > 0 && trait != current_trait()) {
      # Update current trait
      current_trait(trait)
      
      # Clear error message and hide the container
      output$error_message <- renderText("")
      shinyjs::hide("error_message_container")
      
      # Clear clicked point data when a new trait is searched
      clicked_data(NULL)
      
      # Try to load the scan data for the specific trait
      tryCatch({
        message("Searching for trait: ", trait, " in dataset: ", input$selected_dataset)
        data <- trait_scan(file_directory, input$selected_dataset, trait)
        scan_data(data)
      }, error = function(e) {
        # Show error message
        output$error_message <- renderText({
          paste("Error:", e$message)
        })
        shinyjs::show("error_message_container")
        scan_data(NULL)
      })
    }
  })
  
  # Clear caches when dataset changes
  observeEvent(input$selected_dataset, {
    # Clear all caches to ensure we get fresh data for the new dataset
    rm(list = ls(envir = trait_cache), envir = trait_cache)
    rm(list = ls(envir = peaks_cache), envir = peaks_cache)
    
    # Remember the current trait
    current_trait_value <- input$which_trait
    
    # Clear scan data and clicked data
    scan_data(NULL)
    clicked_data(NULL)
    
    # Clear error message and hide the container
    output$error_message <- renderText("")
    shinyjs::hide("error_message_container")
    
    # Clear peak selection
    updateSelectizeInput(session, "which_peak", choices = NULL, selected = NULL)
    
    message("Cleared caches for dataset change to: ", input$selected_dataset)
    
    # If there was a trait already entered, automatically search for it in the new dataset
    if(nchar(current_trait_value) > 0) {
      message("Automatically searching for trait: ", current_trait_value, " in new dataset: ", input$selected_dataset)
      
      # Try to load the scan data for the same trait in the new dataset
      tryCatch({
        data <- trait_scan(file_directory, input$selected_dataset, current_trait_value)
        current_trait(current_trait_value)  # Update the current trait reactive
        scan_data(data)
        message("Successfully loaded trait data for new dataset")
      }, error = function(e) {
        # Show error message
        output$error_message <- renderText({
          paste("Error in new dataset:", e$message)
        })
        shinyjs::show("error_message_container")
        # Don't clear the trait input even if it fails - let the user see what they searched for
      })
    }
  })
  
  # Keep the search button functionality for backward compatibility
  observeEvent(input$search_trait, {
    req(input$selected_dataset, input$which_trait)
    
    # Clear error message
    output$error_message <- renderText("")
    
    # Clear clicked point data when a new trait is searched
    clicked_data(NULL)
    
    # Update current trait
    current_trait(input$which_trait)
    
    # Try to load the scan data for the specific trait
    tryCatch({
      data <- trait_scan(file_directory, input$selected_dataset, input$which_trait)
      scan_data(data)
    }, error = function(e) {
      output$error_message <- renderText({
        paste("Error:", e$message)
      })
      scan_data(NULL)
    })
  })
  
  # Also clear clicked data when trait input changes
  observeEvent(input$which_trait, {
    clicked_data(NULL)
  })
 
  # Main reactive expression for plot data
  plot_obj <- reactive({
      req(scan_data())
      QTL_plot_visualizer(scan_data(), current_trait(), input$LOD_thr, markers)
  })
 
  # Add reactive value to store the official gene symbol
  official_gene_symbol <- reactiveVal("")
  
  # Add reactive value to track color mode
  use_alternating_colors <- reactiveVal(TRUE)
  
  # Add observer for color toggle button
  observeEvent(input$toggle_colors, {
    use_alternating_colors(!use_alternating_colors())
    # Toggle the active class on the lever switch based on the new state
    runjs(sprintf("
      const toggle = document.getElementById('color_toggle');
      if (%s) {
        toggle.classList.add('active');
      } else {
        toggle.classList.remove('active');
      }
    ", tolower(use_alternating_colors())))
  })
  
  # Add observers for preset buttons
  observeEvent(input$preset_1to1, {
    updateNumericInput(session, "plot_width", value = 800)
    updateNumericInput(session, "plot_height", value = 800)
  })
  
  observeEvent(input$preset_3to2, {
    updateNumericInput(session, "plot_width", value = 900)
    updateNumericInput(session, "plot_height", value = 600)
  })
  
  observeEvent(input$preset_16to9, {
    updateNumericInput(session, "plot_width", value = 1600)
    updateNumericInput(session, "plot_height", value = 900)
  })
  
  # Update official gene symbol and automatically select highest peak when scan data changes
  observeEvent(scan_data(), {
    if (!is.null(scan_data()) && "Phenotype" %in% colnames(scan_data())) {
      # Set official gene symbol
      official_gene_symbol(unique(scan_data()$Phenotype)[1])
      
      # Automatically find and select the highest peak when a trait is loaded
      if (nrow(scan_data()) > 0) {
        message("Finding highest peak for trait: ", official_gene_symbol())
        
        # Get peaks data for the selected trait
        peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait)
        
        # Make sure we have data and the required columns
        if (nrow(peaks_data) > 0 && all(c("marker", "lod") %in% colnames(peaks_data))) {
          # Filter by LOD threshold and sort to find highest
          filtered_peaks <- peaks_data %>% 
            filter(lod >= input$LOD_thr) %>%
            arrange(desc(lod))
          
          if (nrow(filtered_peaks) > 0) {
            # Get the highest peak
            highest_peak <- filtered_peaks$marker[1]
            message("Found highest peak: ", highest_peak, " with LOD: ", round(filtered_peaks$lod[1], 2))
            
            # Update dropdown to select this peak
            updateSelectizeInput(session, "which_peak", selected = highest_peak)
            
            # Find this peak in the plot data
            plot_data <- plot_obj()[[2]]
            peak_point <- plot_data %>% filter(markers == highest_peak)
            
            if (nrow(peak_point) > 0) {
              # Create a fake click event at this peak's position
              fake_click <- list(
                x = if(input$selected_chr == "All") peak_point$BPcum[1] else peak_point$position[1],
                y = peak_point$LOD[1],
                curveNumber = 0,
                pointNumber = which(plot_data$markers == highest_peak)[1]
              )
              
              # Update clicked data to automatically show point info for highest peak
              clicked_data(fake_click)
              message("Auto-selected highest peak with LOD: ", round(filtered_peaks$lod[1], 2))
            } else {
              message("Could not find highest peak marker in plot data")
            }
          } else {
            message("No peaks found above threshold (", input$LOD_thr, ") for this trait")
          }
        } else {
          if (nrow(peaks_data) == 0) {
            message("No peaks found for this trait")
          } else {
            message("Required columns missing from peaks data. Have: ", 
                    paste(colnames(peaks_data), collapse=", "))
          }
        }
      }
    }
  })
  
  # Handle clicked points display
  observeEvent(event_data("plotly_click", source = "scan_plot"), {
      clicked_data(event_data("plotly_click", source = "scan_plot"))
  })
  
  # Also update clicked_data when a peak is selected in the dropdown
  observeEvent(input$which_peak, {
      req(input$which_peak, plot_obj())
      
      # Get the peak data
      peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait)
      selected_peak <- peaks_data %>% filter(marker == input$which_peak)
      
      if (nrow(selected_peak) > 0) {
          # Find this peak in the plot data
          plot_data <- plot_obj()[[2]]
          peak_point <- plot_data %>% filter(markers == selected_peak$marker[1])
          
          if (nrow(peak_point) > 0) {
              # Create a fake click event at this peak's position
              fake_click <- list(
                  x = if(input$selected_chr == "All") peak_point$BPcum[1] else peak_point$position[1],
                  y = peak_point$LOD[1],
                  curveNumber = 0,
                  pointNumber = which(plot_data$markers == selected_peak$marker[1])[1]
              )
              
              # Update the clicked data as if the plot was clicked
              clicked_data(fake_click)
              
              # Log that we're updating the clicked data
              message("Updated clicked data for peak: ", input$which_peak)
          }
      }
  })
  
  # Main reactive expression for plot data
  plot_base <- reactive({
      req(plot_obj())
      plot_data <- plot_obj()[[2]]
      
      # Filter data based on selected chromosome
      if (input$selected_chr != "All") {
          chr_num <- switch(input$selected_chr,
              "X" = 20,
              "Y" = 21,
              "M" = 22,
              as.numeric(input$selected_chr)
          )
          plot_data <- plot_data %>% filter(chr == chr_num)
      }
      
      # Calculate chromosome axis positions
      if (input$selected_chr == "All") {
          # For all chromosomes view
          axisdf <- plot_data %>%
              group_by(chr) %>%
              summarise(center = mean(BPcum))
          
          # Convert chromosome numbers to proper labels
          chr_labels <- as.character(axisdf$chr)
          chr_labels[chr_labels == "20"] <- "X"
          chr_labels[chr_labels == "21"] <- "Y"
          chr_labels[chr_labels == "22"] <- "M"
          
          # Create the base ggplot object for all chromosomes
          p <- ggplot(plot_data, aes(x = BPcum, y = LOD)) +
              geom_line(aes(color = if(use_alternating_colors()) as.factor(chr) else NULL), size = 0.75) +
              geom_hline(yintercept = input$LOD_thr, color = "#e74c3c",
                        linetype = "dashed", size = 0.8) +
              scale_color_manual(values = if(use_alternating_colors()) rep(c("#3498db", "#2c3e50"), 22) else "#3498db") +
              scale_x_continuous(
                  breaks = axisdf$center,
                  labels = chr_labels,
                  expand = expansion(mult = c(0.01, 0.01))
              )
      } else {
          # For single chromosome view
          # Create the base ggplot object for a single chromosome
          p <- ggplot(plot_data, aes(x = position, y = LOD)) +
              geom_line(color = "#3498db", size = 0.75) +
              geom_hline(yintercept = input$LOD_thr, color = "#e74c3c",
                        linetype = "dashed", size = 0.8) +
              scale_x_continuous(
                  expand = expansion(mult = c(0.02, 0.02))
              )
      }
      
      # Common theme and y-axis settings for both views
      p <- p +
          scale_y_continuous(
              expand = expansion(mult = c(0.02, 0.1))
          ) +
          theme_minimal() +
          theme(
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "#ecf0f1", size = 0.2),
              axis.line = element_line(color = "#2c3e50", size = 0.5),
              axis.text = element_text(size = 11, color = "#2c3e50"),
              axis.title = element_text(size = 12, face = "bold", color = "#2c3e50"),
              axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
              legend.position = "none",
              plot.title = element_text(size = 14, face = "bold", color = "#2c3e50"),
              plot.subtitle = element_text(size = 11, color = "#7f8c8d"),
              plot.background = element_rect(fill = "white", color = NA),
              panel.background = element_rect(fill = "white", color = NA),
              plot.margin = margin(t = 20, r = 20, b = 40, l = 40)
          ) +
          labs(
              x = if(input$selected_chr == "All") "Chromosome" else paste("Position on Chromosome", input$selected_chr, "(Mb)"),
              y = "LOD Score"
          )
      
      return(list(p = p, data = plot_data))
  })
 
  # Create the interactive plotly plot
  output$scan_plot <- renderPlotly({
      req(plot_base(), input$which_trait)
      
      # Get the base plot and data
      plot_result <- plot_base()
      p <- plot_result$p
      plot_data <- plot_result$data
      
      # Show the highest LOD peak
      peaks_info <- peak_finder(file_directory, input$selected_dataset, input$which_trait)
      if (nrow(peaks_info) > 0) {
          peaks_info <- peaks_info %>%
              arrange(desc(lod)) %>%
              slice(1)
          
          peak_point <- plot_data %>%
              filter(markers == peaks_info$marker)
          
          if (nrow(peak_point) > 0) {
              # Add the red diamond at the peak
              if (input$selected_chr == "All") {
                  p <- p + geom_point(data = peak_point,
                                     aes(x = BPcum, y = LOD),
                                     color = "red",
                                     size = 3,
                                     shape = 20)
              } else {
                  p <- p + geom_point(data = peak_point,
                                     aes(x = position, y = LOD),
                                     color = "red",
                                     size = 3,
                                     shape = 20)
              }
          }
      }
      
      # Create formatted trait text for plot title using the official gene symbol
      trait_text <- paste0("<b style='font-size: 24px;'>", official_gene_symbol(), "</b>")
      
      # Create subtitle with peak information
      subtitle <- if(nrow(peaks_info) > 0) {
          chr_label <- if(peaks_info$chr %in% c(20,21,22)) {
              c("X","Y","M")[peaks_info$chr-19]
          } else {
              peaks_info$chr
          }
          paste0(
              "<span style='font-size: 16px;'>",
              "<b>Peak Marker:</b> ", peaks_info$marker,
              " (Chr", chr_label, ":", round(peaks_info$pos, 2), " Mb) | ",
              "<b>LOD:</b> ", round(peaks_info$lod, 2),
              "</span>"
          )
      } else {
          "<span style='font-size: 16px; color: #7f8c8d;'>No significant peaks</span>"
      }
      
      # Convert ggplot to plotly with custom dimensions and removed features
      plt <- ggplotly(p, source = "scan_plot", width = input$plot_width, height = input$plot_height, 
                     tooltip = c("x", "y", "chr")) %>%
          layout(
              title = list(
                  text = paste0(trait_text, '<br>', subtitle),
                  font = list(family = "Arial"),
                  x = 0,
                  xanchor = "left",
                  y = 0.95,  
                  yanchor = "top",
                  pad = list(b = 20)  
              ),
              margin = list(t = 80),  
              hoverlabel = list(
                  bgcolor = "white",
                  font = list(family = "Arial", size = 12, color = "#2c3e50"),
                  bordercolor = "#95a5a6"
              ),
              hovermode = "closest",
              # Add double click event to reset to full view
              doubleclick = if(input$selected_chr != "All") TRUE else FALSE
          )
      
      # Remove unwanted modebar buttons
      plt <- plt %>% config(
          displaylogo = FALSE,
          modeBarButtonsToRemove = c(
              "select2d", "lasso2d", "autoScale2d", "hoverClosestCartesian",
              "hoverCompareCartesian", "toggleSpikelines"
          )
      )
      
      return(plt)
  })
 
  # Add observer for plotly double click event
  observeEvent(event_data("plotly_doubleclick", source = "scan_plot"), {
      if(input$selected_chr != "All") {
          updateSelectInput(session, "selected_chr", selected = "All")
      }
  })
 
  # Update which_peak dropdown when trait is found
  observe({
      req(input$selected_dataset, input$which_trait, input$LOD_thr)
     
      # Get peaks data for the selected trait from the peaks CSV files
      peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait) %>%
          filter(lod >= input$LOD_thr) %>%
          arrange(desc(lod))  # Sort by LOD score descending
     
      if (nrow(peaks_data) > 0) {
          # Create named vector for dropdown
          marker_choices <- peaks_data$marker
          names(marker_choices) <- paste0(peaks_data$marker,
                                       " (Chr", peaks_data$chr, ": ", round(peaks_data$pos, 2),
                                       " Mb, LOD: ", round(peaks_data$lod, 2), ")")
         
          # Select the maximum LOD peak by default
          updateSelectizeInput(session, "which_peak",
              choices = marker_choices,
              selected = marker_choices[1]  # First one is max LOD due to arrange(desc(lod))
          )
      } else {
          updateSelectizeInput(session, "which_peak",
              choices = character(0),
              selected = NULL
          )
      }
  })
  
  output$clicked_point_info <- renderDT({
      event_data <- clicked_data()
      
      # If no click data, try to get information from the selected peak instead
      if (is.null(event_data) && !is.null(input$which_peak)) {
          message("No click data, using selected peak: ", input$which_peak)
          peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait)
          
          # Debug the structure of peaks_data
          message("Peaks data columns: ", paste(colnames(peaks_data), collapse=", "))
          message("Number of rows in peaks_data: ", nrow(peaks_data))
          
          # Find peak by marker, considering possible column name differences
          peak_match <- NULL
          if ("marker" %in% colnames(peaks_data)) {
              peak_match <- peaks_data %>% filter(marker == input$which_peak)
          } else if ("markers" %in% colnames(peaks_data)) {
              peak_match <- peaks_data %>% filter(markers == input$which_peak)
          }
          
          if (!is.null(peak_match) && nrow(peak_match) > 0) {
              message("Found matching peak information row")
              
              # Debug the matched peak data
              message("Peak match columns: ", paste(colnames(peak_match), collapse=", "))
              
              # Create peak info table with safer column access
              peak_info <- data.frame(
                  Marker = ifelse("marker" %in% colnames(peak_match), 
                                 peak_match$marker[1], 
                                 ifelse("markers" %in% colnames(peak_match), 
                                       peak_match$markers[1], 
                                       NA)),
                  Chromosome = if("chr" %in% colnames(peak_match)) {
                      chr_val <- peak_match$chr[1]
                      if(is.numeric(chr_val) && chr_val %in% c(20,21,22)) {
                          c("X","Y","M")[chr_val-19]
                      } else {
                          chr_val
                      }
                  } else {
                      NA
                  },
                  Position = if("pos" %in% colnames(peak_match)) round(peak_match$pos[1], 3) else NA,
                  LOD = if("lod" %in% colnames(peak_match)) {
                      round(peak_match$lod[1], 3)
                  } else if("LOD" %in% colnames(peak_match)) {
                      round(peak_match$LOD[1], 3)
                  } else {
                      NA
                  }
              )
              
              # Add trait column if available
              if("trait" %in% colnames(peak_match)) {
                  peak_info$Trait <- peak_match$trait[1]
              } else if("lodcolumn" %in% colnames(peak_match)) {
                  peak_info$Trait <- peak_match$lodcolumn[1]
              }
              
              # Add Cis column if available
              if("cis" %in% colnames(peak_match)) {
                  peak_info$Cis <- peak_match$cis[1]
              }
              
              # Add confidence interval information if available
              if (all(c("ci_lo", "ci_hi") %in% colnames(peak_match))) {
                  peak_info$CI_Low <- round(peak_match$ci_lo[1], 3)
                  peak_info$CI_High <- round(peak_match$ci_hi[1], 3)
              }
              
              # Add strain effects if available
              if (all(c("A", "B", "C", "D", "E", "F", "G", "H") %in% colnames(peak_match))) {
                  # Add strain effects as columns, with safety checks
                  peak_info$AJ <- round(peak_match$A[1], 3)
                  peak_info$B6 <- round(peak_match$B[1], 3)
                  peak_info$`129` <- round(peak_match$C[1], 3)
                  peak_info$NOD <- round(peak_match$D[1], 3)
                  peak_info$NZO <- round(peak_match$E[1], 3)
                  peak_info$CAST <- round(peak_match$F[1], 3)
                  peak_info$PWK <- round(peak_match$G[1], 3)
                  peak_info$WSB <- round(peak_match$H[1], 3)
              }
              
              # Return the table
              return(datatable(
                  peak_info,
                  options = list(
                      dom = 't',
                      ordering = FALSE,
                      pageLength = 1
                  ),
                  rownames = FALSE,
                  class = 'compact hover',
                  caption = htmltools::tags$caption(
                      style = 'caption-side: top; text-align: left; color: #2c3e50; font-weight: bold; font-size: 14px;',
                      'Peak Information and Strain Effects'
                  )
              ))
          } else {
              message("No matching peak found in peaks data")
          }
      }
      
      if (is.null(event_data)) return(NULL)
     
      plot_data <- plot_obj()[[2]]
      if (is.null(plot_data)) return(NULL)
     
      # Find nearest point
      x_clicked <- event_data$x
      y_clicked <- event_data$y
     
      distances <- sqrt((plot_data$BPcum - x_clicked)^2 + (plot_data$LOD - y_clicked)^2)
      nearest_point <- plot_data[which.min(distances), ]
     
      if (nrow(nearest_point) == 0) return(NULL)
      
      # Check if this point is in the peaks file to get additional information
      peaks_data <- peak_finder(file_directory, input$selected_dataset, input$which_trait)
      
      # Debug
      message("Looking for marker: ", nearest_point$markers, " in peaks data")
      message("Peaks data columns: ", paste(colnames(peaks_data), collapse=", "))
      
      # Find if the clicked marker is in the peaks data (check both marker and markers)
      peak_match <- NULL
      if ("marker" %in% colnames(peaks_data) && nrow(peaks_data) > 0) {
          peak_match <- peaks_data %>% filter(marker == nearest_point$markers)
      } else if ("markers" %in% colnames(peaks_data) && nrow(peaks_data) > 0) {
          peak_match <- peaks_data %>% filter(markers == nearest_point$markers)
      }
      
      if (!is.null(peak_match) && nrow(peak_match) > 0) {
          message("Found matching peak in peaks data")
          
          # Create base peak info with safer column access
          peak_info <- data.frame(
              Marker = ifelse("marker" %in% colnames(peak_match), 
                             peak_match$marker[1], 
                             ifelse("markers" %in% colnames(peak_match), 
                                   peak_match$markers[1], 
                                   nearest_point$markers)),
              Chromosome = if("chr" %in% colnames(peak_match)) {
                  chr_val <- peak_match$chr[1]
                  if(is.numeric(chr_val) && chr_val %in% c(20,21,22)) {
                      c("X","Y","M")[chr_val-19]
                  } else {
                      chr_val
                  }
              } else {
                  if(nearest_point$chr %in% c(20,21,22)) 
                      c("X","Y","M")[nearest_point$chr-19] 
                  else 
                      nearest_point$chr
              },
              Position = if("pos" %in% colnames(peak_match)) {
                  round(peak_match$pos[1], 3)
              } else {
                  round(nearest_point$position, 3)
              },
              LOD = if("lod" %in% colnames(peak_match)) {
                  round(peak_match$lod[1], 3)
              } else if("LOD" %in% colnames(peak_match)) {
                  round(peak_match$LOD[1], 3)
              } else {
                  round(nearest_point$LOD, 3)
              }
          )
          
          # Add trait column if available
          if("trait" %in% colnames(peak_match)) {
              peak_info$Trait <- peak_match$trait[1]
          } else if("lodcolumn" %in% colnames(peak_match)) {
              peak_info$Trait <- peak_match$lodcolumn[1]
          }
          
          # Add Cis column if available
          if("cis" %in% colnames(peak_match)) {
              peak_info$Cis <- peak_match$cis[1]
          }
          
          # Add confidence interval information if available
          if (all(c("ci_lo", "ci_hi") %in% colnames(peak_match))) {
              peak_info$CI_Low <- round(peak_match$ci_lo[1], 3)
              peak_info$CI_High <- round(peak_match$ci_hi[1], 3)
          }
          
          # Add strain effects if available
          if (all(c("A", "B", "C", "D", "E", "F", "G", "H") %in% colnames(peak_match))) {
              # Add strain effects as columns
              peak_info$AJ <- round(peak_match$A[1], 3)
              peak_info$B6 <- round(peak_match$B[1], 3)
              peak_info$`129` <- round(peak_match$C[1], 3)
              peak_info$NOD <- round(peak_match$D[1], 3)
              peak_info$NZO <- round(peak_match$E[1], 3)
              peak_info$CAST <- round(peak_match$F[1], 3)
              peak_info$PWK <- round(peak_match$G[1], 3)
              peak_info$WSB <- round(peak_match$H[1], 3)
          }
          
          # Return the wide format table
          return(datatable(
              peak_info,
              options = list(
                  dom = 't',
                  ordering = FALSE,
                  pageLength = 1
              ),
              rownames = FALSE,
              class = 'compact hover',
              caption = htmltools::tags$caption(
                  style = 'caption-side: top; text-align: left; color: #2c3e50; font-weight: bold; font-size: 14px;',
                  'Peak Information and Strain Effects'
              )
          ))
      } else {
          message("No matching peak found in peaks data for clicked point")
      }
     
      # If no peak match found, return basic point info
      point_info <- data.frame(
          Marker = nearest_point$markers,
          Chromosome = if(nearest_point$chr %in% c(20,21,22))
                        c("X","Y","M")[nearest_point$chr-19] else nearest_point$chr,
          Position = round(nearest_point$position, 3),
          LOD = round(nearest_point$LOD, 3)
      )
      
      datatable(
          point_info,
          options = list(
              dom = 't',
              ordering = FALSE,
              pageLength = 1
          ),
          rownames = FALSE,
          class = 'compact hover'
      )
  })
 
  # Automatically select first dataset on startup
  observe({
      req(input$selected_dataset == "")
      first_dataset <- unique(file_directory$group)[1]
      updateSelectizeInput(session, "selected_dataset", selected = first_dataset)
  })
 
  # Handle allele effects plot
  observeEvent(c(input$which_peak, input$which_trait), {
      req(input$selected_dataset, input$which_trait, input$which_peak)
     
      # Get peaks data from CSV files
      peaks <- peak_finder(file_directory, input$selected_dataset, input$which_trait) %>%
          filter(marker == input$which_peak)
     
      if (nrow(peaks) > 0 && all(c("A","B","C","D","E","F","G","H") %in% colnames(peaks))) {
          peaks_plot <- peaks[1, c("marker","A","B","C","D","E","F","G","H")]
          colnames(peaks_plot)[2:9] <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
          peaks_plot <- reshape2::melt(peaks_plot, id.vars = "marker")
         
          newClrs <- c(
              "AJ" = "#000000", "B6" = "#96989A", "129" = "#E69F00",
              "NOD" = "#0072B2", "NZO" = "#619BFF", "CAST" = "#009E73",
              "PWK" = "#D55E00", "WSB" = "#CC79A7"
          )
          
          # Get the trait name - use official symbol if available, otherwise use input
          trait_name <- if(nchar(official_gene_symbol()) > 0) {
              official_gene_symbol()
          } else {
              input$which_trait
          }
         
          plot_alleles <- ggplot(data = peaks_plot, aes(x = marker, y = value, color = variable)) +
              geom_point(size = 4, alpha = 0.8) +
              scale_color_manual(values = newClrs) +
              theme_minimal() +
              theme(
                  legend.position = "right",
                  legend.title = element_text(size = 12, face = "bold"),
                  axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)  # Make x-axis text horizontal
              ) +
              labs(x = "Marker ID", y = "Founder Allele Effect", color = "Strain",
                   title = paste0("Strain Effects at ", input$which_peak),
                   subtitle = paste0("Trait: ", trait_name))
             
          output$allele_effects <- renderPlot({
              plot_alleles
          })
         
          # Store the plot in a reactive value for download handlers
          strain_plot(plot_alleles)
         
      } else {
          # If no strain effects data available
          output$allele_effects <- renderPlot({
              plot.new()
              text(0.5, 0.5, "No strain effects data available for this peak", cex = 1.2)
          })
          strain_plot(NULL)
      }
  })


  # Add reactive value to store the strain effects plot
  strain_plot <- reactiveVal(NULL)


  # Download handler for strain effects PNG
  output$download_effects_plot_png <- downloadHandler(
      filename = function() {
          # Use official gene symbol instead of user input for consistency
          actual_trait_name <- official_gene_symbol()
          paste0("strain_effects_", actual_trait_name, "_", input$which_peak, "_", format(Sys.time(), "%Y%m%d"), ".png")
      },
      content = function(file) {
          req(strain_plot())
          ggsave(file, strain_plot(), width = 10, height = 7, dpi = 300)
      }
  )

  # Download handler for strain effects PDF
  output$download_effects_plot_pdf <- downloadHandler(
      filename = function() {
          # Use official gene symbol instead of user input for consistency
          actual_trait_name <- official_gene_symbol()
          paste0("strain_effects_", actual_trait_name, "_", input$which_peak, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
          req(strain_plot())
          ggsave(file, strain_plot(), width = 10, height = 7, device = cairo_pdf)
      }
  )

  # Download handler for PNG
  output$download_qtl_plot_png <- downloadHandler(
      filename = function() {
          # Include chromosome info in filename if specific chromosome is selected
          chr_suffix <- if(input$selected_chr != "All") paste0("_chr", input$selected_chr) else ""
          paste0("lod_plot_", input$which_trait, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".png")
      },
      content = function(file) {
          # Get the base plot - this already has the chromosome filtering applied
          plot_result <- plot_base()
          p <- plot_result$p
          plot_data <- plot_result$data
          
          # Get peak information for the current trait
          peaks_info <- peak_finder(file_directory, input$selected_dataset, input$which_trait) %>%
              filter(lod >= input$LOD_thr) %>%
              arrange(desc(lod))
          
          # Filter peaks to only show those in the selected chromosome if applicable
          if(input$selected_chr != "All") {
              chr_num <- switch(input$selected_chr,
                  "X" = 20,
                  "Y" = 21,
                  "M" = 22,
                  as.numeric(input$selected_chr)
              )
              peaks_info <- peaks_info %>% filter(chr == chr_num)
          }
          
          # Get the highest peak (if any)
          if(nrow(peaks_info) > 0) {
              peaks_info <- peaks_info %>% slice(1)
              
              # Find the peak point in the filtered data
              peak_point <- plot_data %>%
                  filter(markers == peaks_info$marker)
              
              if (nrow(peak_point) > 0) {
                  # Add the red diamond at the peak
                  if (input$selected_chr == "All") {
                      p <- p + geom_point(data = peak_point,
                                         aes(x = BPcum, y = LOD),
                                         color = "red",
                                         size = 3,
                                         shape = 18)
                  } else {
                      p <- p + geom_point(data = peak_point,
                                         aes(x = position, y = LOD),
                                         color = "red",
                                         size = 3,
                                         shape = 18)
                  }
              }
          }
          
          # Create plot title and subtitle
          trait_title <- input$which_trait
          chr_info <- if(input$selected_chr != "All") paste0(" (Chromosome ", input$selected_chr, ")") else ""
          
          subtitle <- if(nrow(peaks_info) > 0) {
              chr_label <- if(peaks_info$chr %in% c(20,21,22)) {
                  c("X","Y","M")[peaks_info$chr-19]
              } else {
                  peaks_info$chr
              }
              paste0(
                  "Peak Marker: ", peaks_info$marker,
                  " (Chr", chr_label, ":", round(peaks_info$pos, 2), " Mb) | ",
                  "LOD: ", round(peaks_info$lod, 2)
              )
          } else {
              if(input$selected_chr != "All") {
                  "No significant peaks in this chromosome"
              } else {
                  "No significant peaks"
              }
          }
          
          # Add title and subtitle to the plot
          p <- p + ggtitle(
              label = paste0(trait_title, chr_info),
              subtitle = subtitle
          )
         
          # Save the plot with high resolution
          ggsave(file, p, width = input$plot_width/72, height = input$plot_height/72, dpi = 300)
      }
  )

  # Download handler for PDF
  output$download_qtl_plot_pdf <- downloadHandler(
      filename = function() {
          # Include chromosome info in filename if specific chromosome is selected
          chr_suffix <- if(input$selected_chr != "All") paste0("_chr", input$selected_chr) else ""
          paste0("lod_plot_", input$which_trait, chr_suffix, "_", format(Sys.time(), "%Y%m%d"), ".pdf")
      },
      content = function(file) {
          # Get the base plot - this already has the chromosome filtering applied
          plot_result <- plot_base()
          p <- plot_result$p
          plot_data <- plot_result$data
          
          # Get peak information for the current trait
          peaks_info <- peak_finder(file_directory, input$selected_dataset, input$which_trait) %>%
              filter(lod >= input$LOD_thr) %>%
              arrange(desc(lod))
          
          # Filter peaks to only show those in the selected chromosome if applicable
          if(input$selected_chr != "All") {
              chr_num <- switch(input$selected_chr,
                  "X" = 20,
                  "Y" = 21,
                  "M" = 22,
                  as.numeric(input$selected_chr)
              )
              peaks_info <- peaks_info %>% filter(chr == chr_num)
          }
          
          # Get the highest peak (if any)
          if(nrow(peaks_info) > 0) {
              peaks_info <- peaks_info %>% slice(1)
              
              # Find the peak point in the filtered data
              peak_point <- plot_data %>%
                  filter(markers == peaks_info$marker)
              
              if (nrow(peak_point) > 0) {
                  # Add the red diamond at the peak
                  if (input$selected_chr == "All") {
                      p <- p + geom_point(data = peak_point,
                                         aes(x = BPcum, y = LOD),
                                         color = "red",
                                         size = 3,
                                         shape = 18)
                  } else {
                      p <- p + geom_point(data = peak_point,
                                         aes(x = position, y = LOD),
                                         color = "red",
                                         size = 3,
                                         shape = 18)
                  }
              }
          }
          
          # Create plot title and subtitle
          trait_title <- input$which_trait
          chr_info <- if(input$selected_chr != "All") paste0(" (Chromosome ", input$selected_chr, ")") else ""
          
          subtitle <- if(nrow(peaks_info) > 0) {
              chr_label <- if(peaks_info$chr %in% c(20,21,22)) {
                  c("X","Y","M")[peaks_info$chr-19]
              } else {
                  peaks_info$chr
              }
              paste0(
                  "Peak Marker: ", peaks_info$marker,
                  " (Chr", chr_label, ":", round(peaks_info$pos, 2), " Mb) | ",
                  "LOD: ", round(peaks_info$lod, 2)
              )
          } else {
              if(input$selected_chr != "All") {
                  "No significant peaks in this chromosome"
              } else {
                  "No significant peaks"
              }
          }
          
          # Add title and subtitle to the plot
          p <- p + ggtitle(
              label = paste0(trait_title, chr_info),
              subtitle = subtitle
          )
         
          # Save the plot as PDF
          ggsave(file, p, width = input$plot_width/72, height = input$plot_height/72, device = cairo_pdf)
      }
  )
}