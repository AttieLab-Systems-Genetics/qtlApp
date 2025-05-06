#' Peak App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory
#'
#' @importFrom DT datatable DTOutput renderDT
#' @importFrom shiny moduleServer NS observeEvent plotOutput reactive renderPlot
#'             renderText req selectizeInput setProgress shinyApp textOutput
#'             updateSelectizeInput withProgress
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
#' @importFrom stringr str_split
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes element_blank element_line element_text geom_hline
#'             geom_point ggplot labs scale_color_manual theme theme_bw
#' 
#' @export
peakApp <- function() {
  ui <- bslib::page_sidebar(
    title = "Test Peak",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"),    # "selected_dataset", "LOD_thr"
      mainParUI("main_par"),       # "which_trait"
      peakInput("peak"),           # "which_peak"
      downloadInput("download"),   # downloadButton, filename
      downloadOutput("download")), # plot_table, inputs for Plots or Tables
    bslib::card(
      bslib::card_header("Peaks"),
      peakOutput("peak")),
    bslib::card(
      bslib::card_header("Strain effects"),
      peakUI("peak"))
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      peak_list <- peakServer("peak", main_par, import)
      downloadServer("download", peak_list)
  }
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname peakApp
#' @export
peakServer <- function(id, main_par, import) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ** new version uses marker not trait **

    chosen_trait <- shiny::reactive({
      shiny::req(import(), main_par$which_trait, main_par$selected_dataset)
      get_selected_trait(import(),
                         main_par$which_trait, main_par$selected_dataset)
    })
    peak_table <- shiny::reactive({
      shiny::req(main_par$selected_dataset, chosen_trait())
      shiny::withProgress(
        message = paste("peaks of", chosen_trait(), "in progress"),
        value = 0, {
          shiny::setProgress(1)
          subset(
            suppressMessages(peak_finder(import()$file_directory,
                                         main_par$selected_dataset)),
            trait == chosen_trait())
        })
    })

    # --- Debug: Log the dimensions of the peak_table before returning --- 
    observeEvent(peak_table(), {
        pt <- peak_table()
        message("--- Debug peakServer output ---")
        message("peak_table class: ", class(pt)[1])
        message("peak_table dimensions: ", paste(dim(pt), collapse="x"))
        message("-----------------------------")
    }, ignoreNULL = FALSE) # Log even if NULL
    # ------------------------------------------------------------------

    # Table ov peaks info ------------------------------------------------------
    output$peak_table <- DT::renderDT({DT::datatable(
      peak_table(),
      options = list(paging = TRUE,    ## paginate the output
                     pageLength = 5,   ## number of rows to output for each page
                     scrollX = TRUE,   ## enable scrolling on X axis
                     scrollY = TRUE,   ## enable scrolling on Y axis
                     autoWidth = TRUE, ## use smart column width handling
                     server = TRUE,    ## use client-side processing
                     dom = 'Bfrtip',
                     buttons = c('csv', 'excel'),
                     columnDefs = list(
                       list(targets = '_all', className = 'dt-center'),
                       list(targets = c(0, 8, 9), visible = FALSE))
      ),
      extensions = 'Buttons',
      selection = 'single',            ## enable selection of a single row
      filter = 'bottom',               ## include column filters at the bottom
      rownames = TRUE)                 ##  show row numbers/names
    })
    # Update peak selection-------------------------------------------
    # TEMPORARY DEBUG: Remove peak_table() from req() to see if it stops the '==' error
    shiny::observeEvent(shiny::req(main_par$LOD_thr()), { 
      # This will likely error later if peak_table() isn't ready, but tests the req() call
      ordered_markers <- highest_peaks(peak_table(), main_par$LOD_thr())$marker
      if(!is.null(ordered_markers)) {
        shiny::updateSelectizeInput(session, "which_peak",
          choices = ordered_markers, selected = ordered_markers[1],
          options = list(maxItems = 1, maxOptions = 5), server = TRUE)
      }
    })
    # Show allele effects.------------------------------------------------------
    output$allele_effects <- renderUI({
      shiny::renderPlot({
        # ** Error: `ui_element` must be a Shiny tag. **
        print(shiny::req(allele_plot()))# |>
          #shinycssloaders::withSpinner(color="#0dc5c1")
      })
    })
    allele_plot <- shiny::reactive({
      shiny::req(peak_table(), input$which_peak)
      peak <- pivot_peaks(peak_table(), input$which_peak)
      ggplot_alleles(peak)
    })
    file_name <- shiny::reactive({
      # Evaluate which_trait reactive
      trait_val <- shiny::req(main_par$which_trait())
      peak_val <- shiny::req(input$which_peak)
      instanceID <- paste(trait_val, peak_val, sep = "_") 
      
      # Debug before comparison
      sel_chr_val <- main_par$selected_chr()
      
      # Defensive check before comparison
      if (is.atomic(sel_chr_val) && is.character(sel_chr_val) && length(sel_chr_val) == 1) {
        # Evaluate selected_chr reactive and compare                        
        if(shiny::req(sel_chr_val) != "All") { 
          # Evaluate selected_chr reactive for paste
          instanceID <- paste0(instanceID, "_chr", sel_chr_val)
        }
      } else {
        warning("peakServer file_name: sel_chr_val is not a single character string!")
      }
      
      paste("peak", instanceID, sep = "_")
    })
    # Return `peak_list` = reactiveValues containing elements `filename`, `tables` and `plots`.
    # --- MODIFICATION: Return ONLY the peak_table reactive directly --- 
    # shiny::reactiveValues(
    #   filename = file_name,
    #   tables = shiny::reactiveValues(
    #     peak = peak_table),
    #   plots  = shiny::reactiveValues(
    #     alleles = allele_plot)
    # )
    return(peak_table) # Return the reactive expression
    # ------------------------------------------------------------------

  })
}
#' @rdname peakApp
#' @export
peakInput <- function(id) {
  ns <- shiny::NS(id)
  list(
    shiny::helpText("Choose a peak to see the strain effects.",
      "This only applies to the additive scans."),
    shiny::selectizeInput(ns("which_peak"),
      label = "Choose peak",
      choices = NULL,
      multiple = TRUE,
      options = list(
        placeholder = 'Search...'
      )))
}
#' @rdname peakApp
#' @export
peakUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("allele_effects"))
}
#' @rdname peakApp
#' @export
peakOutput <- function(id) {
  ns <- shiny::NS(id)
  DT::DTOutput(ns("peak_table"))
}
