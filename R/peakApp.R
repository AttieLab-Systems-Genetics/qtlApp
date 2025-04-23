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
    title = "Test Scan",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"),
      mainParUI("main_par"),
      peakInput("peak")
    ),
    peakOutput("peak")
  )
  server <- function(input, output, session) {
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)
      peakServer("peak", main_par, import)
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
    shiny::observeEvent(shiny::req(peak_table(), main_par$LOD_thr), {
      ordered_markers <- highest_peaks(peak_table(), main_par$LOD_thr)$marker
      if(!is.null(ordered_markers)) {
        shiny::updateSelectizeInput(session, "which_peak",
          choices = ordered_markers, selected = ordered_markers[1],
          options = list(maxItems = 1, maxOptions = 5), server = TRUE)
      }
    })
    # Show allele effects.------------------------------------------------------
    output$allele_effects <- renderUI({
      shiny::req(peak_table(), input$which_peak)
      peak <- pivot_peaks(peak_table(), input$which_peak)
      plot_alleles <- ggplot_alleles(peak)
      shiny::renderPlot({
        print(plot_alleles)
      })
    })
    # Return
    peak_table
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
peakOutput <- function(id) {
  ns <- shiny::NS(id)
  list(
    bslib::card(
      bslib::card_header("Peaks"),
      DT::DTOutput(ns("peak_table"))
    ),
    bslib::card(
      bslib::card_header("Strain effects"),
      shiny::uiOutput(ns("allele_effects")) |>
        shinycssloaders::withSpinner(color="#0dc5c1"))
  )
}
