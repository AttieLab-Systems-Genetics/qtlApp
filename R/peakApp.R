#' Peak App Module
#'
#' @param id shiny identifier
#' @param import reactive list with file_directory
#'
#' @importFrom DT datatable DTOutput renderDT
#' @importFrom shiny moduleServer NS observeEvent plotOutput reactive renderPlot renderText
#'             req selectizeInput setProgress shinyApp textOutput updateSelectizeInput withProgress
#' @importFrom bslib card card_header page_sidebar sidebar
#' @importFrom shinycssloaders withSpinner
#' @importFrom stringr str_split
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes element_blank element_line element_text geom_hline geom_point ggplot
#'             labs scale_color_manual theme theme_bw
#' 
#' @export
peakApp <- function() {
  source(system.file("shinyApp/qtlSetup.R", package = "qtlApp"))
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

    chosen_trait <- shiny::reactive({
      shiny::req(main_par$which_trait)
      stringr::str_split(
        stringr::str_split(main_par$which_trait, pattern=" [(]")[[1]][2],
        pattern="[)]")[[1]][1]
    })
    peaks <- shiny::reactive({
      shiny::req(main_par$group, chosen_trait())
      shiny::withProgress(
        message = paste("peaks of", chosen_trait(), "in progress"),
        value = 0, {
          shiny::setProgress(1)
          subset(peak_finder(import()$file_directory, main_par$group), trait == chosen_trait())
        })
    })

    # get the peaks info ------------------------------------------------------------------
    output$peaks <- DT::renderDT({DT::datatable(
      peaks(),
      options = list(paging = TRUE,    ## paginate the output
                     pageLength = 5,  ## number of rows to output for each page
                     scrollX = TRUE,   ## enable scrolling on X axis
                     scrollY = TRUE,   ## enable scrolling on Y axis
                     autoWidth = TRUE, ## use smart column width handling
                     server = TRUE,   ## use client-side processing
                     dom = 'Bfrtip',
                     buttons = c('csv', 'excel'),
                     columnDefs = list(list(targets = '_all', className = 'dt-center'),
                                       list(targets = c(0, 8, 9), visible = FALSE))
      ),
      extensions = 'Buttons',
      selection = 'single', ## enable selection of a single row
      filter = 'bottom',              ## include column filters at the bottom
      rownames = TRUE)                ##  show row numbers/names
    })
        
    # Update peak selection-------------------------------------------
    shiny::observeEvent(peaks(), {
      shiny::updateSelectizeInput(session, "which_peak",
        choices = peaks()$marker, options = list(maxItems = 1, maxOptions = 5), server = TRUE)
    })
    # Show allele effects.-----------------------------------------------------------
    output$allele_effects <- renderUI({
      shiny::req(peaks(), input$which_peak)
      # Check if scan data exists and is additive
      if (peaks()$intcovar[1] == "none") {
        # set peaks
        peak <- subset(peaks(), marker == input$which_peak)  # Changed from marker.id to marker
        
        # Check if we have data after filtering
        if (nrow(peak) > 0) {
          # Select and rename columns
          peak <- peak[c("marker","A","B","C","D","E","F","G","H")]
          colnames(peak)[2:9] <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
          # Reshape data
          peak <- reshape2::melt(peak, id.vars = "marker")
            
          # Define colors
          newClrs <- c(
                "AJ" = "#000000",
                "B6" = "#96989A",
                "129" = "#E69F00",
                "NOD" = "#0072B2",
                "NZO" = "#619BFF",
                "CAST" = "#009E73",
                "PWK" = "#D55E00",
                "WSB" = "#CC79A7")
          # Create plot
          plot_alleles <- ggplot2::ggplot(data = peak,
            ggplot2::aes(x = marker, y = value, color = variable)) +
            ggplot2::geom_point(size = 10) +
            ggplot2::scale_color_manual(values = newClrs) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              legend.text = ggplot2::element_text(size = 18),
              panel.border = ggplot2::element_blank(),
              panel.grid.major.x = ggplot2::element_blank(),
              panel.grid.minor.x = ggplot2::element_blank(),
              axis.line = ggplot2::element_line(colour = "black"),
              axis.text = ggplot2::element_text(size = 18),
              axis.title = ggplot2::element_text(size = 20)) +
            ggplot2::labs(x = "Marker ID", y = "Founder allele effect", color = "Strain") +
            ggplot2::geom_hline(yintercept = 0, color = "black")
            
          shiny::renderPlot({
            print(plot_alleles)
          })
        } else {
          shiny::renderText({"No data available for selected peak"})
        }
      } else {
        shiny::renderText({"Strain effects only available for additive scans"})
      }
    })
  })
}
#' @rdname peakApp
#' @export
peakInput <- function(id) {
  ns <- shiny::NS(id)
  list(
    shiny::helpText("Choose a peak to see the strain effects. This only applies to the additive scans."),
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
      DT::DTOutput(ns("peaks"))
    ),
    bslib::card(
      bslib::card_header("Strain effects"),
      shiny::uiOutput(ns("allele_effects")) |>
        shinycssloaders::withSpinner(color="#0dc5c1"))
  )
}
