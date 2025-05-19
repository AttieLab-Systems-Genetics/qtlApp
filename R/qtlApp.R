#' QTL App
#' 
#' @param id shiny identifier
#' @param import reactive list with file_directory, annotation_list and markers
#' 
#' @importFrom shiny helpText moduleServer NS shinyApp tags
#' @importFrom bslib page_sidebar sidebar bs_theme
#' @importFrom shinyjs useShinyjs
#' 
#' @export
qtlApp <- function() {
  if (!exists("create_fluid_page", mode = "function")) {
    source("R/ui_styles.R")
  }
  
  if (exists("create_fluid_page", mode = "function") && 
      exists("create_title_panel", mode = "function")) {
    
    ui <- create_fluid_page(
      shinyjs::useShinyjs(),
      
      tags$head(tags$style(custom_css)),
      
      create_title_panel(
        "Pre-scanned QTL Visualizer for Diet DO Study",
        "Interactive visualization tool for QTL analysis"
      ),
      
      create_fluid_row(
        create_column(3,
          create_well_panel(
            qtlInput("qtl")
          )
        ),
        
        create_column(9,
          qtlOutput("qtl")
        )
      )
    )
    
  } else {
    ui <- bslib::page_sidebar(
      title = "Pre-scanned QTL Visualizer for Diet DO Study (qtlApp)",
      shinyjs::useShinyjs(),
      sidebar = bslib::sidebar("side_panel",
      qtlInput("qtl")),
    bslib::navset_tab(
      id = "main_qtl_tabs",
      bslib::nav_panel("LOD Plot",
        bslib::card(
          bslib::card_header("LOD Plot Controls"),
          bslib::card_body(
            shiny::div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px; flex-wrap: wrap;",
              shiny::h4("LOD Score Plot", style = "margin: 0 15px 0 0; color: #2c3e50; font-weight: 600;"),
              shiny::div(style = "display: flex; align-items: center; gap: 10px; flex-grow: 1; justify-content: flex-end;",
                shiny::div(style = "display: flex; gap: 5px; margin-right: 15px;",
                  if (exists("create_button", mode = "function")) {
                    shiny::tagList(
                      create_button(shiny::NS("qtl", shiny::NS("scan_list", "preset_1to1")), "1:1", class = "btn-sm btn-light"),
                      create_button(shiny::NS("qtl", shiny::NS("scan_list", "preset_3to2")), "3:2", class = "btn-sm btn-light"),
                      create_button(shiny::NS("qtl", shiny::NS("scan_list", "preset_16to9")), "16:9", class = "btn-sm btn-light")
                    )
                  } else {
                    shiny::tagList(
                      shiny::actionButton(shiny::NS("qtl", shiny::NS("scan_list", "preset_1to1")), "1:1", class = "btn-sm"),
                      shiny::actionButton(shiny::NS("qtl", shiny::NS("scan_list", "preset_3to2")), "3:2", class = "btn-sm"),
                      shiny::actionButton(shiny::NS("qtl", shiny::NS("scan_list", "preset_16to9")), "16:9", class = "btn-sm")
                    )
                  }
                ),
                if (exists("create_download_button", mode = "function")) {
                  shiny::tagList(
                    create_download_button(shiny::NS("qtl", shiny::NS("scan_list", "download_qtl_plot_png")), "PNG", class = "btn-sm"),
                    create_download_button(shiny::NS("qtl", shiny::NS("scan_list", "download_qtl_plot_pdf")), "PDF", class = "btn-sm")
                  )
                } else {
                  shiny::tagList(
                    shiny::downloadButton(shiny::NS("qtl", shiny::NS("scan_list", "download_qtl_plot_png")), "PNG", class = "btn-sm"),
                    shiny::downloadButton(shiny::NS("qtl", shiny::NS("scan_list", "download_qtl_plot_pdf")), "PDF", class = "btn-sm")
                  )
                }
              )
            ),
            shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 5px; flex-wrap: wrap;",
              shiny::div(style = "display: flex; align-items: center; gap: 10px;",
                if (exists("create_numeric_input", mode = "function")) {
                  create_numeric_input(shiny::NS("qtl", shiny::NS("scan_list", "plot_width")), "Width:", value = 1200, min = 400, max = 2000, step = 50, width = "100px")
                } else {
                  shiny::numericInput(shiny::NS("qtl", shiny::NS("scan_list", "plot_width")), "Width:", value = 1200, min = 400, max = 2000, step = 50, width = "100px")
                },
                if (exists("create_numeric_input", mode = "function")) {
                  create_numeric_input(shiny::NS("qtl", shiny::NS("scan_list", "plot_height")), "Height:", value = 600, min = 300, max = 1200, step = 50, width = "100px")
                } else {
                  shiny::numericInput(shiny::NS("qtl", shiny::NS("scan_list", "plot_height")), "Height:", value = 600, min = 300, max = 1200, step = 50, width = "100px")
                }
              ),
              if (exists("create_lever_switch", mode = "function")) {
                create_lever_switch(shiny::NS("qtl", shiny::NS("scan_list", "color_toggle")))
              } else {
                shiny::checkboxInput(shiny::NS("qtl", shiny::NS("scan_list", "color_toggle")), "Alt Colors", value = TRUE)
              }
            )
          )
        ),
        scanlyOutput(shiny::NS("qtl", "scanly")),
        bslib::card(
          bslib::card_header("Clicked Peak on Scan"),
          scanlyUI(shiny::NS("qtl", "scanly")),
          max_height = "150px" 
        )
      ),
      bslib::nav_panel("Peaks Table",
        bslib::card(
          peakOutput(shiny::NS("qtl", "peak_list"))
        )
      ),
      bslib::nav_panel("Cis/Trans Plot",
        cisTransPlotUI(shiny::NS("qtl", "cistrans"))
      )
    )
  )
  }
  
  server <- function(input, output, session) {
    qtlServer("qtl")
  }
  
  shiny::shinyApp(ui = ui, server = server)
}
#' @rdname qtlApp
#' @export
qtlServer <- function(id) {
    shiny::moduleServer(id, function(input, output, session) {
      ns <- session$ns
      import <- importServer("import")
      main_par <- mainParServer("main_par", import)

      # Create a shared cache environment for peak data
      peaks_shared_cache <- new.env(parent = emptyenv())
      # Create a shared cache environment for trait scan data
      trait_shared_cache <- new.env(parent = emptyenv())

      # Clear caches when dataset changes
      shiny::observeEvent(main_par$selected_dataset(), {
          rm(list = ls(envir = trait_shared_cache), envir = trait_shared_cache)
          rm(list = ls(envir = peaks_shared_cache), envir = peaks_shared_cache)
      })

      # Pass the trait_cache to scanServer
      scan_list <- scanServer("scan_list", main_par, import, trait_cache = trait_shared_cache)
      # Pass the shared peaks_cache to peakServer
      peak_list <- peakServer("peak_list", main_par, import, peaks_cache = peaks_shared_cache)
      # Add cisTransPlotServer and pass the shared cache and import reactives
      cisTransPlotServer("cistrans", import_reactives = import, peaks_cache = peaks_shared_cache)

      scanlyServer("scanly", main_par, scan_list, peak_list)
  })
}
#' @rdname qtlApp
#' @export
qtlInput <- function(id) {
  ns <- shiny::NS(id)
  list(
    shiny::helpText("Select your dataset, trait to show, and other options"),
    mainParInput(ns("main_par")), # "group", "LOD_thr"
    mainParUI(ns("main_par")),    # "which_trait", "selected_chr"
    # Add input for the cis/trans plot
    cisTransPlotInput(ns("cistrans")),
    peakInput(ns("peak_list")))  # "which_peak", "alleles" actionButton
}
#' @rdname qtlApp
#' @export
qtlOutput <- function(id) {
  ns <- shiny::NS(id)
  list(
    cisTransPlotUI(ns("cistrans")),
    scanlyOutput(ns("scanly")),
    scanlyUI(ns("scanly")),
    peakUI(ns("peak_list")),
    peakOutput(ns("peak_list"))
  )
}
