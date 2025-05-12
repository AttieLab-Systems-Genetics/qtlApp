#' Merge App
#'
#' @param id identifier for shiny reactive
#' @param primary_list,secondary_list reactiveValues objects
#'
#' @importFrom shiny isolate isTruthy moduleServer NS
#'             observeEvent reactive reactiveValues shinyApp textAreaInput
#' @importFrom bslib card card_header page_sidebar sidebar
#' @export
mergeApp <- function(id) {
  ui <- bslib::page_sidebar(
    title = "Test Merge",
    sidebar = bslib::sidebar("side_panel", width = 400,
      mergeUI("merged_list"),    # filename (optional)
      downloadInput("download"), # plot_table, inputs for Plots or Tables
      downloadOutput("download") # downloadButton, filename
    ),
    bslib::card(
      bslib::card_header("Merge Preview"),
      downloadUI("download")
    )
  )
  server <- function(input, output, session) { 
    # Test sets
    primary_list <- shiny::reactiveValues(
      filename = "panelID_instanceID",
      no_download = c("plots_skip"),
      plots = shiny::reactiveValues(
        skip = shiny::reactive(plot_null("skip")),
        none = shiny::reactive(plot_null("nothing")),
        some = shiny::reactive(plot_null("something"))),
      tables = shiny::reactiveValues(
        twelve = shiny::reactive(matrix(1:12,nrow=3)),
        twenty = shiny::reactive(matrix(1:20,nrow=4)))
    )
    secondary_list <- shiny::reactiveValues(
      filename = "secondary",
      no_download = c("tables_ten"),
      plots = shiny::reactiveValues(
        a = shiny::reactive(plot_null("aaaa")),
        b = shiny::reactive(plot_null("bbbb"))),
      tables = shiny::reactiveValues(
        ten = shiny::reactive(matrix(1:10,nrow=2)),
        eight = shiny::reactive(matrix(1:8,nrow=4)))
    )
    merged_list <- mergeServer("merged_list", primary_list, secondary_list)
    downloadServer("download", merged_list)
  }
  shiny::shinyApp(ui, server)
}
#' @rdname mergeApp
#' @export
mergeServer <- function(id, primary_list, secondary_list) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Create new reactiveValues for `plots` and `tables`.
    merged_plots <- shiny::reactiveValues()
    merged_tables <- shiny::reactiveValues()
    shiny::isolate({
      # Populate with `primary_list` `plots` and `tables`.
      for(i in names(primary_list$plots)) {
        merged_plots[[i]] <- primary_list$plots[[i]]
      }
      for(i in names(primary_list$tables)) {
        merged_tables[[i]] <- primary_list$tables[[i]]
      }
      # Merge in data from `secondary_list`.
      for(i in names(secondary_list$plots)) {
        if(!(i %in% names(merged_plots)))
          merged_plots[[i]] <- secondary_list$plots[[i]]
      }
      for(i in names(secondary_list$tables)) {
        if(!(i %in% names(merged_tables)))
           merged_tables[[i]] <- secondary_list$tables[[i]]
      }
    })
    output$filename <- renderUI({
      filename <- primary_list$filename
      shiny::textAreaInput(ns("filename"), "Test File Basename:", filename)
    })
    # Create new reactiveValues `merged_list`.
    merged_list <- shiny::reactiveValues(
      filename = shiny::isolate({
        if(shiny::isTruthy(input$filename))
          input$filename
        else
          primary_list$filename
      }),
      no_download = shiny::isolate({
        unique(c(primary_list$no_download, secondary_list$no_download))
      }),
      plots = merged_plots,
      tables = merged_tables)
    # Return
    merged_list
  })
}
#' @rdname mergeApp
#' @export
mergeUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("filename"))
}
