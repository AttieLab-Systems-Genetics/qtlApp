#' Merge App
#'
#' @param id identifier for shiny reactive
#' @param prefix static prefix for filename
#' @param main_par input parameters from calling routine
#' @param download_list reactiveValues with  postfix,plot,table
#' @return nothing 
#'
#' @importFrom shiny isolate moduleServer NS
#'             observeEvent reactive reactiveValues shinyApp  
#' @importFrom bslib card card_header page_sidebar sidebar
#' @export
mergeApp <- function(id) {
  ui <- bslib::page_sidebar(
    title = "Test Merge",
    sidebar = bslib::sidebar("side_panel", width = 400,
      downloadInput("download"), 
      downloadOutput("download") 
    ),
    bslib::card(
      bslib::card_header("Merge Preview"),
      downloadUI("download")
    )
  )
  server <- function(input, output, session) { 
    primary_list <- shiny::reactiveValues(
      filename = shiny::reactive("panelID_instanceID"),
      plots = shiny::reactiveValues(
        none = shiny::reactive(plot_null("nothing")),
        some = shiny::reactive(plot_null("something"))),
      tables = shiny::reactiveValues(
        twelve = shiny::reactive(matrix(1:12,nrow=3)),
        twenty = shiny::reactive(matrix(1:20,nrow=4)))
    )
    secondary_list <- shiny::reactiveValues(
      filename = shiny::reactive("secondary"),
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

    merged_plots <- shiny::isolate({
      out <- primary_list$plots   
      
      s_plots <- secondary_list$plots 
      if (!is.null(s_plots) && inherits(s_plots, "reactivevalues")) {
        for(i in names(s_plots)) {
          out[[i]] <- s_plots[[i]] 
        }
      }
      out
    })
    
    merged_tables <- shiny::isolate({
      out <- primary_list$tables 
      s_tables <- secondary_list$tables 
      if (!is.null(s_tables) && inherits(s_tables, "reactivevalues")) {
        for(i in names(s_tables)) {
          out[[i]] <- s_tables[[i]] 
        }
      }
      out
    })
    
    merged_list <- shiny::reactiveValues(
      filename = shiny::reactive({
        if (shiny::is.reactive(primary_list$filename)) {
          primary_list$filename()
        } else {
          primary_list$filename
        }
      }),
      plots = merged_plots,
      tables = merged_tables
    )
    
    # Return
    merged_list
  })
}
