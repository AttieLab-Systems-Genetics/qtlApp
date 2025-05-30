 # Main Application Module Server
mainAppServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      
      message("Main app server function started.")
      
      
      import_reactives <- importServer("import") # Correct call
      # -----------------------------------------------------------------
      
      # Reactive for error message
      error_message_reactive <- reactiveVal("") 
      
      
      main_par_reactives <- mainParServer("main_par", import_reactives) # Use plain ID
      # -------------------------------------------------------
      
     
      scan_list_reactives <- scanServer("scan_list", main_par_reactives, import_reactives) # Use plain ID
      # -----------------------------------------
                                          
      
      peak_table_direct_reactive <- peakServer("peak", main_par_reactives, import_reactives) # Use plain ID
      # -----------------------------------------
      
      
      scanlyServer("scanly", main_par_reactives, scan_list_reactives, peak_table_direct_reactive) # Use plain ID
      # --------------------------------------------
                   
      
      filename_reactive <- reactive({
          shiny::req(scan_list_reactives$filename()) 
          scan_list_reactives$filename() 
      })
      
      download_items <- reactiveValues(
          filename = filename_reactive, 
          plots = reactiveValues(
              scan_plot = reactive(scan_list_reactives$scan_plot()) 
          ),
          tables = reactiveValues(
              scan_table = reactive(scan_list_reactives$scan_table()) 
          )
      )

      # Call downloadServer instances
      downloadServer("download_qtl", download_items) # Use plain ID
      downloadServer("download_allele", download_items) # Use plain ID
      # ---------------------------------------------

      # Render the error message text
      output$error_message <- renderText({ error_message_reactive() })
      
      message("Main app server function setup complete.")
    }
  )
}
