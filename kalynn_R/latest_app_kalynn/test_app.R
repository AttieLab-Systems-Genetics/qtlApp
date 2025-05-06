library(shiny)
library(bslib)

# Ensure all module files are sourced if not handled by Docker's app.R automatically
# (Assuming Dockerfile's main app.R sources these from the R/ directory)

ui <- bslib::page_sidebar(
  title = "Test Import and MainPar Modules",
  sidebar = bslib::sidebar("side_panel",
    mainParInput("main_par") # UI for main parameters
  ),
  mainParUI("main_par")     # Additional UI for main parameters (like trait selection)
)

server <- function(input, output, session) {
  import_data <- importServer("import") # Assuming importServer is self-contained or uses a default path
  main_par_outputs <- mainParServer("main_par", import_data)
  
  # You can add observers here to print reactive values from main_par_outputs for debugging
  # observe({
  #   print(paste("Selected Dataset:", main_par_outputs$selected_dataset()))
  #   print(paste("Selected Trait:", main_par_outputs$which_trait()))
  # })
}

shinyApp(ui = ui, server = server) 



