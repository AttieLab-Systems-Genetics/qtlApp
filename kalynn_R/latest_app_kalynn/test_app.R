library(shiny)
library(bslib)

ui <- 
    bslib::page_sidebar(
    title = "Test Trait Module",
    sidebar = bslib::sidebar("side_panel",
      mainParInput("main_par"),
      traitUI("trait"),
      traitOutput("trait"),
      mainParUI("main_par"),       # "which_trait"
      peakInput("peak"),           # "which_peak"
      downloadInput("download"),   # downloadButton, filename
      downloadOutput("download")
    ),
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
 #traitServer("trait", main_par, import)
 peak_list <- peakServer("peak", main_par, import)
 downloadServer("download", peak_list)
}

shinyApp(ui = ui, server = server) 



