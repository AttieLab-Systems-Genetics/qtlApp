mainParServer <- function(id, file_directory, annotation_list) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Set data type.
    set_type <- shiny::reactive({
      shiny::req(input$selected_dataset)
      type <- subset(file_directory, group==input$selected_dataset)$trait_type[1]
      assign("set_test", type, envir = .GlobalEnv)
      type
    })
    trait_list <- traitServer("trait", file_directory, annotation_list, set_type)
    trait_id <- shiny::reactive({
      switch(shiny::req(set_type()),
        Genes    = "gene.id",
        Isoforms = "transcript.id",
        Clinical = "data_name")
    })

    # Update trait choices.
    shiny::observeEvent(shiny::req(set_type(), trait_list(), trait_id()), {
      choices = paste0("(", trait_list()[[trait_id()]], ")")
      if(set_type() == "Clinical") {
        choices = paste(annotation_list$clinical$ID, choices)
      } else {
        choices = paste(annotation_list$genes$symbol, choices)
      }
      shiny::updateSelectizeInput(session, "which_trait",
        choices = choices, options = list(maxItems = 1, maxOptions = 5), server = TRUE)
    })

    output$returns <- shiny::renderPrint({
      cat("selected_dataset =", input$selected_dataset,
          "\nwhich_trait =", input$which_trait)
    })
    # Return.
    input
  })
}
mainParInput <- function(id) {
  ns <- shiny::NS(id)
  shiny::selectizeInput(ns("selected_dataset"),
    label = "Choose a dataset to display",
    choices = unique(file_directory$group),
    multiple = FALSE,
    options = list(placeholder = 'Search...'))
}
mainParUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::selectizeInput(ns("which_trait"),
    label = "Choose the trait",
    choices = NULL,
    multiple = FALSE,
    options = list(placeholder = 'Search...'))
}
mainParOutput <- function(id) {
  ns <- shiny::NS(id)
  shiny::verbatimTextOutput(ns("returns"))
}
mainParApp <- function(id) {
    source("qtlSetup.R")
    source("traitServer.R")
    ui <- bslib::page(
        title = "Test Main Par",
        mainParInput("main_par"),
        mainParUI("main_par"),
        mainParOutput("main_par")
    )
    server <- function(input, output, session) {
        mainParServer("main_par", file_directory, annotation_list)
    }
    shiny::shinyApp(ui = ui, server = server)
}