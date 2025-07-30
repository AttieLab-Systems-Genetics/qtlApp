#' Main Application UI
#'
#' Defines the main user interface for the QTL Scan Visualizer application.
#'
#' @return A Shiny UI object.
#' @importFrom shiny NS tagList uiOutput h3 h4 h5 p div icon hr actionButton selectizeInput selectInput conditionalPanel
#' @importFrom bslib page_sidebar sidebar navset_pill nav_panel card card_header card_body
#' @importFrom shinyjs useShinyjs
#' @importFrom htmltools tags
#' @importFrom shinycssloaders withSpinner
#' @importFrom plotly plotlyOutput
#' @export
mainUI <- function() {
    bslib::page_sidebar(
        shinyjs::useShinyjs(),
        tags$head(
            # Viewport meta tag for responsive design
            tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0, user-scalable=yes"),

            # Autosizing CSS and JavaScript
            tags$link(rel = "stylesheet", type = "text/css", href = "autosize.css"),
            tags$script(src = "autosize.js"),

            # Custom app styles
            tags$style(custom_css)
        ),

        # Top navigation bar
        div(
            style = "background: linear-gradient(135deg, #2c3e50, #3498db); padding: 15px; margin-bottom: 20px; border-radius: 8px; text-align: center;",
            h3("QTL Scan Visualizer",
                style = "color: white; margin: 0; font-weight: bold;"
            )
        ),
        sidebar = bslib::sidebar(
            width = 600, # Increased sidebar width a bit more for better screen coverage

            # Dataset Category Selection - above all tabs
            div(
                style = "padding: 15px; margin-bottom: 20px; background: linear-gradient(135deg, #f8f9fa, #e9ecef); border-radius: 8px; border: 2px solid #3498db;",
                h4("Dataset Category", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold; text-align: center;"),
                shiny::selectInput(shiny::NS("app_controller", "dataset_category_selector"),
                    NULL,
                    choices = c("Loading..." = ""),
                    width = "100%"
                ),
                p("Select the type of biological data to analyze",
                    style = "font-size: 12px; color: #6c757d; margin: 5px 0 0 0; text-align: center;"
                )
            ),

            # Tabbed sidebar content
            bslib::navset_pill(
                id = "sidebar_tabs",
                selected = "LOD peaks", # Set default tab to LOD peaks

                # Tab 1: Overview Plot (Manhattan/Cis-Trans) - now the default tab
                bslib::nav_panel(
                    "LOD peaks",
                    div(
                        style = "padding: 10px;",

                        # LOD Threshold Control
                        h5("Peak Filtering", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"),
                        # Dynamic LOD threshold slider that updates based on scan type
                        uiOutput(shiny::NS("app_controller", "lod_threshold_slider")),
                        p("Filters peaks shown in the plot below",
                            style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 15px 0;"
                        ),

                        # Peak Selection Dropdown for future peak differences analysis
                        hr(style = "border-top: 1px solid #bdc3c7; margin: 15px 0;"),
                        h5("🎯 Peak Analysis", style = "color: #2c3e50; margin-bottom: 10px; font-weight: bold;"),
                        shiny::uiOutput(shiny::NS("app_controller", "peak_selection_sidebar")),
                        p("Control interaction analysis for the overview plots below",
                            style = "font-size: 11px; color: #7f8c8d; margin: 5px 0 15px 0;"
                        ),
                        hr(style = "border-top: 1px solid #bdc3c7; margin: 15px 0;"),
                        div(
                            id = "overview-plot-container",
                            class = "overview-plot-container",
                            style = "height: 65vh; min-height: 400px; max-height: 800px; border: 1px solid #bdc3c7; border-radius: 5px;",
                            shiny::uiOutput(shiny::NS("app_controller", "conditional_plot_ui"))
                        ),
                        p("Click on points to view detailed LOD scans. Plot titles show dataset and analysis type.",
                            style = "font-size: 11px; color: #7f8c8d; margin: 10px 0 0 0; text-align: center;"
                        )
                    )
                ),

                # Tab 2: Data Search and Controls - now secondary
                bslib::nav_panel(
                    "Data Search",
                    div(
                        style = "padding: 10px;",

                        # Trait search section
                        h5("🔍 Trait Search", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
                        selectizeInput(shiny::NS("app_controller", "trait_search_input"),
                            "Search for traits:",
                            choices = NULL,
                            selected = NULL,
                            multiple = FALSE,
                            options = list(
                                placeholder = "Type to search (e.g., Gapdh, Insulin, PI_38_3)",
                                maxItems = 1,
                                maxOptions = 10,
                                create = FALSE
                            ),
                            width = "100%"
                        ),
                        div(
                            style = "text-align: center; margin-top: 10px;",
                        )
                    )
                ),
            ),

            # Horizontal separator
            hr(style = "border-top: 2px solid #3498db; margin: 20px 0;"),

            # Additional Analyses section below the main tabs
            h5("📈 Additional Analyses", style = "color: #2c3e50; margin-bottom: 15px; font-weight: bold;"),
            bslib::navset_pill(
                id = "additional_analyses_tabs",

                # Profile Plot tab
                bslib::nav_panel(
                    "Profile Plot",
                    profilePlotUI(shiny::NS("app_controller", "profile_plot_module"))
                ),

                # Correlation tab
                bslib::nav_panel(
                    "Correlation",
                    div(
                        style = "padding: 10px;",
                        div(
                            id = "correlation-plot-container",
                            class = "sidebar-plot-container",
                            style = "height: 50vh; min-height: 300px; max-height: 500px; border: 1px solid #bdc3c7; border-radius: 5px; overflow: hidden;",
                            shinycssloaders::withSpinner(
                                plotly::plotlyOutput(shiny::NS("app_controller", "correlation_plot_output"),
                                    height = "100%", width = "100%"
                                )
                            )
                        ),
                        p("Correlation analysis visualization coming soon",
                            style = "font-size: 11px; color: #7f8c8d; margin: 10px 0 0 0; text-align: center;"
                        )
                    )
                )
            )
        ),

        # Simplified main area - now with a dynamic header
        bslib::card(
            id = "lod_scan_card",
            bslib::card_header(shiny::uiOutput(shiny::NS("app_controller", "main_plot_title"))),
            bslib::card_body(
                shiny::uiOutput(shiny::NS("app_controller", "lod_scan_plot_ui_placeholder"))
            )
        )
    )
}
