#' LOD Threshold Controls Module (Overview/Peaks)
#'
#' Centralized threshold helpers:
#' - get_scan_threshold(): for main LOD scans (reference line)
#' - get_diff_threshold(): for difference plots (reference lines)
#'
#' Finalized values:
#' - Additive: 7.5
#' - Diet interactive: 10.5
#' - Sex interactive: 10.5
#' - Sexbydiet interactive: 15.7
#' - Diet difference: 4.1
#' - Sex difference: 4.1
#' - Sexbydiet difference: 9.5
#'
#' Provides a dynamic LOD threshold slider for overview/peaks UI slider
#' for additive vs interactive difference scans. Exposes a debounced
#' `lod_threshold()` reactive.
#'
#' @param id Module ID
#' @param interaction_type_reactive Reactive that returns 'none' | 'sex' | 'diet' | 'sex_diet'
#'
#' @return list(lod_threshold = reactive)
#' @importFrom shiny moduleServer NS uiOutput renderUI sliderInput reactive
#' @export
lodThresholdUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("lod_threshold_slider"))
}

#' @rdname lodThresholdUI
#' @export
lodThresholdServer <- function(id, interaction_type_reactive) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$lod_threshold_slider <- shiny::renderUI({
      itype <- if (!is.null(interaction_type_reactive)) interaction_type_reactive() else "none"
      scan_info <- switch(itype,
        sex = list(type = "Sex Difference", min = 4.1),
        diet = list(type = "Diet Difference", min = 4.1),
        sex_diet = list(type = "Sex x Diet Difference", min = 9.5),
        list(type = "Additive", min = 7.5)
      )
      current_val <- input$LOD_thr
      default_val <- if (!is.null(current_val) && current_val >= scan_info$min) current_val else scan_info$min
      shiny::sliderInput(ns("LOD_thr"),
        label = paste0("LOD Threshold (", scan_info$type, " scan):"),
        min = scan_info$min, max = 20, value = default_val, step = 0.5,
        width = "100%"
      )
    })

    lod_threshold <- shiny::reactive({
      itype <- if (!is.null(interaction_type_reactive)) interaction_type_reactive() else "none"
      fallback <- if (itype == "none") 7.5 else 4.1
      input$LOD_thr %||% fallback
    }) %>% shiny::debounce(300)

    return(list(lod_threshold = lod_threshold))
  })
}

#' Get main LOD scan threshold by interaction type
#' @param interaction_type one of 'none','sex','diet','sex_diet'
#' @return numeric threshold
#' @export
get_scan_threshold <- function(interaction_type) {
  if (is.null(interaction_type) || identical(interaction_type, "none")) return(7.5)
  if (identical(interaction_type, "sex")) return(10.5)
  if (identical(interaction_type, "diet")) return(10.5)
  if (identical(interaction_type, "sex_diet")) return(15.7)
  7.5
}

#' Get difference plot threshold by interaction type
#' @param interaction_type one of 'sex','diet','sex_diet'
#' @return numeric threshold
#' @export
get_diff_threshold <- function(interaction_type) {
  if (is.null(interaction_type)) return(NA_real_)
  if (identical(interaction_type, "sex")) return(4.1)
  if (identical(interaction_type, "diet")) return(4.1)
  if (identical(interaction_type, "sex_diet")) return(9.5)
  NA_real_
}


