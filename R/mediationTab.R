#' Mediation Tab UI
#'
#' Creates a single dropdown listing mediator datasets (the SECOND token in
#' filenames that match: "<first>_additive_<second>_mediator_all_mediators.csv"),
#' grouped by interaction type inferred from the FIRST token:
#' Additive (All mice), Sex (Female), Sex (Male), Diet (HC), Diet (HF).
#'
#' @return Shiny UI for the Mediation tab
#' @importFrom shiny tags div selectInput h6 conditionalPanel
#' @export
mediation_tab_ui <- function() {
    # Discover mediation files
    candidate_dirs <- c("/data/dev/miniViewer_3.0", file.path(getwd(), "data"))
    files <- unlist(lapply(candidate_dirs, function(d) {
        if (!dir.exists(d)) {
            return(character(0))
        }
        list.files(d, pattern = "_mediator_all_mediators\\.csv$", full.names = TRUE, recursive = TRUE)
    }))

    if (length(files) == 0) {
        return(
            shiny::div(
                style = "padding: 12px; color: #7f8c8d;",
                shiny::tags$em("No mediation datasets detected.")
            )
        )
    }

    parsed <- lapply(files, function(f) {
        nm <- basename(f)
        m <- regexec("^(.+?)_additive_(.+?)_mediator_all_mediators\\.csv$", nm)
        grp <- regmatches(nm, m)[[1]]
        if (length(grp) == 3) list(first = grp[2], second = grp[3], file = f) else NULL
    })
    parsed <- Filter(Negate(is.null), parsed)

    if (length(parsed) == 0) {
        return(
            shiny::div(
                style = "padding: 12px; color: #7f8c8d;",
                shiny::tags$em("No mediation datasets parsed.")
            )
        )
    }

    # Build a data.frame of parsed entries for easier processing
    first_vals <- vapply(parsed, function(x) x$first, character(1))
    second_vals <- vapply(parsed, function(x) x$second, character(1))
    file_vals <- vapply(parsed, function(x) x$file, character(1))

    df <- data.frame(
        first = first_vals,
        second = second_vals,
        file = file_vals,
        stringsAsFactors = FALSE
    )

    # Build grouped choices using SECOND token labels, grouped by interaction from FIRST
    add_mask <- grepl("_all_mice$", df$first, ignore.case = TRUE)
    sex_f_mask <- grepl("qtlxsex_in_female_mice$", df$first, ignore.case = TRUE)
    sex_m_mask <- grepl("qtlxsex_in_male_mice$", df$first, ignore.case = TRUE)
    diet_hc_mask <- grepl("qtlxdiet_in_HC_mice$", df$first, ignore.case = TRUE)
    diet_hf_mask <- grepl("qtlxdiet_in_HF_mice$", df$first, ignore.case = TRUE)

    # Build grouped choices separated by interaction type
    additive_choices <- list()
    sex_choices <- list()
    diet_choices <- list()

    add_seconds <- sort(unique(df$second[add_mask]))
    if (length(add_seconds) > 0) {
        additive_choices[["Additive (All mice)"]] <- stats::setNames(
            paste0("add_all_mice:", add_seconds),
            add_seconds
        )
    }

    sex_f_seconds <- sort(unique(df$second[sex_f_mask]))
    if (length(sex_f_seconds) > 0) {
        sex_choices[["Interactive - Sex (Female)"]] <- stats::setNames(
            paste0("sex_female:", sex_f_seconds),
            sex_f_seconds
        )
    }

    sex_m_seconds <- sort(unique(df$second[sex_m_mask]))
    if (length(sex_m_seconds) > 0) {
        sex_choices[["Interactive - Sex (Male)"]] <- stats::setNames(
            paste0("sex_male:", sex_m_seconds),
            sex_m_seconds
        )
    }

    diet_hc_seconds <- sort(unique(df$second[diet_hc_mask]))
    if (length(diet_hc_seconds) > 0) {
        diet_choices[["Interactive - Diet (HC)"]] <- stats::setNames(
            paste0("diet_hc:", diet_hc_seconds),
            diet_hc_seconds
        )
    }

    diet_hf_seconds <- sort(unique(df$second[diet_hf_mask]))
    if (length(diet_hf_seconds) > 0) {
        diet_choices[["Interactive - Diet (HF)"]] <- stats::setNames(
            paste0("diet_hf:", diet_hf_seconds),
            diet_hf_seconds
        )
    }

    if (length(additive_choices) == 0 && length(sex_choices) == 0 && length(diet_choices) == 0) {
        return(
            shiny::div(
                style = "padding: 12px; color: #7f8c8d;",
                shiny::tags$em("No mediation datasets available.")
            )
        )
    }

    # Defaults per panel
    default_add <- {
        if (!is.null(additive_choices[["Additive (All mice)"]])) additive_choices[["Additive (All mice)"]][[1]] else NULL
    }
    default_sex <- {
        avail <- unlist(sex_choices, use.names = FALSE)
        if (length(avail) > 0) avail[[1]] else NULL
    }
    default_diet <- {
        avail <- unlist(diet_choices, use.names = FALSE)
        if (length(avail) > 0) avail[[1]] else NULL
    }

    {
        panels <- list()
        # Robust JS conditions that work even if the input is namespaced or missing
        cond_find_interaction <- function(target) {
            if (target == "none") {
                return("(function(){var k=Object.keys(input||{}),v=null;for(var i=0;i<k.length;i++){if(k[i].endsWith('interaction_type_selector')){v=input[k[i]];break;}}return(v===null||v===undefined||v==='none');})()")
            } else if (target == "sex") {
                return("(function(){var k=Object.keys(input||{}),v=null;for(var i=0;i<k.length;i++){if(k[i].endsWith('interaction_type_selector')){v=input[k[i]];break;}}return(v==='sex');})()")
            } else if (target == "diet") {
                return("(function(){var k=Object.keys(input||{}),v=null;for(var i=0;i<k.length;i++){if(k[i].endsWith('interaction_type_selector')){v=input[k[i]];break;}}return(v==='diet');})()")
            } else {
                return("false")
            }
        }

        if (length(additive_choices) > 0) {
            panels[[length(panels) + 1]] <- shiny::conditionalPanel(
                condition = cond_find_interaction("none"),
                shiny::div(
                    style = "padding: 12px;",
                    shiny::selectInput(
                        inputId = "mediation_dataset_selector_additive",
                        label = "Mediator dataset (Additive)",
                        choices = additive_choices,
                        selected = default_add,
                        width = "100%"
                    )
                )
            )
        }
        if (length(sex_choices) > 0) {
            panels[[length(panels) + 1]] <- shiny::conditionalPanel(
                condition = cond_find_interaction("sex"),
                shiny::div(
                    style = "padding: 12px;",
                    shiny::selectInput(
                        inputId = "mediation_dataset_selector_sex",
                        label = "Mediator dataset (Sex interaction)",
                        choices = sex_choices,
                        selected = default_sex,
                        width = "100%"
                    )
                )
            )
        }
        if (length(diet_choices) > 0) {
            panels[[length(panels) + 1]] <- shiny::conditionalPanel(
                condition = cond_find_interaction("diet"),
                shiny::div(
                    style = "padding: 12px;",
                    shiny::selectInput(
                        inputId = "mediation_dataset_selector_diet",
                        label = "Mediator dataset (Diet interaction)",
                        choices = diet_choices,
                        selected = default_diet,
                        width = "100%"
                    )
                )
            )
        }

        if (length(panels) == 0) {
            shiny::div(
                style = "padding: 12px; color: #7f8c8d;",
                shiny::tags$em("No mediation datasets available.")
            )
        } else {
            do.call(shiny::tagList, panels)
        }
    }
}
