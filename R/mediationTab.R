#' Mediation Tab UI
#'
#' Creates a single dropdown listing mediator datasets (the SECOND token in
#' filenames that match: "<first>_additive_<second>_mediator_all_mediators.csv"),
#' grouped by interaction type inferred from the FIRST token:
#' Additive (All mice), Sex (Female), Sex (Male), Diet (HC), Diet (HF).
#'
#' @param current_category Optional string of the current dataset category
#'   (e.g., "Liver Genes", "Liver Isoforms", "Liver Lipids", "Clinical Traits",
#'   "Plasma Metabolites"). When provided, the dropdown will be filtered to
#'   only show mediation datasets whose FIRST token family matches the category
#'   (e.g., for "Liver Genes" show files starting with "liver_genes_" and
#'   "liver_isoforms_" so users can mediate genes against isoforms too).
#' @return Shiny UI for the Mediation tab
#' @importFrom shiny tags div selectInput h6 conditionalPanel
#' @export
mediation_tab_ui <- function(current_category = NULL) {
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

    # Optionally filter by current dataset category to avoid showing unrelated families
    if (!is.null(current_category) && nzchar(current_category)) {
        # Map category to allowed FIRST-family prefixes
        allowed_prefixes <- switch(current_category,
            "Liver Genes" = c("liver_genes_", "liver_isoforms_"),
            "Liver Isoforms" = c("liver_isoforms_", "liver_genes_"),
            "Liver Lipids" = c("liver_lipids_", "liver_genes_", "liver_isoforms_"),
            "Clinical Traits" = c("clinical_traits_", "liver_genes_", "liver_isoforms_"),
            "Plasma Metabolites" = c("plasma_metabolites_", "liver_genes_", "liver_isoforms_"),
            NULL
        )
        if (!is.null(allowed_prefixes)) {
            keep <- Reduce(`|`, lapply(allowed_prefixes, function(pfx) startsWith(df$first, pfx)))
            df <- df[keep, , drop = FALSE]
        }
    }

    # Build grouped choices using SECOND token labels, grouped by interaction from FIRST
    add_mask <- grepl("_all_mice$", df$first, ignore.case = TRUE)
    sex_f_mask <- grepl("qtlxsex_in_female_mice$", df$first, ignore.case = TRUE)
    sex_m_mask <- grepl("qtlxsex_in_male_mice$", df$first, ignore.case = TRUE)
    diet_hc_mask <- grepl("qtlxdiet_in_HC_mice$", df$first, ignore.case = TRUE)
    diet_hf_mask <- grepl("qtlxdiet_in_HF_mice$", df$first, ignore.case = TRUE)

    # Humanize SECOND token labels (e.g., "liver_genes" -> "Liver Genes")
    prettify_second <- function(x) {
        labs <- c(
            liver_genes = "Liver Genes",
            liver_isoforms = "Liver Isoforms",
            liver_lipids = "Liver Lipids",
            clinical_traits = "Clinical Traits",
            plasma_metabolites = "Plasma Metabolites"
        )
        out <- labs[tolower(x)]
        ifelse(is.na(out), x, out)
    }

    # Build grouped choices separated by interaction type
    additive_choices <- list()
    sex_choices <- list()
    diet_choices <- list()

    add_rows <- df[add_mask, , drop = FALSE]
    if (nrow(add_rows) > 0) {
        add_rows <- add_rows[!duplicated(add_rows$second), , drop = FALSE]
        additive_choices[["Additive (All mice)"]] <- stats::setNames(
            add_rows$file,
            prettify_second(add_rows$second)
        )
    }

    sex_f_rows <- df[sex_f_mask, , drop = FALSE]
    if (nrow(sex_f_rows) > 0) {
        sex_f_rows <- sex_f_rows[!duplicated(sex_f_rows$second), , drop = FALSE]
        sex_choices[["Interactive - Sex (Female)"]] <- stats::setNames(
            sex_f_rows$file,
            prettify_second(sex_f_rows$second)
        )
    }

    sex_m_rows <- df[sex_m_mask, , drop = FALSE]
    if (nrow(sex_m_rows) > 0) {
        sex_m_rows <- sex_m_rows[!duplicated(sex_m_rows$second), , drop = FALSE]
        sex_choices[["Interactive - Sex (Male)"]] <- stats::setNames(
            sex_m_rows$file,
            prettify_second(sex_m_rows$second)
        )
    }

    diet_hc_rows <- df[diet_hc_mask, , drop = FALSE]
    if (nrow(diet_hc_rows) > 0) {
        diet_hc_rows <- diet_hc_rows[!duplicated(diet_hc_rows$second), , drop = FALSE]
        diet_choices[["Interactive - Diet (HC)"]] <- stats::setNames(
            diet_hc_rows$file,
            prettify_second(diet_hc_rows$second)
        )
    }

    diet_hf_rows <- df[diet_hf_mask, , drop = FALSE]
    if (nrow(diet_hf_rows) > 0) {
        diet_hf_rows <- diet_hf_rows[!duplicated(diet_hf_rows$second), , drop = FALSE]
        diet_choices[["Interactive - Diet (HF)"]] <- stats::setNames(
            diet_hf_rows$file,
            prettify_second(diet_hf_rows$second)
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
        # Use the known app namespace id used across the app
        cond_add <- "input['app_controller-interaction_type_selector'] == 'none'"
        cond_sex <- "input['app_controller-interaction_type_selector'] == 'sex'"
        cond_diet <- "input['app_controller-interaction_type_selector'] == 'diet'"

        if (length(additive_choices) > 0) {
            panels[[length(panels) + 1]] <- shiny::conditionalPanel(
                condition = cond_add,
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
                condition = cond_sex,
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
                condition = cond_diet,
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

        content <- list()
        if (length(panels) == 0) {
            content[[length(content) + 1]] <- shiny::div(
                style = "padding: 12px; color: #7f8c8d;",
                shiny::tags$em("No mediation datasets available.")
            )
        } else {
            content[[length(content) + 1]] <- do.call(shiny::tagList, panels)
        }
        # Mediation plot placeholder (rendered in app server)
        content[[length(content) + 1]] <- shiny::div(
            style = "padding: 8px 12px;",
            plotly::plotlyOutput(shiny::NS("app_controller", "mediation_plot"), height = "380px") %>%
                shinycssloaders::withSpinner(type = 8, color = "#8e44ad")
        )
        do.call(shiny::tagList, content)
    }
}
