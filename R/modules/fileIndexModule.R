#' File Index Module (dataset/category utilities)
#'
#' Centralizes access to file_index.csv (import_reactives()$file_directory),
#' validates required columns, and provides a stable selected dataset group
#' with HC_HF additive auto-selection.
#'
#' @param id Module ID
#' @param import_reactives Reactive list that includes $file_directory
#' @param selected_category_reactive Reactive that returns the current dataset category string
#'
#' @return list(file_index, dataset_category, selected_dataset_group)
#' @importFrom shiny moduleServer NS reactive req
#' @importFrom data.table as.data.table
#' @export
fileIndexServer <- function(id, import_reactives, selected_category_reactive) {
  shiny::moduleServer(id, function(input, output, session) {
    # ns <- session$ns  # not needed

    file_index_dt <- shiny::reactive({
      shiny::req(import_reactives()$file_directory)
      dt <- data.table::as.data.table(import_reactives()$file_directory)
      shiny::validate(
        shiny::need("dataset_category" %in% names(dt), "Error: 'dataset_category' column missing in file_index.csv."),
        shiny::need("group" %in% names(dt), "Error: 'group' column missing in file_index.csv.")
      )
      dt
    })

    dataset_category <- shiny::reactive({
      selected_category_reactive()
    })

    selected_dataset_group <- shiny::reactive({
      selected_cat <- dataset_category()
      shiny::req(selected_cat, file_index_dt())

      datasets_in_category <- file_index_dt()[dataset_category == selected_cat, ]
      specific_datasets_choices <- unique(datasets_in_category$group)

      # Find the appropriate HC_HF dataset (additive) or the first available dataset
      hc_hf_dataset <- NULL

      if (selected_cat %in% c("Liver Genes", "Liver Lipids", "Clinical Traits", "Plasma Metabolites", "Liver Isoforms")) {
        pattern <- switch(selected_cat,
          "Liver Genes" = "^HC_HF.*Liver.*Genes",
          "Liver Lipids" = "^HC_HF.*Liver.*Lipid",
          "Clinical Traits" = "^HC_HF.*Clinical",
          "Plasma Metabolites" = "^HC_HF.*Plasma.*Metabol",
          "Liver Isoforms" = "^HC_HF.*Liver.*Isoform"
        )
        hc_hf_dataset <- specific_datasets_choices[grepl(pattern, specific_datasets_choices, ignore.case = TRUE) &
          !grepl("interactive", specific_datasets_choices, ignore.case = TRUE)]
      }

      if (!is.null(hc_hf_dataset) && length(hc_hf_dataset) > 0) {
        return(hc_hf_dataset[1])
      }

      if (length(specific_datasets_choices) > 0) {
        return(specific_datasets_choices[1])
      }

      return(NULL)
    })

    return(list(
      file_index = file_index_dt,
      dataset_category = dataset_category,
      selected_dataset_group = selected_dataset_group
    ))
  })
}


