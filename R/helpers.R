# Get trait type
get_trait_type <- function(import, selected_group) {
    file_directory <- import$file_directory
    file_directory <- subset(file_directory, group == selected_group)
    trait_type <- tolower(file_directory$trait_type[1])
    trait_type
}
get_trait_list <- function(import, trait_type) {
    annotation_list <- import$annotation_list
    annotation_list[[trait_type]]
}
