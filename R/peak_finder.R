# find the peaks
peak_finder <- function(file_dir, selected_dataset){
  file_dir <- subset(file_dir, group == selected_dataset)
  file_dir <- subset(file_dir, file_type == "peaks")
  peaks <- read.csv(file_dir$File_path)
  peaks <- peaks |> dplyr::relocate("marker","lodcolumn","chr","pos","lod")
  #peaks <- peaks |> dplyr::relocate("lodcolumn","chr","pos","lod")
  colnames(peaks)[which(colnames(peaks) == "lodcolumn")]<-"trait"
  return(peaks)
}