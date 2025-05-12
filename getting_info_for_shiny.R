#' for getting necessary information about the new data
#' 
dataset <- readRDS(
  "/Volumes/adattie/General/Diet DO study/Brian_Yandell_check_here/dataset.liver.hc.genes.rds")
markers <- readRDS(
  "/Volumes/adattie/General/Users/mkeller3/Projects/DO_diet_mapping_files/true_input_files/gigamuga_markers_v2.RDS")

# Tunable paramters: 
#   shiny::selectInput with multiple = TRUE
qtl_chrs <- as.character(c(1:19),"X")
gene_chrs <- as.character(c(1:19),"X")
#   shiny::sliderInput with step = 5, min = 0, max = max(peaks$lod), round = TRUE
min_lod <- 50

# for additve scans=======================================================================
peaks <- dataset$lod.peaks$additive

# if not using the isoforms, add the gene information
gene_info <- dataset$annot.mrna
colnames(gene_info)[which(colnames(gene_info)=="chr")]<-"gene_chr"
colnames(gene_info)[which(colnames(gene_info)=="start")]<-"gene_start"
colnames(gene_info)[which(colnames(gene_info)=="end")]<-"gene_end"
peaks <- dplyr::inner_join(peaks, gene_info, by=c("gene_id"="gene_id")) |>
  dplyr::filter(gene_chr %in% c(1:19,"X"))

# Find min pos by chr.
markermin <- markers |>
  dplyr::filter(!is.na(bp_mm10),
                chr %in% c(1:19,"X")) |>
  dplyr::group_by(chr) |>
  dplyr::summarise(pos = min(bp_mm10),
                   .groups = "drop") |>
  dplyr::ungroup() |>
  dplyr::mutate(chr = factor(chr, c(1:19,"X"))) |>
  dplyr::arrange(chr)
markermin <- array(markermin$pos, dimnames = list(markermin$chr))

# Subtract min pos from gene_start by chr
peaks <- 
  dplyr::mutate(peaks,
                gene_start = pmax(0, gene_start - markermin[gene_chr]))

m <- match(peaks$marker.id, markers$marker)
peaks <- dplyr::mutate(peaks,
  qtl_chr = markers$chr[m],
  qtl_chr = factor(qtl_chr, c("X", 19:1)),
  qtl_pos = pmax(0, markers$bp_mm10[m] - markermin[qtl_chr]),
  gene_chr = factor(gene_chr, c(1:19, "X")))

ggplot2::ggplot(
  peaks |> 
    dplyr::filter(qtl_chr %in% qtl_chrs,
                  gene_chr %in% gene_chrs,
                  lod >= min_lod) |>
    dplyr::mutate(gene_start = gene_start * 1e-6,
                  qtl_pos = qtl_pos * 1e-6)) +
  ggplot2::aes(x = gene_start, y = qtl_pos, col = lod) +
  ggplot2::facet_grid(qtl_chr ~ gene_chr,
                      scales = "free") +
  ggplot2::geom_point(size = 0.25) +
  ggplot2::theme(panel.spacing = ggplot2::unit(0,'lines'),
                 axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
