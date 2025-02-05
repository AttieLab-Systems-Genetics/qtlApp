#' This is the viewer script to view previously made qtl2 scans for various traits
#' it's designed to look at additive, interactive, and difference (interactive-additive) plots
#' all the scan traces are pre-computed, so this just reads the table of relevant values
#' and plots them
#'
#' @author Chris Emfinger, PhD. Couldn't really view anything well
#'
#' In addition to the relevant package citations, parts of the datatable code
#' followed formats shown on
#' https://clarewest.github.io/blog/post/making-tables-shiny/x
#'
#' some help with shiny
#' https://shiny.rstudio.com/gallery
#' https://shiny.posit.co/r/gallery
#' posit.cloud/spaces/298214/join




# load libraries=============================================================
require("rstudioapi")
require("dplyr")
require("stringr")
require("tidyverse")
require("BiocManager")
require("ggplot2")
require("qtl2")
require("grid")
require("ggrepel")
require("gridGraphics")
require("ggpubr")
require("shiny")
require("shinyFiles")
require("bslib")
require("spsComps")
require("DT")
require("shinyjs")
require("shinycssloaders")
require("data.table")
require("reshape2")


# set shiny options==========================================================================
options(shiny.maxRequestSize = 20000*1024^2)  # Increase to 20GB needed for Genoprobs, etc


# load data==================================================================================
# load the file directory
#setwd(selectDirectory())


# load data==================================================================================
# load the file directory
#setwd(selectDirectory())
file_directory <- read.csv("/data/dev/miniViewer_3.0/file_index.csv")


# load the chromosomal breaks
chr_breaks <- read.csv("/data/dev/miniViewer_3.0/chromosomal_sep_mm11.csv")


# load the annotation data for the different traits
# for making the object
# annotation_list <- list()
# annotation_list$isoforms <- mediation_isoforms$RNA_info$Liver[c(1,3,ncol(mediation_isoforms$RNA_info$Liver))]
# annotation_list$genes <- mediation_genes$RNA_info$Liver[c(1,6)]
# saveRDS(annotation_list, file="annotation_list.rds")
annotation_list <- readRDS("/data/dev/miniViewer_3.0/annotation_list.rds")


# load the markers
markers <- readRDS(file.path("/data/dev/miniViewer_3.0/CHTC_dietDO_markers_RDSgrcm39.rds"))


file_directory$group <- paste0(file_directory$diet," ",file_directory$trait_compartment, " ", file_directory$trait_type, ", ", file_directory$scan_type)


# set microfunctions============================================================================
# plot the QTL if doing a new scan
QTL_plot_visualizer <- function(qtl.temp, phenotype, LOD_thr, mrkrs) {
  # The LOD values are already in the "LOD" column
  qtl.temp <- qtl.temp[, c("marker", "LOD")]
  colnames(qtl.temp) <- c("markers", "LOD")
  
  mrkrs2 <- mrkrs[c("marker", "chr", "bp_grcm39")]
  colnames(mrkrs2) <- c("markers", "chr", "position")
  mrkrs2$position <- as.numeric(mrkrs2$position)/(10^6)
  
  # Join the data
  qtl.temp <- merge(qtl.temp, mrkrs2, by.x="markers", by.y="markers", all.x=FALSE)
  qtl.temp <- na.omit(qtl.temp)
  
  # Convert chr to numeric
  qtl.temp$chr <- as.character(qtl.temp$chr)
  qtl.temp$order <- sapply(qtl.temp$chr, function(x) {
    if (x == "X") return(20)
    else if (x == "Y") return(21)
    else if (x == "M") return(22)
    else return(as.numeric(x))
  })
  
  # Remove any rows where chromosome conversion failed
  qtl.temp <- qtl.temp[!is.na(qtl.temp$order), ]
  qtl.temp$chr <- qtl.temp$order
  
  # Create plot object
  qtl_plot_obj <- qtl.temp %>%
    group_by(chr) %>%
    summarise(chr_len = max(position)) %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
    select(-chr_len) %>%
    left_join(qtl.temp, ., by = c("chr" = "chr")) %>%
    arrange(order, position) %>%
    mutate(BPcum = position + tot)
  
  # Create axis labels
  axisdf = qtl_plot_obj %>%
    group_by(order) %>%
    summarize(center = (max(BPcum) + min(BPcum))/2)
  
  # Convert chromosome labels
  axisdf$order[axisdf$order == 20] <- "X"
  axisdf$order[axisdf$order == 21] <- "Y"
  axisdf$order[axisdf$order == 22] <- "M"
  axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19),"X","Y","M"))
  
  # Create plot with x-axis
  plot_QTL <- ggplot(qtl_plot_obj, aes(x=BPcum, y=LOD)) +
    geom_line(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
    scale_color_manual(values = rep(c("black", "darkgrey"), 22)) +
    scale_x_continuous(
      label = axisdf$order, 
      breaks = axisdf$center,
      expand = expansion(mult = 0.15)  # Increased padding
    ) +
    ylim(0, max(qtl_plot_obj$LOD, na.rm=TRUE) * 1.25) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      axis.text.x = element_text(
        angle = 45,           # Rotate labels 45 degrees
        hjust = 1,           # Adjust horizontal position
        vjust = 1,           # Adjust vertical position
        margin = margin(t = 10)  # Add margin at top of labels
      ),
      # Add more bottom margin to accommodate rotated labels
      plot.margin = margin(b = 40, l = 20, r = 20, t = 20, unit = "pt")
    ) +
    xlab("Chromosome") +
    ylab("LOD") +
    geom_hline(aes(yintercept=LOD_thr), color="black", linetype="dashed")
  
  plots <- list(plot_QTL, qtl_plot_obj)
  return(plots)
}
# find the trait
trait_scan <- function(file_dir, selected_dataset, selected_trait) {
  # Subset the data
  file_dir <- subset(file_dir, group == selected_dataset)
  file_dir <- subset(file_dir, file_type == "scans")
  
  if (nrow(file_dir) == 0) {
    stop("No matching files found for the selected dataset")
  }
  
  if (!file.exists(file_dir$File_path[1])) {
    stop("File does not exist: ", file_dir$File_path[1])
  }
  
  scan_data <- tryCatch({
    # Read all data first
    data <- fread(file_dir$File_path[1], drop = "Which_mice")
    
    # Filter for rows where Phenotype matches selected_trait
    data <- data[data$Phenotype == selected_trait, ]
    
    if (nrow(data) == 0) {
      stop("No data found for phenotype: ", selected_trait)
    }
    
    data
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
  
  return(scan_data)
}


# find the peaks
peak_finder <- function(file_dir, selected_dataset){
  file_dir <- subset(file_dir, group == selected_dataset)
  file_dir <- subset(file_dir, file_type == "peaks")
  peaks <- read.csv(file_dir$File_path)
  peaks <- peaks %>% relocate("marker","lodcolumn","chr","pos","lod")
  #peaks <- peaks %>% relocate("lodcolumn","chr","pos","lod")
  colnames(peaks)[which(colnames(peaks)=="lodcolumn")]<-"trait"
  return(peaks)
}


# set UI======================================================================================================
ui <- page_sidebar(
  # sets the page
  # sets the title
  titlePanel("Pre-scanned QTL visualizer, implemented for Diet DO study"),
  #start shinyjs
  useShinyjs(),
  # sets sidebar
  sidebar = sidebar(
    id = "side_panel",
    helpText(
      "Select your dataset, trait to show, and other options"
    ),
    
    selectizeInput(
      inputId = "selected_dataset",
      label = "Choose a dataset to display",
      choices = unique(file_directory$group),
      multiple = FALSE,
      options = list(
        placeholder = 'Search...'
      )
    ),
    sliderInput(
      "LOD_thr",
      label = "LOD threshold for evaluation",
      min = 4,
      max = 20,
      value = 7.5,
      round = TRUE
    ),
    #helpText("Display table of trait ID information (e.g. gene symbols)? "),
    #actionButton("search_ID", "Show IDs"),
    selectizeInput(
      inputId = "which_trait",
      label = "Choose the trait",
      choices = NULL,
      multiple = FALSE,
      options = list(
        placeholder = 'Search...'
      )
    ),
    actionButton("scan", "Show the LOD scan"),
    helpText("Choose a peak to see the strain effects. This only applies to the additive scans."),
    selectizeInput(
      inputId = "which_peak",
      label = "Choose peak",
      choices = NULL,
      multiple = TRUE,
      options = list(
        placeholder = 'Search...'
      )),
    actionButton("alleles", "Show Effects")
  ),
  # sets main panel
  mainPanel(
    # set scrollbar
    div(style = "overflow-y: scroll;"),
    position = "right",
    #set starting layout
    #card(
    #  card_header("Trait info"),
    #  DT::dataTableOutput('search_results')
    #),
    card(
      card_header("Datasets available"),
      DT::dataTableOutput('available_data')
    ),
    card(
      card_header("LOD profile"),
      plotOutput("scan_plot", click = "plot_click") %>% withSpinner(color="#0dc5c1"),
      verbatimTextOutput("scan_points")
    ),
    card(
      card_header("Peaks"),
      DT::dataTableOutput("peaks")
    ),
    card(
      card_header("Strain effects"),
      plotOutput("allele_effects") %>% withSpinner(color="#0dc5c1")
    ),
  ))


# set server==================================================================================================
server <- function(input, output, session) {
  # set the display data for all of the options--------------------------------------------------------
  output$available_data <- DT::renderDT({DT::datatable(
    file_directory,
    options = list(paging = TRUE,    ## paginate the output
                   pageLength = 5,  ## number of rows to output for each page
                   scrollX = TRUE,   ## enable scrolling on X axis
                   scrollY = TRUE,   ## enable scrolling on Y axis
                   autoWidth = TRUE, ## use smart column width handling
                   server = TRUE,   ## use client-side processing
                   dom = 'Bfrtip',
                   buttons = c('csv', 'excel'),
                   columnDefs = list(list(targets = '_all', className = 'dt-center'),
                                     list(targets = c(0, 8, 9), visible = FALSE))
    ),
    extensions = 'Buttons',
    selection = 'single', ## enable selection of a single row
    filter = 'bottom',              ## include column filters at the bottom
    rownames = TRUE                ##  show row numbers/names
  )
  })
  # create the reactive objects---------------------------------------------------------------------
  # covariate reactive values
  annot_reactive <- reactiveValues(
    annot_obj = NULL,
    trait_list = NULL
  )
  # scan reactive
  scan_reactive <- reactiveValues(
    scan_data = NULL
  )
  # create the trait list---------------------------------------------------------------------
  observeEvent(input$selected_dataset,{
    # update expression
    req(input$selected_dataset)
    file_directory <- subset(file_directory, group==input$selected_dataset)
    set_type <- file_directory$trait_type[1]
    assign("set_test",set_type, envir = .GlobalEnv)
    if (set_type == "Genes"){
      updateSelectizeInput(session,
                           "which_trait",
                           choices = paste0(annotation_list$genes$symbol, " (",annotation_list$genes$gene.id,")"),
                           options = list(maxItems = 1,
                                          maxOptions = 5),
                           server = TRUE
      )
      annot_reactive$trait_list <- data.frame(annotation_list$genes)
    }
    if (set_type == "Isoforms"){
      updateSelectizeInput(session,
                           "which_trait",
                           choices = paste0(annotation_list$genes$symbol, " (",annotation_list$isoforms$transcript.id,")"),
                           options = list(maxItems = 1,
                                          maxOptions = 5),
                           server = TRUE
      )
      annot_reactive$trait_list <- data.frame(annotation_list$isoforms)
    }
    if (set_type == "Clinical"){
      updateSelectizeInput(session,
                           "which_trait",
                           choices = paste0(annotation_list$clinical$ID, " (", annotation_list$clinical$data_name,")"),
                           options = list(maxItems = 1,
                                          maxOptions = 5),
                           server = TRUE
      )
      annot_reactive$trait_list <- data.frame(annotation_list$clinical)
    }
    assign("test_annot",annot_reactive$trait_list, envir = .GlobalEnv)
  })
  #observeEvent(input$search_ID, {
  #  req(input$selected_dataset)
  #  test_annot <- annot_reactive$trait_list
  ##  output$search_results <- DT::renderDT({DT::datatable(
  #    data.frame(annot_reactive$trait_list),
  #    options = list(paging = TRUE,    ## paginate the output
  #                   pageLength = 5,  ## number of rows to output for each page
  #                   scrollX = TRUE,   ## enable scrolling on X axis
  #                   scrollY = TRUE,   ## enable scrolling on Y axis
  #                   autoWidth = TRUE, ## use smart column width handling
  #                   server = TRUE,   ## use client-side processing
  #                   dom = 'Bfrtip',
  #                   buttons = c('csv', 'excel'),
  #                   columnDefs = list(list(targets = '_all', className = 'dt-center'),
  #                                     list(targets = c(0, 8, 9), visible = FALSE))
  #    ),
  #    extensions = 'Buttons',
  #    selection = 'single', ## enable selection of a single row
  #    filter = 'bottom',              ## include column filters at the bottom
  #    rownames = TRUE                ##  show row numbers/names
  #  )
  #  })
  #})
  
  # create the scan list---------------------------------------------------------------------
  observeEvent(input$scan, {
    # update expression
    req(input$selected_dataset)
    req(input$which_trait)
    chosen_trait <- str_split(input$which_trait, pattern=" [(]")[[1]][2]
    chosen_trait <- str_split(chosen_trait, pattern="[)]")[[1]][1]
    scans <- trait_scan(file_directory, input$selected_dataset, chosen_trait)
    scan_reactive$scan_data <- scans
    scan_plot <- QTL_plot_visualizer(scans, input$which_trait, input$LOD_thr, markers)
    output$scan_plot <- renderPlot({scan_plot[[1]]})
    output$scan_points <-  renderPrint({
      nearPoints(scan_plot[[2]], input$plot_click, xvar = "BPcum", yvar = "LOD", threshold = 10, maxpoints = 1,
                 addDist = TRUE) })
  })
  
  # get the peaks info ------------------------------------------------------------------
  observeEvent(input$selected_dataset, {
    peaks <- peak_finder(file_directory, input$selected_dataset)
    output$peaks <- DT::renderDT({DT::datatable(
      peaks,
      options = list(paging = TRUE,    ## paginate the output
                     pageLength = 5,  ## number of rows to output for each page
                     scrollX = TRUE,   ## enable scrolling on X axis
                     scrollY = TRUE,   ## enable scrolling on Y axis
                     autoWidth = TRUE, ## use smart column width handling
                     server = TRUE,   ## use client-side processing
                     dom = 'Bfrtip',
                     buttons = c('csv', 'excel'),
                     columnDefs = list(list(targets = '_all', className = 'dt-center'),
                                       list(targets = c(0, 8, 9), visible = FALSE))
      ),
      extensions = 'Buttons',
      selection = 'single', ## enable selection of a single row
      filter = 'bottom',              ## include column filters at the bottom
      rownames = TRUE                ##  show row numbers/names
    )
    })
  })
  # get the peak selected for strain effects-------------------------------------------
  observeEvent(input$which_trait, {
    req(input$selected_dataset)
    req(input$which_trait)
    chosen_trait <- str_split(input$which_trait, pattern=" [(]")[[1]][2]
    chosen_trait <- str_split(chosen_trait, pattern="[)]")[[1]][1]
    peaks <- peak_finder(file_directory, input$selected_dataset)
    peaks <- subset(peaks, trait == chosen_trait)
    updateSelectizeInput(session,
                         "which_peak",
                         choices = peaks$marker,
                         options = list(maxItems = 1,
                                        maxOptions = 5),
                         server = TRUE
    )
  })


#MAKE THIS REACTIVE INSTEAD ??????

  # show the strain effects-----------------------------------------------------------
  observeEvent(c(input$alleles, input$which_peak), {  # Now observing both the button and peak selection
    req(input$selected_dataset)
    req(input$which_trait)
    req(input$which_peak)
    
    peaks <- peak_finder(file_directory, input$selected_dataset)
    
    # Check if scan data exists and is additive
    if (peaks$intcovar[1] == "none") {
        # set trait
        chosen_trait <- str_split(input$which_trait, pattern=" [(]")[[1]][2]
        chosen_trait <- str_split(chosen_trait, pattern="[)]")[[1]][1]
        
        # set peaks
        peaks <- subset(peaks, trait == chosen_trait)
        peaks <- subset(peaks, marker == input$which_peak)  # Changed from marker.id to marker
        
        # Check if we have data after filtering
        if (nrow(peaks) > 0) {
            # Select and rename columns
            peaks <- peaks[c("marker","A","B","C","D","E","F","G","H")]
            colnames(peaks)[2:9] <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
            
            # Reshape data
            peaks <- reshape2::melt(peaks, id.vars = "marker")
            
            # Define colors
            newClrs <- c(
                "AJ" = "#000000",
                "B6" = "#96989A",
                "129" = "#E69F00",
                "NOD" = "#0072B2",
                "NZO" = "#619BFF",
                "CAST" = "#009E73",
                "PWK" = "#D55E00",
                "WSB" = "#CC79A7"
            )
            
            # Create plot
            plot_alleles <- ggplot(data = peaks, aes(x = marker, y = value, color = variable)) +
                geom_point(size = 10) +
                scale_color_manual(values = newClrs) +
                theme_bw() +
                theme(
                    legend.text = element_text(size = 18),
                    panel.border = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.text = element_text(size = 18),
                    axis.title = element_text(size = 20)
                ) +
                labs(x = "Marker ID", y = "Founder allele effect", color = "Strain") +
                geom_hline(yintercept = 0, color = "black")
            
            output$allele_effects <- renderPlot({
                print(plot_alleles)
            })
        } else {
            output$allele_effects <- renderText({"No data available for selected peak"})
        }
    } else {
        output$allele_effects <- renderText({"Strain effects only available for additive scans"})
    }
  })
}


# set app ====================================================================================================
shiny::shinyApp(ui = ui, server = server)
