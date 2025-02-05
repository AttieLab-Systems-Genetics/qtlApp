#' This is the viewer script to view previously made qtl2 scans for various traits
#' it's designed to look at additive, interactive, and difference (interactive-additive) plots
#' all the scan traces are pre-computed, so this just reads the table of relevant values
#' and plots them
#'
#' @author Chris Emfinger, PhD. Couldn't really view anything well
#' 
#' In addition to the relevant package citations, parts of the datatable code 
#' followed formats shown on 
#' https://clarewest.github.io/blog/post/making-tables-shiny/
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
file_directory <- read.csv("C:\\Users\\chris\\Downloads\\file_directory\\file_index.csv")

# load the chromosomal breaks
chr_breaks <- read.csv("C:\\Users\\chris\\Downloads\\file_directory\\chromosomal_sep_mm11.csv")

# load the annotation data for the different traits
# for making the object
# annotation_list <- list()
# annotation_list$isoforms <- mediation_isoforms$RNA_info$Liver[c(1,3,ncol(mediation_isoforms$RNA_info$Liver))]
# annotation_list$genes <- mediation_genes$RNA_info$Liver[c(1,6)]
# saveRDS(annotation_list, file="annotation_list.rds")
annotation_list <- readRDS("C:\\Users\\chris\\Downloads\\file_directory\\annotation_list.rds")

# load the markers
markers <- readRDS(file.path("\\\\research.drive.wisc.edu\\mkeller3\\General\\mapping_files_for_mark\\CHTC_dietDO_markers_RDSgrcm39.rds"))

# make the dataset list-------------------------------------------------------------------------
file_directory$group <- paste0(file_directory$diet," ",file_directory$trait_compartment, " ", file_directory$trait_type, ", ", file_directory$scan_type)

# set microfunctions============================================================================
# plot the QTL if doing a new scan
QTL_plot_visualizer <- function(qtl.temp, phenotype, LOD_thr, mrkrs){
  qtl_info <- unique(qtl.temp[c("ID_code", # the ID code to match it with everything else,
                                "sex", # the sexes
                                "diet", # the diet
                                "covars_additive", # the additive covariates
                                "covars_interactive", # the interactive covariates
                                "trait_compartment", # the trait compartment
                                "transform_type", # the transform type
                                "trait_type", # the trait type
                                "scan_type", # the scan type
                                "scan_date")])
  #remove the metadata
  qtl.temp <- qtl.temp[-which(colnames(qtl.temp)==c("ID_code", # the ID code to match it with everything else,
                          "sex", # the sexes
                          "diet", # the diet
                          "covars_additive", # the additive covariates
                          "covars_interactive", # the interactive covariates
                          "trait_compartment", # the trait compartment
                          "transform_type", # the transform type
                          "trait_type", # the trait type
                          "scan_type", # the scan type
                          "scan_date"))]
  
  #get the rowname/colname info
  clnms <- colnames(qtl.temp)
  markers_ID <- rownames(qtl.temp)
  qtl.temp <- data.frame(qtl.temp)
  
  #derive the chromosomal location info from the markers
  #add markers
  qtl.temp$markers <- markers_ID
  clnms2 <- colnames(mrkrs)
  which_mrkr_clmn <-grep(pattern="marker", clnms2)
  which_pos_clnm <- grep(pattern="bp", clnms2)
  clnms2[which_mrkr_clmn]<-"markers"
  clnms2[which_pos_clnm]<-"position"
  colnames(mrkrs)<-clnms2
  mrkrs2 <- mrkrs[c("markers","chr","position")]
  mrkrs2$position <- mrkrs2$position/(10^6)
  
  qtl.temp <- qtl.temp %>%
    left_join(., mrkrs2, by=c("markers"="markers"))
  qtl.temp$order <- as.numeric(qtl.temp$chr)
  
  qtl.temp[which(qtl.temp$chr=="X"),"order"]<-20
  qtl.temp[which(qtl.temp$chr=="Y"),"order"]<-21
  qtl.temp[which(qtl.temp$chr=="M"),"order"]<-22
  qtl.temp$chr <- qtl.temp$order
  
  qtl_plot_obj <- qtl.temp %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(position)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    
    # Add this info to the initial dataset
    left_join(qtl.temp, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(order, position) %>%
    mutate( BPcum=position+tot)
  
  #x axis mods: we do not want to display the cumulative 
  #position of SNP in bp, but just show the chromosome name instead.
  axisdf = qtl_plot_obj %>% group_by(order) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  axisdf$order[which(axisdf$order==20)]<-"X"
  axisdf$order[which(axisdf$order==21)]<-"Y"
  axisdf$order[which(axisdf$order==22)]<-"M"
  axisdf$order <- factor(axisdf$order, levels=c(as.character(1:19),"X","Y","M"))
  
  colnames(qtl_plot_obj)[1]<-"LOD"
  qtl_plot_obj$chr[which(qtl_plot_obj$chr==20)]<-"X"
  qtl_plot_obj$chr[which(qtl_plot_obj$chr==21)]<-"Y"
  qtl_plot_obj$chr[which(qtl_plot_obj$chr==22)]<-"M"
  
  #create the ggplot object
  #add label
  grob <- grobTree(textGrob(paste0(phenotype), x=0.1,  y=.95, hjust=0,
                            gp=gpar(col="red", fontsize=10, fontface="italic")))
  grob2 <- grobTree(textGrob(paste0("Additive covariates: ", paste0(qtl_info$covars_additive[1])), x=0.1,  y=.9, hjust=0,
                             gp=gpar(col="blue", fontsize=10, fontface="italic")))
  grob3 <- grobTree(textGrob(paste0("Interactive covariates: ", paste0(qtl_info$covars_interactive[1])), x=0.1,  y=.85, hjust=0,
                             gp=gpar(col="orange", fontsize=10, fontface="italic")))
  grob4 <- grobTree(textGrob(paste0("Diet: ", paste0(qtl_info$diet[1])), x=0.1,  y=.8, hjust=0,
                             gp=gpar(col="darkblue", fontsize=10, fontface="italic")))
  grob5 <- grobTree(textGrob(paste0("Sexes: ", paste0(qtl_info$sex[1])), x=0.1,  y=.75, hjust=0,
                             gp=gpar(col="darkblue", fontsize=10, fontface="italic")))
  grob6 <- grobTree(textGrob(paste0(paste0(qtl_info$com[1])," ",qtl_info$trait_compartment[1], " ",qtl_info$transform_type[1]," ",qtl_info$trait_type[1]), x=0.1,  y=.7, hjust=0,
                             gp=gpar(col="darkblue", fontsize=10, fontface="italic")))
  plot_QTL <- ggplot(qtl_plot_obj, aes(x=BPcum, y=LOD))+
    # Show all points
    geom_line(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
    scale_color_manual(values = rep(c("black", "darkgrey"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$order, breaks= axisdf$center ) +
    #scale_y_continuous(expand = c(0, 0) ) +     
    # remove space between plot area and x axis
    ylim(0,c(max(qtl_plot_obj$LOD)+.25*max(qtl_plot_obj$LOD)))+
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    
    #change the axis lines
    ggplot2::theme(axis.line = element_line(colour = "black"))+
    
    #change the axis labels
    xlab("Chromosomal position")+
    ylab("LOD")+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20))+
    
    #add significant LOD
    geom_hline(aes(yintercept=LOD_thr,
                   linetype="Chosen LOD"), color="black")+
    
    #add annotations
    annotation_custom(grob)+
    annotation_custom(grob2)+
    annotation_custom(grob3)+
    annotation_custom(grob4)+
    annotation_custom(grob5)+
    annotation_custom(grob6)
  
  #show(plot_QTL)
  plots <- list(plot_QTL, qtl_plot_obj)
  return(plots)
}

# find the trait
trait_scan <- function(file_dir, selected_dataset, selected_trait){
  file_dir <- subset(file_dir, group == selected_dataset)
  file_dir <- subset(file_dir, file_type == "scans")
  scan_data <- fread(file_dir$File_path[1], 
                     select = c("V1", #the row names
                                "ID_code", # the ID code to match it with everything else,
                                "sex", # the sexes
                                "diet", # the diet
                                "covars_additive", # the additive covariates
                                "covars_interactive", # the interactive covariates
                                "trait_compartment", # the trait compartment
                                "transform_type", # the transform type
                                "trait_type", # the trait type
                                "scan_type", # the scan type
                                "scan_date", # the scan date
                                selected_trait)) # the trait
  
  scan_data <- data.frame(scan_data)
  rownames(scan_data)<-scan_data$V1
  scan_data <- scan_data[-which(colnames(scan_data)=="V1")]
  return(scan_data)
}

# find the peaks
peak_finder <- function(file_dir, selected_dataset){
  file_dir <- subset(file_dir, group == selected_dataset)
  file_dir <- subset(file_dir, file_type == "peaks")
  peaks <- read.csv(file_dir$File_path)
  peaks <- peaks %>% relocate("marker.id","lodcolumn","chr","pos","lod")
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
                         choices = peaks$marker.id,
                         options = list(maxItems = 1,
                                        maxOptions = 5),
                         server = TRUE
    )
  })
  
  # show the strain effects-----------------------------------------------------------
  observeEvent(input$alleles, {
    req(input$selected_dataset)
    req(input$which_trait)
    req(input$which_peak)
    if (scan_reactive$scan_data$scan_type[1] == "additive"){
      # set trait
      chosen_trait <- str_split(input$which_trait, pattern=" [(]")[[1]][2]
      chosen_trait <- str_split(chosen_trait, pattern="[)]")[[1]][1]
      # set peaks
      peaks <- peak_finder(file_directory, input$selected_dataset)
      peaks <- subset(peaks, trait == chosen_trait)
      peaks <- subset(peaks, marker.id == input$which_peak)
      peaks <- peaks[c("marker.id","A","B","C","D","E","F","G","H")]
      colnames(peaks)[2:9]<-c(c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB"))
      peaks <- reshape2::melt(peaks, id=c("marker.id"))
      # set colors
      newClrs <- c(1,2,3,4,5,6,7,8)
      newClrs[1]<-"#000000"
      newClrs[2]<-"#96989A"
      newClrs[3]<-"#E69F00"
      newClrs[4]<-"#0072B2"
      newClrs[5]<-"#619BFF"
      newClrs[6]<-"#009E73"
      newClrs[7]<-"#D55E00"
      newClrs[8]<-"#CC79A7"
      names(newClrs)<-c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
      # set plot
      plot_alleles <- ggplot(data=peaks, aes(x=marker.id, y=value, color = variable))
      # set as point plot- this worked and bar didn't
      plot_alleles <- plot_alleles + geom_point(size=10) +
        # add custom colors
        scale_color_manual(values = c(newClrs))+
        # Custom the theme:
        theme_bw() +
        theme( 
          legend.text=element_text(size=18),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )+
        
        #change the axis lines
        ggplot2::theme(axis.line = element_line(colour = "black"))+
        
        #change the axis labels
        xlab("Marker ID")+
        ylab("Founder allele effect")+
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=20))+
        
        #add significant LOD
        geom_hline(aes(yintercept=0), color="black")
  
      
      output$allele_effects <- renderPlot({plot_alleles})
    }
    if (scan_reactive$scan_data$scan_type[1] != "additive"){
      output$allele_effects <- renderText({"Not additive"})
    }
  })
  
  
}

# set app ====================================================================================================
shiny::shinyApp(ui = ui, server = server)