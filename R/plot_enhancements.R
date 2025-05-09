# Enhanced plot styling and functionality module
# This module contains functions for creating modern, interactive plots with enhanced styling

# Function to create a modern plot theme
create_modern_theme <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#f0f0f0", size = 0.5),
      panel.grid.minor = element_line(color = "#f0f0f0", size = 0.25),
      axis.line = element_line(color = "#2c3e50", size = 0.5),
      axis.text = element_text(color = "#2c3e50", size = 12),
      axis.title = element_text(color = "#2c3e50", size = 14, face = "bold"),
      plot.title = element_text(color = "#2c3e50", size = 16, face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "#7f8c8d", size = 14, hjust = 0),
      legend.title = element_text(color = "#2c3e50", size = 12, face = "bold"),
      legend.text = element_text(color = "#2c3e50", size = 11),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "#f8f9fa", color = NA),
      strip.text = element_text(color = "#2c3e50", size = 12, face = "bold")
    )
}

# Function to create a modern color palette
create_modern_palette <- function(n) {
  colors <- c(
    "#3498db", # Blue
    "#2c3e50", # Dark Blue
    "#e74c3c", # Red
    "#f1c40f", # Yellow
    "#9b59b6", # Purple
    "#1abc9c", # Turquoise
    "#e67e22", # Orange
    "#34495e"  # Navy
  )
  
  if (n <= length(colors)) {
    return(colors[1:n])
  } else {
    return(colorRampPalette(colors)(n))
  }
}

# Function to create a modern hover template
create_hover_template <- function() {
  paste(
    "<b>%{x}</b><br>",
    "Value: %{y:.2f}<br>",
    "<extra></extra>"
  )
}

# Function to create a modern plot layout
create_plot_layout <- function(title = NULL, subtitle = NULL, xaxis_title = NULL, yaxis_title = NULL) {
  list(
    title = list(
      text = title,
      font = list(size = 24, color = "#2c3e50"),
      x = 0.05
    ),
    xaxis = list(
      title = list(
        text = xaxis_title,
        font = list(size = 16, color = "#2c3e50")
      ),
      showgrid = TRUE,
      gridcolor = "#f0f0f0",
      zeroline = TRUE,
      zerolinecolor = "#2c3e50",
      zerolinewidth = 1
    ),
    yaxis = list(
      title = list(
        text = yaxis_title,
        font = list(size = 16, color = "#2c3e50")
      ),
      showgrid = TRUE,
      gridcolor = "#f0f0f0",
      zeroline = TRUE,
      zerolinecolor = "#2c3e50",
      zerolinewidth = 1
    ),
    plot_bgcolor = "white",
    paper_bgcolor = "white",
    margin = list(t = 80, b = 60, l = 60, r = 40),
    showlegend = TRUE,
    legend = list(
      bgcolor = "white",
      bordercolor = "#f0f0f0",
      borderwidth = 1,
      x = 1,
      y = 1,
      xanchor = "right",
      yanchor = "top"
    ),
    hovermode = "closest",
    hoverlabel = list(
      bgcolor = "white",
      bordercolor = "#2c3e50",
      font = list(size = 14, color = "#2c3e50")
    )
  )
}

# Function to create a modern plot configuration
create_plot_config <- function() {
  list(
    displayModeBar = TRUE,
    displaylogo = FALSE,
    modeBarButtonsToRemove = c("lasso2d", "select2d"),
    toImageButtonOptions = list(
      format = "png",
      filename = "plot",
      height = 800,
      width = 1200,
      scale = 2
    )
  )
}

# Function to create a modern plot annotation
create_plot_annotation <- function(text, x, y, showarrow = TRUE) {
  list(
    text = text,
    x = x,
    y = y,
    showarrow = showarrow,
    font = list(size = 14, color = "#2c3e50"),
    bgcolor = "white",
    bordercolor = "#2c3e50",
    borderwidth = 1,
    borderpad = 4,
    xanchor = "left",
    yanchor = "bottom"
  )
}

# Function to create a modern plot shape
create_plot_shape <- function(type, x0, x1, y0, y1, line = list(color = "#e74c3c", width = 2, dash = "dot")) {
  list(
    type = type,
    x0 = x0,
    x1 = x1,
    y0 = y0,
    y1 = y1,
    line = line,
    fillcolor = "rgba(231, 76, 60, 0.1)"
  )
}

# Function to create a modern plot marker
create_plot_marker <- function(size = 8, color = "#3498db", symbol = "circle") {
  list(
    size = size,
    color = color,
    symbol = symbol,
    line = list(
      color = "white",
      width = 1
    )
  )
}

# Function to create a modern plot line
create_plot_line <- function(width = 2, color = "#3498db", dash = "solid") {
  list(
    width = width,
    color = color,
    dash = dash
  )
}

# Function to create a modern plot area
create_plot_area <- function(color = "#3498db", opacity = 0.2) {
  list(
    color = color,
    opacity = opacity
  )
}

# Function to create a modern plot error bar
create_plot_error_bar <- function(color = "#2c3e50", width = 1, thickness = 1) {
  list(
    color = color,
    width = width,
    thickness = thickness
  )
}


# Function to customize chromosome plot appearance
customize_chromosome_plot <- function(p, selected_chr = "All") {
  # Base customizations
  p <- p + theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
  
  # Set appropriate x-ßaxis label based on view
  x_label <- if(selected_chr == "All") "Chromosome" else paste("Position on Chr", selected_chr, "(Mb)")
  p <- p + labs(x = x_label, y = "LOD Score")
  
  return(p)
}

# Function to add peak markers to a plot
add_peak_markers <- function(p, peak_data, selected_chr = "All", color = "#e74c3c", size = 3, shape = 18) {
  if (nrow(peak_data) == 0) return(p)
  
  xvar <- if(selected_chr == "All") "BPcum" else "position"
  
  p + geom_point(
    data = peak_data,
    aes_string(x = xvar, y = "LOD"),
    color = color,
    size = size,
    shape = shape
  )
}

# Function to create highlight bands for QTL regions
create_highlight_regions <- function(p, regions, selected_chr = "All", fill = "rgba(231, 76, 60, 0.1)", 
                                     outline = "#e74c3c", alpha = 0.2) {
  if (length(regions) == 0) return(p)
  
  xvar <- if(selected_chr == "All") "BPcum" else "position"
  
  for (i in seq_along(regions)) {
    region <- regions[[i]]
    p <- p + geom_rect(
      xmin = region$start, 
      xmax = region$end,
      ymin = -Inf, 
      ymax = Inf,
      fill = fill, 
      color = outline,
      alpha = alpha
    )
  }
  
  return(p)
}

# Function to create modern plotly hover text
create_hover_text <- function(data, trait_column = "trait", chr_column = "chr", pos_column = "pos", lod_column = "lod") {
  paste0(
    "<b>", data[[trait_column]], "</b><br>",
    "Chr: ", data[[chr_column]], "<br>",
    "Position: ", round(data[[pos_column]], 2), " Mb<br>",
    "LOD: ", round(data[[lod_column]], 2)
  )
} 