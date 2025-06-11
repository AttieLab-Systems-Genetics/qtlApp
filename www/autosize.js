// QTL App Autosizing JavaScript - Ensures consistent behavior across all browsers and devices

(function() {
  'use strict';

  // Configuration
  const CONFIG = {
    RESIZE_DEBOUNCE_MS: 250,
    PLOT_RESIZE_DELAY_MS: 100,
    MIN_SIDEBAR_WIDTH: 350,
    MAX_SIDEBAR_WIDTH: 600,
    MOBILE_BREAKPOINT: 768
  };

  // Utility functions
  const debounce = (func, wait) => {
    let timeout;
    return function executedFunction(...args) {
      const later = () => {
        clearTimeout(timeout);
        func(...args);
      };
      clearTimeout(timeout);
      timeout = setTimeout(later, wait);
    };
  };

  const isMobile = () => window.innerWidth <= CONFIG.MOBILE_BREAKPOINT;

  // Force viewport meta tag if not present
  const ensureViewportMeta = () => {
    let viewport = document.querySelector('meta[name="viewport"]');
    if (!viewport) {
      viewport = document.createElement('meta');
      viewport.name = 'viewport';
      viewport.content = 'width=device-width, initial-scale=1.0, user-scalable=yes';
      document.head.appendChild(viewport);
    } else {
      viewport.content = 'width=device-width, initial-scale=1.0, user-scalable=yes';
    }
  };

  // Resize plotly plots
  const resizePlotlyPlots = () => {
    const plotlyDivs = document.querySelectorAll('.plotly-graph-div');
    plotlyDivs.forEach(div => {
      if (window.Plotly && div._fullLayout) {
        setTimeout(() => {
          window.Plotly.Plots.resize(div);
        }, CONFIG.PLOT_RESIZE_DELAY_MS);
      }
    });
  };

  // Adjust sidebar width based on screen size
  const adjustSidebarWidth = () => {
    const sidebar = document.querySelector('.sidebar');
    if (!sidebar) return;

    const viewportWidth = window.innerWidth;
    
    if (isMobile()) {
      sidebar.style.width = '100vw';
      sidebar.style.minWidth = '100vw';
      sidebar.style.maxWidth = '100vw';
    } else {
      const calculatedWidth = Math.max(
        CONFIG.MIN_SIDEBAR_WIDTH,
        Math.min(CONFIG.MAX_SIDEBAR_WIDTH, viewportWidth * 0.35)
      );
      sidebar.style.width = `${calculatedWidth}px`;
      sidebar.style.minWidth = `${CONFIG.MIN_SIDEBAR_WIDTH}px`;
      sidebar.style.maxWidth = `${CONFIG.MAX_SIDEBAR_WIDTH}px`;
    }
  };

  // Set consistent plot heights
  const setPlotHeights = () => {
    const plotSelectors = [
      '#cis_trans_plot_output',
      '#manhattan_plot_output'
    ];
    
    const lodScanSelectors = [
      '#lod_scan_plot_output'
    ];

    const sidebarPlotSelectors = [
      '.sidebar .plotly-graph-div'
    ];

    const viewportHeight = window.innerHeight;
    
    // Main plots (65vh with constraints)
    plotSelectors.forEach(selector => {
      const element = document.querySelector(selector);
      if (element) {
        const height = Math.max(400, Math.min(800, viewportHeight * 0.65));
        element.style.height = `${height}px`;
      }
    });

    // LOD scan plots (70vh with constraints)
    lodScanSelectors.forEach(selector => {
      const element = document.querySelector(selector);
      if (element) {
        const height = Math.max(450, Math.min(900, viewportHeight * 0.70));
        element.style.height = `${height}px`;
      }
    });

    // Sidebar plots (50vh with constraints)
    document.querySelectorAll('.sidebar .plotly-graph-div').forEach(element => {
      const height = Math.max(300, Math.min(500, viewportHeight * 0.50));
      element.style.height = `${height}px`;
    });

    // Mobile adjustments
    if (isMobile()) {
      [...plotSelectors, ...lodScanSelectors].forEach(selector => {
        const element = document.querySelector(selector);
        if (element) {
          const height = Math.max(300, viewportHeight * 0.50);
          element.style.height = `${height}px`;
        }
      });
    }
  };

  // Handle orientation change on mobile devices
  const handleOrientationChange = () => {
    setTimeout(() => {
      adjustSidebarWidth();
      setPlotHeights();
      resizePlotlyPlots();
    }, 500); // Delay to allow orientation change to complete
  };

  // Main resize handler
  const handleResize = debounce(() => {
    adjustSidebarWidth();
    setPlotHeights();
    resizePlotlyPlots();
  }, CONFIG.RESIZE_DEBOUNCE_MS);

  // Initialize autosizing
  const initializeAutosizing = () => {
    ensureViewportMeta();
    adjustSidebarWidth();
    setPlotHeights();
    
    // Set up event listeners
    window.addEventListener('resize', handleResize);
    window.addEventListener('orientationchange', handleOrientationChange);
    
    // Handle Shiny connection
    $(document).on('shiny:connected', () => {
      setTimeout(() => {
        adjustSidebarWidth();
        setPlotHeights();
        resizePlotlyPlots();
      }, 1000);
    });

    // Handle Shiny output updates
    $(document).on('shiny:value', (event) => {
      if (event.target.classList.contains('plotly-graph-div')) {
        setTimeout(() => {
          resizePlotlyPlots();
        }, CONFIG.PLOT_RESIZE_DELAY_MS);
      }
    });

    // Handle tab changes
    $(document).on('shown.bs.tab', () => {
      setTimeout(() => {
        resizePlotlyPlots();
      }, CONFIG.PLOT_RESIZE_DELAY_MS);
    });

    // Handle sidebar toggle if present
    $(document).on('click', '[data-bs-toggle="sidebar"]', () => {
      setTimeout(() => {
        adjustSidebarWidth();
        resizePlotlyPlots();
      }, 300);
    });

    // Periodic check for new plotly plots
    setInterval(() => {
      const plotlyDivs = document.querySelectorAll('.plotly-graph-div');
      plotlyDivs.forEach(div => {
        if (window.Plotly && div._fullLayout && !div._autoSizeInitialized) {
          div._autoSizeInitialized = true;
          window.Plotly.Plots.resize(div);
        }
      });
    }, 2000);
  };

  // Enhanced plotly configuration
  const enhancePlotlyConfig = () => {
    if (window.Plotly) {
      // Set default config for all plotly plots
      window.Plotly.setPlotConfig({
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: [
          'sendDataToCloud',
          'editInChartStudio'
        ]
      });
    }
  };

  // Initialize when DOM is ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initializeAutosizing);
  } else {
    initializeAutosizing();
  }

  // Initialize plotly enhancements when plotly is available
  if (window.Plotly) {
    enhancePlotlyConfig();
  } else {
    // Wait for plotly to load
    const checkPlotly = setInterval(() => {
      if (window.Plotly) {
        enhancePlotlyConfig();
        clearInterval(checkPlotly);
      }
    }, 100);
  }

  // Expose utilities globally for debugging
  window.qtlAppAutosize = {
    resizePlots: resizePlotlyPlots,
    adjustSidebar: adjustSidebarWidth,
    setHeights: setPlotHeights,
    isMobile: isMobile,
    config: CONFIG
  };

})(); 