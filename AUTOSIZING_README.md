# QTL App Autosizing Implementation

This document describes the comprehensive autosizing solution implemented for the QTL App to ensure consistent layout and user experience across all devices and screen sizes.

## Overview

The autosizing system ensures that:

- All users see the same proportional layout regardless of screen size
- Plots maintain optimal visibility and interaction
- The interface adapts gracefully to mobile, tablet, and desktop screens
- Performance remains optimal across all device types

## Components

### 1. CSS Framework (`www/autosize.css`)

**Key Features:**

- Responsive viewport units (vh, vw, clamp)
- Flexible sidebar sizing (350px - 600px range)
- Consistent plot height constraints
- Mobile-first responsive breakpoints
- High DPI display optimizations

**Plot Sizing Strategy:**

- Main plots: 65vh (400px min, 800px max)
- LOD scan plots: 70vh (450px min, 900px max)
- Sidebar plots: 50vh (300px min, 500px max)
- Mobile plots: 50vh (300px min)

### 2. JavaScript Controller (`www/autosize.js`)

**Functionality:**

- Dynamic viewport meta tag enforcement
- Debounced resize handling (250ms)
- Plotly plot auto-resizing
- Orientation change handling
- Shiny integration hooks

**Key Functions:**

- `adjustSidebarWidth()`: Responsive sidebar sizing
- `setPlotHeights()`: Consistent plot height management
- `resizePlotlyPlots()`: Force plotly plot resizing
- `handleOrientationChange()`: Mobile orientation support

### 3. R Helper Functions (`R/plot_enhancements.R`)

**New Functions:**

- `make_plotly_responsive()`: Applies consistent responsive config to plotly plots
- `responsive_plotly_output()`: Creates responsive plotly UI elements

### 4. App Integration

**Modified Files:**

- `R/scanApp_monolithic_backup.R`: Added autosizing resources to UI head
- `R/cisTransPlotApp.R`: Updated to use responsive sizing
- Plot containers now use CSS classes and IDs for targeting

## Device Breakpoints

### Mobile (≤768px)

- Sidebar: 100vw (full width)
- Plots: 50vh height
- Stacked layout (sidebar above main content)

### Tablet (769px - 1024px)

- Sidebar: 30vw (300px - 500px)
- Standard plot heights
- Side-by-side layout

### Desktop (≥1025px)

- Sidebar: 35vw (350px - 600px)
- Full plot heights
- Optimized for large screens

### Large Screens (≥1400px)

- Container max-width: 1400px (centered)
- Sidebar max-width: 500px
- Prevents excessive stretching

## Implementation Details

### CSS Classes Added

- `.overview-plot-container`: Main overview plots
- `.sidebar-plot-container`: Sidebar plots
- `.qtl-app-container`: Main app container

### JavaScript Events Handled

- Window resize (debounced)
- Orientation change
- Shiny connection/updates
- Tab changes
- Sidebar toggles

### Plotly Configuration

- `responsive: true`
- `autosize: true`
- Consistent margins and config
- Hardware acceleration via WebGL

## Browser Support

**Fully Supported:**

- Chrome 60+
- Firefox 55+
- Safari 12+
- Edge 79+

**Graceful Degradation:**

- Internet Explorer 11 (basic responsive features)
- Older mobile browsers (fallback sizing)

## Performance Optimizations

1. **Debounced Resize**: Prevents excessive recalculations
2. **CSS Hardware Acceleration**: Uses transform3d and will-change
3. **Efficient Selectors**: Minimal DOM queries
4. **Lazy Plotly Resize**: Only resizes when necessary
5. **Memory Management**: Cleanup of event listeners

## Testing Recommendations

### Device Testing

- iPhone (various sizes)
- iPad (portrait/landscape)
- Android tablets
- Desktop monitors (1080p, 1440p, 4K)
- Ultrawide monitors

### Browser Testing

- Chrome DevTools device simulation
- Firefox responsive design mode
- Safari Web Inspector
- Real device testing

### Functionality Testing

- Plot interaction (zoom, pan, hover)
- Sidebar responsiveness
- Orientation changes
- Window resizing
- Tab switching

## Debugging

The autosizing system exposes debugging utilities:

```javascript
// In browser console
window.qtlAppAutosize.resizePlots(); // Force plot resize
window.qtlAppAutosize.adjustSidebar(); // Adjust sidebar width
window.qtlAppAutosize.setHeights(); // Set plot heights
window.qtlAppAutosize.isMobile(); // Check if mobile
window.qtlAppAutosize.config; // View configuration
```

## Customization

### Adding New Plot Types

1. Add CSS selector to `autosize.css`
2. Add element targeting in `autosize.js`
3. Use `responsive_plotly_output()` for UI
4. Apply `make_plotly_responsive()` to plotly objects

### Modifying Breakpoints

Update `CONFIG.MOBILE_BREAKPOINT` in `autosize.js` and corresponding CSS media queries.

### Adjusting Plot Heights

Modify the height calculations in `setPlotHeights()` function and corresponding CSS.

## Maintenance

### Regular Tasks

- Test on new browser versions
- Update breakpoints for new device sizes
- Monitor performance on low-end devices
- Update CSS for new Shiny/bslib versions

### Known Issues

- Some older browsers may not support CSS clamp()
- iOS Safari viewport height quirks (handled with JS)
- Plotly resize timing on slow devices (debounced)

## Future Enhancements

1. **Container Queries**: When browser support improves
2. **Intersection Observer**: For lazy plot loading
3. **Web Workers**: For heavy resize calculations
4. **CSS Grid**: Enhanced layout flexibility
5. **Progressive Enhancement**: Better fallbacks

---

For questions or issues with the autosizing system, refer to this documentation or check the browser console for debugging information.
