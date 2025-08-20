## qtlApp

### Summary

- **What it is**: Modular R/Shiny package for QTL (Quantitative Trait Loci) analysis and visualization.
- **What it does**: Fast LOD scans, Manhattan and cis/trans plots, peak exploration with allele effects, interactive sex/diet analyses with on‑the‑fly difference plots, trait search, correlation and profile plots.
- **How it’s built**: Small, focused Shiny modules; ggplot2/plotly visualizations; efficient data access (fst/data.table); clear reactive patterns with debouncing; consistent styling via `bslib`.

### Key Changes (Summer)

- **Plotting & UX**
  - Responsive layout so pages fit across screen sizes; horizontal zoom and chromosome-aware zoom for scans and Manhattan/peaks.
  - Smaller, carded plots; consistent titles; colored threshold lines; thinner overlays; white divider at thresholds.
  - Allele effects plot switched to horizontal layout; peaks plot squared; overlay thickness tuned (additive thinner).
  - Additive vs Interactive difference profile displayed and kept in sync with zoom state; LOD line colors harmonized with threshold/overlay; distinct color for subtraction.
  - Axis/hover polish: added x‑axis for Manhattan/peaks; hover uses position (not `BPcum`); replaced `N/A` symbols with gene IDs; normalized specific x‑axis labels (e.g., `'_'` over `'x'` where applicable).
  - Removed clicked-point side info in favor of cleaner plot interactions; color-by control simplified to single-choice and correctly updates point colors.
- **Analysis Features**
  - Sex/Diet interactive scans with on‑the‑fly subtraction from additive.
  - Manhattan plot difference shown only for interactive analyses.
  - Cis/trans plotting added; sex×diet Manhattan and cis plots added.
  - Split-by LOD overlays (e.g., Female vs Male; HC vs HF) using additive scans.
  - Dual LOD thresholds exposed and tuned (additive 7.5, interactive 10.5), with colored lines and overlays.
  - Peaks tied to strain effects (default to highest peak; removed dropdown friction).
  - Restored legacy peaks behavior and enriched hover with peak metadata.
  - Added 4 Mb context when locating peaks; improved manipulability of peaks plot.
  - Permutation-based p-value support (e.g., 0.1) integrated into thresholding displays.
- **Data & Content**
  - Loaded metabolites dataset (small) and ensured gene symbol fixes.
  - Preserved selected gene/trait when switching datasets; ensured count of mice scanned is displayed (e.g., N/1157) and phenotype-aware.
  - Added Profile Plot and Correlation sections.
- **Performance, State & Quality**
  - Debounced expensive reactives; guarded against double-firing; improved additive/interactive state sync.
  - General directory cleanup and consistency work across modules and helpers.

### In Progress / TODO

- **Export & Sharing**: Add PNG/PDF buttons; Download as CSV.
- **SNP/Variants**: Mouse SNP Wizard backend (clone, strip query logic), integrate `founder_variants_all_csqs.sqlite`; ingest `annotated_peak_summaries/qtl_by_covar` (split by sex/diet); SNP association.
- **Peaks & Windows**: Consolidate 4 Mb search window logic in code paths.
- **Interactions & Clicks**: Fix click behavior on additive scan.
- **Data Sources & Annotations**: Biomart integration.
- **Usability Defaults**: Make split-by scans/peaks the default; show informative blurb when none found.
- **Performance & Stability**: Address memory failure conditions.
- **Correlation**: Add number of mice to correlation (requires re-run).
- **Trait/ID Gaps**: Fix cases where `ENSMUSG` traits don’t show peaks.
- **Domain additions**: Consider liver splice junctions.

### System Structure

- **Entry Points**
  - `app.R`: Main Shiny app wiring, sourcing order, global options.
  - `mainUI.R`: High-level UI composition and layout components.
- **Core Modules (Server/UI)**
  - `scanPlotModule.R`: LOD scans (additive, interactive, overlays, difference logic).
  - `splitByLodOverlayModule.R`: Split-by additive overlays (Female vs Male; HC vs HF), follows main zoom.
  - `alleleEffectsModule.R`: Founder allele effects visualization tied to selected peak.
  - `peaksTableModule.R`: Peaks table with hover/click metadata; ties into scans and effects.
  - `datasetSelectionModule.R`, `interactiveAnalysisModule.R`: Dataset and interaction type controls, state sync.
  - `traitSearchModule.R`, `traitProcessingModule.R`: Large-scale trait search and preprocessing.
  - `profilePlotApp.R`, `correlationApp.R`: Profile and correlation analysis sections.
  - `manhattanPlotApp.R`, `cisTransPlotApp.R`: Genome-wide Manhattan and single-panel cis/trans.
  - `downloadApp.R`: Download-related UI/actions (extend for CSV/PNG/PDF).
- **Computation & Data Access**
  - `trait_scan.R`: Per-trait scan retrieval; bridges to visualization.
  - `peak_finder.R`, `peak_info.R`: Peak detection and metadata extraction.
  - `fst_rows.R`, `data_handling.R`, `import_data.R`: Efficient I/O and dataset wiring.
  - `QTL_plot_visualizer.R`: Normalizes scan data for plotting (chr/position/BPcum joins with markers).
- **Visualization Helpers**
  - `ggplot_qtl_scan.R`, `ggplotly_qtl_scan.R`: Scan plotting (ggplot2/plotly variants).
  - `ggplot_alleles.R`: Allele effects plotting.
  - `plot_enhancements.R`, `ui_styles.R`, `plot_null.R`: Themes, styles, and fallbacks.
- **Apps/Utilities**
  - `scanApp.R`, `cisTransPlotApp.R`, `manhattanPlotApp.R`, `interactiveSubtractionApp.R`: Focused app shells for specific views.
  - `helpers.R`: Shared utilities (validation, formatting, conversions, reactive helpers).
  - `preprocess_pheno_data.R`, `update_annotation_list_from_csvs.R`: Data prep and annotation tooling.
  - `scripts/`, `kalynn_R/`: Processing pipelines and project-specific scripts.
- **Package Files**
  - `DESCRIPTION`, `NAMESPACE`, `man/`: Package metadata and documentation.
  - `inst/`, `data/`, `docs/`: Installed assets, datasets, and docs.

### Notes

- Interaction analysis supports additive and sex/diet modes with automatic on‑the‑fly subtraction for difference plots; split-by overlays render additive-only comparisons for interpretability.
- Default LOD thresholds: additive 7.5; interactive 10.5. Colored threshold lines are used consistently across plots and overlays.
- The system uses debounced reactives and stable state patterns to prevent circular dependencies and double-firing during UI updates.
