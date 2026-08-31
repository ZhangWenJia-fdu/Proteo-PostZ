# Reusable Shiny UI helpers and analysis-card layout
size_inputs <- function(prefix, width_pt = default_plot_width_pt, height_pt = default_plot_height_pt, show_note = TRUE) {
  controls <- tagList(
    numericInput(paste0(prefix, "_w_pt"), "PDF width (pt)", value = width_pt, min = 100),
    numericInput(paste0(prefix, "_h_pt"), "PDF height (pt)", value = height_pt, min = 100)
  )
  if (isTRUE(show_note)) controls <- tagList(controls, div(class = "small-note", "Vector PDF uses points: 250 pt = 3.47 in. Default text: axis title 12 pt, tick labels 8 pt, legend/title 8 pt."))
  controls
}

heatmap_size_inputs <- function(prefix) {
  if (prefix == "cor") return(size_inputs(prefix, default_cor_heatmap_width_pt, default_cor_heatmap_height_pt))
  if (prefix %in% c("exprhm", "feature_hm")) return(size_inputs(prefix, default_expr_heatmap_width_pt, default_expr_heatmap_height_pt))
  size_inputs(prefix)
}

export_input <- function(id) checkboxInput(paste0(id, "_export_csv"), "Export corresponding plot data CSV", TRUE)
palette_input <- function(id) selectInput(paste0(id, "_palette"), "Palette", choices = palette_choices)

select_ml_feature_ids <- function(source, ml_results, n) {
  source <- match.arg(source, c("rf", "l1", "rfl1", "union"))
  final_ids <- function(x) {
    x <- unique(trimws(as.character(x %||% character())))
    x[nzchar(x) & !is.na(x)]
  }
  available <- ml_results[vapply(ml_results, function(x) length(final_ids(x)) > 0, logical(1))]
  if (source != "union") {
    ids <- final_ids(ml_results[[source]])
    if (length(ids) == 0) stop(source, " selected proteins have not been generated in the current session. Please run that machine-learning analysis first.")
    return(head(ids, n))
  }
  if (length(available) == 0) stop("No ML selected proteins are available in the current session. Please run RF, L1, or RF + L1 first.")
  if (length(available) == 1) return(head(final_ids(available[[1]]), n))
  head(unique(unlist(lapply(available, final_ids), use.names = FALSE)), n)
}

heatmap_display_controls <- function(prefix) {
  tagList(
    selectInput(paste0(prefix, "_row_cluster"), "Row clustering", choices = c("Hierarchical" = "hclust", "K-means" = "kmeans", "None" = "none")),
    conditionalPanel(sprintf("input.%s_row_cluster == 'hclust'", prefix), numericInput(paste0(prefix, "_row_hclust_k"), "Row hierarchical clustering k", 1, min = 1, step = 1), div(class = "small-note", "Cuts the row dendrogram into k clusters.")),
    conditionalPanel(sprintf("input.%s_row_cluster == 'kmeans'", prefix), numericInput(paste0(prefix, "_row_k"), "Row K-means k", 1, min = 1, step = 1)),
    selectInput(paste0(prefix, "_col_cluster"), "Column clustering", choices = c("Hierarchical" = "hclust", "K-means" = "kmeans", "None" = "none")),
    conditionalPanel(sprintf("input.%s_col_cluster == 'hclust'", prefix), numericInput(paste0(prefix, "_col_hclust_k"), "Column hierarchical clustering k", 1, min = 1, step = 1), div(class = "small-note", "Cuts the column dendrogram into k clusters; cluster annotation is shown above the heatmap.")),
    conditionalPanel(sprintf("input.%s_col_cluster == 'kmeans'", prefix), numericInput(paste0(prefix, "_col_k"), "Column K-means k", 1, min = 1, step = 1)),
    selectInput(paste0(prefix, "_row_labels"), "Row labels", choices = c("Accession" = "accession", "Protein name" = "protein_name", "Gene name" = "gene_name", "None" = "none"), selected = "none"),
    selectInput(paste0(prefix, "_color_scheme"), "Heatmap color scheme", choices = heatmap_color_scheme_choices, selected = "blue_white_red"),
    selectInput(paste0(prefix, "_group_palette"), "Group annotation colors", choices = annotation_palette_choices, selected = "npg"),
    conditionalPanel(sprintf("input.%s_group_palette == 'custom'", prefix), uiOutput(paste0(prefix, "_group_custom_colors"))),
    selectInput(paste0(prefix, "_cluster_palette"), "Cluster annotation colors", choices = annotation_palette_choices, selected = "lancet"),
    conditionalPanel(sprintf("input.%s_col_cluster != 'none' && input.%s_cluster_palette == 'custom'", prefix, prefix), uiOutput(paste0(prefix, "_cluster_custom_colors")))
  )
}

clean_group_suffix <- function(x) {
  x <- trimws(as.character(x %||% "1"))
  x <- sub("^Group", "", x, ignore.case = TRUE)
  x[x == ""] <- "1"
  x
}

venn_note <- uiOutput("venn_note")
upset_note <- uiOutput("upset_note")

ml_common_inputs <- function(prefix) {
  tagList(
    numericInput(paste0(prefix, "_seed"), "Random seed", 123, min = 1, step = 1),
    selectInput(paste0(prefix, "_split_mode"), "Train/test split mode", choices = split_mode_choices, selected = "auto"),
    numericInput(paste0(prefix, "_train_prop"), "Training set proportion", 0.7, min = 0.1, max = 0.95, step = 0.05),
    checkboxInput(paste0(prefix, "_small_sample"), "Allow small-sample exploratory ML", FALSE),
    div(class = "small-note", "Default strict mode requires at least 6 samples per group for ML; train/test split requires at least 8 per group. Small-sample mode is exploratory only.")
  )
}

stability_inputs <- function(prefix, mode = c("rf", "l1", "rfl1")) {
  mode <- match.arg(mode)
  controls <- list(
    numericInput(paste0(prefix, "_stability_repeats"), "Stability selection repeats", 50, min = 1, max = 1000, step = 1),
    numericInput(paste0(prefix, "_stability_top_var_n"), "Stability prefilter: top variance features", 200, min = 10, step = 10)
  )
  if (mode == "rf") controls <- c(controls, list(
    numericInput(paste0(prefix, "_stability_top20_weight"), "Stability weight: RF top 20 frequency", 0.35, min = 0, max = 1, step = 0.05),
    numericInput(paste0(prefix, "_stability_top50_weight"), "Stability weight: RF top 50 frequency", 0.25, min = 0, max = 1, step = 0.05),
    numericInput(paste0(prefix, "_stability_gini_weight"), "Stability weight: mean Gini", 0.10, min = 0, max = 1, step = 0.05)
  ))
  if (mode == "rfl1") controls <- c(controls, list(
    numericInput(paste0(prefix, "_stability_rf_top20_weight"), "Stability weight: RF top 20 frequency", 0.35, min = 0, max = 1, step = 0.05),
    numericInput(paste0(prefix, "_stability_rf_top50_weight"), "Stability weight: RF top 50 frequency", 0.25, min = 0, max = 1, step = 0.05),
    numericInput(paste0(prefix, "_stability_l1_weight"), "Stability weight: L1 selection frequency", 0.30, min = 0, max = 1, step = 0.05),
    numericInput(paste0(prefix, "_stability_gini_weight"), "Stability weight: mean Gini", 0.10, min = 0, max = 1, step = 0.05)
  ))
  note <- switch(mode,
    rf = "Each repeat uses the selected training proportion for a stratified split. RF features are ranked by selection frequency and mean Gini stability.",
    l1 = "Each repeat uses the selected training proportion for a stratified split. L1 features are ranked by how often they are selected.",
    rfl1 = "Each repeat uses the selected training proportion for a stratified split. The four weights determine how RF and L1 stability measures contribute to the combined feature score; their total should be 1."
  )
  tagList(controls, div(class = "small-note", note))
}

make_short_sample_names <- function(group_suffix) {
  suffix <- clean_group_suffix(group_suffix)
  counters <- list()
  vapply(suffix, function(s) {
    counters[[s]] <<- (counters[[s]] %||% 0) + 1
    paste0(s, "-", counters[[s]])
  }, character(1))
}

analysis_card <- function(id, title, controls) {
  card(class = "analysis-card", full_screen = TRUE,
    card_header(title),
    layout_columns(col_widths = c(4, 8),
      div(class = "analysis-controls", controls, div(class = "control-actions", actionButton(paste0("reset_", id), "Restore defaults", class = "btn-outline-secondary btn-sm"))),
      div(class = "analysis-results",
        div(class = "preview-wrap", imageOutput(paste0(id, "_preview"), height = "auto")),
        uiOutput(paste0(id, "_preview_controls")),
        br(),
        verbatimTextOutput(paste0(id, "_status"))
      )
    )
  )
}
quantitative_qc_card <- function(controls) {
  card(class = "analysis-card", full_screen = TRUE,
    card_header("Quantitative QC"),
    layout_columns(col_widths = c(4, 8),
      div(class = "analysis-controls", controls,
        div(class = "control-actions", actionButton("reset_qc", "Restore defaults", class = "btn-outline-secondary btn-sm"))),
      div(class = "analysis-results",
        navset_tab(
          nav_panel("Sample missingness", div(class = "preview-wrap", imageOutput("qc_sample_missingness_preview", height = "auto"))),
          nav_panel("Protein missingness", div(class = "preview-wrap", imageOutput("qc_protein_missingness_preview", height = "auto"))),
          nav_panel("Missingness heatmap", div(class = "preview-wrap", imageOutput("qc_missingness_heatmap_preview", height = "auto"))),
          nav_panel("Missingness vs abundance", div(class = "preview-wrap", imageOutput("qc_missingness_abundance_preview", height = "auto"))),
          nav_panel("Sample abundance", div(class = "preview-wrap", imageOutput("qc_sample_abundance_preview", height = "auto")))
        ),
        br(),
        verbatimTextOutput("qc_status")
      )
    )
  )
}
quant_page <- function(...) div(class = "quant-page-viewport", layout_columns(col_widths = c(6, 6), ...))
