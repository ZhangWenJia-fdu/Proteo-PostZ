# ProteoPostZ Formal V2.0.2
# Developed by Wenjia Zhang

options(shiny.maxRequestSize = 1024^3)
library(shiny)
library(bslib)
library(DT)
library(dplyr)

source(file.path("R", "analysis_core.R"), encoding = "UTF-8")
`%||%` <- function(a, b) if (!is.null(a)) a else b

app_version <- "Formal V2.0.2"
app_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
package_root <- normalizePath(file.path(app_root, ".."), winslash = "/", mustWork = FALSE)
annotation_dir <- file.path(app_root, "annotations")
default_output <- file.path(package_root, "outputs")
default_annotation_file <- file.path(annotation_dir, "uniprot_reviewed_human_9606_annotations.csv")
default_plot_width_pt <- 250
default_plot_height_pt <- 250
default_cor_heatmap_width_pt <- 300
default_cor_heatmap_height_pt <- 200
default_expr_heatmap_width_pt <- 600
default_expr_heatmap_height_pt <- 800

palette_choices <- c(
  "NPG / Nature Publishing: red-cyan-green-navy" = "npg",
  "Lancet: blue-red-green-purple" = "lancet",
  "JAMA: muted teal-orange-blue" = "jama",
  "NEJM: red-blue-orange-green" = "nejm",
  "UChicago: maroon-grey-gold" = "uchicago"
)
annotation_palette_choices <- c(palette_choices, "Custom #RRGGBB colors" = "custom")
heatmap_color_scheme_choices <- c(
  "Blue-White-Red" = "blue_white_red",
  "Red-White-Blue" = "red_white_blue",
  "Purple-White-Orange" = "purple_white_orange",
  "Green-Black-Red" = "green_black_red"
)

split_mode_choices <- c("Auto" = "auto", "Cross-validation only" = "cross_validation_only", "Train/test split" = "train_test_split")
lambda_selection_choices <- c("lambda.1se" = "lambda.1se", "lambda.min" = "lambda.min")
feature_source_choices <- c(
  "Random forest selected proteins" = "rf",
  "L1 selected proteins" = "l1",
  "RF + L1 selected proteins" = "rfl1",
  "All available ML selected proteins / union" = "union"
)
input_family_choices <- c("DIA result" = "dia", "DDA result" = "dda", "Standard quantitative matrix" = "standard_matrix")
dia_format_choices <- c("Auto-detect DIA format" = "auto", "DIA-NN" = "diann", "Spectronaut" = "spectronaut")
dda_format_choices <- c("Auto-detect DDA format" = "auto", "FragPipe / MSFragger" = "fragpipe", "PEAKS" = "peaks", "MaxQuant" = "maxquant")
standard_zero_mode_choices <- c(
  "0 is a valid quantitative value" = "zero_is_value",
  "0 represents missing or unquantified" = "zero_is_missing"
)

app_css <- "
body, .form-control, .selectize-input, .btn, table { font-family: 'Microsoft YaHei', Arial, sans-serif; }
.en, .brand { font-family: Arial, sans-serif; }
.sidebar { overflow: hidden !important; padding-bottom: 0; }
.sidebar .sidebar-content { height: 100%; display: flex; flex-direction: column; min-height: 0; }
.sidebar-scroll { flex: 1 1 auto; min-height: 0; overflow-y: auto; padding-right: 4px; }
.sidebar-fixed { flex: 0 0 auto; padding-top: 12px; border-top: 1px solid #d2d6da; background: transparent; }
.sidebar-fixed .form-group { margin-bottom: 10px; }
.input-control-grid { display: grid; grid-template-columns: minmax(220px, 1fr) minmax(260px, 1fr); gap: 10px 14px; align-items: start; }
.input-control-span-2 { grid-column: 1 / -1; }
.input-control-grid > div { min-width: 0; }
.input-control-grid .form-group { margin-bottom: 0; }
.input-control-grid .form-control, .input-control-grid .selectize-control { max-width: 100%; }
.brand { padding-top: 10px; margin-top: 12px; border-top: 1px solid #d2d6da; background: transparent; font-style: italic; color: #555; font-size: 13px; display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1fr); gap: 0 18px; }
.brand > div { min-width: 0; }
.brand > div:first-child { grid-column: 1 / -1; }
.brand .developer-name { font-family: 'Segoe Script', 'Brush Script MT', 'Lucida Handwriting', cursive; font-size: 16px; color: #3b4b54; }
.brand .institution-cn, .brand .lab-line { font-family: 'Microsoft YaHei', Arial, sans-serif; font-style: normal; color: #263942; font-size: 13px; font-weight: 650; margin-top: 6px; }
.brand .institution-en, .brand .lab-line-en { font-family: Arial, sans-serif; font-style: normal; color: #40515c; font-size: 12px; margin-top: 2px; }
.brand .institution-cn { grid-column: 1; grid-row: 2; }
.brand .institution-en { grid-column: 1; grid-row: 3; }
.brand .lab-line { grid-column: 2; grid-row: 2; }
.brand .lab-line-en { grid-column: 2; grid-row: 3; }
.brand-grid { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1fr); gap: 8px 18px; margin-top: 8px; }
.brand-grid > div { min-width: 0; }
.card { border-radius: 6px; box-shadow: 0 1px 8px rgba(20, 35, 50, 0.06); }
.small-note { color: #666; font-size: 12px; line-height: 1.35; }
.detect-note { color: #2f5f73; font-size: 12px; line-height: 1.35; margin: 4px 0 8px 0; white-space: pre-wrap; }
.sample-warning { color: #B2182B; font-size: 13px; font-weight: 600; margin: 6px 0; }
.input-error { color: #B2182B; font-size: 13px; font-weight: 600; white-space: pre-wrap; margin-top: 8px; }
.zero-mode-detail { margin: 4px 0 8px 22px; color: #555; font-size: 12px; line-height: 1.35; }
.analysis-card .card-header { font-weight: 650; }
.quant-page-viewport { height: calc(100vh - 190px); min-height: 0; min-width: 0; overflow: hidden; }
.quant-page-viewport > .bslib-grid { height: 100%; min-height: 0; grid-auto-rows: minmax(0, 1fr); align-content: stretch; overflow: hidden; }
.quant-page-viewport .analysis-card { height: 100%; min-height: 0; overflow: hidden; }
.quant-page-viewport .analysis-card > .card-body { min-height: 0; height: calc(100% - 44px); overflow: hidden; }
.quant-page-viewport .analysis-card > .card-body > .bslib-grid { height: 100%; min-height: 0; overflow: hidden; }
.quant-page-viewport .analysis-card .analysis-controls, .quant-page-viewport .analysis-card .analysis-results { min-height: 0; height: 100%; overflow-y: auto; overflow-x: hidden; padding-right: 8px; }
.quant-page-viewport .analysis-card .preview-wrap { min-height: 180px; }
.input-main-stack { height: calc(100vh - 150px); min-height: 620px; min-width: 380px; display: flex; flex-direction: column; gap: 12px; overflow: hidden; }
.input-half-card { flex: 1 1 0; min-height: 0; }
.input-half-card > .card-body { overflow: auto; }
.control-actions { margin-top: 12px; display: flex; gap: 8px; flex-wrap: wrap; }
.preview-wrap { min-height: 280px; display: flex; align-items: center; justify-content: center; background: #fafafa; border: 1px solid #eee; border-radius: 6px; padding: 8px; }
.preview-wrap img { max-width: 100%; height: auto; }
.preview-nav { margin-top: 8px; display: flex; align-items: center; justify-content: center; gap: 10px; }
@media (max-width: 980px) {
  .input-control-grid { grid-template-columns: 1fr; }
  .input-control-span-2 { grid-column: 1; }
  .input-main-stack { min-width: 320px; }
  .quant-page-viewport { height: calc(100vh - 170px); min-height: 0; }
}
@media (max-width: 1250px), (max-height: 760px) {
  .input-control-grid { grid-template-columns: 1fr; }
  .input-control-span-2 { grid-column: 1; }
  .brand { grid-template-columns: 1fr; gap: 0; }
  .brand > div:first-child, .brand .institution-cn, .brand .institution-en, .brand .lab-line, .brand .lab-line-en { grid-column: 1; grid-row: auto; }
  .brand-grid { grid-template-columns: 1fr; gap: 4px; }
}
"

size_inputs <- function(prefix, width_pt = default_plot_width_pt, height_pt = default_plot_height_pt) {
  tagList(
    numericInput(paste0(prefix, "_w_pt"), "PDF width (pt)", value = width_pt, min = 100),
    numericInput(paste0(prefix, "_h_pt"), "PDF height (pt)", value = height_pt, min = 100),
    div(class = "small-note", "Vector PDF uses points: 250 pt = 3.47 in. Default text: axis title 12 pt, tick labels 8 pt, legend/title 8 pt.")
  )
}

heatmap_size_inputs <- function(prefix) {
  if (prefix == "cor") return(size_inputs(prefix, default_cor_heatmap_width_pt, default_cor_heatmap_height_pt))
  if (prefix %in% c("exprhm", "feature_hm")) return(size_inputs(prefix, default_expr_heatmap_width_pt, default_expr_heatmap_height_pt))
  size_inputs(prefix)
}

export_input <- function(id) checkboxInput(paste0(id, "_export_csv"), "Export corresponding plot data CSV", TRUE)
palette_input <- function(id) selectInput(paste0(id, "_palette"), "Palette", choices = palette_choices)

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

quant_page <- function(...) div(class = "quant-page-viewport", layout_columns(col_widths = c(6, 6), ...))

ui <- page_navbar(
  title = paste("ProteoPostZ", app_version),
  theme = bs_theme(version = 5, primary = "#155F83", bootswatch = "flatly"),
  header = tags$head(tags$style(HTML(app_css))),
  nav_panel("Input",
    layout_sidebar(
      sidebar = sidebar(width = "50%",
        div(class = "sidebar-scroll",
          div(class = "input-control-grid",
            div(class = "input-control-span-2", textInput("file_path", "Result file path (.csv/.tsv)", value = "")),
            radioButtons("input_family", "Input category", choices = input_family_choices, selected = "dia", inline = FALSE),
            conditionalPanel("input.input_family == 'dia'",
              selectInput("dia_format", "Software format", choices = dia_format_choices, selected = "auto")
            ),
            conditionalPanel("input.input_family == 'dda'",
              selectInput("dda_format", "Software format", choices = dda_format_choices, selected = "auto")
            ),
            div(class = "input-control-span-2",
              conditionalPanel("input.input_family != 'standard_matrix'",
                selectInput("row_id", "Feature identifier used as matrix row name", choices = c("Accessions" = "accession", "Protein name" = "protein_name", "Gene name" = "gene_name"), selected = "accession"),
                div(class = "small-note", "Accessions are recommended for stable annotation mapping. Gene names may be missing for some accessions, especially when the search FASTA is not one-gene-one-protein."),
                uiOutput("detected_format_message")
              )
            ),
            div(class = "input-control-span-2",
              conditionalPanel("input.input_family == 'standard_matrix'",
                div(class = "small-note", "Different software packages and processing workflows may use zero differently. Select the mode that matches the data source; the software will not infer this automatically."),
                div(class = "small-note", "For physicochemical property analysis, standard matrix feature identifiers are matched to the annotation Accession column as-is. Use UniProt accessions in the first column if you want built-in annotation matching."),
                selectInput("standard_zero_mode", "0 values and missing values - Please select how zero should be interpreted", choices = c("Please select" = "", standard_zero_mode_choices), selected = ""),
                div(class = "zero-mode-detail", strong("0 is a valid quantitative value"), br(), "Zero is retained as a valid quantitative value. Only explicit missing entries, such as blank cells, whitespace-only cells, NA, and NaN, are treated as missing values."),
                div(class = "zero-mode-detail", strong("0 represents missing or unquantified"), br(), "Zero represents a missing or unquantified entry. Numeric zeros, blank cells, whitespace-only cells, NA, and NaN are all treated as missing values.")
              )
            )
          )
        ),
        div(class = "sidebar-fixed",
          textInput("outdir", "Output directory", value = default_output),
          actionButton("load_data", "Load data", class = "btn-primary"),
          div(class = "input-error", textOutput("input_error")),
          div(class = "brand", div("Developed by ", span(class = "developer-name", "Wenjia Zhang")), div(class = "institution-cn", "复旦大学化学系"), div(class = "institution-en", "Department of Chemistry, Fudan University"), div(class = "lab-line", "现代色谱分离分析实验室"), div(class = "lab-line-en", "Modern Chromatography Separation and Analysis"))
        )
      ),
      div(class = "input-main-stack",
        card(class = "input-half-card", full_screen = TRUE,
          card_header(uiOutput("input_overview_title")),
          layout_columns(col_widths = c(4, 4, 4), h3(textOutput("n_proteins")), h3(textOutput("n_samples")), h3(textOutput("n_groups"))),
          DTOutput("input_summary"),
          DTOutput("counts_table")),
        card(class = "input-half-card", full_screen = TRUE,
          card_header("Sample groups"),
          div(class = "small-note", "Enter group suffixes only, for example 1, 2, A, or B. Internal group names become Group1, Group2, GroupA, or GroupB. Edit sample names to control names shown in downstream plots."),
          div(class = "control-actions", actionButton("auto_short_names", "Auto-generate short sample names", class = "btn-outline-primary btn-sm")),
          div(class = "sample-warning", textOutput("sample_name_warning")),
          uiOutput("group_inputs")
        )
      )
    )
  ),
  nav_panel("Qualitative plots",
    layout_columns(col_widths = c(6, 6),
      analysis_card("idbar", "Sample protein count barplot", tagList(uiOutput("idbar_count_control"), uiOutput("idbar_count_note"), size_inputs("idbar"), palette_input("idbar"), export_input("idbar"), actionButton("run_idbar", "Generate barplot", class = "btn-primary"))),
      analysis_card("venn", "Venn diagram", tagList(venn_note, numericInput("venn_min_reps", "Minimum replicates detected in group", 1, min = 1, step = 1), size_inputs("venn"), palette_input("venn"), export_input("venn"), actionButton("run_venn", "Generate Venn", class = "btn-primary"))),
      analysis_card("upset", "UpSet plot", tagList(upset_note, numericInput("upset_min_reps", "Minimum replicates detected in group", 1, min = 1, step = 1), size_inputs("upset"), export_input("upset"), actionButton("run_upset", "Generate UpSet", class = "btn-primary"))),
      analysis_card("phys", "Physicochemical property distributions", tagList(numericInput("phys_min_reps", "Minimum replicates detected in group", 1, min = 1, step = 1), selectInput("annotation_preset", "Built-in annotation", choices = c("Human reviewed Swiss-Prot" = "human", "Mouse reviewed Swiss-Prot" = "mouse", "C. elegans UniProtKB Swiss-Prot + TrEMBL" = "celegans", "Custom annotation path below" = "custom"), selected = "human"), textInput("annotation_file", "Custom annotation table", value = default_annotation_file), size_inputs("phys"), palette_input("phys"), export_input("phys"), actionButton("run_phys", "Generate property plots", class = "btn-primary")))
    )
  ),
  nav_panel("Quantitative plots",
    navset_tab(
      nav_panel("Page 1",
        quant_page(
      analysis_card("cor", "Sample correlation heatmap", tagList(selectInput("cor_method", "Correlation method", choices = c("Pearson" = "pearson", "Spearman" = "spearman")), selectInput("cor_order", "Sample order", choices = c("Original" = "original", "By group" = "group", "Hierarchical clustering" = "hclust"), selected = "original"), conditionalPanel("input.cor_order == 'hclust'", selectInput("cor_cluster_distance", "Clustering distance", choices = c("1 - correlation" = "one_minus_correlation", "Euclidean distance of correlation profiles" = "correlation_euclidean"), selected = "one_minus_correlation"), selectInput("cor_cluster_linkage", "Linkage", choices = c("Complete" = "complete", "Average" = "average"), selected = "complete"), numericInput("cor_cluster_k", "Clusters (k)", 1, min = 1, step = 1), div(class = "small-note", "The dendrogram is cut into k clusters; gaps are shown between clusters.")), numericInput("cor_digits", "Displayed decimals", 2, min = 0, max = 5), numericInput("cor_font", "Number font size", 8, min = 4), selectInput("cor_color", "Heatmap color", choices = c("Blue-white-red" = "blue_white_red", "Purple-white-orange" = "purple_white_orange")), numericInput("cor_min", "Legend min", -1), numericInput("cor_max", "Legend max", 1), heatmap_size_inputs("cor"), export_input("cor"), actionButton("run_cor", "Generate correlation heatmap", class = "btn-primary"))),
      analysis_card("rank", "Rank-abundance plot", tagList(size_inputs("rank"), palette_input("rank"), export_input("rank"), actionButton("run_rank", "Generate rank-abundance", class = "btn-primary"))),
      analysis_card("cv", "Within-group CV ridgeline", tagList(numericInput("cv_max", "CV x-axis max (%)", 60, min = 1), size_inputs("cv"), palette_input("cv"), export_input("cv"), actionButton("run_cv", "Generate CV ridgeline", class = "btn-primary")))
        )
      ),
      nav_panel("Page 2",
        quant_page(
      analysis_card("pca", "PCA plot", tagList(numericInput("pca_min_valid", "Minimum valid fraction", 0.5, min = 0.1, max = 1, step = 0.05), size_inputs("pca"), palette_input("pca"), export_input("pca"), actionButton("run_pca", "Generate PCA", class = "btn-primary"))),
      analysis_card("umap", "UMAP plot", tagList(numericInput("umap_min_valid", "Minimum valid fraction", 0.5, min = 0.1, max = 1, step = 0.05), numericInput("umap_neighbors", "n_neighbors", 10, min = 2), numericInput("umap_min_dist", "min_dist", 0.1, min = 0, max = 1, step = 0.05), size_inputs("umap"), palette_input("umap"), export_input("umap"), actionButton("run_umap", "Generate UMAP", class = "btn-primary"))),
      analysis_card("tsne", "t-SNE plot", tagList(numericInput("tsne_min_valid", "Minimum valid fraction", 0.5, min = 0.1, max = 1, step = 0.05), textInput("tsne_perplexity", "Perplexity", "Auto"), numericInput("tsne_max_iter", "Max iterations", 1000, min = 250, step = 250), size_inputs("tsne"), palette_input("tsne"), export_input("tsne"), actionButton("run_tsne", "Generate t-SNE", class = "btn-primary")))
        )
      ),
      nav_panel("Page 3",
        quant_page(
      analysis_card("volcano", "Volcano plot", tagList(uiOutput("volcano_groups"), selectInput("volcano_fc_method", "log2FC calculation", choices = c("Mean difference after log2(x+1)" = "log2_then_diff", "log2 of raw mean ratio" = "ratio_then_log2"), selected = "log2_then_diff"), selectInput("volcano_test_method", "Statistical test", choices = c("limma" = "limma", "t test" = "ttest"), selected = "limma"), selectInput("volcano_sig_metric", "P-value threshold type", choices = c("BH-adjusted p value" = "adj_p", "P value" = "raw_p"), selected = "adj_p"), numericInput("log2fc", "log2FC cutoff", 1, min = 0), numericInput("adj_p_cutoff", "BH-adjusted p value cutoff", 0.05, min = 0, max = 1), numericInput("raw_p_cutoff", "P value cutoff", 0.05, min = 0, max = 1), size_inputs("volcano"), export_input("volcano"), actionButton("run_volcano", "Generate volcano", class = "btn-primary"))),
      analysis_card("exprhm", "Expression heatmap", tagList(numericInput("hm_top", "Top variable proteins", 100, min = 5), heatmap_display_controls("hm"), heatmap_size_inputs("exprhm"), export_input("exprhm"), actionButton("run_exprhm", "Generate heatmap", class = "btn-primary")))
        )
      ),
      nav_panel("Page 4",
        quant_page(
      analysis_card("rf", "Random forest feature selection", tagList(ml_common_inputs("rf"), numericInput("rf_top", "Top RF selected features", 30, min = 5), numericInput("rf_ntree", "Random forest ntree", 500, min = 50), textInput("rf_mtry", "Random forest mtry", "Auto"), div(class = "small-note", "Feature ranking uses MeanDecreaseGini when available."), size_inputs("rf"), palette_input("rf"), export_input("rf"), actionButton("run_rf", "Generate RF selection", class = "btn-primary"))),
      analysis_card("l1", "L1 feature selection", tagList(ml_common_inputs("l1"), numericInput("l1_alpha", "L1 alpha (1 = LASSO; 0-1 = elastic net)", 1, min = 0, max = 1), selectInput("l1_lambda", "Lambda selection", choices = lambda_selection_choices, selected = "lambda.1se"), textInput("l1_folds", "Cross-validation folds", "Auto"), numericInput("l1_top", "Top L1 features", 50, min = 5), size_inputs("l1"), palette_input("l1"), export_input("l1"), actionButton("run_l1", "Generate L1 feature plot", class = "btn-primary"))),
      analysis_card("rfl1", "RF + L1 combined feature selection", tagList(ml_common_inputs("rfl1"), numericInput("rfl1_top", "Top combined features", 50, min = 5), numericInput("rfl1_ntree", "Random forest ntree", 500, min = 50), textInput("rfl1_mtry", "Random forest mtry", "Auto"), numericInput("rfl1_alpha", "L1 alpha (1 = LASSO; 0-1 = elastic net)", 1, min = 0, max = 1), selectInput("rfl1_lambda", "Lambda selection", choices = lambda_selection_choices, selected = "lambda.1se"), textInput("rfl1_folds", "Cross-validation folds", "Auto"), size_inputs("rfl1"), palette_input("rfl1"), export_input("rfl1"), actionButton("run_rfl1", "Generate combined features", class = "btn-primary"))),
      analysis_card("feature_umap", "Feature-protein UMAP", tagList(selectInput("feature_umap_source", "Selected feature source", choices = feature_source_choices, selected = "rf"), numericInput("feature_top_umap", "Use top N selected features", 50, min = 5), numericInput("feature_umap_neighbors", "n_neighbors", 10, min = 2), numericInput("feature_umap_min_dist", "min_dist", 0.1, min = 0, max = 1, step = 0.05), size_inputs("feature_umap"), palette_input("feature_umap"), export_input("feature_umap"), actionButton("run_feature_umap", "Generate feature UMAP", class = "btn-primary"))),
      analysis_card("feature_hm", "Feature-protein heatmap", tagList(selectInput("feature_hm_source", "Selected feature source", choices = feature_source_choices, selected = "rf"), numericInput("feature_top_hm", "Use top N selected features", 50, min = 5), heatmap_display_controls("feature_hm"), heatmap_size_inputs("feature_hm"), export_input("feature_hm"), actionButton("run_feature_hm", "Generate feature heatmap", class = "btn-primary")))
        )
      ),
      nav_panel("Page 5",
        quant_page(
      analysis_card("sling", "Slingshot pseudotime", tagList(selectInput("sling_reduction", "Reduction", choices = c("PCA", "UMAP")), div(class = "small-note", "If groups have time/order: Start = earliest/control, End = final group. For exploratory analysis: Start = control, End = None."), uiOutput("slingshot_groups"), numericInput("sling_top", "Top pseudotime proteins", 50, min = 2), size_inputs("sling"), numericInput("sling_heatmap_w_pt", "Top heatmap PDF width (pt)", default_expr_heatmap_width_pt, min = 100), numericInput("sling_heatmap_h_pt", "Top heatmap PDF height (pt)", default_expr_heatmap_height_pt, min = 100), palette_input("sling"), export_input("sling"), actionButton("run_sling", "Generate Slingshot", class = "btn-primary")))
        )
      )
    )
  ),
  nav_panel("Exported files", card(card_header("Current output directory"), textOutput("outdir_text"), DTOutput("file_table"))),
  nav_spacer(),
  nav_item(actionButton("exit_app", "Exit app / 退出程序", class = "btn-danger btn-sm"))
)

server <- function(input, output, session) {
  rv <- reactiveValues(data = NULL, groups = NULL, preview = list(), preview_files = list(), preview_index = list(), status = list(), rf_features = NULL, load_error = NULL)
  ids <- c("idbar","venn","upset","phys","cor","rank","cv","pca","umap","tsne","volcano","exprhm","rf","l1","rfl1","feature_umap","feature_hm","sling")

  observeEvent(input$exit_app, {
    showNotification("Exiting app and stopping the Shiny session...", type = "message", duration = 2)
    session$close()
    session$onFlushed(function() shiny::stopApp(), once = TRUE)
  }, ignoreInit = TRUE)

  for (id in ids) {
    local({
      .id <- id
      output[[paste0(.id, "_preview")]] <- renderImage({
        p <- rv$preview[[.id]]
        if (is.null(p) || !file.exists(p)) return(NULL)
        list(src = p, contentType = "image/png", alt = paste(.id, "preview"), width = "100%", `data-refresh` = rv$preview_token[[.id]] %||% 0)
      }, deleteFile = FALSE)
      output[[paste0(.id, "_preview_controls")]] <- renderUI({
        files <- rv$preview_files[[.id]]
        if (is.null(files) || length(files) <= 1) return(NULL)
        idx <- rv$preview_index[[.id]] %||% 1
        div(class = "preview-nav",
          actionButton(paste0(.id, "_preview_prev"), "\u524d\u4e00\u5f20", class = "btn-outline-secondary btn-sm"),
          span(class = "small-note", paste(idx, "/", length(files))),
          actionButton(paste0(.id, "_preview_next"), "\u540e\u4e00\u5f20", class = "btn-outline-secondary btn-sm")
        )
      })
      observeEvent(input[[paste0(.id, "_preview_prev")]], {
        files <- rv$preview_files[[.id]]
        if (is.null(files) || length(files) <= 1) return(NULL)
        idx <- rv$preview_index[[.id]] %||% 1
        idx <- if (idx <= 1) length(files) else idx - 1
        rv$preview_index[[.id]] <- idx
        rv$preview[[.id]] <- files[[idx]]
      }, ignoreInit = TRUE)
      observeEvent(input[[paste0(.id, "_preview_next")]], {
        files <- rv$preview_files[[.id]]
        if (is.null(files) || length(files) <= 1) return(NULL)
        idx <- rv$preview_index[[.id]] %||% 1
        idx <- if (idx >= length(files)) 1 else idx + 1
        rv$preview_index[[.id]] <- idx
        rv$preview[[.id]] <- files[[idx]]
      }, ignoreInit = TRUE)
      output[[paste0(.id, "_status")]] <- renderText({ rv$status[[.id]] %||% "No plot generated yet." })
    })
  }

  is_standard_matrix_data <- function(dat = rv$data) identical(dat$input_family %||% dat$input_source %||% dat$software %||% "", "standard_matrix")
  is_dda_data <- function(dat = rv$data) identical(dat$input_family %||% "", "dda")
  count_metric_choices <- function(dat = rv$data) {
    if (is_standard_matrix_data(dat)) {
      c("Quantified protein count" = "Quantified_Protein_Count")
    } else if (identical(dat$input_source %||% "", "peaks") && identical(dat$peaks_subtype %||% "", "lfq")) {
      c(
        "Identification count (unavailable)" = "Identified_Protein_Count",
        "Quantified protein count" = "Quantified_Protein_Count"
      )
    } else {
      c(
        "Identification count" = "Identified_Protein_Count",
        "Quantified protein count" = "Quantified_Protein_Count"
      )
    }
  }
  active_count_column <- reactive({
    req(rv$data)
    choices <- unname(count_metric_choices(rv$data))
    requested <- input$count_metric %||% if (is_standard_matrix_data() || (identical(rv$data$input_source %||% "", "peaks") && identical(rv$data$peaks_subtype %||% "", "lfq"))) "Quantified_Protein_Count" else "Identified_Protein_Count"
    if (requested %in% choices) requested else choices[[1]]
  })
  count_metric_text <- function(count_col = active_count_column()) {
    switch(
      count_col,
      Identified_Protein_Count = list(
        title = "Identification count",
        y_label = "Identified protein count",
        barplot_title = "Sample identification count"
      ),
      Quantified_Protein_Count = list(
        title = "Quantified protein count",
        y_label = "Quantified protein count",
        barplot_title = if (is_standard_matrix_data()) {
          "Available quantitative values"
        } else {
          "Sample quantified protein count"
        }
      ),
      list(
        title = count_col,
        y_label = count_col,
        barplot_title = count_col
      )
    )
  }
  clear_loaded_input <- function() {
    old_pngs <- unlist(rv$preview_files, use.names = FALSE)
    if (length(old_pngs) > 0) unlink(old_pngs[file.exists(old_pngs)], force = TRUE)
    rv$data <- NULL
    rv$groups <- NULL
    rv$rf_features <- NULL
    rv$preview <- list()
    rv$preview_files <- list()
    rv$preview_index <- list()
    rv$status <- list()
    rv$load_error <- NULL
  }
  observeEvent(input$input_family, {
    clear_loaded_input()
    if (!identical(input$input_family, "standard_matrix")) updateSelectInput(session, "standard_zero_mode", selected = "")
  }, ignoreInit = TRUE)
  observeEvent(input$dia_format, { if (!is.null(rv$data)) clear_loaded_input() }, ignoreInit = TRUE)
  observeEvent(input$dda_format, { if (!is.null(rv$data)) clear_loaded_input() }, ignoreInit = TRUE)
  observeEvent(input$file_path, { clear_loaded_input() }, ignoreInit = TRUE)
  observeEvent(input$standard_zero_mode, {
    if (is_standard_matrix_data()) clear_loaded_input()
  }, ignoreInit = TRUE)
  output$detected_format_message <- renderUI({
    path <- trimws(input$file_path %||% "")
    family <- input$input_family %||% "dia"
    if (!nzchar(path) || identical(family, "standard_matrix")) return(NULL)
    if (!file.exists(path)) return(div(class = "detect-note", "Format detection: waiting for an existing file path."))
    selected <- if (family == "dia") input$dia_format %||% "auto" else input$dda_format %||% "auto"
    detected <- tryCatch(
      detect_input_format(path, family),
      error = function(e) list(id = NA_character_, label = "Detection failed", evidence = conditionMessage(e), ambiguous = FALSE)
    )
    selected_note <- if (identical(selected, "auto")) "Loader will use the detected format." else paste0("Manual override selected: ", selected, ".")
    msg <- paste0("Detected format: ", detected$label, "\n", detected$evidence, "\n", selected_note)
    if (isTRUE(detected$ambiguous)) msg <- paste0(msg, "\nMultiple signatures matched: ", paste(detected$all_hits, collapse = ", "), ". The first registered match is shown.")
    div(class = "detect-note", msg)
  })
  standard_input_summary_df <- function(dat) {
    total_cells <- length(dat$analysis_quant_matrix)
    total_missing <- sum(dat$missingness_matrix)
    available <- dat$available_quantitative_value_count
    data.frame(
      Metric = c(
        "File type",
        "First column name",
        "Feature count",
        "Sample count",
        "Zero handling mode",
        "Original zero cell count",
        "Explicit missing value count",
        "Total missing value count after selected mode",
        "Number of available quantitative values",
        "Available quantitative value proportion",
        "Duplicate features present",
        "All-missing sample columns present"
      ),
      Value = c(
        toupper(tools::file_ext(input$file_path)),
        dat$feature_column_name,
        nrow(dat$quantity),
        ncol(dat$quantity),
        paste0(dat$missingness_mode, " - ", dat$missingness_mode_label),
        dat$raw_zero_cell_count,
        dat$explicit_missing_value_count,
        total_missing,
        available,
        if (total_cells > 0) sprintf("%.2f%%", 100 * available / total_cells) else "NA",
        "No",
        "No"
      ),
      stringsAsFactors = FALSE
    )
  }
  standard_mode_note <- function() {
    if (is_standard_matrix_data()) paste0("Zero handling mode: ", rv$data$missingness_mode, " - ", rv$data$missingness_mode_label) else NULL
  }
  input_format_summary_df <- function(dat) {
    data.frame(
      Metric = c("Input category", "Detected/selected format", "Data level", "Feature count", "Sample count", "Format evidence"),
      Value = c(
        dat$input_family %||% dat$input_source %||% dat$software %||% "unknown",
        dat$format_label %||% dat$input_source %||% dat$software %||% "unknown",
        dat$data_level %||% "protein",
        nrow(dat$quantity),
        ncol(dat$quantity),
        dat$format_evidence %||% "Manual input format."
      ),
      stringsAsFactors = FALSE
    )
  }
  outdir <- reactive({ normalizePath(input$outdir, winslash = "/", mustWork = FALSE) })
  observe({ dir.create(outdir(), recursive = TRUE, showWarnings = FALSE) })
  pts <- function(prefix) c(input[[paste0(prefix, "_w_pt")]] / 72, input[[paste0(prefix, "_h_pt")]] / 72)
  analysis_dir <- function(id) { d <- file.path(outdir(), id); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
  reset_common <- function(id, palette = TRUE, export = TRUE, width_pt = default_plot_width_pt, height_pt = default_plot_height_pt) {
    updateNumericInput(session, paste0(id, "_w_pt"), value = width_pt)
    updateNumericInput(session, paste0(id, "_h_pt"), value = height_pt)
    if (palette) updateSelectInput(session, paste0(id, "_palette"), selected = "npg")
    if (export) updateCheckboxInput(session, paste0(id, "_export_csv"), value = TRUE)
  }
  observeEvent(input$reset_idbar, { reset_common("idbar") })
  observeEvent(input$reset_venn, { updateNumericInput(session, "venn_min_reps", value = 1); reset_common("venn") })
  observeEvent(input$reset_upset, { updateNumericInput(session, "upset_min_reps", value = 1); reset_common("upset", palette = FALSE) })
  observeEvent(input$reset_phys, { updateNumericInput(session, "phys_min_reps", value = 1); updateSelectInput(session, "annotation_preset", selected = "human"); updateTextInput(session, "annotation_file", value = default_annotation_file); reset_common("phys") })
  observeEvent(input$reset_cor, { cluster_k <- if (is.null(rv$data)) 1 else length(levels(group_info()$Group)); cor_auto_cluster_k(cluster_k); updateSelectInput(session, "cor_method", selected = "pearson"); updateSelectInput(session, "cor_order", selected = "original"); updateSelectInput(session, "cor_cluster_distance", selected = "one_minus_correlation"); updateSelectInput(session, "cor_cluster_linkage", selected = "complete"); updateNumericInput(session, "cor_cluster_k", value = cluster_k); updateNumericInput(session, "cor_digits", value = 2); updateNumericInput(session, "cor_font", value = 8); updateSelectInput(session, "cor_color", selected = "blue_white_red"); updateNumericInput(session, "cor_min", value = -1); updateNumericInput(session, "cor_max", value = 1); reset_common("cor", palette = FALSE, width_pt = default_cor_heatmap_width_pt, height_pt = default_cor_heatmap_height_pt) })
  observeEvent(input$reset_rank, { reset_common("rank") })
  observeEvent(input$reset_cv, { updateNumericInput(session, "cv_max", value = 60); reset_common("cv") })
  observeEvent(input$reset_pca, { updateNumericInput(session, "pca_min_valid", value = 0.5); reset_common("pca") })
  observeEvent(input$reset_umap, { updateNumericInput(session, "umap_min_valid", value = 0.5); updateNumericInput(session, "umap_neighbors", value = 10); updateNumericInput(session, "umap_min_dist", value = 0.1); reset_common("umap") })
  observeEvent(input$reset_tsne, { updateNumericInput(session, "tsne_min_valid", value = 0.5); updateTextInput(session, "tsne_perplexity", value = "Auto"); updateNumericInput(session, "tsne_max_iter", value = 1000); reset_common("tsne") })
  observeEvent(input$reset_volcano, { updateSelectInput(session, "volcano_fc_method", selected = "log2_then_diff"); updateSelectInput(session, "volcano_test_method", selected = "limma"); updateSelectInput(session, "volcano_sig_metric", selected = "adj_p"); updateNumericInput(session, "log2fc", value = 1); updateNumericInput(session, "adj_p_cutoff", value = 0.05); updateNumericInput(session, "raw_p_cutoff", value = 0.05); reset_common("volcano", palette = FALSE) })
  reset_heatmap_controls <- function(prefix) {
    cluster_k <- if (is.null(rv$data)) 1 else length(levels(group_info()$Group))
    heatmap_auto_cluster_k[[prefix]] <- cluster_k
    updateSelectInput(session, paste0(prefix, "_row_cluster"), selected = "hclust")
    updateNumericInput(session, paste0(prefix, "_row_hclust_k"), value = 1)
    updateNumericInput(session, paste0(prefix, "_row_k"), value = 1)
    updateSelectInput(session, paste0(prefix, "_col_cluster"), selected = "hclust")
    updateNumericInput(session, paste0(prefix, "_col_hclust_k"), value = cluster_k)
    updateNumericInput(session, paste0(prefix, "_col_k"), value = cluster_k)
    updateSelectInput(session, paste0(prefix, "_row_labels"), selected = "none")
    updateSelectInput(session, paste0(prefix, "_color_scheme"), selected = "blue_white_red")
    updateSelectInput(session, paste0(prefix, "_group_palette"), selected = "npg")
    updateSelectInput(session, paste0(prefix, "_cluster_palette"), selected = "lancet")
  }
  observeEvent(input$reset_exprhm, { updateNumericInput(session, "hm_top", value = 100); reset_heatmap_controls("hm"); reset_common("exprhm", palette = FALSE, width_pt = default_expr_heatmap_width_pt, height_pt = default_expr_heatmap_height_pt) })
  observeEvent(input$reset_rf, { updateNumericInput(session, "rf_seed", value = 123); updateSelectInput(session, "rf_split_mode", selected = "auto"); updateNumericInput(session, "rf_train_prop", value = 0.7); updateCheckboxInput(session, "rf_small_sample", value = FALSE); updateNumericInput(session, "rf_top", value = 30); updateNumericInput(session, "rf_ntree", value = 500); updateTextInput(session, "rf_mtry", value = "Auto"); reset_common("rf") })
  observeEvent(input$reset_l1, { updateNumericInput(session, "l1_seed", value = 123); updateSelectInput(session, "l1_split_mode", selected = "auto"); updateNumericInput(session, "l1_train_prop", value = 0.7); updateCheckboxInput(session, "l1_small_sample", value = FALSE); updateNumericInput(session, "l1_alpha", value = 1); updateSelectInput(session, "l1_lambda", selected = "lambda.1se"); updateTextInput(session, "l1_folds", value = "Auto"); updateNumericInput(session, "l1_top", value = 50); reset_common("l1") })
  observeEvent(input$reset_rfl1, { updateNumericInput(session, "rfl1_seed", value = 123); updateSelectInput(session, "rfl1_split_mode", selected = "auto"); updateNumericInput(session, "rfl1_train_prop", value = 0.7); updateCheckboxInput(session, "rfl1_small_sample", value = FALSE); updateNumericInput(session, "rfl1_top", value = 50); updateNumericInput(session, "rfl1_ntree", value = 500); updateTextInput(session, "rfl1_mtry", value = "Auto"); updateNumericInput(session, "rfl1_alpha", value = 1); updateSelectInput(session, "rfl1_lambda", selected = "lambda.1se"); updateTextInput(session, "rfl1_folds", value = "Auto"); reset_common("rfl1") })
  observeEvent(input$reset_feature_umap, { updateSelectInput(session, "feature_umap_source", selected = "rf"); updateNumericInput(session, "feature_top_umap", value = 50); updateNumericInput(session, "feature_umap_neighbors", value = 10); updateNumericInput(session, "feature_umap_min_dist", value = 0.1); reset_common("feature_umap") })
  observeEvent(input$reset_feature_hm, { updateSelectInput(session, "feature_hm_source", selected = "rf"); updateNumericInput(session, "feature_top_hm", value = 50); reset_heatmap_controls("feature_hm"); reset_common("feature_hm", palette = FALSE, width_pt = default_expr_heatmap_width_pt, height_pt = default_expr_heatmap_height_pt) })
  observeEvent(input$reset_sling, { updateSelectInput(session, "sling_reduction", selected = "PCA"); updateNumericInput(session, "sling_top", value = 50); updateNumericInput(session, "sling_heatmap_w_pt", value = default_expr_heatmap_width_pt); updateNumericInput(session, "sling_heatmap_h_pt", value = default_expr_heatmap_height_pt); reset_common("sling") })
  annotation_path <- reactive({
    preset <- input$annotation_preset %||% "human"
    if (preset == "human") return(file.path(annotation_dir, "uniprot_reviewed_human_9606_annotations.csv"))
    if (preset == "mouse") return(file.path(annotation_dir, "uniprot_reviewed_mouse_10090_annotations.csv"))
    if (preset == "celegans") return(file.path(annotation_dir, "uniprot_all_celegans_6239_annotations.csv"))
    input$annotation_file
  })

  preview_pdf <- function(pdf) {
    png_base <- file.path(dirname(pdf), paste0(tools::file_path_sans_ext(basename(pdf)), "_preview_", as.integer(Sys.time()), "_", sample.int(1000000, 1)))
    png_template <- paste0(png_base, "_%d.%s")
    png <- paste0(png_base, "_1.png")
    ok <- FALSE
    if (requireNamespace("pdftools", quietly = TRUE) && file.exists(pdf)) {
      try({ pdftools::pdf_convert(pdf, format = "png", pages = 1, dpi = 144, filenames = png_template, verbose = FALSE); ok <- file.exists(png) }, silent = TRUE)
    }
    if (ok) png else NA_character_
  }
  json_escape <- function(x) {
    x <- as.character(x %||% "")
    x <- gsub("\\\\", "\\\\\\\\", x)
    x <- gsub('"', '\\"', x)
    x <- gsub("\n", "\\\\n", x, fixed = TRUE)
    x
  }
  analysis_parameters <- function(id) {
    switch(id,
      pca = list(min_valid_fraction = input$pca_min_valid),
      umap = list(min_valid_fraction = input$umap_min_valid, n_neighbors = input$umap_neighbors, min_dist = input$umap_min_dist),
      tsne = list(min_valid_fraction = input$tsne_min_valid, perplexity = input$tsne_perplexity, max_iter = input$tsne_max_iter),
      rf = list(seed = input$rf_seed, split_mode = input$rf_split_mode, train_prop = input$rf_train_prop, small_sample = isTRUE(input$rf_small_sample), top_n = input$rf_top, ntree = input$rf_ntree, mtry = input$rf_mtry),
      l1 = list(seed = input$l1_seed, split_mode = input$l1_split_mode, train_prop = input$l1_train_prop, small_sample = isTRUE(input$l1_small_sample), alpha = input$l1_alpha, lambda = input$l1_lambda, cv_folds = input$l1_folds, top_n = input$l1_top),
      rfl1 = list(seed = input$rfl1_seed, split_mode = input$rfl1_split_mode, train_prop = input$rfl1_train_prop, small_sample = isTRUE(input$rfl1_small_sample), top_n = input$rfl1_top, ntree = input$rfl1_ntree, mtry = input$rfl1_mtry, alpha = input$rfl1_alpha, lambda = input$rfl1_lambda, cv_folds = input$rfl1_folds),
      volcano = list(reference_group = input$volcano_a, comparison_group = input$volcano_b, fc_method = input$volcano_fc_method, test_method = input$volcano_test_method, significance_metric = input$volcano_sig_metric, log2fc_cutoff = input$log2fc, adjusted_p_cutoff = input$adj_p_cutoff, raw_p_cutoff = input$raw_p_cutoff),
      cor = list(method = input$cor_method, sample_order = input$cor_order, clustering_distance = input$cor_cluster_distance, clustering_linkage = input$cor_cluster_linkage, clustering_k = input$cor_cluster_k),
      exprhm = list(top_n = input$hm_top, row_clustering = input$hm_row_cluster, row_hierarchical_k = input$hm_row_hclust_k, row_kmeans_k = input$hm_row_k, column_clustering = input$hm_col_cluster, column_hierarchical_k = input$hm_col_hclust_k, column_kmeans_k = input$hm_col_k, row_labels = input$hm_row_labels, color_scheme = input$hm_color_scheme, group_palette = input$hm_group_palette, group_annotation_colors = heatmap_annotation_colors("hm", "group", levels(group_info()$Group)), cluster_palette = input$hm_cluster_palette, cluster_annotation_colors = heatmap_annotation_colors("hm", "cluster", heatmap_cluster_labels("hm"))),
      feature_hm = list(source = input$feature_hm_source, top_n = input$feature_top_hm, row_clustering = input$feature_hm_row_cluster, row_hierarchical_k = input$feature_hm_row_hclust_k, row_kmeans_k = input$feature_hm_row_k, column_clustering = input$feature_hm_col_cluster, column_hierarchical_k = input$feature_hm_col_hclust_k, column_kmeans_k = input$feature_hm_col_k, row_labels = input$feature_hm_row_labels, color_scheme = input$feature_hm_color_scheme, group_palette = input$feature_hm_group_palette, group_annotation_colors = heatmap_annotation_colors("feature_hm", "group", levels(group_info()$Group)), cluster_palette = input$feature_hm_cluster_palette, cluster_annotation_colors = heatmap_annotation_colors("feature_hm", "cluster", heatmap_cluster_labels("feature_hm"))),
      list()
    )
  }
  write_analysis_manifest <- function(id, pdfs, csvs, note = NULL) {
    if (!requireNamespace("jsonlite", quietly = TRUE) || is.null(rv$data)) return(invisible(NULL))
    dat <- data_for_analysis()
    files <- unique(normalizePath(c(pdfs[file.exists(pdfs)], csvs[file.exists(csvs)]), winslash = "/", mustWork = FALSE))
    manifest <- list(
      app = list(name = "ProteoPostZ", version = app_version, generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
      input = list(
        file = normalizePath(input$file_path, winslash = "/", mustWork = FALSE),
        category = dat$input_family %||% dat$input_source %||% dat$software %||% "",
        format = dat$format_label %||% dat$input_source %||% dat$software %||% "",
        detection_evidence = dat$format_evidence %||% "",
        data_level = dat$data_level %||% "protein",
        features = nrow(dat$quantity),
        samples = ncol(dat$quantity),
        zero_handling = if (is_standard_matrix_data(dat)) paste0(dat$missingness_mode, " - ", dat$missingness_mode_label) else "not applicable"
      ),
      grouping = sample_metadata(),
      analysis = list(id = id, parameters = analysis_parameters(id), note = note %||% ""),
      preprocessing = list(default_matrix_preprocessing = "Most quantitative plots and ML use log2(x + 1), valid-feature filtering, and median imputation inside preprocess_expr where applicable."),
      outputs = files
    )
    json_path <- file.path(outdir(), "analysis_manifest.json")
    jsonlite::write_json(manifest, json_path, pretty = TRUE, auto_unbox = TRUE, na = "null")
    html_path <- file.path(outdir(), "analysis_summary.html")
    rows <- paste0("<li>", json_escape(files), "</li>", collapse = "\n")
    html <- paste0(
      "<!doctype html><html><head><meta charset=\"utf-8\"><title>ProteoPostZ Analysis Summary</title>",
      "<style>body{font-family:Arial,'Microsoft YaHei',sans-serif;margin:24px;line-height:1.45;color:#263942}code{background:#f4f6f7;padding:2px 4px;border-radius:3px}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:5px 8px}</style></head><body>",
      "<h1>ProteoPostZ Analysis Summary</h1>",
      "<p><strong>Version:</strong> ", json_escape(app_version), "</p>",
      "<p><strong>Input:</strong> <code>", json_escape(manifest$input$file), "</code></p>",
      "<p><strong>Format:</strong> ", json_escape(manifest$input$category), " / ", json_escape(manifest$input$format), "</p>",
      "<p><strong>Detection evidence:</strong> ", json_escape(manifest$input$detection_evidence), "</p>",
      "<p><strong>Latest analysis:</strong> ", json_escape(id), "</p>",
      "<p><strong>Manifest:</strong> <code>", json_escape(json_path), "</code></p>",
      "<h2>Output files</h2><ul>", rows, "</ul>",
      "</body></html>"
    )
    writeLines(html, html_path, useBytes = TRUE)
    invisible(list(json = json_path, html = html_path))
  }
  finish <- function(id, pdfs, csvs = character(), export_csv = TRUE, note = NULL) {
    if (!export_csv && length(csvs) > 0) unlink(csvs[file.exists(csvs)], force = TRUE)
    valid_pdfs <- pdfs[file.exists(pdfs)]
    pngs <- character()
    old_pngs <- rv$preview_files[[id]]
    if (!is.null(old_pngs)) unlink(old_pngs[file.exists(old_pngs)], force = TRUE)
    if (length(valid_pdfs) > 0) {
      pngs <- vapply(valid_pdfs, preview_pdf, character(1))
      pngs <- pngs[!is.na(pngs) & file.exists(pngs)]
    }
    rv$preview_files[[id]] <- pngs
    rv$preview_index[[id]] <- if (length(pngs) > 0) 1 else NULL
    rv$preview[[id]] <- if (length(pngs) > 0) pngs[[1]] else NULL
    rv$preview_token[[id]] <- as.numeric(Sys.time())
    csv_msg <- if (export_csv && any(file.exists(csvs))) paste("CSV:", paste(normalizePath(csvs[file.exists(csvs)], winslash = "/"), collapse = "\n     ")) else "CSV: not generated"
    pdf_msg <- paste("PDF:", paste(normalizePath(pdfs[file.exists(pdfs)], winslash = "/"), collapse = "\n     "))
    if (!is.null(rv$data)) try(data.table::fwrite(sample_metadata(), file.path(outdir(), "sample_metadata.csv")), silent = TRUE)
    manifest_paths <- tryCatch(write_analysis_manifest(id, pdfs, if (export_csv) csvs else character(), note), error = function(e) NULL)
    manifest_msg <- if (!is.null(manifest_paths)) paste("Summary:", normalizePath(manifest_paths$html, winslash = "/", mustWork = FALSE), "\nManifest:", normalizePath(manifest_paths$json, winslash = "/", mustWork = FALSE)) else NULL
    rv$status[[id]] <- paste(c(note, pdf_msg, csv_msg, manifest_msg), collapse = "\n")
  }

  fail_analysis <- function(id, label, err) {
    while (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
    detail <- conditionMessage(err)
    if (!nzchar(detail)) detail <- "Current input is incomplete or not suitable for this analysis."
    msg <- paste0(label, " failed: ", detail)
    rv$status[[id]] <- msg
    showNotification(msg, type = "error", duration = 10)
    invisible(NULL)
  }
  run_analysis <- function(id, label, expr) {
    rv$status[[id]] <- paste0(label, " running...")
    withProgress(message = paste("Running", label), value = 0.1, {
      tryCatch({
        force(expr)
        incProgress(0.9)
      }, error = function(e) fail_analysis(id, label, e))
    })
  }
  sample_input_id <- function(prefix, i) paste0(prefix, "_", i)
  sample_metadata <- reactive({
    req(rv$data)
    originals <- rv$data$samples
    idx <- seq_along(originals)
    display <- vapply(idx, function(i) input[[sample_input_id("sample_name", i)]] %||% originals[[i]], character(1))
    display <- trimws(display)
    suffix <- vapply(idx, function(i) input[[sample_input_id("grp_suffix", i)]] %||% rv$groups[[originals[[i]]]] %||% "1", character(1))
    suffix <- clean_group_suffix(suffix)
    data.frame(
      original_name = originals,
      display_name = display,
      group_suffix = suffix,
      group = paste0("Group", suffix),
      stringsAsFactors = FALSE
    )
  })

  sample_name_problem <- reactive({
    req(rv$data)
    meta <- sample_metadata()
    if (any(!nzchar(meta$display_name))) return("Sample names cannot be empty.")
    dup <- unique(meta$display_name[duplicated(meta$display_name)])
    if (length(dup) > 0) return(paste("Duplicate sample names:", paste(dup, collapse = ", ")))
    ""
  })

  valid_sample_names <- reactive({
    req(rv$data)
    identical(sample_name_problem(), "")
  })

  data_for_analysis <- reactive({
    req(rv$data)
    meta <- sample_metadata()
    validate(
      need(valid_sample_names(), sample_name_problem())
    )
    dat <- rv$data
    colnames(dat$quantity) <- meta$display_name
    colnames(dat$qualitative) <- meta$display_name
    if (!is.null(dat$analysis_quant_matrix)) colnames(dat$analysis_quant_matrix) <- meta$display_name
    if (!is.null(dat$missingness_matrix)) colnames(dat$missingness_matrix) <- meta$display_name
    if (!is.null(dat$explicit_missing_matrix)) colnames(dat$explicit_missing_matrix) <- meta$display_name
    if (!is.null(dat$ibaq)) colnames(dat$ibaq) <- meta$display_name
    dat$counts$OriginalSample <- dat$counts$Sample
    dat$counts$Sample <- meta$display_name[match(dat$counts$OriginalSample, meta$original_name)]
    dat$samples <- meta$display_name
    dat$sample_metadata <- meta
    dat
  })

  observeEvent(input$load_data, {
    req(input$file_path)
    clear_loaded_input()
    withProgress(message = "Loading data", value = 0.2, {
      tryCatch({
        source_type <- input$input_family %||% "dia"
        if (source_type == "standard_matrix") {
          mode <- input$standard_zero_mode %||% ""
          if (!nzchar(mode)) stop("Please select how zero should be interpreted before loading a standard quantitative matrix.")
          dat <- load_input_dataset(input$file_path, "standard_matrix", "standard_matrix", "d", input$row_id, mode)
        } else if (source_type == "dia") {
          dat <- load_input_dataset(input$file_path, "dia", input$dia_format %||% "auto", "d", input$row_id)
        } else if (source_type == "dda") {
          dat <- load_input_dataset(input$file_path, "dda", input$dda_format %||% "auto", "d", input$row_id)
        } else {
          stop("Unknown input category: ", source_type)
        }
        rv$data <- dat
        rv$groups <- rep("1", length(dat$samples)); names(rv$groups) <- dat$samples
        if (is_standard_matrix_data(dat)) {
          write_matrix_csv(dat$quantity, file.path(outdir(), "standard_matrix_analysis_quantity_matrix.csv"), id_col = "FeatureID")
          write_matrix_csv(dat$raw_quant_matrix, file.path(outdir(), "standard_matrix_raw_quantity_matrix.csv"), id_col = "FeatureID")
          write_matrix_csv(dat$missingness_matrix, file.path(outdir(), "standard_matrix_missingness_matrix.csv"), id_col = "FeatureID")
          data.table::fwrite(dat$counts, file.path(outdir(), "standard_matrix_available_quantitative_value_counts.csv"))
          data.table::fwrite(standard_input_summary_df(dat), file.path(outdir(), "standard_matrix_input_summary.csv"))
        } else {
          write_matrix_csv(dat$quantity, file.path(outdir(), "protein_sample_quantity_matrix.csv"))
          write_matrix_csv(dat$qualitative, file.path(outdir(), "protein_sample_qualitative_matrix.csv"))
          if (!is.null(dat$ibaq)) write_matrix_csv(dat$ibaq, file.path(outdir(), "protein_sample_ibaq_matrix.csv"))
          data.table::fwrite(dat$counts, file.path(outdir(), "identified_protein_counts.csv"))
        }
        data.table::fwrite(dat$meta, file.path(outdir(), if (is_standard_matrix_data(dat)) "feature_metadata.csv" else "protein_metadata.csv"))
        data.table::fwrite(data.frame(original_name = dat$samples, display_name = dat$samples, group_suffix = "1", group = "Group1", stringsAsFactors = FALSE), file.path(outdir(), "sample_metadata.csv"))
      }, error = function(e) {
        rv$load_error <- conditionMessage(e)
        showNotification(rv$load_error, type = "error", duration = 10)
      })
    })
  })
  observeEvent(input$auto_short_names, {
    req(rv$data)
    meta <- sample_metadata()
    generated <- make_short_sample_names(meta$group_suffix)
    for (i in seq_len(nrow(meta))) updateTextInput(session, sample_input_id("sample_name", i), value = generated[[i]])
  })

  group_info <- reactive({ req(rv$data); meta <- sample_metadata(); make_group_info(meta$display_name, meta$group) })
  cor_auto_cluster_k <- reactiveVal(1)
  heatmap_auto_cluster_k <- reactiveValues(hm = 1, feature_hm = 1)
  heatmap_cluster_labels <- function(prefix) {
    req(rv$data)
    mode <- input[[paste0(prefix, "_col_cluster")]] %||% "hclust"
    if (mode == "none") return(character())
    k_id <- if (mode == "hclust") paste0(prefix, "_col_hclust_k") else paste0(prefix, "_col_k")
    k <- suppressWarnings(as.integer(input[[k_id]] %||% length(levels(group_info()$Group))))
    if (!is.finite(k) || k < 1) k <- 1
    paste0("Cluster ", seq_len(min(k, ncol(data_for_analysis()$quantity))))
  }
  heatmap_custom_color_ui <- function(prefix, type) {
    req(rv$data)
    labels <- if (type == "group") levels(group_info()$Group) else heatmap_cluster_labels(prefix)
    palette <- if (type == "group") "npg" else "lancet"
    defaults <- stats::setNames(sci_palette(length(labels), palette), labels)
    tagList(lapply(seq_along(labels), function(i) textInput(paste0(prefix, "_", type, "_color_", i), labels[[i]], value = defaults[[i]], placeholder = "#RRGGBB")))
  }
  output$hm_group_custom_colors <- renderUI({ heatmap_custom_color_ui("hm", "group") })
  output$hm_cluster_custom_colors <- renderUI({ heatmap_custom_color_ui("hm", "cluster") })
  output$feature_hm_group_custom_colors <- renderUI({ heatmap_custom_color_ui("feature_hm", "group") })
  output$feature_hm_cluster_custom_colors <- renderUI({ heatmap_custom_color_ui("feature_hm", "cluster") })
  heatmap_annotation_colors <- function(prefix, type, labels) {
    if (length(labels) == 0) return(NULL)
    palette <- input[[paste0(prefix, "_", type, "_palette")]] %||% if (type == "group") "npg" else "lancet"
    if (palette != "custom") return(stats::setNames(sci_palette(length(labels), palette), labels))
    values <- vapply(seq_along(labels), function(i) trimws(input[[paste0(prefix, "_", type, "_color_", i)]] %||% ""), character(1))
    stats::setNames(values, labels)
  }
  heatmap_row_labels <- function(mode) {
    if (identical(mode, "none")) return(NULL)
    dat <- data_for_analysis()
    meta <- dat$meta
    value_col <- switch(mode, accession = "Accession", protein_name = "ProteinName", gene_name = "GeneName", NULL)
    if (is.null(value_col) || !value_col %in% colnames(meta)) {
      labels <- rep("", nrow(dat$quantity))
    } else {
      key_col <- if ("AnalysisID" %in% colnames(meta)) "AnalysisID" else "RowID"
      labels <- as.character(meta[[value_col]][match(rownames(dat$quantity), meta[[key_col]])])
      labels[is.na(labels)] <- ""
    }
    stats::setNames(labels, rownames(dat$quantity))
  }
  sync_heatmap_cluster_k <- function(prefix, group_k) {
    ids <- c(paste0(prefix, "_col_hclust_k"), paste0(prefix, "_col_k"))
    prior_k <- heatmap_auto_cluster_k[[prefix]]
    for (id in ids) {
      current_k <- suppressWarnings(as.integer(input[[id]] %||% prior_k))
      if (is.na(current_k) || current_k == prior_k) updateNumericInput(session, id, value = group_k)
    }
    heatmap_auto_cluster_k[[prefix]] <- group_k
  }
  observeEvent(group_info(), {
    group_k <- length(levels(group_info()$Group))
    current_k <- suppressWarnings(as.integer(input$cor_cluster_k %||% cor_auto_cluster_k()))
    if (is.na(current_k) || current_k == cor_auto_cluster_k()) {
      updateNumericInput(session, "cor_cluster_k", value = group_k)
    }
    cor_auto_cluster_k(group_k)
    sync_heatmap_cluster_k("hm", group_k)
    sync_heatmap_cluster_k("feature_hm", group_k)
  }, ignoreInit = FALSE)
  output$group_inputs <- renderUI({
    req(rv$data)
    originals <- rv$data$samples
    tagList(lapply(seq_along(originals), function(i) {
      layout_columns(col_widths = c(5, 3, 4),
        div(class = "small-note", originals[[i]]),
        textInput(sample_input_id("grp_suffix", i), "Group suffix", value = rv$groups[[originals[[i]]]] %||% "1", placeholder = "1"),
        textInput(sample_input_id("sample_name", i), "Display sample name", value = originals[[i]])
      )
    }))
  })
  output$sample_name_warning <- renderText({
    req(rv$data)
    sample_name_problem()
  })
  output$input_error <- renderText({ rv$load_error %||% "" })
  output$input_overview_title <- renderUI({
    if (is_standard_matrix_data()) {
      "Available quantitative values overview"
    } else {
      "Sample count"
    }
  })
  output$n_proteins <- renderText({ req(rv$data); paste0(nrow(rv$data$quantity), if (is_standard_matrix_data()) " features" else " proteins") })
  output$n_samples <- renderText({ req(rv$data); paste0(ncol(rv$data$quantity), " samples") })
  output$n_groups <- renderText({ req(rv$data); paste0(length(unique(group_info()$Group)), " groups") })
  output$input_summary <- renderDT({
    req(rv$data)
    df <- if (is_standard_matrix_data()) standard_input_summary_df(rv$data) else input_format_summary_df(rv$data)
    datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })
  output$idbar_count_control <- renderUI({
    req(rv$data)
    selectInput(
      "count_metric",
      "Count type",
      choices = count_metric_choices(rv$data),
      selected = if (is_standard_matrix_data() || (identical(rv$data$input_source %||% "", "peaks") && identical(rv$data$peaks_subtype %||% "", "lfq"))) "Quantified_Protein_Count" else "Identified_Protein_Count"
    )
  })
  output$idbar_count_note <- renderUI({
    req(rv$data)
    if (is_standard_matrix_data()) {
      return(div(
        class = "small-note",
        "Standard quantitative matrix does not provide independent sample-level identification evidence. Counts shown here are available quantitative values only."
      ))
    }
    if (identical(rv$data$input_source %||% "", "maxquant")) {
      return(div(
        class = "small-note",
        rv$data$count_approximation_note_en %||% ""
      ))
    }
    if (identical(rv$data$input_source %||% "", "peaks") && identical(rv$data$peaks_subtype %||% "", "lfq")) {
      return(div(
        class = "small-note",
        rv$data$count_approximation_note_en %||% "Sample-level identification counts are unavailable for this PEAKS LFQ result."
      ))
    }
    NULL
  })
  output$counts_table <- renderDT({
    req(rv$data)
    datatable(dplyr::left_join(data_for_analysis()$counts, sample_metadata(), by = c("Sample" = "display_name", "OriginalSample" = "original_name")), options = list(pageLength = 8))
  })
  venn_upset_note_ui <- function() {
    if (is_standard_matrix_data()) {
      return(div(class = "small-note", "Venn/UpSet use group-level quantified-feature sets based on non-missing available quantitative values after the selected zero handling mode. Recommended: Venn for 2-4 groups; UpSet for 5 or more groups."))
    } else {
      return(div(class = "small-note", "Venn/UpSet use group-level protein sets, not sample-level sets. Each group set is decided by Minimum replicates detected in group. Recommended: Venn for 2-4 groups; UpSet for 5 or more groups."))
    }
  }
  output$venn_note <- renderUI({ venn_upset_note_ui() })
  output$upset_note <- renderUI({ venn_upset_note_ui() })
  output$volcano_groups <- renderUI({ req(rv$data); gs <- levels(group_info()$Group); tagList(selectInput("volcano_a", "Reference group", choices = gs), selectInput("volcano_b", "Comparison group", choices = gs, selected = gs[min(2, length(gs))])) })
  output$slingshot_groups <- renderUI({ req(rv$data); gs <- levels(group_info()$Group); tagList(selectInput("sling_start", "Start group", choices = gs, selected = gs[1]), selectInput("sling_end", "Optional end group", choices = c("None", gs), selected = "None")) })

  observeEvent(input$run_idbar, {
    run_analysis("idbar", count_metric_text()$barplot_title, {
      req(rv$data)
      d <- analysis_dir("idbar")
      wh <- pts("idbar")
      count_col <- active_count_column()
      labels <- count_metric_text(count_col)
      if (identical(count_col, "Identified_Protein_Count") && all(is.na(data_for_analysis()$counts[[count_col]]))) {
        stop(rv$data$count_approximation_note_en %||% "Sample-level identification counts are unavailable for this input.", call. = FALSE)
      }
      if (is_standard_matrix_data()) {
        pdf <- file.path(d, "available_quantitative_values_barplot.pdf")
        csv <- file.path(d, "available_quantitative_values_group_summary.csv")
      } else if (identical(count_col, "Identified_Protein_Count")) {
        pdf <- file.path(d, "identification_barplot.pdf")
        csv <- file.path(d, "identification_group_summary.csv")
      } else {
        pdf <- file.path(d, "quantified_protein_barplot.pdf")
        csv <- file.path(d, "quantified_protein_group_summary.csv")
      }
      plot_sample_count_bar(
        data_for_analysis()$counts,
        group_info(),
        count_col,
        pdf,
        csv,
        wh[1],
        wh[2],
        input$idbar_palette,
        y_label = labels$y_label,
        plot_title = labels$barplot_title
      )
      finish(
        "idbar",
        pdf,
        c(csv, sub("\\.csv$", "_sample_counts.csv", csv)),
        input$idbar_export_csv,
        paste(c(labels$title, standard_mode_note()), collapse = "\n")
      )
    })
  })

  make_sets <- function(min_reps) identified_by_group(data_for_analysis()$qualitative, group_info(), min_reps)
  validate_current_sets <- function(sets, min_reps, analysis) {
    tryCatch(validate_group_sets(sets, min_reps, analysis), error = function(e) {
      msg <- conditionMessage(e)
      if (is_standard_matrix_data()) msg <- gsub("protein sets", "quantified-feature sets", msg, fixed = TRUE)
      stop(msg, call. = FALSE)
    })
  }
  membership_for_current_input <- function(sets) {
    set_df <- make_set_membership(sets)
    if (is_standard_matrix_data() && "ProteinID" %in% colnames(set_df)) colnames(set_df)[colnames(set_df) == "ProteinID"] <- "FeatureID"
    set_df
  }
  observeEvent(input$run_venn, { run_analysis("venn", "Venn diagram", { req(rv$data); d <- analysis_dir("venn"); wh <- pts("venn"); sets <- validate_current_sets(make_sets(input$venn_min_reps), input$venn_min_reps, "venn"); pdf <- file.path(d, "venn.pdf"); csv <- file.path(d, "venn_membership.csv"); set_df <- membership_for_current_input(sets); if (input$venn_export_csv) data.table::fwrite(set_df, csv); grDevices::pdf(pdf, width = wh[1], height = wh[2]); grid::grid.draw(VennDiagram::venn.diagram(sets, filename = NULL, fill = sci_palette(length(sets), input$venn_palette), alpha = 0.45, cex = 0.8, cat.cex = 0.8, margin = 0.08)); grDevices::dev.off(); finish("venn", pdf, csv, input$venn_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_upset, { run_analysis("upset", "UpSet plot", { req(rv$data); d <- analysis_dir("upset"); wh <- pts("upset"); sets <- validate_current_sets(make_sets(input$upset_min_reps), input$upset_min_reps, "upset"); pdf <- file.path(d, "upset.pdf"); csv <- file.path(d, "upset_membership.csv"); set_df <- membership_for_current_input(sets); if (nrow(set_df) == 0 || ncol(set_df) - 1 < 2) stop("UpSet plot requires a membership table with at least one ", if (is_standard_matrix_data()) "feature" else "protein", " and at least 2 non-empty group columns."); if (input$upset_export_csv) data.table::fwrite(set_df, csv); grDevices::pdf(pdf, width = max(wh[1], 5), height = max(wh[2], 4)); UpSetR::upset(as.data.frame(set_df[, -1, drop = FALSE]), nsets = length(sets), order.by = "freq"); grDevices::dev.off(); finish("upset", pdf, csv, input$upset_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_phys, { run_analysis("phys", "Physicochemical property distributions", {
    req(rv$data)
    d <- analysis_dir("phys")
    wh <- pts("phys")
    sets <- make_sets(input$phys_min_reps)
    meta <- data_for_analysis()$meta
    if (is_standard_matrix_data()) {
      sets_acc <- lapply(sets, function(ids) unique(ids[nzchar(ids) & !is.na(ids)]))
      annotation_note <- "Standard matrix feature identifiers are matched to the annotation Accession column as-is; if few entries match, the first column may not be UniProt accessions for the selected species."
    } else {
      sets_acc <- lapply(sets, function(ids) unique(analysis_ids_to_accessions(meta, ids)))
      annotation_note <- NULL
    }
    run_physicochemical(sets_acc, annotation_path(), d, wh[1], wh[2], input$phys_palette)
    phys_pdfs <- list.files(d, pattern = "\\.pdf$", full.names = TRUE)
    phys_csvs <- list.files(d, pattern = "\\.csv$", full.names = TRUE)
    finish("phys", phys_pdfs, phys_csvs, input$phys_export_csv, paste(c(standard_mode_note(), annotation_note, paste("Annotation:", annotation_path())), collapse = "\n"))
  }) })

  observeEvent(input$run_cor, { run_analysis("cor", "Correlation heatmap", { req(rv$data); d <- analysis_dir("cor"); wh <- pts("cor"); pdf <- file.path(d, paste0(input$cor_method, "_correlation_heatmap.pdf")); csv <- file.path(d, paste0(input$cor_method, "_correlation_ordered_matrix.csv")); plot_correlation_heatmap(data_for_analysis()$quantity, group_info(), pdf, csv, input$cor_method, input$cor_order, input$cor_cluster_distance, input$cor_cluster_linkage, input$cor_cluster_k, input$cor_digits, input$cor_font, input$cor_color, input$cor_min, input$cor_max, wh[1], wh[2]); finish("cor", pdf, list.files(d, pattern = "\\.csv$", full.names = TRUE), input$cor_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_rank, { run_analysis("rank", "Rank-abundance plot", { req(rv$data); d <- analysis_dir("rank"); wh <- pts("rank"); pdf <- file.path(d, "rank_abundance.pdf"); csv <- file.path(d, "rank_abundance_data.csv"); plot_rank_abundance(data_for_analysis()$quantity, group_info(), pdf, csv, wh[1], wh[2], input$rank_palette); finish("rank", pdf, csv, input$rank_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_cv, { run_analysis("cv", "CV ridgeline", { req(rv$data); d <- analysis_dir("cv"); wh <- pts("cv"); pdf <- file.path(d, "cv_ridgeline.pdf"); csv <- file.path(d, "cv_values.csv"); plot_cv_ridges(data_for_analysis()$quantity, group_info(), pdf, csv, input$cv_max, wh[1], wh[2], input$cv_palette); finish("cv", pdf, c(csv, sub("\\.csv$", "_median.csv", csv)), input$cv_export_csv, standard_mode_note()) }) })

  pca_data <- function(min_valid) { used <- preprocess_expr(data_for_analysis()$quantity, TRUE, min_valid); sample_mat <- t(used); pca <- prcomp(sample_mat, center = TRUE, scale. = TRUE); var <- summary(pca)$importance[2, 1:2] * 100; df <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2]) |> left_join(group_info(), by = "Sample"); list(df = df, var = var) }
  observeEvent(input$run_pca, { run_analysis("pca", "PCA plot", { req(rv$data); d <- analysis_dir("pca"); wh <- pts("pca"); res <- pca_data(input$pca_min_valid); pdf <- file.path(d, "PCA_plot.pdf"); csv <- file.path(d, "PCA_coordinates.csv"); if (input$pca_export_csv) data.table::fwrite(res$df, csv); p <- ggplot2::ggplot(res$df, ggplot2::aes(PC1, PC2, color = Group)) + ggplot2::geom_point(size = 2.4) + ggplot2::scale_color_manual(values = sci_palette(length(levels(group_info()$Group)), input$pca_palette)) + theme_sci() + ggplot2::labs(x = sprintf("PC1 (%.2f%%)", res$var[1]), y = sprintf("PC2 (%.2f%%)", res$var[2])); ggplot2::ggsave(pdf, p, width = wh[1], height = wh[2]); finish("pca", pdf, csv, input$pca_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_umap, { run_analysis("umap", "UMAP plot", { req(rv$data); d <- analysis_dir("umap"); wh <- pts("umap"); used <- preprocess_expr(data_for_analysis()$quantity, TRUE, input$umap_min_valid); sample_mat <- t(used); set.seed(123); nn <- min(input$umap_neighbors, max(2, nrow(sample_mat) - 1)); um <- uwot::umap(sample_mat, n_neighbors = nn, min_dist = input$umap_min_dist, metric = "euclidean", verbose = FALSE); df <- data.frame(Sample = rownames(sample_mat), UMAP1 = um[,1], UMAP2 = um[,2]) |> left_join(group_info(), by = "Sample"); pdf <- file.path(d, "UMAP_plot.pdf"); csv <- file.path(d, "UMAP_coordinates.csv"); if (input$umap_export_csv) data.table::fwrite(df, csv); p <- ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, color = Group)) + ggplot2::geom_point(size = 2.4) + ggplot2::scale_color_manual(values = sci_palette(length(levels(group_info()$Group)), input$umap_palette)) + theme_sci(); ggplot2::ggsave(pdf, p, width = wh[1], height = wh[2]); finish("umap", pdf, csv, input$umap_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_tsne, { run_analysis("tsne", "t-SNE plot", { req(rv$data); if (!requireNamespace("Rtsne", quietly = TRUE)) stop("t-SNE requires the Rtsne package, which is not available in the current runtime."); d <- analysis_dir("tsne"); wh <- pts("tsne"); used <- preprocess_expr(data_for_analysis()$quantity, TRUE, input$tsne_min_valid); sample_mat <- t(used); if (nrow(sample_mat) < 4) stop("t-SNE requires at least four samples after filtering."); max_perplexity <- max(1, floor((nrow(sample_mat) - 1) / 3)); requested <- parse_auto_integer(input$tsne_perplexity, NA_integer_); perplexity <- if (is.na(requested)) min(30, max_perplexity) else min(requested, max_perplexity); if (perplexity < 1) stop("t-SNE perplexity must be at least 1 and less than one-third of the sample count."); set.seed(123); ts <- Rtsne::Rtsne(sample_mat, dims = 2, perplexity = perplexity, max_iter = input$tsne_max_iter, check_duplicates = FALSE, pca = TRUE, verbose = FALSE); df <- data.frame(Sample = rownames(sample_mat), tSNE1 = ts$Y[,1], tSNE2 = ts$Y[,2], Perplexity = perplexity, MaxIter = input$tsne_max_iter) |> left_join(group_info(), by = "Sample"); pdf <- file.path(d, "tSNE_plot.pdf"); csv <- file.path(d, "tSNE_coordinates.csv"); if (input$tsne_export_csv) data.table::fwrite(df, csv); p <- ggplot2::ggplot(df, ggplot2::aes(tSNE1, tSNE2, color = Group)) + ggplot2::geom_point(size = 2.4) + ggplot2::scale_color_manual(values = sci_palette(length(levels(group_info()$Group)), input$tsne_palette)) + theme_sci() + ggplot2::labs(x = "t-SNE 1", y = "t-SNE 2"); ggplot2::ggsave(pdf, p, width = wh[1], height = wh[2]); finish("tsne", pdf, csv, input$tsne_export_csv, standard_mode_note()) }) })

  observeEvent(input$run_volcano, { run_analysis("volcano", "Volcano plot", { req(rv$data); d <- analysis_dir("volcano"); wh <- pts("volcano"); pdf <- file.path(d, paste0("volcano_", input$volcano_b, "_vs_", input$volcano_a, ".pdf")); csv <- file.path(d, paste0("volcano_", input$volcano_b, "_vs_", input$volcano_a, ".csv")); run_volcano(data_for_analysis()$quantity, group_info(), input$volcano_a, input$volcano_b, pdf, csv, input$log2fc, input$adj_p_cutoff, input$raw_p_cutoff, input$volcano_fc_method, input$volcano_test_method, input$volcano_sig_metric, wh[1], wh[2]); finish("volcano", pdf, csv, input$volcano_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_exprhm, { run_analysis("exprhm", "Expression heatmap", { req(rv$data); d <- analysis_dir("exprhm"); wh <- pts("exprhm"); pdf <- file.path(d, "expression_heatmap.pdf"); csv <- file.path(d, "expression_heatmap_values.csv"); plot_expression_heatmap(data_for_analysis()$quantity, group_info(), pdf, csv, input$hm_top, input$hm_row_cluster, input$hm_col_cluster, input$hm_row_k, input$hm_col_k, input$hm_row_hclust_k, input$hm_col_hclust_k, input$hm_color_scheme, heatmap_annotation_colors("hm", "group", levels(group_info()$Group)), heatmap_annotation_colors("hm", "cluster", heatmap_cluster_labels("hm")), heatmap_row_labels(input$hm_row_labels), max(wh[1], 3), max(wh[2], 3)); finish("exprhm", pdf, list.files(d, pattern = "\\.csv$", full.names = TRUE), input$exprhm_export_csv, standard_mode_note()) }) })

  observeEvent(input$run_rf, { run_analysis("rf", "Random forest feature selection", { req(rv$data); d <- analysis_dir("rf"); wh <- pts("rf"); top <- run_random_forest_selection(data_for_analysis()$quantity, group_info(), d, input$rf_top, input$rf_ntree, input$rf_seed, input$rf_split_mode, input$rf_train_prop, input$rf_mtry, isTRUE(input$rf_small_sample)); rv$rf_features <- top; imp <- data.table::fread(file.path(d, "random_forest_importance.csv"), data.table = FALSE); plot_df <- head(imp, input$rf_top); plot_df$ProteinID <- factor(plot_df$ProteinID, levels = rev(plot_df$ProteinID)); pdf <- file.path(d, "random_forest_top_importance.pdf"); csvs <- list.files(d, pattern = "\\.csv$", full.names = TRUE); p <- ggplot2::ggplot(plot_df, ggplot2::aes(ProteinID, RFImportance, fill = RFImportance)) + ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::scale_fill_gradient(low = "#DCEAF7", high = sci_palette(1, input$rf_palette)) + theme_sci() + ggplot2::labs(x = NULL, y = "MeanDecreaseGini"); ggplot2::ggsave(pdf, p, width = wh[1], height = max(wh[2], 4)); finish("rf", pdf, csvs, input$rf_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_l1, { run_analysis("l1", "L1 feature selection", { req(rv$data); d <- analysis_dir("l1"); wh <- pts("l1"); top <- run_l1_selection(data_for_analysis()$quantity, group_info(), d, input$l1_top, input$l1_alpha, input$l1_seed, input$l1_split_mode, input$l1_train_prop, input$l1_lambda, input$l1_folds, isTRUE(input$l1_small_sample)); scores <- data.table::fread(file.path(d, "l1_feature_scores.csv"), data.table = FALSE); plot_df <- head(scores, input$l1_top); plot_df$ProteinID <- factor(plot_df$ProteinID, levels = rev(plot_df$ProteinID)); pdf <- file.path(d, "l1_top_coefficients.pdf"); csvs <- list.files(d, pattern = "\\.csv$", full.names = TRUE); if (nrow(plot_df) > 0) { p <- ggplot2::ggplot(plot_df, ggplot2::aes(ProteinID, L1Score, fill = L1Score)) + ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::scale_fill_gradient(low = "#E9EDF3", high = sci_palette(1, input$l1_palette)) + theme_sci() + ggplot2::labs(x = NULL, y = "Sum absolute L1 coefficient"); ggplot2::ggsave(pdf, p, width = wh[1], height = max(wh[2], 4)) }; finish("l1", pdf, csvs, input$l1_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_rfl1, { run_analysis("rfl1", "RF + L1 feature selection", { req(rv$data); d <- analysis_dir("rfl1"); wh <- pts("rfl1"); top <- run_feature_selection(data_for_analysis()$quantity, group_info(), d, input$rfl1_top, input$rfl1_ntree, input$rfl1_alpha, input$rfl1_seed, input$rfl1_split_mode, input$rfl1_train_prop, input$rfl1_mtry, input$rfl1_lambda, input$rfl1_folds, isTRUE(input$rfl1_small_sample)); rv$rf_features <- top; rf <- data.table::fread(file.path(d, "random_forest_importance.csv"), data.table = FALSE); l1 <- data.table::fread(file.path(d, "l1_feature_scores.csv"), data.table = FALSE); summary_df <- data.frame(ProteinID = top, CombinedRank = seq_along(top), stringsAsFactors = FALSE) |> dplyr::left_join(rf, by = "ProteinID") |> dplyr::left_join(l1, by = "ProteinID") |> dplyr::mutate(Source = dplyr::case_when(!is.na(RFImportance) & !is.na(L1Score) ~ "RF + L1", !is.na(RFImportance) ~ "RF only", !is.na(L1Score) ~ "L1 only", TRUE ~ "Other")); summary_csv <- file.path(d, "combined_feature_summary.csv"); data.table::fwrite(summary_df, summary_csv); plot_df <- summary_df; plot_df$ProteinID <- factor(plot_df$ProteinID, levels = rev(plot_df$ProteinID)); pdf <- file.path(d, "rf_l1_combined_features.pdf"); p <- ggplot2::ggplot(plot_df, ggplot2::aes(ProteinID, CombinedRank, fill = Source)) + ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = c("RF + L1" = sci_palette(1, input$rfl1_palette), "RF only" = "#8FB9D9", "L1 only" = "#E8A35D", "Other" = "grey70")) + ggplot2::scale_y_reverse() + theme_sci() + ggplot2::labs(x = NULL, y = "Combined rank"); ggplot2::ggsave(pdf, p, width = wh[1], height = max(wh[2], 4)); finish("rfl1", pdf, list.files(d, pattern = "\\.csv$", full.names = TRUE), input$rfl1_export_csv, standard_mode_note()) }) })
  read_feature_ids <- function(files, label) {
    files <- files[file.exists(files)]
    if (length(files) == 0) stop(label, " selected proteins have not been generated yet. Please run the corresponding machine-learning analysis first.")
    f <- files[[length(files)]]
    dat <- data.table::fread(f, data.table = FALSE)
    id_col <- intersect(c("ProteinID", "ProteinName", "RowID"), colnames(dat))
    ids <- if (length(id_col) > 0) dat[[id_col[[1]]]] else dat[[1]]
    ids <- unique(trimws(as.character(ids)))
    ids[nzchar(ids) & !is.na(ids)]
  }
  ml_feature_files <- function(source) {
    switch(source,
      rf = list(files = c(list.files(analysis_dir("rf"), pattern = "^top.*_rf_features\\.csv$", full.names = TRUE), file.path(analysis_dir("rf"), "random_forest_importance.csv")), label = "Random forest"),
      l1 = list(files = c(list.files(analysis_dir("l1"), pattern = "^top.*_l1_features\\.csv$", full.names = TRUE), file.path(analysis_dir("l1"), "l1_feature_scores.csv")), label = "L1"),
      rfl1 = list(files = c(list.files(analysis_dir("rfl1"), pattern = "^top.*_rf_l1_union_features\\.csv$", full.names = TRUE), file.path(analysis_dir("rfl1"), "combined_feature_summary.csv")), label = "RF + L1"),
      union = list(files = character(), label = "All available ML"),
      stop("Unknown selected feature source: ", source)
    )
  }
  get_ml_features <- function(source, n) {
    if (source == "union") {
      feats <- unique(unlist(lapply(c("rf", "l1", "rfl1"), function(src) {
        info <- ml_feature_files(src)
        tryCatch(read_feature_ids(info$files, info$label), error = function(e) character())
      }), use.names = FALSE))
      if (length(feats) == 0) stop("No ML selected proteins are available yet. Please run RF, L1, or RF + L1 first.")
      return(head(feats, n))
    }
    info <- ml_feature_files(source)
    head(read_feature_ids(info$files, info$label), n)
  }
  observeEvent(input$run_feature_umap, { run_analysis("feature_umap", "Feature UMAP", { req(rv$data); d <- analysis_dir("feature_umap"); wh <- pts("feature_umap"); feats <- get_ml_features(input$feature_umap_source, input$feature_top_umap); feats <- intersect(feats, rownames(data_for_analysis()$quantity)); if (length(feats) < 2) stop("Feature UMAP requires at least two selected proteins present in the current quantity matrix."); used <- preprocess_expr(data_for_analysis()$quantity[feats, , drop = FALSE], TRUE, 0.5); if (nrow(used) < 2) stop("Feature UMAP requires at least two selected proteins after preprocessing."); sample_mat <- t(used); set.seed(123); nn <- min(input$feature_umap_neighbors, max(2, nrow(sample_mat)-1)); um <- uwot::umap(sample_mat, n_neighbors = nn, min_dist = input$feature_umap_min_dist, metric = "euclidean", verbose = FALSE); df <- data.frame(Sample = rownames(sample_mat), UMAP1 = um[,1], UMAP2 = um[,2]) |> left_join(group_info(), by = "Sample"); pdf <- file.path(d, "feature_UMAP.pdf"); csv <- file.path(d, "feature_UMAP_coordinates.csv"); source_csv <- file.path(d, "feature_UMAP_selected_proteins.csv"); if (input$feature_umap_export_csv) { data.table::fwrite(df, csv); data.table::fwrite(data.frame(ProteinID = rownames(used), Source = input$feature_umap_source), source_csv) }; p <- ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, color = Group)) + ggplot2::geom_point(size = 2.4) + ggplot2::scale_color_manual(values = sci_palette(length(levels(group_info()$Group)), input$feature_umap_palette)) + theme_sci(); ggplot2::ggsave(pdf, p, width = wh[1], height = wh[2]); finish("feature_umap", pdf, c(csv, source_csv), input$feature_umap_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_feature_hm, { run_analysis("feature_hm", "ML selected expression heatmap", { req(rv$data); d <- analysis_dir("feature_hm"); wh <- pts("feature_hm"); feats <- get_ml_features(input$feature_hm_source, input$feature_top_hm); feats <- intersect(feats, rownames(data_for_analysis()$quantity)); if (length(feats) < 2) stop("Feature heatmap requires at least two selected proteins present in the current quantity matrix."); pdf <- file.path(d, "feature_heatmap.pdf"); csv <- file.path(d, "feature_heatmap_values.csv"); source_csv <- file.path(d, "feature_heatmap_selected_proteins.csv"); data.table::fwrite(data.frame(ProteinID = feats, Source = input$feature_hm_source), source_csv); plot_expression_heatmap(data_for_analysis()$quantity[feats, , drop = FALSE], group_info(), pdf, csv, input$feature_top_hm, input$feature_hm_row_cluster, input$feature_hm_col_cluster, input$feature_hm_row_k, input$feature_hm_col_k, input$feature_hm_row_hclust_k, input$feature_hm_col_hclust_k, input$feature_hm_color_scheme, heatmap_annotation_colors("feature_hm", "group", levels(group_info()$Group)), heatmap_annotation_colors("feature_hm", "cluster", heatmap_cluster_labels("feature_hm")), heatmap_row_labels(input$feature_hm_row_labels), max(wh[1], 3), max(wh[2], 3)); finish("feature_hm", pdf, list.files(d, pattern = "\\.csv$", full.names = TRUE), input$feature_hm_export_csv, standard_mode_note()) }) })
  observeEvent(input$run_sling, { run_analysis("sling", "Slingshot", { req(rv$data); d <- analysis_dir("sling"); wh <- pts("sling"); hm_wh <- c(input$sling_heatmap_w_pt / 72, input$sling_heatmap_h_pt / 72); run_slingshot_pseudotime(data_for_analysis()$quantity, group_info(), d, input$sling_reduction, input$sling_start, input$sling_end, wh[1], wh[2], input$sling_palette, top_n = input$sling_top, heatmap_width = hm_wh[1], heatmap_height = hm_wh[2]); finish("sling", list.files(d, pattern = "\\.pdf$", full.names = TRUE), list.files(d, pattern = "\\.(csv|txt)$", full.names = TRUE), input$sling_export_csv, paste(c(standard_mode_note(), "Outputs are under the selected output directory / sling."), collapse = "\n")) }) })

  output$outdir_text <- renderText(outdir())
  output$file_table <- renderDT({ dir.create(outdir(), recursive = TRUE, showWarnings = FALSE); files <- list.files(outdir(), recursive = TRUE, full.names = FALSE); info <- file.info(file.path(outdir(), files)); datatable(data.frame(File = files, SizeKB = round(info$size / 1024, 1), Modified = info$mtime), options = list(pageLength = 12)) })
}

shinyApp(ui, server)
