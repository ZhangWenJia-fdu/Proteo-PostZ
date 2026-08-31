# ProteoPostZ dependency manifest and package-loading helpers
# ProteoPostZ core functions
# Developed by Wenjia Zhang

`%||%` <- function(a, b) if (!is.null(a)) a else b

dependency_manifest <- list(
  cran_runtime = c(
    "data.table", "dplyr", "tidyr", "ggplot2", "pheatmap", "RColorBrewer",
    "VennDiagram", "UpSetR", "ggridges", "uwot", "randomForest", "glmnet",
    "Peptides", "pdftools", "tibble", "shiny", "bslib", "DT"
  ),
  bioconductor_runtime = c(
    "slingshot", "SingleCellExperiment", "S4Vectors", "limma"
  ),
  r_base_runtime = c(
    "grid", "grDevices", "graphics", "stats", "utils", "tools"
  ),
  cran_annotation_update = c(
    "httr"
  ),
  project_runtime = c(
    "app/R/analysis_core.R",
    "app/annotations/uniprot_all_celegans_6239_annotations.csv",
    "app/annotations/uniprot_reviewed_human_9606_annotations.csv",
    "app/annotations/uniprot_reviewed_mouse_10090_annotations.csv",
    "portable/R-4.5.1",
    "portable/Rlibs"
  ),
  system_runtime = c(
    "Windows shell for launcher scripts",
    "Microsoft .NET runtime for launcher executable",
    "local loopback TCP access to 127.0.0.1:3840",
    "native DLL dependencies bundled with portable R and compiled R packages"
  ),
  uncertain_system_runtime = c(
    "pdftools native PDF rendering dependencies when not fully bundled by the Windows binary package"
  )
)

# Quantitative QC is maintained as an independent module. It is sourced after
# this file's shared plotting helpers are defined (see the footer below).

required_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "pheatmap", "RColorBrewer",
  "VennDiagram", "UpSetR", "ggridges", "uwot", "randomForest", "glmnet",
  "Peptides", "grid", "pdftools", "slingshot", "SingleCellExperiment", "S4Vectors",
  "tibble"
)

load_required_packages <- function() {
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "), ". Please install them before running the app.")
  }
}
