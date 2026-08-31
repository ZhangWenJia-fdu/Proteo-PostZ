# ProteoPostZ core compatibility source entry point
# Each implementation is defined in exactly one module.

core_dir <- if (file.exists(file.path("R", "dependencies.R"))) {
  normalizePath("R", winslash = "/", mustWork = TRUE)
} else if (file.exists(file.path("app", "R", "dependencies.R"))) {
  normalizePath(file.path("app", "R"), winslash = "/", mustWork = TRUE)
} else {
  dirname(normalizePath(sys.frame(1)$ofile %||% "R/analysis_core.R", winslash = "/", mustWork = FALSE))
}
source_order <- c(
  "dependencies.R", "import_io.R", "qualitative_analysis.R", "preprocessing.R",
  "quantitative_qc.R", "quantitative_analysis.R", "dimension_reduction.R",
  "differential_analysis.R", "heatmap_analysis.R", "ml_analysis.R",
  "annotation.R", "trajectory_analysis.R"
)
invisible(lapply(source_order, function(f) source(file.path(core_dir, f), encoding = "UTF-8")))
