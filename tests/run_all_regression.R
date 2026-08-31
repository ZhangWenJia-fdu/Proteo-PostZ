script_file <- NULL
if (sys.nframe() > 0L && "ofile" %in% names(parent.frame())) script_file <- get("ofile", envir = parent.frame())
if (is.null(script_file) || !nzchar(script_file)) script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1L])
tests_dir <- normalizePath(dirname(sub("^--file=", "", script_file)), winslash = "/", mustWork = TRUE)
package_root <- normalizePath(file.path(tests_dir, ".."), winslash = "/", mustWork = TRUE)
setwd(package_root)
source(file.path("app", "R", "analysis_core.R"), encoding = "UTF-8")
library(shiny)
source(file.path("app", "R", "ui_helpers.R"), encoding = "UTF-8")
source(file.path("tests", "helpers", "regression_helpers.R"), encoding = "UTF-8")
source(file.path("tests", "helpers", "synthetic_canonical_data.R"), encoding = "UTF-8")
source(file.path("tests", "helpers", "synthetic_format_writers.R"), encoding = "UTF-8")
regression_files <- list.files(file.path(tests_dir, "regression"), pattern = "[.]R$", full.names = TRUE)
for (path in regression_files) source(path, encoding = "UTF-8")

run_group <- function(id, fn, ctx) {
  record_test(id, fn(ctx))
}

root_tmp <- file.path(tempdir(), "proteopostz_regression")
ctx <- list(
  package_root = package_root,
  tests_dir = tests_dir,
  fixtures = file.path(root_tmp, "fixtures"),
  outputs = file.path(root_tmp, "outputs"),
  canonical = make_canonical_fixture(seed = 101, n_features = 80, n_samples = 12)
)
dir.create(ctx$fixtures, recursive = TRUE, showWarnings = FALSE)
dir.create(ctx$outputs, recursive = TRUE, showWarnings = FALSE)

run_group("importers", test_importers, ctx)
run_group("quantitative_qc", test_quantitative_qc, ctx)
run_group("dimension_reduction", test_dimension_reduction, ctx)
run_group("analysis_outputs", test_analysis_outputs, ctx)
run_group("machine_learning", test_ml, ctx)
run_group("trajectory", test_trajectory, ctx)
run_group("gui_acceptance_followup", test_gui_acceptance_followup, ctx)
if (tolower(Sys.getenv("PROTEOPOSTZ_HIGH_FIDELITY", "")) %in% c("1", "true", "yes")) run_group("high_fidelity_ml", test_ml_high_fidelity, ctx)
run_group("app_source_static", function(z) {
  paths <- list.files(file.path(z$package_root, "app"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
  for (p in paths) parse(file = p)
  expect_true(file.exists(file.path(z$package_root, "app", "R", "analysis_core.R")), "analysis_core.R missing")
  app_text <- paste(readLines(file.path(z$package_root, "app", "app.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("draw_upset_membership\\(set_df\\)", app_text), "GUI UpSet observer must call the production membership renderer")
  expect_false(grepl("UpSetR::upset\\(as.data.frame\\(set_df", app_text), "GUI UpSet observer must not bypass the production renderer")
  expect_false(grepl("phys_pdfs <- list.files\\(d", app_text), "GUI physicochemical preview must not scan stale directory PDFs")
  expect_false(grepl("finish\\(\"sling\", list.files\\(d", app_text), "GUI Slingshot preview must not scan stale directory outputs")
  expect_true(grepl("select_ml_feature_ids\\(source, rv\\$ml_results, n\\)", app_text), "Feature modules must use current-session ML result state")
  expect_true(grepl("ml_results = list\\(rf = NULL, l1 = NULL, rfl1 = NULL\\)", app_text), "ML result state must be initialized")
  expect_false(grepl("read_feature_ids <- function|ml_feature_files <- function", app_text), "Feature source selection must not use stale output-directory scans")
}, ctx)

res <- get(".regression_results", envir = .GlobalEnv)
cat("ProteoPostZ regression test\n")
cat("Package root: ", package_root, "\n", sep = "")
cat("High-fidelity mode: ", tolower(Sys.getenv("PROTEOPOSTZ_HIGH_FIDELITY", "")) %in% c("1", "true", "yes"), "\n", sep = "")
cat("PASS: ", sum(res$Status == "PASS"), "  FAIL: ", sum(res$Status == "FAIL"), "  SKIP: ", sum(res$Status == "SKIP"), "\n", sep = "")
if (nrow(res)) print(res[, c("ID", "Status", "Message")], row.names = FALSE)
if (length(.regression_warnings)) {
  cat("Warnings captured: ", length(.regression_warnings), "\n", sep = "")
  print(.regression_warnings)
}
if (keep_test_artifacts()) cat("Artifacts kept at: ", root_tmp, "\n", sep = "")
if (any(res$Status == "FAIL")) quit(status = 1L, save = "no")
