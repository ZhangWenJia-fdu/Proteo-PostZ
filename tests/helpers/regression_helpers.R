# Long-lived, version-agnostic regression harness helpers.
.regression_results <- data.frame(ID = character(), Status = character(), Message = character(), stringsAsFactors = FALSE)
.regression_warnings <- character()
keep_test_artifacts <- function() tolower(Sys.getenv("PROTEOPOSTZ_KEEP_TEST_ARTIFACTS", "")) %in% c("1", "true", "yes")
expect_true <- function(value, message) { if (!isTRUE(value)) stop(message, call. = FALSE); invisible(TRUE) }
expect_false <- function(value, message) expect_true(!isTRUE(value), message)
expect_equal <- function(actual, expected, message, tolerance = NULL) {
  ok <- if (is.null(tolerance)) isTRUE(all.equal(actual, expected, check.attributes = FALSE)) else isTRUE(all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE))
  if (!ok) stop(message, call. = FALSE)
  invisible(TRUE)
}
expect_near <- function(actual, expected, tolerance = 1e-8, message = "Values are not within tolerance") {
  expect_true(length(actual) == length(expected) && all(is.finite(actual)) && all(is.finite(expected)) && max(abs(actual - expected)) <= tolerance, paste0(message, "; max absolute difference = ", max(abs(actual - expected))))
}
expect_file_nonempty <- function(path, message = NULL) {
  info <- if (length(path) == 1 && file.exists(path)) file.info(path) else NULL
  expect_true(!is.null(info) && is.finite(info[["size"]]) && info[["size"]] > 0, message %||% paste("Expected non-empty file:", path))
}
expect_columns_include <- function(dat, columns, message = NULL) expect_true(all(columns %in% colnames(dat)), message %||% paste("Missing columns:", paste(setdiff(columns, colnames(dat)), collapse = ", ")))
expect_columns_exclude <- function(dat, columns, message = NULL) expect_true(!any(columns %in% colnames(dat)), message %||% paste("Unexpected columns:", paste(intersect(columns, colnames(dat)), collapse = ", ")))
expect_no_duplicates <- function(x, message = "Unexpected duplicate identifiers") expect_true(!anyDuplicated(x), message)
expect_set_equal <- function(actual, expected, message = "Sets differ") expect_true(setequal(actual, expected), paste0(message, "; missing=", paste(setdiff(expected, actual), collapse = ","), "; extra=", paste(setdiff(actual, expected), collapse = ",")))
record_test <- function(id, expr, required = TRUE) {
  result <- tryCatch({
    withCallingHandlers(force(expr), warning = function(w) { .regression_warnings <<- unique(c(.regression_warnings, conditionMessage(w))); invokeRestart("muffleWarning") })
    data.frame(ID = id, Status = "PASS", Message = "", stringsAsFactors = FALSE)
  }, error = function(e) data.frame(ID = id, Status = if (required) "FAIL" else "SKIP", Message = conditionMessage(e), stringsAsFactors = FALSE))
  .regression_results <<- rbind(.regression_results, result)
  invisible(result)
}
require_namespace_or_skip <- function(pkg, id, expr = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    record_test(id, stop("Optional package unavailable: ", pkg, call. = FALSE), required = FALSE)
    return(FALSE)
  }
  if (is.null(expr)) return(TRUE)
  record_test(id, expr)
  TRUE
}
