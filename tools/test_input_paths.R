resolve_external_test_file <- function(env_name) {
  path <- Sys.getenv(env_name, unset = "")
  if (!nzchar(path)) {
    stop("Set environment variable ", env_name, " to the local regression input file.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Regression input file from ", env_name, " does not exist.", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
