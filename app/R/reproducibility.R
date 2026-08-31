# Runtime metadata helpers used by analysis_manifest.json.
# Optional packages are recorded as not installed rather than failing analysis.
runtime_package_metadata <- function() {
  packages <- unique(unlist(dependency_manifest[c("cran_runtime", "bioconductor_runtime", "r_base_runtime")], use.names = FALSE))
  out <- lapply(packages, function(pkg) {
    installed <- requireNamespace(pkg, quietly = TRUE)
    list(version = if (installed) as.character(utils::packageVersion(pkg)) else NULL,
         source = if (pkg %in% dependency_manifest$bioconductor_runtime) "Bioconductor" else if (pkg %in% dependency_manifest$cran_runtime) "CRAN" else "base/recommended R",
         status = if (installed) "installed" else "not installed")
  })
  names(out) <- packages
  out
}
