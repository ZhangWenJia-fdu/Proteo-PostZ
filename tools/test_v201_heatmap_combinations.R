args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE) else normalizePath("tools/test_v201_heatmap_combinations.R", winslash = "/", mustWork = TRUE)
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(file.path(package_root, "app"))
source("R/analysis_core.R", encoding = "UTF-8")

out <- file.path(package_root, "outputs", "v201_heatmap_combinations")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
quantity <- matrix(c(1, 2, 4, 8, 2, 3, 5, 9, 3, 5, 6, 10, 4, 7, 8, 12), nrow = 4, byrow = TRUE)
rownames(quantity) <- paste0("P", seq_len(nrow(quantity)))
colnames(quantity) <- paste0("S", seq_len(ncol(quantity)))
groups <- make_group_info(colnames(quantity), c("A", "A", "B", "B"))
modes <- c("hclust", "kmeans", "none")
run_heatmap_set <- function(mat, prefix) {
  for (row_method in modes) {
    for (col_method in modes) {
      id <- paste(prefix, row_method, col_method, sep = "_")
      plot_expression_heatmap(
        mat,
        groups,
        file.path(out, paste0(id, ".pdf")),
        file.path(out, paste0(id, ".csv")),
        top_n = nrow(mat),
        row_cluster = row_method,
        col_cluster = col_method,
        row_kmeans_k = 2,
        col_kmeans_k = 2,
        width = 3,
        height = 3
      )
    }
  }
}
run_heatmap_set(quantity, "expression")
run_heatmap_set(quantity[1:3, , drop = FALSE], "feature")
expected <- length(modes)^2
expression_actual <- length(list.files(out, pattern = "^expression_.*\\.pdf$"))
feature_actual <- length(list.files(out, pattern = "^feature_.*\\.pdf$"))
if (expression_actual != expected || feature_actual != expected) stop("Heatmap combination test did not generate all expected PDFs")
cat("HEATMAP_COMBINATIONS_OK expression=", expression_actual, " feature=", feature_actual, "\n", sep = "")
