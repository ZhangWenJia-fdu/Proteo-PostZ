# Shared expression preprocessing; behavior retained from V2.0.3
preprocess_expr <- function(mat, log2_transform = TRUE, min_valid_fraction = 0.5) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  if (log2_transform) mat <- log2(mat + 1)
  keep_n <- max(1, ceiling(ncol(mat) * min_valid_fraction))
  mat <- mat[rowSums(!is.na(mat)) >= keep_n, , drop = FALSE]
  imp <- t(apply(mat, 1, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))
  imp[complete.cases(imp), , drop = FALSE]
}
