# Pure numerical helpers for sample-level dimension reduction.
# Preprocessing, grouping, plotting, and file output remain in app.R.

prepare_pca_input <- function(sample_matrix) {
  sample_matrix <- as.matrix(sample_matrix)
  feature_variance <- apply(sample_matrix, 2, stats::var, na.rm = TRUE)
  keep <- is.finite(feature_variance) & feature_variance > 0
  if (sum(keep) < 2) stop("PCA requires at least two non-constant features after preprocessing.")
  sample_matrix[, keep, drop = FALSE]
}

run_pca_core <- function(sample_matrix) {
  sample_matrix <- prepare_pca_input(sample_matrix)
  pca <- stats::prcomp(sample_matrix, center = TRUE, scale. = TRUE)
  variance_percent <- summary(pca)$importance[2, 1:2] * 100
  coordinates <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  list(
    model = pca,
    coordinates = coordinates,
    variance_percent = variance_percent,
    sample_matrix = sample_matrix
  )
}

run_umap_core <- function(sample_matrix, n_neighbors, min_dist, seed = 123) {
  set.seed(seed)
  nn <- min(n_neighbors, max(2, nrow(sample_matrix) - 1))
  umap_result <- uwot::umap(
    sample_matrix,
    n_neighbors = nn,
    min_dist = min_dist,
    metric = "euclidean",
    verbose = FALSE
  )
  coordinates <- data.frame(
    Sample = rownames(sample_matrix),
    UMAP1 = umap_result[, 1],
    UMAP2 = umap_result[, 2],
    stringsAsFactors = FALSE
  )
  list(
    model = umap_result,
    coordinates = coordinates,
    n_neighbors = nn
  )
}

resolve_tsne_perplexity <- function(value, n_samples) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])) {
    requested <- NA_integer_
  } else {
    value_text <- trimws(as.character(value[[1]]))
    if (!nzchar(value_text) || tolower(value_text) == "auto") {
      requested <- NA_integer_
    } else {
      requested <- suppressWarnings(as.integer(value_text))
      if (is.na(requested) || requested < 1) {
        stop("Expected Auto or a positive integer, got: ", value_text)
      }
    }
  }

  max_perplexity <- max(1, floor((n_samples - 1) / 3))
  perplexity <- if (is.na(requested)) {
    min(30, max_perplexity)
  } else {
    min(requested, max_perplexity)
  }
  if (perplexity < 1) {
    stop("t-SNE perplexity must be at least 1 and less than one-third of the sample count.")
  }
  perplexity
}

run_tsne_core <- function(sample_matrix, perplexity, max_iter, seed = 123) {
  set.seed(seed)
  tsne_result <- Rtsne::Rtsne(
    sample_matrix,
    dims = 2,
    perplexity = perplexity,
    max_iter = max_iter,
    check_duplicates = FALSE,
    pca = TRUE,
    verbose = FALSE
  )
  coordinates <- data.frame(
    Sample = rownames(sample_matrix),
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2],
    stringsAsFactors = FALSE
  )
  list(
    model = tsne_result,
    coordinates = coordinates,
    perplexity = perplexity,
    max_iter = max_iter
  )
}
