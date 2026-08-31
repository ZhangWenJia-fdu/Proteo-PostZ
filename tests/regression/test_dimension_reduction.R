test_dimension_reduction <- function(ctx) {
  f <- make_binary_ml_fixture(seed = 606, n_features = 120, n_per_group = 12)
  used <- preprocess_expr(f$quantity, TRUE, 0.5)
  record_test("pca_basic", { r <- run_pca_core(t(used)); expect_equal(nrow(r$coordinates), 24, "PCA sample count differs"); expect_true(all(is.finite(as.matrix(r$coordinates[, c("PC1", "PC2")]))), "PCA coordinates are not finite"); expect_true(all(r$variance_percent >= 0 & r$variance_percent <= 100), "PCA variance is invalid") })
  record_test("pca_constant_feature", { x <- cbind(zero = rep(0, 6), constant = rep(5, 6), v1 = 1:6, v2 = c(2, 4, 1, 5, 3, 6)); rownames(x) <- paste0("S", 1:6); r <- run_pca_core(x); expect_set_equal(colnames(r$sample_matrix), c("v1", "v2"), "PCA constant-feature filter differs") })
  record_test("pca_constant_feature_failure", { x <- matrix(0, nrow = 6, ncol = 3); rownames(x) <- paste0("S", 1:6); msg <- tryCatch({ run_pca_core(x); "" }, error = function(e) conditionMessage(e)); expect_true(grepl("at least two non-constant features", msg, fixed = TRUE), "PCA constant-feature failure message is unclear") })
  record_test("umap_basic", { r <- run_umap_core(t(used), 10, 0.1, seed = 123); expect_equal(dim(r$coordinates), c(24, 3), "UMAP coordinate dimensions differ"); expect_equal(r$n_neighbors, 10, "UMAP neighbor parameter differs") })
  record_test("tsne_basic", { r <- run_tsne_core(t(used), resolve_tsne_perplexity("Auto", ncol(used)), 250, seed = 123); expect_equal(nrow(r$coordinates), 24, "t-SNE sample count differs"); expect_equal(r$perplexity, 7, "t-SNE Auto perplexity differs") })
}
