# Slingshot trajectory and pseudotime analysis
run_slingshot_pseudotime <- function(mat, group_info, outdir, reduction = c("PCA", "UMAP"), start_group = NULL, end_group = NULL, width = 3.3, height = 3.3, palette = "npg", n_neighbors = 10, min_dist = 0.1, seed = 123, top_n = 50, heatmap_width = 4.17, heatmap_height = 5.56) {
  reduction <- match.arg(reduction)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  missing <- c(
    if (!requireNamespace("uwot", quietly = TRUE)) "uwot" else character(),
    if (!requireNamespace("slingshot", quietly = TRUE)) "slingshot" else character(),
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) "SingleCellExperiment" else character(),
    if (!requireNamespace("pheatmap", quietly = TRUE)) "pheatmap" else character()
  )
  if (length(missing) > 0) stop("Slingshot pseudotime requires missing package(s): ", paste(missing, collapse = ", "), ".")
  used <- preprocess_expr(mat, TRUE, 0.5)
  if (nrow(used) < 3) stop("Slingshot requires at least three proteins after preprocessing.")
  sample_mat <- t(used)
  sample_groups <- group_info$Group[match(rownames(sample_mat), group_info$Sample)]
  if (any(is.na(sample_groups))) stop("Some samples do not have group labels.")
  if (length(unique(sample_groups)) < 2) stop("Slingshot requires at least two groups/clusters.")
  if (nrow(sample_mat) < 4) stop("Slingshot requires at least four samples for a meaningful trajectory.")
  if (!is.null(start_group) && nzchar(start_group) && !start_group %in% as.character(sample_groups)) stop("Start group is not present in the current group labels: ", start_group)
  if (!is.null(end_group) && nzchar(end_group) && end_group != "None" && !end_group %in% as.character(sample_groups)) stop("End group is not present in the current group labels: ", end_group)
  pca_result <- run_pca_core(sample_mat)
  sample_mat <- pca_result$sample_matrix
  used <- t(sample_mat)
  pca <- pca_result$model
  pca_coords <- pca$x[, 1:2, drop = FALSE]
  colnames(pca_coords) <- c("PC1", "PC2")
  set.seed(seed)
  nn <- min(n_neighbors, max(2, nrow(sample_mat) - 1))
  umap_coords <- uwot::umap(pca_coords, n_neighbors = nn, min_dist = min_dist, metric = "euclidean", verbose = FALSE)
  rownames(umap_coords) <- rownames(sample_mat)
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  pca_df <- data.frame(Sample = rownames(pca_coords), Group = as.character(sample_groups), PC1 = pca_coords[, 1], PC2 = pca_coords[, 2], check.names = FALSE)
  umap_df <- data.frame(Sample = rownames(umap_coords), Group = as.character(sample_groups), UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2], check.names = FALSE)
  data.table::fwrite(pca_df, file.path(outdir, "slingshot_input_PCA_coordinates.csv"))
  data.table::fwrite(umap_df, file.path(outdir, "slingshot_input_UMAP_coordinates.csv"))
  coords <- if (reduction == "PCA") pca_coords else umap_coords
  colnames(coords) <- c("Dim1", "Dim2")
  start_group <- if (is.null(start_group) || !nzchar(start_group)) as.character(sample_groups[1]) else start_group
  end_group <- if (is.null(end_group) || !nzchar(end_group) || end_group == "None") NULL else end_group
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = used),
    colData = data.frame(Sample = rownames(sample_mat), Group = factor(sample_groups, levels = levels(factor(sample_groups))), Cluster = factor(sample_groups), row.names = rownames(sample_mat))
  )
  SingleCellExperiment::reducedDims(sce)$PCA <- pca_coords
  SingleCellExperiment::reducedDims(sce)$UMAP <- umap_coords
  sce <- tryCatch(
    slingshot::slingshot(sce, clusterLabels = "Cluster", reducedDim = reduction, start.clus = start_group, end.clus = end_group),
    error = function(e) {
      stop("Slingshot could not infer a trajectory from the current group labels and ", reduction, " coordinates. Check that stages form a clear progression, or try the other reduction. Original error: ", conditionMessage(e), call. = FALSE)
    }
  )
  pst <- slingshot::slingPseudotime(sce)
  pst_df <- data.frame(Sample = colnames(sce), Group = as.character(sce$Group), Cluster = as.character(sce$Cluster), pst, check.names = FALSE)
  data.table::fwrite(pst_df, file.path(outdir, "slingshot_pseudotime_by_sample.csv"))
  data.table::fwrite(pst_df, file.path(outdir, "slingshot_sample_pseudotime.csv"))
  data.table::fwrite(data.frame(Group = levels(factor(sample_groups))), file.path(outdir, "slingshot_group_clusters.csv"))
  capture.output(slingshot::slingLineages(sce), file = file.path(outdir, "slingshot_lineages.txt"))
  cols <- sci_palette(length(levels(factor(sample_groups))), palette)
  names(cols) <- levels(factor(sample_groups))
  plot_coords <- if (reduction == "PCA") pca_coords else umap_coords
  plot_xlab <- if (reduction == "PCA") "PC1" else "UMAP1"
  plot_ylab <- if (reduction == "PCA") "PC2" else "UMAP2"
  grDevices::pdf(file.path(outdir, "slingshot_pseudotime_trajectory.pdf"), width = width, height = height)
  graphics::plot(plot_coords[, 1], plot_coords[, 2], col = cols[as.character(sample_groups)], pch = 16, xlab = plot_xlab, ylab = plot_ylab, main = paste("Slingshot trajectory on", reduction))
  graphics::text(plot_coords[, 1], plot_coords[, 2], labels = rownames(plot_coords), pos = 3, cex = 0.55)
  graphics::legend("topright", legend = names(cols), col = cols, pch = 16, cex = 0.75, bty = "n")
  graphics::lines(slingshot::SlingshotDataSet(sce), lwd = 2, col = "black")
  grDevices::dev.off()
  p_stage <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = Group)) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(values = cols) +
    theme_sci() +
    ggplot2::labs(title = "UMAP by stage")
  ggplot2::ggsave(file.path(outdir, "UMAP_by_stage.pdf"), p_stage, width = width, height = height)
  pt_col <- colnames(pst)[1]
  umap_pt_df <- dplyr::left_join(umap_df, pst_df[, c("Sample", pt_col), drop = FALSE], by = "Sample")
  colnames(umap_pt_df)[colnames(umap_pt_df) == pt_col] <- "Pseudotime"
  p_pt <- ggplot2::ggplot(umap_pt_df, ggplot2::aes(UMAP1, UMAP2, color = Pseudotime)) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_gradient(low = "#3B4CC0", high = "#B40426", na.value = "grey80") +
    theme_sci() +
    ggplot2::labs(title = "UMAP by pseudotime")
  ggplot2::ggsave(file.path(outdir, "UMAP_by_pseudotime.pdf"), p_pt, width = width, height = height)
  pt_vec <- umap_pt_df$Pseudotime
  names(pt_vec) <- umap_pt_df$Sample
  valid_samples <- names(pt_vec)[!is.na(pt_vec)]
  if (length(valid_samples) < 3) stop("Slingshot produced too few samples with finite pseudotime for protein association analysis.")
  expr_for_cor <- used[, valid_samples, drop = FALSE]
  pt_vec <- pt_vec[valid_samples]
  cor_res <- data.frame(ProteinID = rownames(expr_for_cor), SpearmanR = NA_real_, Pvalue = NA_real_)
  for (i in seq_len(nrow(expr_for_cor))) {
    x <- expr_for_cor[i, ]
    if (length(unique(x)) > 2 && length(unique(pt_vec)) > 2) {
      ct <- suppressWarnings(stats::cor.test(x, pt_vec, method = "spearman"))
      cor_res$SpearmanR[i] <- unname(ct$estimate)
      cor_res$Pvalue[i] <- ct$p.value
    }
  }
  cor_res$FDR <- p.adjust(cor_res$Pvalue, method = "BH")
  cor_res <- cor_res[order(cor_res$FDR, -abs(cor_res$SpearmanR)), , drop = FALSE]
  data.table::fwrite(cor_res, file.path(outdir, "pseudotime_associated_proteins_spearman.csv"))
  top_proteins <- head(cor_res$ProteinID[!is.na(cor_res$FDR)], top_n)
  if (length(top_proteins) >= 2) {
    ordered_samples <- names(sort(pt_vec))
    heat_mat <- expr_for_cor[top_proteins, ordered_samples, drop = FALSE]
    heat_z <- t(scale(t(heat_mat)))
    heat_z[!is.finite(heat_z)] <- 0
    ann_col <- data.frame(Group = umap_pt_df$Group[match(colnames(heat_z), umap_pt_df$Sample)], Pseudotime = pt_vec[colnames(heat_z)])
    rownames(ann_col) <- colnames(heat_z)
    grDevices::pdf(file.path(outdir, "Top_pseudotime_proteins_heatmap.pdf"), width = heatmap_width, height = heatmap_height)
    pheatmap::pheatmap(heat_z, color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100), breaks = seq(-2, 2, length.out = 101), cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = ann_col, show_colnames = FALSE, fontsize_row = 7, border_color = NA, main = paste0("Top ", length(top_proteins), " pseudotime proteins"))
    grDevices::dev.off()
  }
  invisible(pst_df)
}
