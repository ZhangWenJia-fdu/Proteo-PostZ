# Standard quantitative analysis
plot_correlation_heatmap <- function(mat, group_info, out_pdf, out_csv, method = "pearson", order_mode = "original", cluster_distance = "one_minus_correlation", cluster_linkage = "complete", cluster_k = NULL, digits = 2, fontsize_number = 8, color_scheme = "blue_white_red", min_cor = -1, max_cor = 1, width = 3.3, height = 3.3) {
  used <- preprocess_expr(mat, TRUE, 0.3)
  cor_mat <- cor(used, method = method, use = "pairwise.complete.obs")
  if (!is.finite(min_cor) || !is.finite(max_cor) || min_cor >= max_cor) {
    stop("Correlation legend min must be smaller than max.")
  }
  if (ncol(cor_mat) < 2) stop("Correlation heatmap requires at least two samples.")

  block_gaps <- function(labels) {
    if (length(labels) < 2) return(NULL)
    gaps <- which(labels[-1] != labels[-length(labels)])
    if (length(gaps) == 0) NULL else gaps
  }

  plot_mat <- cor_mat
  plot_cluster_rows <- FALSE
  plot_cluster_cols <- FALSE
  plot_gaps <- NULL
  if (order_mode == "group") {
    ord <- unlist(lapply(levels(group_info$Group), function(g) {
      smp <- group_info$Sample[group_info$Group == g]
      intersect(smp, colnames(cor_mat))
    }))
    cor_mat <- cor_mat[ord, ord, drop = FALSE]
    plot_mat <- cor_mat
    plot_gaps <- block_gaps(as.character(group_info$Group[match(colnames(cor_mat), group_info$Sample)]))
  } else if (order_mode == "hclust") {
    if (is.null(cluster_k)) cluster_k <- length(levels(group_info$Group))
    cluster_k <- as.integer(cluster_k)
    if (!is.finite(cluster_k) || cluster_k < 1 || cluster_k > ncol(cor_mat)) {
      stop("Hierarchical clustering k must be between 1 and the number of samples.")
    }
    if (cluster_distance == "one_minus_correlation") {
      distance_matrix <- 1 - cor_mat
      diag(distance_matrix) <- 0
      if (any(!is.finite(distance_matrix[upper.tri(distance_matrix)]))) {
        stop("Hierarchical clustering requires finite pairwise correlations. Check for samples with insufficient valid quantitative values.")
      }
      distance_object <- as.dist(pmax(distance_matrix, 0))
    } else if (cluster_distance == "correlation_euclidean") {
      if (any(!is.finite(cor_mat))) {
        stop("Hierarchical clustering requires finite correlations. Check for samples with insufficient valid quantitative values.")
      }
      distance_object <- stats::dist(cor_mat, method = "euclidean")
    } else {
      stop("Unsupported hierarchical clustering distance.")
    }
    hc <- stats::hclust(distance_object, method = cluster_linkage)
    cor_mat <- cor_mat[hc$order, hc$order, drop = FALSE]
    plot_cluster_rows <- hc
    plot_cluster_cols <- hc
    plot_gaps <- cluster_k
  }
  sample_group <- as.character(group_info$Group[match(colnames(cor_mat), group_info$Sample)])
  data.table::fwrite(data.frame(Sample = rownames(cor_mat), Group = sample_group, cor_mat, check.names = FALSE), out_csv)
  long_csv <- sub("\\.csv$", "_long.csv", out_csv)
  long_df <- as.data.frame(cor_mat) |>
    tibble::rownames_to_column(".cor_row_sample") |>
    tidyr::pivot_longer(-.cor_row_sample, names_to = ".cor_col_sample", values_to = "Correlation") |>
    dplyr::mutate(
      Group1 = as.character(group_info$Group[match(.cor_row_sample, group_info$Sample)]),
      Group2 = as.character(group_info$Group[match(.cor_col_sample, group_info$Sample)]),
      Method = method,
      OrderedSample1 = match(.cor_row_sample, rownames(cor_mat)),
      OrderedSample2 = match(.cor_col_sample, colnames(cor_mat))
    ) |>
    dplyr::transmute(Method, Sample1 = .cor_row_sample, Group1, OrderedSample1, Sample2 = .cor_col_sample, Group2, OrderedSample2, Correlation)
  data.table::fwrite(long_df, long_csv)
  ann <- data.frame(Group = group_info$Group[match(colnames(plot_mat), group_info$Sample)])
  rownames(ann) <- colnames(plot_mat)
  cols <- if (color_scheme == "purple_white_orange") colorRampPalette(c("#5E3C99", "white", "#E66101"))(100) else colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  breaks <- seq(min_cor, max_cor, length.out = 101)
  nums <- matrix(sprintf(paste0("%.", digits, "f"), plot_mat), nrow(plot_mat), dimnames = dimnames(plot_mat))
  grDevices::pdf(out_pdf, width = width, height = height)
  pheatmap::pheatmap(plot_mat, color = cols, breaks = breaks, cluster_rows = plot_cluster_rows, cluster_cols = plot_cluster_cols, cutree_rows = if (order_mode == "hclust") plot_gaps else NA, cutree_cols = if (order_mode == "hclust") plot_gaps else NA, gaps_row = if (order_mode == "group") plot_gaps else NULL, gaps_col = if (order_mode == "group") plot_gaps else NULL, annotation_col = ann, annotation_row = ann, display_numbers = nums, fontsize_number = fontsize_number, border_color = NA)
  grDevices::dev.off()
  cor_mat
}
plot_rank_abundance <- function(mat, group_info, out_pdf, out_csv, width = 3.3, height = 3.3, palette = "npg", transform = c("log2", "log10")) {
  transform <- match.arg(transform)
  transform_label <- paste0(transform, "(abundance + 1)")
  df <- as.data.frame(mat) |>
    tibble::rownames_to_column("ProteinID") |>
    tidyr::pivot_longer(-ProteinID, names_to = "Sample", values_to = "Raw_Abundance") |>
    dplyr::filter(!is.na(Raw_Abundance), is.finite(Raw_Abundance)) |>
    dplyr::group_by(Sample) |>
    dplyr::arrange(dplyr::desc(Raw_Abundance), .by_group = TRUE) |>
    dplyr::mutate(Rank = dplyr::row_number(), Transformed_Abundance = if (transform == "log2") log2(Raw_Abundance + 1) else log10(Raw_Abundance + 1), Transform = transform_label) |>
    dplyr::ungroup() |>
    dplyr::left_join(group_info, by = "Sample")
  data.table::fwrite(df, out_csv)
  p <- ggplot2::ggplot(df, ggplot2::aes(Rank, Transformed_Abundance, color = Group, group = Sample)) + ggplot2::geom_line(alpha = 0.55, linewidth = 0.35) + theme_sci() + ggplot2::scale_color_manual(values = sci_palette(length(levels(group_info$Group)), palette)) + ggplot2::labs(y = transform_label)
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
  invisible(df)
}
plot_cv_ridges <- function(mat, group_info, out_pdf, out_csv, cv_max = 60, width = 3.3, height = 3.3, palette = "npg") {
  rows <- lapply(levels(group_info$Group), function(g) {
    smp <- intersect(group_info$Sample[group_info$Group == g], colnames(mat))
    if (length(smp) < 2) return(NULL)
    vals <- mat[, smp, drop = FALSE]
    cv <- apply(vals, 1, function(x) if (sum(!is.na(x)) >= 2 && mean(x, na.rm = TRUE) > 0) 100 * sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) else NA_real_)
    data.frame(ProteinID = rownames(mat), Group = g, CV = cv)
  })
  df <- dplyr::bind_rows(rows) |> dplyr::filter(!is.na(CV), CV >= 0, CV <= cv_max)
  med <- df |> dplyr::group_by(Group) |> dplyr::summarise(MedianCV = median(CV), .groups = "drop")
  data.table::fwrite(df, out_csv)
  data.table::fwrite(med, sub("\\.csv$", "_median.csv", out_csv))
  p <- ggplot2::ggplot(df, ggplot2::aes(CV, Group, fill = Group)) + ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) + ggplot2::geom_segment(data = med, ggplot2::aes(x = MedianCV, xend = MedianCV, y = as.numeric(factor(Group)) - 0.35, yend = as.numeric(factor(Group)) + 0.35), inherit.aes = FALSE, linetype = "dashed") + ggplot2::coord_cartesian(xlim = c(0, cv_max)) + ggplot2::scale_fill_manual(values = sci_palette(length(levels(group_info$Group)), palette)) + theme_sci() + ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
}
