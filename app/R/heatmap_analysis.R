# Expression and feature heatmap analysis
heatmap_scheme_colors <- function(scheme = "blue_white_red", n = 100) {
  anchors <- switch(scheme,
    blue_white_red = c("#2166AC", "white", "#B2182B"),
    red_white_blue = c("#B2182B", "white", "#2166AC"),
    purple_white_orange = c("#5E3C99", "white", "#E66101"),
    green_black_red = c("#1B9E77", "#000000", "#D73027"),
    stop("Unsupported heatmap color scheme.")
  )
  grDevices::colorRampPalette(anchors)(n)
}

validate_annotation_colors <- function(colors, labels, annotation_name) {
  if (length(labels) == 0) return(character())
  colors <- as.character(colors)
  if (is.null(names(colors))) names(colors) <- labels
  colors <- colors[labels]
  if (any(is.na(colors)) || any(!grepl("^#[0-9A-Fa-f]{6}$", colors))) {
    stop(annotation_name, " colors must be valid #RRGGBB values for every displayed label.")
  }
  colors
}

plot_expression_heatmap <- function(mat, group_info, out_pdf, out_csv, top_n = 100, row_cluster = "hclust", col_cluster = "hclust", row_kmeans_k = 4, col_kmeans_k = 4, row_hclust_k = 2, col_hclust_k = 2, color_scheme = "blue_white_red", group_annotation_colors = NULL, cluster_annotation_colors = NULL, row_labels = NULL, width = 4.5, height = 4.5) {
  used <- preprocess_expr(mat, TRUE, 0.5)
  if (nrow(used) < 2 || ncol(used) < 2) stop("Expression heatmap requires at least two proteins and two samples after filtering.")
  vars <- apply(used, 1, var, na.rm = TRUE)
  selected <- rownames(used)[order(vars, decreasing = TRUE)[seq_len(min(top_n, length(vars)))]]
  used <- used[rownames(used) %in% selected, , drop = FALSE]
  scaled <- t(scale(t(used)))
  scaled[!is.finite(scaled)] <- 0
  group_labels <- as.character(group_info$Group[match(colnames(scaled), group_info$Sample)])
  if (anyNA(group_labels)) stop("Some heatmap samples do not have group labels.")
  group_levels <- levels(group_info$Group)
  group_colors <- validate_annotation_colors(group_annotation_colors %||% stats::setNames(sci_palette(length(group_levels), "npg"), group_levels), group_levels, "Group annotation")

  safe_k <- function(k, max_k, label) {
    k <- suppressWarnings(as.integer(k))
    if (!is.finite(k) || k < 1) stop(label, " must be at least 1.")
    min(k, max_k)
  }
  cluster_gaps <- function(cluster) {
    if (length(cluster) < 2) return(NULL)
    gaps <- which(cluster[-1] != cluster[-length(cluster)])
    if (length(gaps) == 0) NULL else gaps
  }
  set.seed(123)
  row_km <- if (row_cluster == "kmeans") kmeans(scaled, centers = safe_k(row_kmeans_k, nrow(scaled), "Row K-means k"))$cluster else rep(NA_integer_, nrow(scaled))
  set.seed(123)
  col_km <- if (col_cluster == "kmeans") kmeans(t(scaled), centers = safe_k(col_kmeans_k, ncol(scaled), "Column K-means k"))$cluster else rep(NA_integer_, ncol(scaled))
  row_order <- if (row_cluster == "kmeans") rownames(scaled)[order(row_km, seq_along(row_km))] else rownames(scaled)
  col_order <- if (col_cluster == "kmeans") colnames(scaled)[order(col_km, seq_along(col_km))] else colnames(scaled)
  plot_scaled <- scaled[row_order, col_order, drop = FALSE]
  plot_row_labels <- NULL
  if (!is.null(row_labels)) {
    if (is.null(names(row_labels))) stop("Heatmap row labels must be named by the current protein identifiers.")
    plot_row_labels <- as.character(row_labels[match(rownames(plot_scaled), names(row_labels))])
    plot_row_labels[is.na(plot_row_labels)] <- ""
  }
  row_kmeans_gaps <- if (row_cluster == "kmeans") cluster_gaps(row_km[match(rownames(plot_scaled), rownames(scaled))]) else NULL
  col_kmeans_gaps <- if (col_cluster == "kmeans") cluster_gaps(col_km[match(colnames(plot_scaled), colnames(scaled))]) else NULL

  row_hc <- if (row_cluster == "hclust") stats::hclust(stats::dist(plot_scaled), method = "complete") else FALSE
  col_hc <- if (col_cluster == "hclust") stats::hclust(stats::dist(t(plot_scaled)), method = "complete") else FALSE
  row_hc_k <- if (row_cluster == "hclust") safe_k(row_hclust_k, nrow(plot_scaled), "Row hierarchical clustering k") else NA_integer_
  col_hc_k <- if (col_cluster == "hclust") safe_k(col_hclust_k, ncol(plot_scaled), "Column hierarchical clustering k") else NA_integer_
  row_hc_cluster <- if (row_cluster == "hclust") stats::cutree(row_hc, k = row_hc_k) else rep(NA_integer_, nrow(plot_scaled))
  col_hc_cluster <- if (col_cluster == "hclust") stats::cutree(col_hc, k = col_hc_k) else rep(NA_integer_, ncol(plot_scaled))
  col_cluster_values <- if (col_cluster == "hclust") col_hc_cluster else if (col_cluster == "kmeans") col_km[match(colnames(plot_scaled), colnames(scaled))] else rep(NA_integer_, ncol(plot_scaled))
  cluster_levels <- if (col_cluster == "none") character() else paste0("Cluster ", sort(unique(col_cluster_values)))
  ann <- data.frame(Group = factor(group_labels[match(colnames(plot_scaled), colnames(scaled))], levels = group_levels), stringsAsFactors = FALSE)
  if (length(cluster_levels) > 0) ann$Cluster <- factor(paste0("Cluster ", col_cluster_values), levels = cluster_levels)
  rownames(ann) <- colnames(plot_scaled)
  annotation_colors <- list(Group = group_colors)
  if (length(cluster_levels) > 0) {
    cluster_colors <- validate_annotation_colors(cluster_annotation_colors %||% stats::setNames(sci_palette(length(cluster_levels), "lancet"), cluster_levels), cluster_levels, "Cluster annotation")
    annotation_colors$Cluster <- cluster_colors
  }
  # Reserve vertical space for column annotations and plot margins, then keep each
  # visible label below the available per-row height to avoid overlap.
  row_label_fontsize <- if (is.null(plot_row_labels)) {
    10
  } else {
    available_row_height_pt <- max(1, height * 72 - 115) / nrow(plot_scaled)
    max(1.5, min(10, 0.75 * available_row_height_pt))
  }
  grDevices::pdf(out_pdf, width = width, height = height)
  ph <- pheatmap::pheatmap(plot_scaled, scale = "none", color = heatmap_scheme_colors(color_scheme), cluster_rows = row_hc, cluster_cols = col_hc, cutree_rows = row_hc_k, cutree_cols = col_hc_k, gaps_row = row_kmeans_gaps, gaps_col = col_kmeans_gaps, kmeans_k = NA, annotation_col = ann, annotation_colors = annotation_colors, labels_row = plot_row_labels, show_rownames = !is.null(plot_row_labels), fontsize_row = row_label_fontsize, border_color = NA)
  grDevices::dev.off()
  row_order <- if (row_cluster == "hclust" && !is.null(ph$tree_row)) {
    rownames(plot_scaled)[ph$tree_row$order]
  } else {
    row_order
  }
  col_order <- if (col_cluster == "hclust" && !is.null(ph$tree_col)) colnames(plot_scaled)[ph$tree_col$order] else col_order
  ordered_scaled <- scaled[row_order, col_order, drop = FALSE]
  row_hc_ordered <- if (row_cluster == "hclust") row_hc_cluster[match(row_order, rownames(plot_scaled))] else rep(NA_integer_, length(row_order))
  col_hc_ordered <- if (col_cluster == "hclust") col_hc_cluster[match(col_order, colnames(plot_scaled))] else rep(NA_integer_, length(col_order))
  data.table::fwrite(data.frame(ProteinID = rownames(ordered_scaled), KmeansCluster = row_km[match(rownames(ordered_scaled), rownames(scaled))], HierarchicalCluster = row_hc_ordered, ordered_scaled, check.names = FALSE), out_csv)
  data.table::fwrite(data.frame(RowOrder = seq_along(row_order), ProteinID = row_order, KmeansCluster = row_km[match(row_order, rownames(scaled))], HierarchicalCluster = row_hc_ordered), sub("\\.csv$", "_row_order.csv", out_csv))
  data.table::fwrite(data.frame(ColOrder = seq_along(col_order), Sample = col_order, Group = as.character(group_info$Group[match(col_order, group_info$Sample)]), KmeansCluster = col_km[match(col_order, colnames(scaled))], HierarchicalCluster = col_hc_ordered), sub("\\.csv$", "_col_order.csv", out_csv))
}
