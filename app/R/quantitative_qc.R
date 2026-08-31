# ProteoPostZ V2.1.0 quantitative QC
# This module intentionally reads the importer-produced quantitative matrix
# directly. It must not call preprocess_expr() or modify the input matrix.

prepare_quantitative_qc <- function(mat, group_info) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  if (nrow(mat) < 1L || ncol(mat) < 1L) stop("Quantitative QC requires at least one protein and one sample.")
  if (is.null(rownames(mat)) || any(!nzchar(trimws(rownames(mat))))) stop("Quantitative QC requires protein row names.")
  if (is.null(colnames(mat)) || any(!nzchar(trimws(colnames(mat))))) stop("Quantitative QC requires sample column names.")
  if (!all(c("Sample", "Group") %in% colnames(group_info))) stop("group_info must contain Sample and Group columns.")

  sample_order <- colnames(mat)
  groups <- as.character(group_info$Group[match(sample_order, group_info$Sample)])
  groups[is.na(groups) | !nzchar(trimws(groups))] <- "Unassigned"
  missing_matrix <- is.na(mat)
  log2_mat <- suppressWarnings(log2(mat + 1))
  log2_mat[!is.finite(log2_mat)] <- NA_real_
  observed_median <- function(x) {
    if (!any(is.finite(x))) NA_real_ else stats::median(x[is.finite(x)])
  }
  observed_mean <- function(x) {
    if (!any(is.finite(x))) NA_real_ else mean(x[is.finite(x)])
  }

  sample_summary <- data.frame(
    Sample = sample_order, Group = groups,
    Total_Proteins = nrow(mat),
    Quantified_Proteins = colSums(!missing_matrix),
    Missing_Values = colSums(missing_matrix),
    Missing_Percent = colMeans(missing_matrix) * 100,
    Median_Log2_Abundance = vapply(seq_len(ncol(log2_mat)), function(i) observed_median(log2_mat[, i]), numeric(1)),
    Mean_Log2_Abundance = vapply(seq_len(ncol(log2_mat)), function(i) observed_mean(log2_mat[, i]), numeric(1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  protein_missing <- rowSums(missing_matrix)
  protein_summary <- data.frame(
    ProteinID = rownames(mat), Total_Samples = ncol(mat),
    Quantified_Samples = rowSums(!missing_matrix), Missing_Samples = protein_missing,
    Missing_Percent = protein_missing / ncol(mat) * 100,
    Median_Log2_Abundance = vapply(seq_len(nrow(log2_mat)), function(i) observed_median(log2_mat[i, ]), numeric(1)),
    Mean_Log2_Abundance = vapply(seq_len(nrow(log2_mat)), function(i) observed_mean(log2_mat[i, ]), numeric(1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  group_values <- unique(groups)
  group_keys <- make.names(group_values, unique = TRUE)
  group_column_map <- data.frame(Group = group_values, Column_Key = group_keys, stringsAsFactors = FALSE)
  for (i in seq_along(group_values)) {
    idx <- groups == group_values[i]
    key <- group_keys[i]
    qn <- rowSums(!missing_matrix[, idx, drop = FALSE])
    protein_summary[[paste0(key, "_Quantified_N")]] <- qn
    protein_summary[[paste0(key, "_Total_N")]] <- sum(idx)
    protein_summary[[paste0(key, "_Quantified_Percent")]] <- qn / sum(idx) * 100
  }

  missing_levels <- seq.int(0L, ncol(mat))
  counts <- tabulate(factor(protein_missing, levels = missing_levels), nbins = length(missing_levels))
  distribution <- data.frame(
    Missing_Samples = missing_levels,
    Missing_Percent = missing_levels / ncol(mat) * 100,
    Protein_Count = as.integer(counts), stringsAsFactors = FALSE
  )
  distribution$Axis_Label <- paste0(distribution$Missing_Samples, "\n(", format(round(distribution$Missing_Percent, 1), trim = TRUE, scientific = FALSE), "%)")
  list(mat = mat, log2_mat = log2_mat, missing_matrix = missing_matrix,
       sample_summary = sample_summary, protein_summary = protein_summary,
       protein_missingness_distribution = distribution, group_column_map = group_column_map)
}

qc_group_colors <- function(groups, palette = "npg") {
  groups <- unique(as.character(groups))
  stats::setNames(sci_palette(length(groups), palette), groups)
}

qc_sample_axis_settings <- function(n_samples) {
  if (n_samples <= 10) return(list(size = 8, angle = 45))
  if (n_samples <= 20) return(list(size = 7, angle = 45))
  if (n_samples <= 40) return(list(size = 6, angle = 60))
  list(size = 5, angle = 60)
}

qc_text_theme <- function() {
  theme_sci() + ggplot2::theme(
    axis.title = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 8),
    legend.title = ggplot2::element_text(size = 8, face = "bold"),
    legend.text = ggplot2::element_text(size = 7),
    plot.title = ggplot2::element_text(size = 10, face = "bold")
  )
}

qc_missingness_axis_settings <- function(total_samples) {
  step <- if (total_samples <= 15) 1 else if (total_samples <= 30) 2 else if (total_samples <= 60) 5 else 10
  shown <- sort(unique(c(seq.int(0L, total_samples, by = step), total_samples)))
  size <- if (total_samples <= 15) 8 else if (total_samples <= 30) 7 else if (total_samples <= 60) 6 else 5
  list(shown = shown, size = size, angle = if (total_samples <= 15) 0 else 45)
}

qc_heatmap_font_settings <- function(n_samples, n_proteins) {
  list(
    fontsize_col = if (n_samples <= 15) 8 else if (n_samples <= 30) 7 else if (n_samples <= 50) 6 else 5,
    fontsize_row = if (n_proteins <= 50) 6 else if (n_proteins <= 100) 5 else 5,
    show_rownames = n_proteins <= 100
  )
}

plot_qc_sample_missingness <- function(qc, out_pdf, palette = "npg", width = 5, height = 4) {
  df <- qc$sample_summary
  df$Sample <- factor(df$Sample, levels = df$Sample)
  df$Group <- factor(df$Group, levels = unique(df$Group))
  axis <- qc_sample_axis_settings(nrow(df))
  p <- ggplot2::ggplot(df, ggplot2::aes(Sample, Missing_Percent, fill = Group)) +
    ggplot2::geom_col(width = 0.75) + ggplot2::scale_fill_manual(values = qc_group_colors(df$Group, palette)) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0, .03))) +
    ggplot2::labs(x = NULL, y = "Missing values (%)", fill = "Group") + qc_text_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = axis$size, angle = axis$angle, hjust = 1, vjust = 1))
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
  invisible(p)
}

plot_qc_protein_missingness_distribution <- function(qc, out_pdf, width = 5, height = 4) {
  df <- qc$protein_missingness_distribution
  df$Missing_Samples <- factor(df$Missing_Samples, levels = df$Missing_Samples)
  axis <- qc_missingness_axis_settings(ncol(qc$mat))
  labels <- ifelse(as.integer(levels(df$Missing_Samples)) %in% axis$shown, df$Axis_Label, "")
  p <- ggplot2::ggplot(df, ggplot2::aes(Missing_Samples, Protein_Count)) + ggplot2::geom_col(width = .75) +
    ggplot2::scale_x_discrete(labels = stats::setNames(labels, levels(df$Missing_Samples)), drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .05))) +
    ggplot2::labs(x = "Missing samples per protein\n(missing percentage)", y = "Number of proteins") + qc_text_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = axis$size, angle = axis$angle, hjust = 1, vjust = 1))
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
  invisible(p)
}

plot_qc_missingness_heatmap <- function(qc, out_pdf, heatmap_top_n = 200, palette = "npg", width = 6, height = 8) {
  ps <- qc$protein_summary
  eligible <- which(ps$Missing_Samples > 0 & ps$Quantified_Samples > 0)
  if (!length(eligible)) {
    p <- ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) + ggplot2::geom_blank() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No partially missing proteins were detected.\nAll proteins are either fully quantified or fully missing.") + ggplot2::theme_void()
    ggplot2::ggsave(out_pdf, p, width = width, height = height)
    return(invisible(NULL))
  }
  eligible <- eligible[order(-ps$Missing_Percent[eligible], ps$ProteinID[eligible])]
  selected <- head(eligible, max(1L, as.integer(heatmap_top_n)))
  binary <- 1 - qc$missing_matrix[selected, , drop = FALSE]
  rownames(binary) <- ps$ProteinID[selected]
  annotation_col <- data.frame(Group = qc$sample_summary$Group, row.names = qc$sample_summary$Sample, stringsAsFactors = FALSE)
  annotation_colors <- list(Group = qc_group_colors(annotation_col$Group, palette))
  fonts <- qc_heatmap_font_settings(ncol(binary), nrow(binary))
  pheatmap::pheatmap(binary, cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotation_col,
    annotation_colors = annotation_colors, show_rownames = fonts$show_rownames, show_colnames = TRUE,
    fontsize = 7, fontsize_row = fonts$fontsize_row, fontsize_col = fonts$fontsize_col,
    fontsize_number = 7, border_color = NA, breaks = c(-.5, .5, 1.5), color = c("#D9D9D9", "#1A1A1A"), legend_breaks = c(0, 1),
    legend_labels = c("Missing", "Quantified"), main = paste0("Missingness heatmap (Top ", nrow(binary), " partially missing proteins)"),
    filename = out_pdf, width = width, height = height)
  invisible(binary)
}

plot_qc_missingness_vs_abundance <- function(qc, out_pdf, width = 4.5, height = 4) {
  df <- qc$protein_summary[is.finite(qc$protein_summary$Median_Log2_Abundance), , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(Median_Log2_Abundance, Missing_Percent)) + ggplot2::geom_point(alpha = .45, size = 1.2) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::labs(x = "Median log2(abundance + 1)", y = "Missing values (%)") + qc_text_theme()
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
  invisible(p)
}

plot_qc_sample_abundance_distribution <- function(qc, out_pdf, palette = "npg", width = 5, height = 4) {
  m <- qc$log2_mat
  df <- data.frame(ProteinID = rep(rownames(m), times = ncol(m)), Sample = rep(colnames(m), each = nrow(m)),
                   Log2_Abundance = as.vector(m), stringsAsFactors = FALSE)
  df <- df[is.finite(df$Log2_Abundance), , drop = FALSE]
  df$Group <- qc$sample_summary$Group[match(df$Sample, qc$sample_summary$Sample)]
  df$Sample <- factor(df$Sample, levels = qc$sample_summary$Sample)
  df$Group <- factor(df$Group, levels = unique(qc$sample_summary$Group))
  axis <- qc_sample_axis_settings(nrow(qc$sample_summary))
  p <- ggplot2::ggplot(df, ggplot2::aes(Sample, Log2_Abundance, fill = Group)) + ggplot2::geom_boxplot(width = .65, outlier.size = .25) +
    ggplot2::scale_fill_manual(values = qc_group_colors(df$Group, palette)) + ggplot2::labs(x = NULL, y = "log2(abundance + 1)", fill = "Group") + qc_text_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = axis$size, angle = axis$angle, hjust = 1, vjust = 1))
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
  invisible(p)
}

run_quantitative_qc <- function(mat, group_info, outdir, heatmap_top_n = 200, palette = "npg", width = 5, height = 4, heatmap_width = 6, heatmap_height = 8) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  qc <- prepare_quantitative_qc(mat, group_info)
  files <- c(
    sample_missingness = file.path(outdir, "quantitative_qc_sample_missingness.pdf"),
    protein_missingness_distribution = file.path(outdir, "quantitative_qc_protein_missingness_distribution.pdf"),
    missingness_heatmap = file.path(outdir, "quantitative_qc_missingness_heatmap.pdf"),
    missingness_vs_abundance = file.path(outdir, "quantitative_qc_missingness_vs_abundance.pdf"),
    sample_abundance_distribution = file.path(outdir, "quantitative_qc_sample_abundance_distribution.pdf"),
    sample_summary = file.path(outdir, "quantitative_qc_sample_summary.csv"),
    protein_summary = file.path(outdir, "quantitative_qc_protein_summary.csv")
  )
  plot_qc_sample_missingness(qc, files[["sample_missingness"]], palette, width, height)
  plot_qc_protein_missingness_distribution(qc, files[["protein_missingness_distribution"]], width, height)
  plot_qc_missingness_heatmap(qc, files[["missingness_heatmap"]], heatmap_top_n, palette, heatmap_width, heatmap_height)
  plot_qc_missingness_vs_abundance(qc, files[["missingness_vs_abundance"]], width, height)
  plot_qc_sample_abundance_distribution(qc, files[["sample_abundance_distribution"]], palette, width, height)
  data.table::fwrite(qc$sample_summary, files[["sample_summary"]])
  data.table::fwrite(qc$protein_summary, files[["protein_summary"]])
  invisible(list(qc = qc, files = unname(files)))
}
