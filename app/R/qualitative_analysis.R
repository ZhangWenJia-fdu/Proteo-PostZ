# Qualitative analysis and shared plotting helpers
make_group_info <- function(samples, groups) {
  groups <- trimws(as.character(groups))
  groups[groups == ""] <- "Group1"
  data.frame(Sample = samples, Group = factor(groups, levels = unique(groups)), stringsAsFactors = FALSE)
}

write_matrix_csv <- function(mat, path, id_col = "ProteinID") {
  out <- data.frame(mat, check.names = FALSE)
  out <- cbind(stats::setNames(data.frame(rownames(mat), check.names = FALSE), id_col), out)
  data.table::fwrite(out, path)
}

plot_sample_count_bar <- function(counts, group_info, count_col, out_pdf, out_csv, width = 3.3, height = 3.3, palette = "npg", y_label = "Protein groups", plot_title = NULL) {
  if (!count_col %in% colnames(counts)) stop("Count column not found: ", count_col)
  df <- dplyr::left_join(counts, group_info, by = "Sample")
  summary <- df |>
    dplyr::group_by(Group) |>
    dplyr::summarise(
      Mean = if (all(is.na(.data[[count_col]]))) NA_real_ else mean(.data[[count_col]], na.rm = TRUE),
      SD = if (sum(!is.na(.data[[count_col]])) < 2) NA_real_ else stats::sd(.data[[count_col]], na.rm = TRUE),
      N = sum(!is.na(.data[[count_col]])),
      .groups = "drop"
    )
  summary$Count_Column <- count_col
  data.table::fwrite(df, sub("\\.csv$", "_sample_counts.csv", out_csv))
  data.table::fwrite(summary, out_csv)
  cols <- sci_palette(length(unique(df$Group)), palette)
  p <- ggplot2::ggplot(summary, ggplot2::aes(Group, Mean, fill = Group)) +
    ggplot2::geom_col(width = 0.65, color = "black", linewidth = 0.25) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.18, linewidth = 0.35) +
    ggplot2::geom_jitter(data = df, ggplot2::aes(x = Group, y = .data[[count_col]]), inherit.aes = FALSE, width = 0.09, size = 1.7) +
    ggplot2::scale_fill_manual(values = cols) + theme_sci() + ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = y_label, title = plot_title %||% y_label)
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
  p
}

plot_identification_bar <- function(counts, group_info, out_pdf, out_csv, width = 3.3, height = 3.3, palette = "npg") {
  plot_sample_count_bar(counts, group_info, "Identified_Protein_Count", out_pdf, out_csv, width, height, palette, "Protein groups")
}

plot_available_quantitative_bar <- function(counts, group_info, out_pdf, out_csv, width = 3.3, height = 3.3, palette = "npg") {
  plot_sample_count_bar(counts, group_info, "Available_Quantitative_Value_Count", out_pdf, out_csv, width, height, palette, "Number of available quantitative values")
}

sci_palette <- function(n, palette = "npg") {
  sets <- list(
    npg = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2"),
    lancet = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A"),
    jama = c("#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"),
    nejm = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91"),
    uchicago = c("#800000", "#767676", "#FFA319", "#8A9045", "#155F83", "#C16622", "#58593F")
  )
  pal <- sets[[palette]] %||% sets$npg
  rep(pal, length.out = n)
}

theme_sci <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", size = 8),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      axis.text = ggplot2::element_text(size = 8, color = "black"),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
      strip.text = ggplot2::element_text(size = 8)
    )
}

identified_by_group <- function(qual_mat, group_info, min_reps = 1) {
  groups <- levels(group_info$Group)
  sets <- lapply(groups, function(g) {
    samples <- group_info$Sample[group_info$Group == g]
    samples <- intersect(samples, colnames(qual_mat))
    rownames(qual_mat)[rowSums(!is.na(qual_mat[, samples, drop = FALSE])) >= min_reps]
  })
  names(sets) <- groups
  sets
}

validate_group_sets <- function(sets, min_reps = 1, analysis = c("venn", "upset")) {
  analysis <- match.arg(analysis)
  set_sizes <- vapply(sets, length, integer(1))
  if (length(set_sizes) == 0 || all(set_sizes == 0)) {
    stop("All group-level protein sets are empty. Minimum replicates detected in group may be too high, or no proteins in the current data meet the detection rule.")
  }
  nonempty <- sets[set_sizes > 0]
  if (length(nonempty) < 2) {
    stop("At least 2 non-empty group-level protein sets are required for Venn/UpSet plots.")
  }
  if (analysis == "venn" && length(nonempty) > 4) {
    stop("Venn diagram is recommended for 2-4 non-empty groups. Please use UpSet plot for 5 or more groups.")
  }
  nonempty
}

make_set_membership <- function(sets) {
  membership <- unique(unlist(sets, use.names = FALSE))
  if (length(membership) == 0) {
    return(data.frame(ProteinID = character(), stringsAsFactors = FALSE))
  }
  set_df <- data.frame(ProteinID = membership, stringsAsFactors = FALSE)
  for (nm in names(sets)) set_df[[nm]] <- as.integer(membership %in% sets[[nm]])
  set_df
}
upset_intersection_data <- function(set_df) {
  group_names <- colnames(set_df)[-1]
  combo <- apply(set_df[, group_names, drop = FALSE], 1, paste, collapse = "")
  counts <- sort(table(combo), decreasing = TRUE)
  labels <- names(counts); values <- as.integer(counts)
  combo_mat <- do.call(rbind, strsplit(labels, "", fixed = TRUE)); storage.mode(combo_mat) <- "integer"
  list(labels = labels, values = values, combo_mat = combo_mat, group_names = group_names)
}
draw_upset_membership <- function(set_df) {
  plot_data <- upset_intersection_data(set_df)
  group_names <- plot_data$group_names
  labels <- plot_data$labels
  values <- plot_data$values
  combo_mat <- plot_data$combo_mat
  n <- length(values)
  graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(2.2, 1.4))
  graphics::par(mar = c(0.5, 4, 2, 1))
  label_cex <- if (n <= 10) 0.95 else if (n <= 20) 0.85 else 0.75
  bar_mid <- graphics::barplot(values, names.arg = rep("", n), col = "#3C5488", border = NA, ylim = c(0, max(1, max(values) * 1.14)), ylab = "Intersection size", main = "UpSet intersections")
  graphics::text(bar_mid, values, labels = as.character(values), pos = 3, offset = 0.25, cex = label_cex, xpd = NA)
  graphics::par(mar = c(4, 4, 0.5, 1))
  y <- rev(seq_along(group_names))
  graphics::plot(seq_len(n), rep(NA_real_, n), type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0.5, n + 0.5), ylim = c(0.5, length(group_names) + 0.5))
  graphics::axis(2, at = y, labels = group_names, las = 1, tick = FALSE, cex.axis = 0.8)
  graphics::axis(1, at = seq_len(n), labels = labels, las = 2, cex.axis = 0.7)
  for (j in seq_len(n)) {
    on <- which(combo_mat[j, ] == 1)
    graphics::points(rep(j, length(y)), y, pch = 16, col = ifelse(seq_along(y) %in% on, "#3C5488", "#D9D9D9"))
    if (length(on) > 1) graphics::lines(rep(j, 2), range(y[on]), col = "#3C5488", lwd = 1.5)
  }
}

plot_venn_upset <- function(sets, outdir, width = 3.3, height = 3.3) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  sets <- validate_group_sets(sets, analysis = "upset")
  set_df <- make_set_membership(sets)
  data.table::fwrite(set_df, file.path(outdir, "group_detected_membership.csv"))
  if (length(sets) >= 2 && length(sets) <= 4) {
    grDevices::pdf(file.path(outdir, "venn.pdf"), width = width, height = height)
    grid::grid.draw(VennDiagram::venn.diagram(sets, filename = NULL, disable.logging = TRUE, fill = sci_palette(length(sets)), alpha = 0.45, cex = 0.8, cat.cex = 0.8, margin = 0.08))
    grDevices::dev.off()
  }
  if (nrow(set_df) == 0 || ncol(set_df) - 1 < 2) {
    stop("UpSet plot requires a membership table with at least one protein and at least 2 non-empty group columns.")
  }
  grDevices::pdf(file.path(outdir, "upset.pdf"), width = max(width, 5), height = max(height, 4))
  draw_upset_membership(set_df)
  grDevices::dev.off()
}
