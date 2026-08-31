# Differential analysis and Volcano plotting
run_volcano <- function(mat, group_info, group_a, group_b, out_pdf, out_csv, log2fc_cutoff = 1, adj_p_cutoff = 0.05, raw_p_cutoff = 0.05, fc_method = c("log2_then_diff", "ratio_then_log2"), test_method = c("limma", "ttest"), sig_metric = c("adj_p", "raw_p"), width = 3.3, height = 3.3) {
  fc_method <- match.arg(fc_method)
  test_method <- match.arg(test_method)
  sig_metric <- match.arg(sig_metric)
  used_log2 <- log2(mat + 1)
  a <- intersect(group_info$Sample[group_info$Group == group_a], colnames(used_log2)); b <- intersect(group_info$Sample[group_info$Group == group_b], colnames(used_log2))
  if (length(a) < 2 || length(b) < 2) stop("Volcano requires at least two samples in each selected group.")
  mean_a_log2 <- rowMeans(used_log2[, a, drop = FALSE], na.rm = TRUE)
  mean_b_log2 <- rowMeans(used_log2[, b, drop = FALSE], na.rm = TRUE)
  mean_a_raw <- rowMeans(mat[, a, drop = FALSE], na.rm = TRUE)
  mean_b_raw <- rowMeans(mat[, b, drop = FALSE], na.rm = TRUE)
  log2fc <- if (fc_method == "log2_then_diff") {
    mean_b_log2 - mean_a_log2
  } else {
    suppressWarnings(log2(mean_b_raw / mean_a_raw))
  }
  log2fc[!is.finite(log2fc)] <- NA_real_
  res <- data.frame(
    ProteinID = rownames(used_log2),
    MeanA_Log2 = mean_a_log2,
    MeanB_Log2 = mean_b_log2,
    MeanA_Raw = mean_a_raw,
    MeanB_Raw = mean_b_raw,
    log2FC = log2fc,
    stringsAsFactors = FALSE
  )
  if (test_method == "limma") {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Volcano plot with limma requires the limma package.")
    }
    sample_group <- factor(c(rep(group_a, length(a)), rep(group_b, length(b))), levels = c(group_a, group_b))
    design <- stats::model.matrix(~ 0 + sample_group)
    colnames(design) <- c(group_a, group_b)
    contrast <- limma::makeContrasts(contrasts = paste0("`", group_b, "`-`", group_a, "`"), levels = design)
    fit <- limma::lmFit(used_log2[, c(a, b), drop = FALSE], design)
    fit2 <- limma::contrasts.fit(fit, contrast)
    fit2 <- limma::eBayes(fit2)
    res$P.Value <- fit2$p.value[, 1]
    p_value_method <- "limma moderated p value"
  } else {
    res$P.Value <- apply(used_log2, 1, function(x) tryCatch(t.test(x[b], x[a])$p.value, error = function(e) NA_real_))
    p_value_method <- "two-sample t-test p value"
  }
  res$BH.Adjusted.P.Value <- p.adjust(res$P.Value, method = "BH")
  sig_values <- if (sig_metric == "adj_p") res$BH.Adjusted.P.Value else res$P.Value
  sig_cutoff <- if (sig_metric == "adj_p") adj_p_cutoff else raw_p_cutoff
  sig_label <- if (sig_metric == "adj_p") "BH-adjusted p value" else "p value"
  res$FCMethod <- fc_method
  res$TestMethod <- test_method
  res$PValueMethod <- p_value_method
  res$AdjustedPValueMethod <- "BH-adjusted p value"
  res$SignificanceMetric <- if (sig_metric == "adj_p") "BH-adjusted p value" else "p value"
  res$SignificanceValue <- sig_values
  res$Regulation <- dplyr::case_when(!is.na(sig_values) & sig_values < sig_cutoff & res$log2FC >= log2fc_cutoff ~ "Up", !is.na(sig_values) & sig_values < sig_cutoff & res$log2FC <= -log2fc_cutoff ~ "Down", TRUE ~ "NotSig")
  data.table::fwrite(res, out_csv)
  plot_df <- res
  plot_df$NegLog10Significance <- -log10(plot_df$SignificanceValue)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(log2FC, NegLog10Significance, color = Regulation)) + ggplot2::geom_point(size = 1.2, alpha = 0.75, na.rm = TRUE) + ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") + ggplot2::geom_hline(yintercept = -log10(sig_cutoff), linetype = "dashed") + ggplot2::scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC", NotSig = "grey75")) + theme_sci() + ggplot2::labs(x = "log2FC", y = paste0("-log10(", sig_label, ")"))
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
}
