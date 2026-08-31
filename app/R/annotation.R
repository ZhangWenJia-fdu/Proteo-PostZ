# Offline physicochemical annotation analysis
run_physicochemical <- function(sets, annotation_file, outdir, width = 3.3, height = 3.3, palette = "npg") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ann <- data.table::fread(annotation_file, data.table = FALSE)
  required <- c("Accession", "GRAVY", "TM_helices", "Subcellular_class", "MW", "pI", "Length")
  miss <- setdiff(required, colnames(ann))
  if (length(miss) > 0) stop("Annotation table missing columns: ", paste(miss, collapse = ", "))
  if (is.null(names(sets)) || any(!nzchar(names(sets)))) stop("Physicochemical analysis requires non-empty experimental group labels.")
  names(sets) <- as.character(names(sets))
  long <- dplyr::bind_rows(lapply(names(sets), function(g) data.frame(Group = g, Accession = sets[[g]], stringsAsFactors = FALSE))) |>
    dplyr::left_join(ann, by = "Accession")
  coverage <- dplyr::bind_rows(lapply(names(sets), function(g) {
    ids <- unique(as.character(sets[[g]])); ids <- ids[!is.na(ids) & nzchar(ids)]
    matched <- sum(ids %in% as.character(ann$Accession))
    data.frame(Group = g, Queried_Proteins = length(ids), Matched_Proteins = matched, Unmatched_Proteins = length(ids) - matched, Match_Percent = if (length(ids)) 100 * matched / length(ids) else NA_real_, stringsAsFactors = FALSE)
  }))
  data.table::fwrite(coverage, file.path(outdir, "annotation_mapping_coverage.csv"))
  if (sum(coverage$Matched_Proteins) == 0) stop("No proteins could be matched to the selected annotation table/species. Check annotation species / ProteinID compatibility.", call. = FALSE)
  tm_status <- if ("AnnotationStatus" %in% colnames(long)) as.character(long$AnnotationStatus) else rep(NA_character_, nrow(long))
  if ("TM_helices" %in% colnames(long)) long$TM_helices[!is.na(tm_status) & tm_status == "not_found_in_uniprot"] <- NA
  data.table::fwrite(long, file.path(outdir, "detected_proteins_annotations_long.csv"))
  props <- c("GRAVY", "MW", "pI", "Length")
  for (prop in props) {
    df <- long |> dplyr::filter(!is.na(.data[[prop]]))
    if (nrow(df) == 0) next
    p <- ggplot2::ggplot(df, ggplot2::aes(Group, .data[[prop]], fill = Group)) + ggplot2::geom_violin(trim = TRUE, alpha = 0.65) + ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.3) + ggplot2::scale_fill_manual(values = sci_palette(length(unique(df$Group)), palette)) + theme_sci() + ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(file.path(outdir, paste0(prop, "_distribution.pdf")), p, width = width, height = height)
  }
  tm_df <- long |>
    dplyr::mutate(
      TM_helices_numeric = suppressWarnings(as.numeric(TM_helices)),
      TM_valid_annotation = !is.na(TM_helices_numeric) & (is.na(tm_status) | tm_status != "not_found_in_uniprot"),
      TM_category = dplyr::case_when(
        TM_valid_annotation & TM_helices_numeric == 0 ~ "0",
        TM_valid_annotation & TM_helices_numeric == 1 ~ "1",
        TM_valid_annotation & TM_helices_numeric == 2 ~ "2",
        TM_valid_annotation & TM_helices_numeric >= 3 ~ ">=3",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(TM_valid_annotation, !is.na(TM_category))
  if (nrow(tm_df) > 0) {
    tm_levels <- c("0", "1", "2", ">=3")
    tm_tab <- tm_df |>
      dplyr::count(Group, TM_category, name = "protein_count") |>
      tidyr::complete(Group, TM_category = tm_levels, fill = list(protein_count = 0)) |>
      dplyr::group_by(Group) |>
      dplyr::mutate(
        denominator = sum(protein_count),
        proportion = dplyr::if_else(denominator > 0, protein_count / denominator, NA_real_)
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(TM_category = factor(TM_category, levels = tm_levels)) |>
      dplyr::arrange(Group, TM_category)
    data.table::fwrite(tm_tab, file.path(outdir, "TM_helices_category_proportions.csv"))
    cols <- sci_palette(length(tm_levels), palette)
    names(cols) <- tm_levels
    p <- ggplot2::ggplot(tm_tab, ggplot2::aes(Group, proportion, fill = TM_category)) +
      ggplot2::geom_col(width = 0.75, color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      theme_sci() +
      ggplot2::labs(x = NULL, y = "Proportion", fill = "TM helices")
    ggplot2::ggsave(file.path(outdir, "TM_helices_category_proportions.pdf"), p, width = width, height = height)
  }
  if (any(!is.na(long$Subcellular_class))) {
    tab <- long |>
      dplyr::filter(!is.na(Subcellular_class), Subcellular_class != "") |>
      dplyr::count(Group, Subcellular_class) |>
      dplyr::group_by(Group) |>
      dplyr::mutate(
        denominator = sum(n),
        proportion = dplyr::if_else(denominator > 0, n / denominator, NA_real_),
        label = dplyr::if_else(proportion >= 0.06, sprintf("%.1f%%\nn=%s", proportion * 100, n), NA_character_)
      ) |>
      dplyr::ungroup()
    data.table::fwrite(tab, file.path(outdir, "Subcellular_class_counts.csv"))
    # Labels below 6% are hidden because narrow stacked segments tend to overlap.
    p <- ggplot2::ggplot(tab, ggplot2::aes(Group, proportion, fill = Subcellular_class)) +
      ggplot2::geom_col(width = 0.75, color = "white", linewidth = 0.2) +
      ggplot2::geom_text(ggplot2::aes(label = label), position = ggplot2::position_stack(vjust = 0.5), size = 2.2, na.rm = TRUE) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      theme_sci() +
      ggplot2::labs(y = "Fraction")
    ggplot2::ggsave(file.path(outdir, "Subcellular_class_distribution.pdf"), p, width = width, height = height)
  }
  invisible(list(data = long, coverage = coverage))
}
