test_quantitative_qc <- function(ctx) {
  record_test("qc_group_specific_missingness", {
    q <- run_quantitative_qc(ctx$canonical$quantity, ctx$canonical$group_info, file.path(ctx$outputs, "qc"), heatmap_top_n = 20, width = 4, height = 3, heatmap_width = 5, heatmap_height = 6)
    pdfs <- q$files[grepl("[.]pdf$", q$files)]; csvs <- q$files[grepl("[.]csv$", q$files)]
    expect_equal(length(pdfs), 5, "Quantitative QC must create five PDFs")
    expect_equal(length(csvs), 2, "Quantitative QC must create two CSV files")
    for (p in q$files) expect_file_nonempty(p)
    sample <- data.table::fread(q$files[[6]], data.table = FALSE)
    protein <- data.table::fread(q$files[[7]], data.table = FALSE)
    expect_columns_include(sample, c("Sample", "Group", "Total_Proteins", "Quantified_Proteins", "Missing_Values", "Missing_Percent", "Median_Log2_Abundance", "Mean_Log2_Abundance"))
    expect_columns_include(protein, c("ProteinID", "Total_Samples", "Quantified_Samples", "Missing_Samples", "Missing_Percent", "Median_Log2_Abundance", "Mean_Log2_Abundance"))
    expect_equal(protein$Missing_Samples[match("P001", protein$ProteinID)], 0, "All-zero legitimate feature should not be missing")
    expect_equal(protein$Missing_Samples[match("P006", protein$ProteinID)], 12, "All-missing feature count differs")
    qc_lines <- readLines("app/R/quantitative_qc.R")
    expect_true(!any(grepl("^[[:space:]]*preprocess_expr[[:space:]]*\\(", qc_lines)), "QC must not call preprocess_expr")
  })
  for (mode in c("zero_is_value", "zero_is_missing")) record_test(paste0("standard_zero_", mode), {
    f <- make_canonical_fixture(seed = 909, n_features = 8, n_samples = 12)
    path <- file.path(ctx$fixtures, paste0(mode, ".csv")); write_synthetic_standard(f, path)
    dat <- read_standard_matrix(path, mode)
    q <- run_quantitative_qc(dat$quantity, f$group_info, file.path(ctx$outputs, mode), heatmap_top_n = 20, width = 4, height = 3, heatmap_width = 5, heatmap_height = 6)
    if (mode == "zero_is_value") expect_true(all(dat$quantity[1, ] == 0), "Mode A changed legitimate zeros") else expect_true(all(is.na(dat$quantity[1, ])), "Mode B did not convert zeros to NA")
    expect_true(all(file.info(q$files)[["size"]] > 0), "Standard QC output contains empty files")
  })
}
