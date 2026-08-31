same_frame <- function(a, b, numeric_cols = character()) {
  if (!identical(a$ProteinID, b$ProteinID)) return(FALSE)
  for (nm in setdiff(names(a), c("ProteinID", "Repeats"))) {
    if (!nm %in% names(b)) return(FALSE)
    if (nm %in% numeric_cols) {
      if (!isTRUE(all.equal(as.numeric(a[[nm]]), as.numeric(b[[nm]]), tolerance = 1e-12, check.attributes = FALSE))) return(FALSE)
    } else if (!isTRUE(all.equal(a[[nm]], b[[nm]], check.attributes = FALSE))) return(FALSE)
  }
  TRUE
}

test_ml_high_fidelity <- function(ctx) {
  if (!require_namespace_or_skip("randomForest", "hf_ml_packages")) return(invisible(NULL))
  if (!require_namespace_or_skip("glmnet", "hf_ml_packages")) return(invisible(NULL))
  f <- make_binary_ml_fixture(seed = 202, n_features = 320, n_per_group = 12)
  root <- file.path(ctx$outputs, "high_fidelity")
  common <- list(seed = 2026, split_mode = "cross_validation_only", train_prop = 0.8, allow_small_sample = FALSE, stability_repeats = 50, stability_sample_fraction = 0.8, stability_top_var_n = 200)
  rf_call <- function(d, top) do.call(run_random_forest_selection, c(list(mat = f$quantity, group_info = f$group_info, outdir = d, top_n = top, rf_ntree = 500, rf_mtry = NA), common))
  l1_call <- function(d, folds = "Auto") do.call(run_l1_selection, c(list(mat = f$quantity, group_info = f$group_info, outdir = d, top_n = 50, l1_alpha = 1, lambda_selection = "lambda.1se", cv_folds = folds), common))
  read_rf <- function(d) data.table::fread(file.path(d, "random_forest_stability_selection.csv"), data.table = FALSE)
  read_combined <- function(d) data.table::fread(file.path(d, "combined_stability_scores.csv"), data.table = FALSE)
  numeric_rf <- c("RFTop20Frequency", "RFTop50Frequency", "MeanDecreaseGini", "MeanDecreaseGiniScaled", "StabilityScore")

  record_test("hf_rf_topn_invariance", {
    dirs <- file.path(root, paste0("rf_top", c(20, 50, 100)))
    tops <- Map(rf_call, dirs, c(20, 50, 100))
    frames <- lapply(dirs, read_rf)
    expect_true(all(vapply(frames, function(x) !anyDuplicated(x$ProteinID) && all(is.finite(x$StabilityScore)), logical(1))), "RF stability ranking invalid")
    expect_true(same_frame(frames[[1]], frames[[2]], numeric_rf) && same_frame(frames[[2]], frames[[3]], numeric_rf), "RF stability statistics changed with final Top N")
    expect_true(length(tops[[1]]) <= 20 && length(tops[[2]]) <= 50 && length(tops[[3]]) <= 100, "RF final Top N output exceeds request")
  })
  record_test("hf_rf_reproducibility", {
    d1 <- file.path(root, "rf_top50"); d2 <- file.path(root, "rf_repeat"); rf_call(d2, 50)
    a <- read_rf(d1); b <- read_rf(d2); expect_true(same_frame(a, b, numeric_rf), "Standalone RF repeated high-fidelity result differs")
    s <- data.table::fread(file.path(d1, "random_forest_ml_settings.csv"), data.table = FALSE); expect_true(any(s$Setting == "RandomSeed" & s$Value == 2026), "RF validation seed was not recorded")
  })
  record_test("hf_seed_sensitivity", {
    p1 <- make_stability_split_plan(factor(rep(c("A", "B"), each = 12)), repeats = 5, train_prop = 0.8, seed = 2026)
    p2 <- make_stability_split_plan(factor(rep(c("A", "B"), each = 12)), repeats = 5, train_prop = 0.8, seed = 2027)
    expect_true(!identical(p1, p2), "Stability split plan did not respond to seed")
  })
  record_test("hf_l1_folds_and_reproducibility", {
    d1 <- file.path(root, "l1_auto"); d2 <- file.path(root, "l1_repeat"); l1_call(d1, "Auto"); l1_call(d2, "Auto")
    a <- data.table::fread(file.path(d1, "l1_stability_selection.csv"), data.table = FALSE); b <- data.table::fread(file.path(d2, "l1_stability_selection.csv"), data.table = FALSE)
    expect_equal(a, b, "Standalone L1 repeated high-fidelity result differs")
    s <- data.table::fread(file.path(d1, "l1_feature_selection_ml_settings.csv"), data.table = FALSE); expect_true(any(s$Setting == "L1Family" & s$Value == "binomial"), "Binary L1 is not binomial")
    p <- data.table::fread(file.path(d1, "l1_stability_resample_performance.csv"), data.table = FALSE); expect_columns_include(p, c("RequestedCVFolds", "EffectiveCVFolds"), "High-fidelity L1 folds"); expect_true(all(p$EffectiveCVFolds == 5), "High-fidelity Auto folds did not resolve to 5")
  })
  rfl1_call <- function(d, top, folds = "Auto") do.call(run_feature_selection, c(list(mat = f$quantity, group_info = f$group_info, outdir = d, top_n = top, rf_ntree = 500, rf_mtry = NA, l1_alpha = 1, lambda_selection = "lambda.1se", cv_folds = folds, stability_rf_top20_weight = 0.35, stability_rf_top50_weight = 0.25, stability_l1_weight = 0.30, stability_gini_weight = 0.10), common))
  record_test("hf_rfl1_topn_and_formula", {
    d50 <- file.path(root, "rfl1_top50"); d100 <- file.path(root, "rfl1_top100"); rfl1_call(d50, 50); rfl1_call(d100, 100); a <- read_combined(d50); b <- read_combined(d100)
    cols <- c("RFTop20Frequency", "RFTop50Frequency", "L1SelectionFrequency", "MeanDecreaseGini", "MeanDecreaseGiniScaled", "StabilityScore"); expect_true(same_frame(a[, c("ProteinID", cols)], b[, c("ProteinID", cols)], cols), "RF+L1 ranking changed with final Top N"); expect_equal(head(b$ProteinID, 50), head(a$ProteinID, 50), "RF+L1 first 50 changed")
    expected <- 0.35 * a$RFTop20Frequency + 0.25 * a$RFTop50Frequency + 0.30 * a$L1SelectionFrequency + 0.10 * a$MeanDecreaseGiniScaled; expect_near(a$StabilityScore, expected, tolerance = 1e-12, message = "RF+L1 score formula differs")
    expect_true(!anyDuplicated(a$ProteinID) && all(is.finite(a$StabilityScore)), "RF+L1 scores invalid")
  })
  record_test("hf_rfl1_reproducibility", {
    d1 <- file.path(root, "rfl1_top50"); d2 <- file.path(root, "rfl1_repeat"); rfl1_call(d2, 50); a <- read_combined(d1); b <- read_combined(d2); cols <- c("RFTop20Frequency", "RFTop50Frequency", "L1SelectionFrequency", "MeanDecreaseGini", "MeanDecreaseGiniScaled", "StabilityScore"); expect_true(same_frame(a[, c("ProteinID", cols)], b[, c("ProteinID", cols)], cols), "RF+L1 repeated high-fidelity result differs")
  })
  record_test("hf_multiclass_rfl1", {
    f3 <- make_multiclass_fixture(seed = 303, n_features = 320, n_per_group = 10); d <- file.path(root, "multiclass_rfl1"); rfl1_call3 <- function(folds) run_feature_selection(f3$quantity, f3$group_info, d, top_n = 50, rf_ntree = 500, l1_alpha = 1, seed = 2026, split_mode = "cross_validation_only", train_prop = 0.8, rf_mtry = NA, lambda_selection = "lambda.1se", cv_folds = folds, allow_small_sample = FALSE, stability_repeats = 50, stability_sample_fraction = 0.8, stability_top_var_n = 200); rfl1_call3("Auto"); x <- read_combined(d); expect_true(nrow(x) > 0 && !anyDuplicated(x$ProteinID) && all(is.finite(x$StabilityScore)), "Multiclass RF+L1 result invalid"); s <- data.table::fread(file.path(d, "rf_l1_combined_feature_selection_ml_settings.csv"), data.table = FALSE); expect_true(any(s$Setting == "L1Family" & s$Value == "multinomial"), "Multiclass L1 family is not multinomial"); expected <- 0.35 * x$RFTop20Frequency + 0.25 * x$RFTop50Frequency + 0.30 * x$L1SelectionFrequency + 0.10 * x$MeanDecreaseGiniScaled; expect_near(x$StabilityScore, expected, tolerance = 1e-12, message = "Multiclass score formula differs")
  })
}
