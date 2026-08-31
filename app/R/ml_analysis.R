# Machine-learning feature selection and evaluation
format_class_counts <- function(class_n) {
  paste(paste0(names(class_n), " (n=", as.integer(class_n), ")"), collapse = ", ")
}

prepare_ml_input <- function(mat, group_info) {
  used <- preprocess_expr(mat, TRUE, 0.5)
  x <- t(used)
  y <- factor(group_info$Group[match(rownames(x), group_info$Sample)])
  keep_samples <- !is.na(y)
  x <- x[keep_samples, , drop = FALSE]
  y <- droplevels(y[keep_samples])
  list(x = x, y = y)
}

check_ml_classes <- function(y, label, min_per_class = 2) {
  class_n <- table(y)
  if (length(class_n) < 2) stop(label, " requires at least two groups.")
  too_small <- names(class_n)[class_n < min_per_class]
  if (length(too_small) > 0) {
    stop(label, " requires at least ", min_per_class, " samples per group. Insufficient group(s): ", paste(paste0(too_small, " (n=", class_n[too_small], ")"), collapse = ", "), ".")
  }
  class_n
}

check_ml_sample_policy <- function(y, label, allow_small_sample = FALSE, strict_min_per_class = 6, exploratory_min_per_class = 2) {
  if (allow_small_sample) {
    return(check_ml_classes(y, label, exploratory_min_per_class))
  }
  class_n <- table(y)
  if (length(class_n) < 2) stop(label, " requires at least two groups.")
  too_small <- names(class_n)[class_n < strict_min_per_class]
  if (length(too_small) > 0) {
    stop(label, " requires at least ", strict_min_per_class, " samples per group in the default strict mode. Insufficient group(s): ", paste(paste0(too_small, " (n=", class_n[too_small], ")"), collapse = ", "), ". Small-sample machine learning is exploratory and can be enabled with 'Allow small-sample exploratory ML'.")
  }
  class_n
}

stratified_train_test_split <- function(y, train_prop = 0.7, seed = 123, min_train_per_class = 2, min_test_per_class = 1) {
  if (!is.finite(train_prop) || train_prop <= 0 || train_prop >= 1) stop("Training set proportion must be between 0 and 1.")
  set.seed(seed)
  train_idx <- integer()
  test_idx <- integer()
  for (cls in levels(y)) {
    idx <- which(y == cls)
    n <- length(idx)
    if (n < min_train_per_class + min_test_per_class) {
      stop("Train/test split is not supported because ", cls, " has n=", n, "; need at least ", min_train_per_class, " training and ", min_test_per_class, " test samples.")
    }
    idx <- sample(idx, n)
    n_train <- ceiling(n * train_prop)
    n_train <- max(min_train_per_class, min(n_train, n - min_test_per_class))
    train_idx <- c(train_idx, idx[seq_len(n_train)])
    test_idx <- c(test_idx, idx[(n_train + 1):n])
  }
  list(train = sort(train_idx), test = sort(test_idx))
}

resolve_ml_training <- function(y, mode = "auto", train_prop = 0.7, seed = 123, min_train_per_class = 2, min_test_per_class = 1, label = "Machine learning", allow_small_sample = FALSE, train_test_min_per_class = 8) {
  mode <- tolower(gsub("[ _-]+", "_", mode %||% "auto"))
  if (!mode %in% c("auto", "cross_validation_only", "train_test_split")) stop("Unknown train/test split mode: ", mode)
  class_n <- table(y)
  if (allow_small_sample) {
    return(list(mode = "Cross-validation only", train = seq_along(y), test = integer(), class_n = class_n, train_class_n = class_n, test_class_n = NULL, auto_note = "Small-sample exploratory ML enabled: independent train/test split disabled; results are for hypothesis generation only."))
  }
  if (mode == "cross_validation_only") {
    return(list(mode = "Cross-validation only", train = seq_along(y), test = integer(), class_n = class_n, train_class_n = class_n, test_class_n = NULL))
  }
  below_train_test <- names(class_n)[class_n < train_test_min_per_class]
  if (length(below_train_test) > 0) {
    msg <- paste0("Train/test split requires at least ", train_test_min_per_class, " samples per group. Insufficient group(s): ", paste(paste0(below_train_test, " (n=", class_n[below_train_test], ")"), collapse = ", "), ".")
    if (mode == "train_test_split") stop(msg)
    return(list(mode = "Cross-validation only", train = seq_along(y), test = integer(), class_n = class_n, train_class_n = class_n, test_class_n = NULL, auto_note = paste("Auto used cross-validation only:", msg)))
  }
  split <- tryCatch(stratified_train_test_split(y, train_prop, seed, min_train_per_class, min_test_per_class), error = function(e) e)
  if (inherits(split, "error")) {
    if (mode == "train_test_split") stop(conditionMessage(split))
    return(list(mode = "Cross-validation only", train = seq_along(y), test = integer(), class_n = class_n, train_class_n = class_n, test_class_n = NULL, auto_note = paste("Auto used cross-validation only:", conditionMessage(split))))
  }
  train_class_n <- table(droplevels(y[split$train]))
  test_class_n <- table(droplevels(y[split$test]))
  bad_train <- names(train_class_n)[train_class_n < min_train_per_class]
  missing_test <- setdiff(levels(y), names(test_class_n))
  bad_test <- names(test_class_n)[test_class_n < min_test_per_class]
  if (length(bad_train) > 0 || length(missing_test) > 0 || length(bad_test) > 0) {
    msg <- paste0("Stratified split cannot provide valid train/test class counts. Train: ", format_class_counts(train_class_n), "; test: ", format_class_counts(test_class_n), ".")
    if (mode == "train_test_split") stop(msg)
    return(list(mode = "Cross-validation only", train = seq_along(y), test = integer(), class_n = class_n, train_class_n = class_n, test_class_n = NULL, auto_note = paste("Auto used cross-validation only:", msg)))
  }
  list(mode = "Train/test split", train = split$train, test = split$test, class_n = class_n, train_class_n = train_class_n, test_class_n = test_class_n)
}

parse_auto_integer <- function(value, default = NA_integer_) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])) return(default)
  value <- trimws(as.character(value[[1]]))
  if (!nzchar(value) || tolower(value) == "auto") return(default)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1) stop("Expected Auto or a positive integer, got: ", value)
  out
}

write_ml_settings <- function(outdir, settings) {
  analysis <- tolower(gsub("[^A-Za-z0-9]+", "_", settings[["Analysis"]] %||% "ml"))
  analysis <- gsub("^_|_$", "", analysis)
  data.table::fwrite(data.frame(Setting = names(settings), Value = unname(unlist(settings)), stringsAsFactors = FALSE), file.path(outdir, paste0(analysis, "_ml_settings.csv")))
}

safe_div <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

classification_eval_tables <- function(actual, predicted, probabilities = NULL, positive_class = NULL) {
  actual <- droplevels(factor(actual))
  predicted <- factor(predicted, levels = levels(actual))
  cm <- table(Actual = actual, Predicted = predicted)
  total <- sum(cm)
  classes <- levels(actual)
  per_class <- dplyr::bind_rows(lapply(classes, function(cls) {
    tp <- cm[cls, cls]
    fn <- sum(cm[cls, ]) - tp
    fp <- sum(cm[, cls]) - tp
    tn <- total - tp - fn - fp
    data.frame(
      Class = cls,
      TP = as.integer(tp),
      FP = as.integer(fp),
      TN = as.integer(tn),
      FN = as.integer(fn),
      Sensitivity = safe_div(tp, tp + fn),
      Specificity = safe_div(tn, tn + fp),
      Precision = safe_div(tp, tp + fp),
      F1 = safe_div(2 * tp, 2 * tp + fp + fn),
      stringsAsFactors = FALSE
    )
  }))
  accuracy <- safe_div(sum(diag(cm)), total)
  auc_value <- NA_real_
  roc_df <- NULL
  if (!is.null(probabilities) && length(classes) == 2 && requireNamespace("pROC", quietly = TRUE)) {
    probs <- as.data.frame(probabilities, check.names = FALSE)
    pos <- positive_class %||% classes[[2]]
    if (pos %in% colnames(probs)) {
      roc_obj <- tryCatch(pROC::roc(response = actual, predictor = probs[[pos]], levels = classes, direction = "<", quiet = TRUE), error = function(e) NULL)
      if (!is.null(roc_obj)) {
        auc_value <- as.numeric(pROC::auc(roc_obj))
        roc_df <- data.frame(
          Specificity = rev(roc_obj$specificities),
          Sensitivity = rev(roc_obj$sensitivities),
          Threshold = rev(roc_obj$thresholds),
          PositiveClass = pos,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  summary <- data.frame(
    Metric = c("Accuracy", "Macro sensitivity", "Macro specificity", "Macro F1", "AUC"),
    Value = c(accuracy, mean(per_class$Sensitivity, na.rm = TRUE), mean(per_class$Specificity, na.rm = TRUE), mean(per_class$F1, na.rm = TRUE), auc_value),
    stringsAsFactors = FALSE
  )
  list(confusion = as.data.frame.matrix(cm), per_class = per_class, summary = summary, roc = roc_df)
}

write_classification_evaluation <- function(actual, predicted, probabilities, outdir, prefix, positive_class = NULL) {
  eval <- classification_eval_tables(actual, predicted, probabilities, positive_class)
  confusion <- cbind(Actual = rownames(eval$confusion), eval$confusion)
  data.table::fwrite(confusion, file.path(outdir, paste0(prefix, "_confusion_matrix.csv")))
  data.table::fwrite(eval$per_class, file.path(outdir, paste0(prefix, "_class_metrics.csv")))
  data.table::fwrite(eval$summary, file.path(outdir, paste0(prefix, "_summary_metrics.csv")))
  if (!is.null(eval$roc)) data.table::fwrite(eval$roc, file.path(outdir, paste0(prefix, "_roc_curve.csv")))
  invisible(eval)
}

make_cv_folds <- function(y, seed = 123, max_folds = 5) {
  class_n <- table(y)
  nfolds <- max(2, min(max_folds, min(class_n)))
  make_stratified_foldid(y, nfolds, seed)
}

resolve_l1_family <- function(y) {
  n_groups <- nlevels(droplevels(y))
  if (n_groups == 2) return("binomial")
  if (n_groups >= 3) return("multinomial")
  stop("L1 classification requires at least two groups.")
}

fit_l1_cv <- function(x, y, l1_alpha, lambda_selection, foldid, grouped = TRUE) {
  family <- resolve_l1_family(y)
  args <- list(
    x = x, y = y, family = family, alpha = l1_alpha,
    type.measure = "class", nfolds = max(foldid), foldid = foldid
  )
  if (identical(family, "multinomial")) {
    args$type.multinomial <- "ungrouped"
    args$grouped <- grouped
  }
  do.call(glmnet::cv.glmnet, args)
}

predict_l1_probabilities <- function(model, x, s, classes) {
  family <- as.character(model$glmnet.fit$call$family %||% resolve_l1_family(factor(classes, levels = classes)))[1]
  if (identical(family, "binomial")) {
    positive <- as.numeric(predict(model, x, s = s, type = "response"))
    probs <- cbind(1 - positive, positive)
    colnames(probs) <- classes
    return(as.data.frame(probs, check.names = FALSE))
  }
  predicted <- predict(model, x, s = s, type = "response")
  probs <- predicted[, , 1, drop = FALSE][, , 1]
  probs <- as.matrix(probs)
  if (is.null(colnames(probs))) colnames(probs) <- classes
  as.data.frame(probs, check.names = FALSE)
}

extract_l1_coefficients <- function(model, s, classes = NULL) {
  family <- as.character(model$glmnet.fit$call$family %||% "multinomial")[1]
  co <- coef(model, s = s)
  if (identical(family, "binomial")) {
    m <- as.matrix(co)
    positive_class <- if (!is.null(classes) && length(classes) >= 2) classes[[2]] else "positive"
    return(data.frame(ProteinID = rownames(m), Class = positive_class, Coefficient = as.numeric(m[, 1]), row.names = NULL))
  }
  dplyr::bind_rows(lapply(names(co), function(cls) {
    m <- as.matrix(co[[cls]])
    data.frame(ProteinID = rownames(m), Class = cls, Coefficient = as.numeric(m[, 1]), row.names = NULL)
  }))
}

rf_predict_prob <- function(model, x) {
  probs <- predict(model, x, type = "prob")
  as.data.frame(probs, check.names = FALSE)
}

write_rf_evaluation <- function(x, y, training, rf_args, rf_model, outdir, seed = 123) {
  if (length(training$test) > 0) {
    test_x <- x[training$test, , drop = FALSE]
    test_y <- droplevels(y[training$test])
    probs <- rf_predict_prob(rf_model, test_x)
    pred <- factor(colnames(probs)[max.col(probs, ties.method = "first")], levels = levels(test_y))
    write_classification_evaluation(test_y, pred, probs, outdir, "random_forest_test")
  }
  foldid <- make_cv_folds(y, seed, max_folds = 5)
  pred <- rep(NA_character_, length(y))
  probs_all <- matrix(NA_real_, nrow = length(y), ncol = length(levels(y)), dimnames = list(names(y), levels(y)))
  for (fold in sort(unique(foldid))) {
    train_idx <- which(foldid != fold)
    test_idx <- which(foldid == fold)
    args <- rf_args
    args$x <- x[train_idx, , drop = FALSE]
    args$y <- droplevels(y[train_idx])
    set.seed(seed + fold)
    fit <- do.call(randomForest::randomForest, args)
    probs <- rf_predict_prob(fit, x[test_idx, , drop = FALSE])
    pred[test_idx] <- colnames(probs)[max.col(probs, ties.method = "first")]
    probs_all[test_idx, colnames(probs)] <- as.matrix(probs)
  }
  write_classification_evaluation(y, factor(pred, levels = levels(y)), as.data.frame(probs_all, check.names = FALSE), outdir, "random_forest_cross_validation")
  data.table::fwrite(data.frame(Sample = rownames(x), Fold = foldid, Actual = as.character(y), Predicted = pred, probs_all, check.names = FALSE), file.path(outdir, "random_forest_cross_validation_predictions.csv"))
}

write_l1_evaluation <- function(x, y, training, l1_alpha, lambda_selection, nfolds, outdir, seed = 123) {
  family <- resolve_l1_family(y)
  if (length(training$test) > 0) {
    x_train <- x[training$train, , drop = FALSE]
    y_train <- droplevels(y[training$train])
    effective <- resolve_l1_folds(y_train, nfolds)$nfolds
    foldid_train <- make_stratified_foldid(y_train, effective, seed)
    fit <- fit_l1_cv(x_train, y_train, l1_alpha, lambda_selection, foldid_train, grouped = TRUE)
    probs <- predict_l1_probabilities(fit, x[training$test, , drop = FALSE], lambda_selection, levels(y_train))
    test_y <- droplevels(y[training$test])
    pred <- factor(colnames(probs)[max.col(probs, ties.method = "first")], levels = levels(test_y))
    write_classification_evaluation(test_y, pred, probs, outdir, "l1_test")
  }
  effective_outer <- resolve_l1_folds(y, nfolds, allow_small_sample = TRUE)$nfolds
  foldid <- make_cv_folds(y, seed, max_folds = effective_outer)
  pred <- rep(NA_character_, length(y))
  probs_all <- matrix(NA_real_, nrow = length(y), ncol = length(levels(y)), dimnames = list(names(y), levels(y)))
  for (fold in sort(unique(foldid))) {
    train_idx <- which(foldid != fold)
    test_idx <- which(foldid == fold)
    y_train <- droplevels(y[train_idx])
    inner_folds <- resolve_l1_folds(y_train, nfolds)$nfolds
    inner_foldid <- make_stratified_foldid(y_train, inner_folds, seed + fold)
    set.seed(seed + fold)
    fit <- fit_l1_cv(x[train_idx, , drop = FALSE], y_train, l1_alpha, lambda_selection, inner_foldid, grouped = TRUE)
    probs <- predict_l1_probabilities(fit, x[test_idx, , drop = FALSE], lambda_selection, levels(y_train))
    pred[test_idx] <- colnames(probs)[max.col(probs, ties.method = "first")]
    probs_all[test_idx, colnames(probs)] <- as.matrix(probs)
  }
  write_classification_evaluation(y, factor(pred, levels = levels(y)), as.data.frame(probs_all, check.names = FALSE), outdir, "l1_cross_validation")
  data.table::fwrite(data.frame(Sample = rownames(x), Fold = foldid, Actual = as.character(y), Predicted = pred, probs_all, check.names = FALSE), file.path(outdir, "l1_cross_validation_predictions.csv"))
}

prepare_stability_input <- function(mat, group_info) {
  used <- as.matrix(mat)
  mode(used) <- "numeric"
  used <- used[rowSums(!is.na(used)) > 0, , drop = FALSE]
  used <- log2(used + 1)
  keep_n <- max(1, ceiling(ncol(used) * 0.5))
  used <- used[rowSums(!is.na(used)) >= keep_n, , drop = FALSE]
  x <- t(used)
  y <- factor(group_info$Group[match(rownames(x), group_info$Sample)])
  keep_samples <- !is.na(y)
  list(x = x[keep_samples, , drop = FALSE], y = droplevels(y[keep_samples]))
}

impute_from_training <- function(train_x, test_x) {
  medians <- apply(train_x, 2, median, na.rm = TRUE)
  medians[!is.finite(medians)] <- 0
  for (j in seq_len(ncol(train_x))) {
    train_x[is.na(train_x[, j]), j] <- medians[j]
    test_x[is.na(test_x[, j]), j] <- medians[j]
  }
  list(train = train_x, test = test_x)
}

make_stability_split_plan <- function(y, repeats = 50, train_prop = 0.8, seed = 123) {
  lapply(seq_len(repeats), function(i) {
    stratified_train_test_split(y, train_prop = train_prop, seed = seed + i - 1, min_train_per_class = 2, min_test_per_class = 1)
  })
}

write_rf_stability <- function(x, y, outdir, top_n = 50, rf_ntree = 500, mtry = NA_integer_, seed = 123, repeats = 50, sample_fraction = 0.8, top_var_n = 200, top20_weight = 0.35, top50_weight = 0.25, gini_weight = 0.10, split_plan = NULL) {
  ranked <- list()
  gini_values <- list()
  performance <- vector("list", repeats)
  split_plan <- split_plan %||% make_stability_split_plan(y, repeats, sample_fraction, seed + 1000)
  for (i in seq_len(repeats)) {
    split <- split_plan[[i]]
    train_x <- x[split$train, , drop = FALSE]
    train_y <- droplevels(y[split$train])
    test_x <- x[split$test, , drop = FALSE]
    test_y <- droplevels(y[split$test])
    train_var <- apply(train_x, 2, var, na.rm = TRUE)
    train_var[!is.finite(train_var)] <- -Inf
    keep <- order(train_var, decreasing = TRUE)[seq_len(min(top_var_n, ncol(train_x)))]
    imputed <- impute_from_training(train_x[, keep, drop = FALSE], test_x[, keep, drop = FALSE])
    args <- list(x = train_x[, keep, drop = FALSE], y = train_y, ntree = rf_ntree, importance = TRUE)
    args$x <- imputed$train
    if (!is.na(mtry)) args$mtry <- mtry
    fit <- do.call(randomForest::randomForest, args)
    pred <- predict(fit, imputed$test)
    recalls <- diag(table(factor(test_y, levels = levels(train_y)), factor(pred, levels = levels(train_y)))) / rowSums(table(factor(test_y, levels = levels(train_y)), factor(pred, levels = levels(train_y))))
    performance[[i]] <- data.frame(Resample = i, TrainSamples = length(split$train), TestSamples = length(split$test), Accuracy = mean(pred == test_y), BalancedAccuracy = mean(recalls, na.rm = TRUE), stringsAsFactors = FALSE)
    imp <- randomForest::importance(fit)
    col <- if ("MeanDecreaseGini" %in% colnames(imp)) "MeanDecreaseGini" else tail(colnames(imp), 1)
    ranked[[i]] <- rownames(imp)[order(imp[, col], decreasing = TRUE)]
    gini_values[[i]] <- setNames(as.numeric(imp[, col]), rownames(imp))
  }
  all_ids <- colnames(x)
  top20_counts <- stats::setNames(numeric(length(all_ids)), all_ids)
  top50_counts <- stats::setNames(numeric(length(all_ids)), all_ids)
  for (ranks in ranked) {
    top20_counts[head(ranks, min(20, length(ranks)))] <- top20_counts[head(ranks, min(20, length(ranks)))] + 1
    top50_counts[head(ranks, min(50, length(ranks)))] <- top50_counts[head(ranks, min(50, length(ranks)))] + 1
  }
  top20_freq <- top20_counts / repeats
  top50_freq <- top50_counts / repeats
  mean_gini <- vapply(all_ids, function(id) sum(vapply(gini_values, function(v) unname(v[id] %||% 0), numeric(1))) / repeats, numeric(1))
  min_gini <- min(mean_gini, na.rm = TRUE)
  max_gini <- max(mean_gini, na.rm = TRUE)
  gini_scaled <- if (!is.finite(min_gini) || !is.finite(max_gini) || max_gini == min_gini) rep(0, length(mean_gini)) else (mean_gini - min_gini) / (max_gini - min_gini)
  score <- top20_weight * top20_freq + top50_weight * top50_freq + gini_weight * gini_scaled
  score[!is.finite(score)] <- 0
  ord <- order(score, top20_freq, top50_freq, mean_gini, all_ids, decreasing = c(TRUE, TRUE, TRUE, TRUE, FALSE), na.last = TRUE, method = "radix")
  out <- data.frame(ProteinID = all_ids[ord], SelectionFrequency = top50_freq[ord], RFTop20Frequency = top20_freq[ord], RFTop50Frequency = top50_freq[ord], MeanDecreaseGini = mean_gini[ord], MeanDecreaseGiniScaled = gini_scaled[ord], StabilityScore = score[ord], Repeats = repeats, stringsAsFactors = FALSE)
  data.table::fwrite(out, file.path(outdir, "random_forest_stability_selection.csv"))
  data.table::fwrite(dplyr::bind_rows(performance), file.path(outdir, "random_forest_stability_resample_performance.csv"))
  invisible(out)
}

write_l1_stability <- function(x, y, outdir, top_n = 50, l1_alpha = 1, lambda_selection = "lambda.1se", seed = 123, repeats = 50, sample_fraction = 0.8, top_var_n = 200, frequency_weight = 1, split_plan = NULL, cv_folds = "Auto") {
  selected <- list()
  performance <- vector("list", repeats)
  split_plan <- split_plan %||% make_stability_split_plan(y, repeats, sample_fraction, seed + 2000)
  for (i in seq_len(repeats)) {
    split <- split_plan[[i]]
    train_x <- x[split$train, , drop = FALSE]
    y_sub <- droplevels(y[split$train])
    test_x <- x[split$test, , drop = FALSE]
    test_y <- droplevels(y[split$test])
    if (min(table(y_sub)) < 2) next
    train_var <- apply(train_x, 2, var, na.rm = TRUE)
    train_var[!is.finite(train_var)] <- -Inf
    keep <- order(train_var, decreasing = TRUE)[seq_len(min(top_var_n, ncol(train_x)))]
    imputed <- impute_from_training(train_x[, keep, drop = FALSE], test_x[, keep, drop = FALSE])
    cv_settings <- resolve_l1_folds(y_sub, cv_folds, allow_small_sample = TRUE)
    foldid <- make_stratified_foldid(y_sub, cv_settings$nfolds, seed + i)
    fit <- fit_l1_cv(imputed$train, y_sub, l1_alpha, lambda_selection, foldid, grouped = FALSE)
    pred <- as.vector(predict(fit, imputed$test, s = lambda_selection, type = "class"))
    recalls <- diag(table(factor(test_y, levels = levels(y_sub)), factor(pred, levels = levels(y_sub)))) / rowSums(table(factor(test_y, levels = levels(y_sub)), factor(pred, levels = levels(y_sub))))
    performance[[i]] <- data.frame(Resample = i, TrainSamples = length(split$train), TestSamples = length(split$test), Accuracy = mean(pred == test_y), BalancedAccuracy = mean(recalls, na.rm = TRUE), RequestedCVFolds = cv_folds, EffectiveCVFolds = cv_settings$nfolds, stringsAsFactors = FALSE)
    co <- extract_l1_coefficients(fit, lambda_selection, levels(y_sub))
    ids <- unique(co$ProteinID[co$Coefficient != 0])
    ids <- setdiff(ids, "(Intercept)")
    selected[[i]] <- ids
  }
  freq <- sort(table(unlist(selected, use.names = FALSE)) / repeats, decreasing = TRUE)
  if (length(freq) == 0) {
    out <- data.frame(ProteinID = character(), SelectionFrequency = numeric(), L1SelectionFrequency = numeric(), StabilityScore = numeric(), Repeats = integer(), stringsAsFactors = FALSE)
  } else {
    out <- data.frame(ProteinID = names(freq), SelectionFrequency = as.numeric(freq), L1SelectionFrequency = as.numeric(freq) * frequency_weight, StabilityScore = as.numeric(freq) * frequency_weight, Repeats = repeats, stringsAsFactors = FALSE)
  }
  data.table::fwrite(out, file.path(outdir, "l1_stability_selection.csv"))
  data.table::fwrite(dplyr::bind_rows(performance), file.path(outdir, "l1_stability_resample_performance.csv"))
  invisible(out)
}

run_random_forest_selection <- function(mat, group_info, outdir, top_n = 50, rf_ntree = 500, seed = 123, split_mode = "auto", train_prop = 0.7, rf_mtry = NA, allow_small_sample = FALSE, stability_repeats = 50, stability_sample_fraction = 0.8, stability_top_var_n = 200, stability_top20_weight = 0.35, stability_top50_weight = 0.25, stability_gini_weight = 0.10, stability_split_plan = NULL) {
  stability <- getOption("proteopostz.stability", list())$rf %||% list()
  stability_repeats <- stability$repeats %||% stability_repeats
  stability_sample_fraction <- stability$sample_fraction %||% stability_sample_fraction
  stability_top_var_n <- stability$top_var_n %||% stability_top_var_n
  stability_top20_weight <- stability$top20_weight %||% stability_top20_weight
  stability_top50_weight <- stability$top50_weight %||% stability_top50_weight
  stability_gini_weight <- stability$gini_weight %||% stability_gini_weight
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ml <- prepare_ml_input(mat, group_info)
  x <- ml$x
  y <- ml$y
  class_n <- check_ml_sample_policy(y, "Random forest", allow_small_sample, strict_min_per_class = 6, exploratory_min_per_class = 2)
  training <- resolve_ml_training(y, split_mode, train_prop, seed, min_train_per_class = 2, min_test_per_class = 1, label = "Random forest", allow_small_sample = allow_small_sample, train_test_min_per_class = 8)
  mtry <- parse_auto_integer(rf_mtry, NA_integer_)
  if (!is.na(mtry)) mtry <- min(mtry, ncol(x))
  set.seed(seed)
  rf_args <- list(x = x[training$train, , drop = FALSE], y = droplevels(y[training$train]), ntree = rf_ntree, importance = TRUE)
  if (!is.na(mtry)) rf_args$mtry <- mtry
  rf <- do.call(randomForest::randomForest, rf_args)
  write_rf_evaluation(x, y, training, rf_args, rf, outdir, seed)
  imp_mat <- randomForest::importance(rf)
  importance_col <- if ("MeanDecreaseGini" %in% colnames(imp_mat)) "MeanDecreaseGini" else tail(colnames(imp_mat), 1)
  imp <- data.frame(ProteinID = rownames(imp_mat), RFImportance = imp_mat[, importance_col], row.names = NULL) |>
    dplyr::arrange(dplyr::desc(RFImportance))
  data.table::fwrite(imp, file.path(outdir, "random_forest_importance.csv"))
  stability_input <- prepare_stability_input(mat, group_info)
  stability <- write_rf_stability(stability_input$x, stability_input$y, outdir, top_n = top_n, rf_ntree = rf_ntree, mtry = mtry, seed = seed, repeats = stability_repeats, sample_fraction = stability_sample_fraction, top_var_n = stability_top_var_n, top20_weight = stability_top20_weight, top50_weight = stability_top50_weight, gini_weight = stability_gini_weight, split_plan = stability_split_plan)
  write_ml_settings(outdir, list(
    Analysis = "Random forest",
    RandomSeed = seed,
    SplitMode = training$mode,
    TrainingSetProportion = train_prop,
    SamplesPerGroup = format_class_counts(class_n),
    TrainingSamplesPerGroup = format_class_counts(training$train_class_n),
    TestSamplesPerGroup = if (is.null(training$test_class_n)) "not used" else format_class_counts(training$test_class_n),
    RandomForestNtree = rf_ntree,
    RandomForestMtry = if (is.na(mtry)) "Auto" else mtry,
    Importance = importance_col,
    StabilitySelectionRepeats = stability_repeats,
    StabilitySampleFraction = stability_sample_fraction,
    StabilityTopVarianceFeatures = stability_top_var_n,
    StabilityTop20Weight = stability_top20_weight,
    StabilityTop50Weight = stability_top50_weight,
    StabilityGiniWeight = stability_gini_weight,
    RFTop20FrequencyDefinition = "Top 20 ranked proteins per stability resample, independent of final output Top N",
    RFTop50FrequencyDefinition = "Top 50 ranked proteins per stability resample, independent of final output Top N",
    MeanGiniScalingUniverse = "All stability-eligible proteins",
    EvaluationOutputs = "confusion_matrix, class_metrics, summary_metrics, ROC/AUC for binary groups, stability_selection",
    SmallSampleExploratoryML = allow_small_sample,
    ReliabilityNote = if (allow_small_sample) "Exploratory only: no independent test set; feature selection may be unstable." else "",
    AutoNote = training$auto_note %||% ""
  ))
  top <- head(stability$ProteinID, min(top_n, nrow(stability)))
  data.table::fwrite(data.frame(ProteinID = top), file.path(outdir, paste0("top", length(top), "_rf_features.csv")))
  write_matrix_csv(mat[top, , drop = FALSE], file.path(outdir, paste0("top", length(top), "_rf_feature_quantity_matrix.csv")))
  top
}

make_stratified_foldid <- function(y, nfolds, seed = 123) {
  set.seed(seed)
  foldid <- integer(length(y))
  for (cls in levels(y)) {
    idx <- which(y == cls)
    idx <- sample(idx, length(idx))
    foldid[idx] <- rep(seq_len(nfolds), length.out = length(idx))
  }
  foldid
}

resolve_l1_folds <- function(y_train, requested = "Auto", allow_small_sample = FALSE) {
  class_n <- table(y_train)
  min_class_n <- min(class_n)
  req <- parse_auto_integer(requested, NA_integer_)
  nfolds <- if (is.na(req)) {
    min(5, min_class_n)
  } else {
    min(req, min_class_n)
  }
  min_supported <- if (allow_small_sample) 4 else 4
  if (nfolds < min_supported) {
    stop("L1 feature selection cross-validation requires at least ", min_supported, " folds and sufficient samples per group in the training data because glmnet requires nfolds > 3. Current training group counts: ", format_class_counts(class_n), ".")
  }
  fold_size <- floor(length(y_train) / nfolds)
  list(nfolds = nfolds, grouped = !allow_small_sample && fold_size >= 3)
}

run_l1_selection <- function(mat, group_info, outdir, top_n = 50, l1_alpha = 1, seed = 123, split_mode = "auto", train_prop = 0.7, lambda_selection = "lambda.1se", cv_folds = "Auto", allow_small_sample = FALSE, stability_repeats = 50, stability_sample_fraction = 0.8, stability_top_var_n = 200, stability_frequency_weight = 1, stability_split_plan = NULL) {
  stability <- getOption("proteopostz.stability", list())$l1 %||% list()
  stability_repeats <- stability$repeats %||% stability_repeats
  stability_sample_fraction <- stability$sample_fraction %||% stability_sample_fraction
  stability_top_var_n <- stability$top_var_n %||% stability_top_var_n
  stability_frequency_weight <- stability$frequency_weight %||% stability_frequency_weight
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  lambda_selection <- match.arg(lambda_selection, c("lambda.1se", "lambda.min"))
  ml <- prepare_ml_input(mat, group_info)
  x <- ml$x
  y <- ml$y
  class_n <- check_ml_sample_policy(y, "L1 feature selection", allow_small_sample, strict_min_per_class = 6, exploratory_min_per_class = 4)
  training <- resolve_ml_training(y, split_mode, train_prop, seed, min_train_per_class = 3, min_test_per_class = 1, label = "L1 feature selection", allow_small_sample = allow_small_sample, train_test_min_per_class = 8)
  x_train <- x[training$train, , drop = FALSE]
  y_train <- droplevels(y[training$train])
  cv_settings <- resolve_l1_folds(y_train, cv_folds, allow_small_sample)
  nfolds <- cv_settings$nfolds
  foldid <- make_stratified_foldid(y_train, nfolds, seed)
  set.seed(seed)
  glmnet_warnings <- character()
  cv_expr <- quote(fit_l1_cv(x_train, y_train, l1_alpha, lambda_selection, foldid, grouped = cv_settings$grouped))
  cv <- if (allow_small_sample) {
    withCallingHandlers(eval(cv_expr), warning = function(w) {
      glmnet_warnings <<- c(glmnet_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  } else {
    eval(cv_expr)
  }
  write_l1_evaluation(x, y, training, l1_alpha, lambda_selection, nfolds, outdir, seed)
  coef_df <- extract_l1_coefficients(cv, lambda_selection, levels(y_train)) |>
    dplyr::filter(ProteinID != "(Intercept)", Coefficient != 0)
  data.table::fwrite(coef_df, file.path(outdir, "l1_nonzero_coefficients.csv"))
  scores <- coef_df |>
    dplyr::group_by(ProteinID) |>
    dplyr::summarise(L1Score = sum(abs(Coefficient)), NonzeroClasses = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(L1Score))
  data.table::fwrite(scores, file.path(outdir, "l1_feature_scores.csv"))
  stability_input <- prepare_stability_input(mat, group_info)
  stability <- write_l1_stability(stability_input$x, stability_input$y, outdir, top_n = top_n, l1_alpha = l1_alpha, lambda_selection = lambda_selection, seed = seed, repeats = stability_repeats, sample_fraction = stability_sample_fraction, top_var_n = stability_top_var_n, frequency_weight = stability_frequency_weight, split_plan = stability_split_plan, cv_folds = cv_folds)
  top <- head(stability$ProteinID, min(top_n, nrow(stability)))
  data.table::fwrite(data.frame(ProteinID = top), file.path(outdir, paste0("top", length(top), "_l1_features.csv")))
  if (length(top) > 0) write_matrix_csv(mat[top, , drop = FALSE], file.path(outdir, paste0("top", length(top), "_l1_feature_quantity_matrix.csv")))
  write_ml_settings(outdir, list(
    Analysis = "L1 feature selection",
    RandomSeed = seed,
    SplitMode = training$mode,
    TrainingSetProportion = train_prop,
    SamplesPerGroup = format_class_counts(class_n),
    TrainingSamplesPerGroup = format_class_counts(training$train_class_n),
    TestSamplesPerGroup = if (is.null(training$test_class_n)) "not used" else format_class_counts(training$test_class_n),
    L1Alpha = l1_alpha,
    LambdaSelection = lambda_selection,
    L1Family = resolve_l1_family(y_train),
    CrossValidationFolds = nfolds,
    EffectiveCVFolds = nfolds,
    RequestedCrossValidationFolds = cv_folds,
    StabilityRequestedCrossValidationFolds = cv_folds,
    StabilityEffectiveCVFolds = "Per resample; see l1_stability_resample_performance.csv",
    CVFoldPolicy = "Auto: min(5, minimum training class size); explicit N: requested N capped by the current training/resample minimum class size, subject to glmnet minimum folds",
    StabilityCVFoldPolicy = "Auto: min(5, minimum resample training class size); explicit N: requested N capped by each resample minimum class size, subject to glmnet minimum folds",
    GroupedCV = cv_settings$grouped,
    StabilitySelectionRepeats = stability_repeats,
    StabilitySampleFraction = stability_sample_fraction,
    StabilityTopVarianceFeatures = stability_top_var_n,
    StabilityFrequencyWeight = stability_frequency_weight,
    EvaluationOutputs = "confusion_matrix, class_metrics, summary_metrics, ROC/AUC for binary groups, stability_selection",
    SmallSampleExploratoryML = allow_small_sample,
    GlmnetWarnings = if (length(glmnet_warnings) > 0) paste(unique(glmnet_warnings), collapse = " | ") else "",
    ReliabilityNote = if (allow_small_sample) "Exploratory only: no independent test set; feature selection may be unstable." else "",
    AutoNote = training$auto_note %||% ""
  ))
  top
}

run_feature_selection <- function(mat, group_info, outdir, top_n = 50, rf_ntree = 500, l1_alpha = 1, seed = 123, split_mode = "auto", train_prop = 0.7, rf_mtry = NA, lambda_selection = "lambda.1se", cv_folds = "Auto", allow_small_sample = FALSE, stability_repeats = 50, stability_sample_fraction = 0.8, stability_top_var_n = 200, stability_rf_top20_weight = 0.35, stability_rf_top50_weight = 0.25, stability_l1_weight = 0.30, stability_gini_weight = 0.10) {
  stability <- getOption("proteopostz.stability", list())$rfl1 %||% list()
  stability_repeats <- stability$repeats %||% stability_repeats
  stability_sample_fraction <- stability$sample_fraction %||% stability_sample_fraction
  stability_top_var_n <- stability$top_var_n %||% stability_top_var_n
  stability_rf_top20_weight <- stability$rf_top20_weight %||% stability_rf_top20_weight
  stability_rf_top50_weight <- stability$rf_top50_weight %||% stability_rf_top50_weight
  stability_l1_weight <- stability$l1_weight %||% stability_l1_weight
  stability_gini_weight <- stability$gini_weight %||% stability_gini_weight
  prior_stability <- getOption("proteopostz.stability", list())
  options(proteopostz.stability = list(
    rf = list(repeats = stability_repeats, sample_fraction = stability_sample_fraction, top_var_n = stability_top_var_n),
    l1 = list(repeats = stability_repeats, sample_fraction = stability_sample_fraction, top_var_n = stability_top_var_n)
  ))
  on.exit(options(proteopostz.stability = prior_stability), add = TRUE)
  stability_input <- prepare_stability_input(mat, group_info)
  shared_splits <- make_stability_split_plan(stability_input$y, stability_repeats, stability_sample_fraction, seed + 1000)
  rf_top <- run_random_forest_selection(mat, group_info, outdir, top_n, rf_ntree, seed, split_mode, train_prop, rf_mtry, allow_small_sample, stability_repeats, stability_sample_fraction, stability_top_var_n, stability_rf_top20_weight, stability_rf_top50_weight, stability_gini_weight, shared_splits)
  if (length(rf_top) < 2) stop("RF + L1 combined requires at least 2 RF-selected candidate proteins before running the L1 stage.")
  l1_top <- run_l1_selection(mat, group_info, outdir, top_n, l1_alpha, seed, split_mode, train_prop, lambda_selection, cv_folds, allow_small_sample, stability_repeats, stability_sample_fraction, stability_top_var_n, stability_l1_weight, shared_splits)
  rf_stability <- data.table::fread(file.path(outdir, "random_forest_stability_selection.csv"), data.table = FALSE)
  l1_stability <- data.table::fread(file.path(outdir, "l1_stability_selection.csv"), data.table = FALSE)
  ids <- unique(c(rf_stability$ProteinID, l1_stability$ProteinID))
  combined <- data.frame(ProteinID = ids, stringsAsFactors = FALSE) |>
    dplyr::left_join(rf_stability[, c("ProteinID", "RFTop20Frequency", "RFTop50Frequency", "MeanDecreaseGini", "MeanDecreaseGiniScaled")], by = "ProteinID") |>
    dplyr::left_join(l1_stability[, c("ProteinID", "L1SelectionFrequency")], by = "ProteinID") |>
    dplyr::mutate(dplyr::across(-ProteinID, ~ tidyr::replace_na(.x, 0)), StabilityScore = stability_rf_top20_weight * RFTop20Frequency + stability_rf_top50_weight * RFTop50Frequency + stability_l1_weight * L1SelectionFrequency + stability_gini_weight * MeanDecreaseGiniScaled) |>
    dplyr::filter(RFTop20Frequency > 0 | RFTop50Frequency > 0 | L1SelectionFrequency > 0 | MeanDecreaseGini > 0) |>
    dplyr::arrange(dplyr::desc(StabilityScore), dplyr::desc(RFTop20Frequency), dplyr::desc(RFTop50Frequency), dplyr::desc(L1SelectionFrequency), dplyr::desc(MeanDecreaseGiniScaled), ProteinID)
  data.table::fwrite(combined, file.path(outdir, "combined_stability_scores.csv"))
  top <- head(combined$ProteinID, min(top_n, nrow(combined)))
  data.table::fwrite(data.frame(ProteinID = top, CombinedRank = seq_along(top), stringsAsFactors = FALSE), file.path(outdir, paste0("top", length(top), "_rf_l1_union_features.csv")))
  data.table::fwrite(data.frame(ProteinID = top, CombinedRank = seq_along(top), StabilityScore = combined$StabilityScore[match(top, combined$ProteinID)], stringsAsFactors = FALSE), file.path(outdir, "combined_feature_summary.csv"))
  if (length(top) > 0) write_matrix_csv(mat[top, , drop = FALSE], file.path(outdir, paste0("top", length(top), "_rf_l1_union_quantity_matrix.csv")))
  write_ml_settings(outdir, list(
    Analysis = "RF + L1 combined feature selection",
    RandomSeed = seed,
    TrainingSetProportion = train_prop,
    StabilitySelectionRepeats = stability_repeats,
    StabilitySampleFraction = stability_sample_fraction,
    StabilityTopVarianceFeatures = stability_top_var_n,
    RFTop20Weight = stability_rf_top20_weight,
    RFTop50Weight = stability_rf_top50_weight,
    L1FrequencyWeight = stability_l1_weight,
    MeanGiniWeight = stability_gini_weight,
    RFNtree = rf_ntree,
    RFMtry = if (is.na(parse_auto_integer(rf_mtry, NA_integer_))) "Auto" else parse_auto_integer(rf_mtry, NA_integer_),
    L1Alpha = l1_alpha,
    LambdaSelection = lambda_selection,
    RequestedCVFolds = cv_folds,
    L1Family = resolve_l1_family(stability_input$y),
    FinalTopN = top_n,
    CandidateRule = "At least one RF Top20, RF Top50, L1 frequency, or positive MeanGini contribution"
  ))
  top
}
