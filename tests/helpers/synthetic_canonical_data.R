# Deterministic synthetic truth used by regression groups.
with_local_seed <- function(seed, code) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  on.exit(if (is.null(old)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  force(code)
}
make_canonical_fixture <- function(seed = 101, n_features = 80, n_samples = 12, groups = NULL) {
  with_local_seed(seed, {
    groups <- groups %||% rep(c("A", "B"), length.out = n_samples)
    samples <- paste0("Sample", seq_len(n_samples)); ids <- paste0("P", sprintf("%03d", seq_len(n_features)))
    mat <- matrix(exp(rnorm(n_features * n_samples, mean = 4, sd = 0.7)), nrow = n_features, dimnames = list(ids, samples))
    mat[1, ] <- 0; mat[2, 1] <- NA_real_; mat[3, groups == "A"] <- NA_real_; mat[4, groups == "B"] <- NA_real_; mat[5, c(1, 2)] <- NA_real_; mat[6, ] <- NA_real_; mat[7, 1] <- 0
    meta <- data.frame(ProteinID = ids, ProteinName = paste0("Protein_", ids), Gene = paste0("Gene_", ids), stringsAsFactors = FALSE)
    list(quantity = mat, qualitative = mat, samples = samples, groups = groups, group_info = make_group_info(samples, groups), meta = meta)
  })
}
make_binary_ml_fixture <- function(seed = 202, n_features = 320, n_per_group = 12) {
  with_local_seed(seed, { groups <- rep(c("A", "B"), each = n_per_group); f <- make_canonical_fixture(seed + 1, n_features, length(groups), groups); signal <- matrix(rnorm(20 * length(groups), sd = 0.15), nrow = 20); signal[1:10, groups == "B"] <- signal[1:10, groups == "B"] + 2.4; signal[11:20, groups == "B"] <- signal[11:20, groups == "B"] + 1; f$quantity[1:20, ] <- exp(signal + 4); f$qualitative <- f$quantity; f })
}
make_multiclass_fixture <- function(seed = 303, n_features = 320, n_per_group = 10) {
  with_local_seed(seed, { groups <- rep(c("A", "B", "C"), each = n_per_group); f <- make_canonical_fixture(seed + 1, n_features, length(groups), groups); for (k in seq_len(3)) { idx <- ((k - 1) * 10 + 1):(k * 10); f$quantity[idx, groups == c("A", "B", "C")[[k]]] <- f$quantity[idx, groups == c("A", "B", "C")[[k]]] * 8 }; f$qualitative <- f$quantity; f })
}
make_trajectory_fixture <- function(seed = 404, n_features = 260, n_per_stage = 6) {
  with_local_seed(seed, { groups <- rep(paste0("Stage", 1:4), each = n_per_stage); f <- make_canonical_fixture(seed + 1, n_features, length(groups), groups); stage <- rep(seq_len(4), each = n_per_stage); f$quantity[1, ] <- 0; f$quantity[2, ] <- 10; f$quantity[3, ] <- exp(3 + stage / 2); f$quantity[4, ] <- exp(5 - stage / 2); f$quantity[5, ] <- exp(4 + ifelse(stage %in% c(2, 3), 1.5, 0)); f$qualitative <- f$quantity; f })
}
