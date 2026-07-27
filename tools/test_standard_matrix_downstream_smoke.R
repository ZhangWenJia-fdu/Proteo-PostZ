coalesce_null <- function(a, b) if (!is.null(a)) a else b
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
script_path <- if (length(file_arg) > 0) normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE) else normalizePath(coalesce_null(sys.frame(1)$ofile, 'tools/test_standard_matrix_downstream_smoke.R'), winslash = '/', mustWork = FALSE)
package_root <- normalizePath(file.path(dirname(script_path), '..'), winslash = '/', mustWork = TRUE)
app_root <- file.path(package_root, 'app')
setwd(app_root)
source('R/analysis_core.R')

tmp <- file.path(package_root, 'outputs', 'standard_matrix_downstream_smoke')
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
write_utf8 <- function(path, lines) writeLines(lines, path, useBytes = TRUE)
expect_true <- function(x, label) if (!isTRUE(x)) stop(label)

std <- file.path(tmp, 'standard_downstream.csv')
write_utf8(std, c(
  'Feature,S1,S2,S3,S4',
  'F1,1,0,2,3',
  'F2,4,5,0,6',
  'F3,7,,8,9',
  'F4,10,11,12,13',
  'F5,14,NA,15,16',
  'F6,17,18,19,20'
))
mode_a <- read_standard_matrix(std, 'zero_is_value')
mode_b <- read_standard_matrix(std, 'zero_is_missing')
expect_true(mode_a$available_quantitative_value_count == 22, 'mode A available count incorrect')
expect_true(mode_b$available_quantitative_value_count == 20, 'mode B available count incorrect')
expect_true(mode_a$quantity['F1', 'S2'] == 0, 'mode A did not retain zero in analysis matrix')
expect_true(is.na(mode_b$quantity['F1', 'S2']), 'mode B did not remove zero from analysis matrix')
expect_true(mode_b$raw_quant_matrix['F1', 'S2'] == 0, 'mode B raw matrix did not preserve zero')

group_info <- make_group_info(mode_b$samples, c('G1', 'G1', 'G2', 'G2'))
plot_available_quantitative_bar(mode_b$counts, group_info, file.path(tmp, 'available_quantitative_values_barplot.pdf'), file.path(tmp, 'available_quantitative_values_group_summary.csv'), 3, 3)
plot_correlation_heatmap(mode_b$quantity, group_info, file.path(tmp, 'correlation.pdf'), file.path(tmp, 'correlation.csv'), width = 3, height = 3)
used <- preprocess_expr(mode_b$quantity, TRUE, 0.5)
expect_true(!is.nan(used['F1', 'S2']), 'PCA preprocessing produced NaN from a zero that should have been NA in mode B')
sample_mat <- t(used)
pca <- prcomp(sample_mat, center = TRUE, scale. = TRUE)
expect_true(nrow(pca$x) == 4, 'PCA did not use all samples')
sets <- identified_by_group(mode_b$quantity, group_info, 1)
set_df <- make_set_membership(sets)
colnames(set_df)[colnames(set_df) == 'ProteinID'] <- 'FeatureID'
data.table::fwrite(set_df, file.path(tmp, 'quantified_feature_membership.csv'))
expect_true('FeatureID' %in% colnames(set_df), 'standard matrix membership export did not use FeatureID')
expect_true(all(file.exists(file.path(tmp, c('available_quantitative_values_barplot.pdf', 'correlation.pdf', 'correlation.csv', 'quantified_feature_membership.csv')))), 'standard matrix downstream smoke outputs missing')

cat('STANDARD_MATRIX_DOWNSTREAM_SMOKE_OK mode_a=22 mode_b=20 zero_mode_b_is_na=', is.na(mode_b$quantity['F1', 'S2']), '\n', sep = '')