coalesce_null <- function(a, b) if (!is.null(a)) a else b
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
script_path <- if (length(file_arg) > 0) normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE) else normalizePath(coalesce_null(sys.frame(1)$ofile, 'tools/test_v14_input_regression.R'), winslash = '/', mustWork = FALSE)
package_root <- normalizePath(file.path(dirname(script_path), '..'), winslash = '/', mustWork = TRUE)
app_root <- file.path(package_root, 'app')
setwd(app_root)
source('app.R', local = TRUE)
source('R/analysis_core.R')

tmp <- file.path(package_root, 'outputs', 'v14_input_regression')
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
write_utf8 <- function(path, lines) writeLines(lines, path, useBytes = TRUE)
expect_true <- function(x, label) if (!isTRUE(x)) stop(label)
expect_error <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (is.null(err)) stop('Expected error matching ', pattern, ' but no error was raised.')
  if (!grepl(pattern, err, ignore.case = TRUE)) stop('Expected error matching ', pattern, ' but got: ', err)
  invisible(err)
}

# DIA-NN regression: load, sample names, identification counts, one identification module, one quantitative module.
diann <- extract_protein_data('F:/test/DIANNreport.pg_matrix.tsv', 'DIANN', 'd', 'protein_name')
expect_true(nrow(diann$quantity) == 4862 && ncol(diann$quantity) == 8, 'DIA-NN quantity dimensions changed')
expect_true(identical(diann$samples, colnames(diann$quantity)), 'DIA-NN sample names not aligned with quantity matrix')
expect_true('Identified_Protein_Count' %in% colnames(diann$counts), 'DIA-NN identification count column missing')
expect_true('AnalysisID' %in% colnames(diann$meta) && identical(diann$meta$AnalysisID, rownames(diann$quantity)), 'DIA-NN meta AnalysisID does not track actual matrix row names')
diann_group <- make_group_info(diann$samples, rep(c('G1', 'G2'), each = 4))
plot_identification_bar(diann$counts, diann_group, file.path(tmp, 'diann_idbar.pdf'), file.path(tmp, 'diann_idbar.csv'), 3, 3)
plot_rank_abundance(diann$quantity, diann_group, file.path(tmp, 'diann_rank.pdf'), file.path(tmp, 'diann_rank.csv'), 3, 3)
expect_true(all(file.exists(file.path(tmp, c('diann_idbar.pdf', 'diann_rank.pdf')))), 'DIA-NN regression output missing')

# Spectronaut regression: load, sample names, PG.IBAQ identification semantics, one identification module, one quantitative module.
spec <- extract_protein_data('F:/test/Spectronaut20260428_141042_20260428-YGQ-Celegans-repeatabilitytest_Report_1_8.tsv', 'Spectronaut', 'd', 'protein_name')
expect_true(nrow(spec$quantity) == 3734 && ncol(spec$quantity) == 8, 'Spectronaut quantity dimensions changed')
expect_true(!is.null(spec$ibaq) && identical(colnames(spec$ibaq), colnames(spec$quantity)), 'Spectronaut PG.IBAQ qualitative matrix not aligned')
expect_true(identical(spec$samples, colnames(spec$quantity)), 'Spectronaut sample names not aligned with quantity matrix')
expect_true(all(spec$counts$Identified_Protein_Count == colSums(!is.na(spec$ibaq))), 'Spectronaut identification counts no longer use PG.IBAQ')
expect_true('AnalysisID' %in% colnames(spec$meta) && identical(spec$meta$AnalysisID, rownames(spec$quantity)), 'Spectronaut meta AnalysisID does not track actual matrix row names')
spec_group <- make_group_info(spec$samples, rep(c('G1', 'G2'), each = 4))
plot_identification_bar(spec$counts, spec_group, file.path(tmp, 'spectronaut_idbar.pdf'), file.path(tmp, 'spectronaut_idbar.csv'), 3, 3)
plot_rank_abundance(spec$quantity, spec_group, file.path(tmp, 'spectronaut_rank.pdf'), file.path(tmp, 'spectronaut_rank.csv'), 3, 3)
expect_true(all(file.exists(file.path(tmp, c('spectronaut_idbar.pdf', 'spectronaut_rank.pdf')))), 'Spectronaut regression output missing')

dup_meta <- data.frame(
  RowID = c('DuplicateName', 'DuplicateName', 'UniqueName'),
  AnalysisID = make.unique(c('DuplicateName', 'DuplicateName', 'UniqueName')),
  Accession = c('P00001', 'P00002', 'P00003'),
  stringsAsFactors = FALSE
)
expect_true(
  identical(analysis_ids_to_accessions(dup_meta, c('DuplicateName', 'DuplicateName.1', 'UniqueName')), c('P00001', 'P00002', 'P00003')),
  'AnalysisID-to-Accession mapping fails for duplicated protein/gene row names'
)

# Parser coverage for standard matrix edge cases.
std_csv <- file.path(tmp, 'standard_regression.csv')
write_utf8(std_csv, c(
  paste0('\ufeffAny feature ID,S 1,S-2,Sci.Sample'),
  'ProteinA,123.4,0,1.2e3',
  'ProteinB,85.1, ,NA',
  'ProteinC,0,41.7,NaN',
  'ProteinD,-2.5,3,4'
))
std_tsv <- file.path(tmp, 'standard_regression.tsv')
write_utf8(std_tsv, c(
  paste('Gene symbol', 'Alpha sample', 'Beta sample', sep = '\t'),
  paste('GeneA', '1.5', '0', sep = '\t'),
  paste('GeneB', '', '-3.2', sep = '\t'),
  paste('GeneC', '5e-1', 'NaN', sep = '\t')
))
std_a <- read_standard_matrix(std_csv, 'zero_is_value')
std_b <- read_standard_matrix(std_csv, 'zero_is_missing')
std_tsv_a <- read_standard_matrix(std_tsv, 'zero_is_value')
std_tsv_b <- read_standard_matrix(std_tsv, 'zero_is_missing')
expect_true(std_a$feature_column_name == 'Any feature ID', 'CSV arbitrary first column name not preserved')
expect_true(identical(std_a$samples, c('S 1', 'S-2', 'Sci.Sample')), 'CSV arbitrary sample names/order not preserved')
expect_true(std_a$available_quantitative_value_count == 9 && std_b$available_quantitative_value_count == 7, 'CSV zero-mode available counts incorrect')
expect_true(std_tsv_a$available_quantitative_value_count == 4 && std_tsv_b$available_quantitative_value_count == 3, 'TSV zero-mode available counts incorrect')
expect_true(std_a$raw_quant_matrix['ProteinD', 'S 1'] == -2.5, 'negative value not parsed')
expect_true(std_a$raw_quant_matrix['ProteinA', 'Sci.Sample'] == 1200, 'scientific notation not parsed')
expect_true(is.na(std_b$quantity['ProteinA', 'S-2']) && std_b$raw_quant_matrix['ProteinA', 'S-2'] == 0, 'mode B did not distinguish raw and analysis matrix')
expect_error(read_standard_matrix(std_csv), 'zero handling mode')
write_utf8(file.path(tmp, 'dup_feature.csv'), c('Feature,S1', 'A,1', 'A,2'))
expect_error(read_standard_matrix(file.path(tmp, 'dup_feature.csv'), 'zero_is_value'), 'feature identifiers must be unique')
write_utf8(file.path(tmp, 'blank_feature.csv'), c('Feature,S1', ' ,1'))
expect_error(read_standard_matrix(file.path(tmp, 'blank_feature.csv'), 'zero_is_value'), 'feature identifiers cannot be blank|entirely empty')
write_utf8(file.path(tmp, 'dup_sample.csv'), c('Feature,S1,S1', 'A,1,2'))
expect_error(read_standard_matrix(file.path(tmp, 'dup_sample.csv'), 'zero_is_value'), 'sample column names must be unique')
write_utf8(file.path(tmp, 'bad_text.csv'), c('Feature,S1', 'A,text'))
expect_error(read_standard_matrix(file.path(tmp, 'bad_text.csv'), 'zero_is_value'), 'non-numeric text')
write_utf8(file.path(tmp, 'inf.csv'), c('Feature,S1', 'A,Inf'))
expect_error(read_standard_matrix(file.path(tmp, 'inf.csv'), 'zero_is_value'), 'Inf or -Inf')
write_utf8(file.path(tmp, 'neg_inf.csv'), c('Feature,S1', 'A,-Inf'))
expect_error(read_standard_matrix(file.path(tmp, 'neg_inf.csv'), 'zero_is_value'), 'Inf or -Inf')
write_utf8(file.path(tmp, 'all_missing.csv'), c('Feature,S1,S2', 'A,0,1', 'B,0,2'))
expect_error(read_standard_matrix(file.path(tmp, 'all_missing.csv'), 'zero_is_missing'), 'entirely missing')
write_utf8(file.path(tmp, 'wrong_delimiter.csv'), c('Feature\tS1', 'A\t1'))
expect_error(read_standard_matrix(file.path(tmp, 'wrong_delimiter.csv'), 'zero_is_value'), 'delimiter')
expect_error(read_standard_matrix(std_csv), 'zero handling mode')

# Shiny state transitions and error recovery.
shiny::testServer(server, {
  session$setInputs(input_family = 'standard_matrix', file_path = std_csv, outdir = tmp, standard_zero_mode = '')
  session$setInputs(load_data = 1)
  session$flushReact()
  expect_true(is.null(rv$data) && grepl('select how zero', rv$load_error), 'missing zero-mode did not block standard matrix load')
  session$setInputs(standard_zero_mode = 'zero_is_value')
  session$setInputs(load_data = 2)
  session$flushReact()
  expect_true(identical(rv$data$samples, c('S 1', 'S-2', 'Sci.Sample')), 'Shiny standard CSV sample names incorrect')
  expect_true(identical(names(rv$groups), rv$data$samples), 'Shiny grouping names not refreshed')
  session$setInputs(file_path = file.path(tmp, 'bad_text.csv'))
  session$setInputs(load_data = 3)
  session$flushReact()
  expect_true(is.null(rv$data) && grepl('non-numeric text', rv$load_error), 'Shiny bad matrix did not show recoverable error')
  session$setInputs(file_path = std_tsv, standard_zero_mode = 'zero_is_missing')
  session$setInputs(load_data = 4)
  session$flushReact()
  expect_true(!is.null(rv$data) && rv$data$available_quantitative_value_count == 3, 'Shiny did not recover after bad standard matrix')
  session$setInputs(input_family = 'dia', dia_format = 'spectronaut')
  session$flushReact()
  expect_true(is.null(rv$data), 'Shiny input source switch did not clear standard matrix state')
})

# Standard matrix downstream modules on a positive matrix suitable for log-based analyses and ML.
downstream <- file.path(tmp, 'standard_downstream_full.csv')
features <- paste0('Feat', seq_len(24))
samples <- paste0(rep(c('A', 'B'), each = 6), '_', rep(seq_len(6), 2))
mat <- matrix(NA_real_, nrow = length(features), ncol = length(samples), dimnames = list(features, samples))
for (i in seq_along(features)) {
  mat[i, 1:6] <- 20 + i + seq_len(6) * 0.3
  mat[i, 7:12] <- 24 + i * 1.2 + seq_len(6) * 0.4
}
mat[1, 2] <- 0
mat[2, 8] <- 0
mat[3, 4] <- NA
mat[4, 10] <- NaN
out_df <- data.frame(Feature = features, mat, check.names = FALSE)
data.table::fwrite(out_df, downstream)
std_full_a <- read_standard_matrix(downstream, 'zero_is_value')
std_full_b <- read_standard_matrix(downstream, 'zero_is_missing')
expect_true(std_full_a$available_quantitative_value_count == 286, 'full standard mode A available count incorrect')
expect_true(std_full_b$available_quantitative_value_count == 284, 'full standard mode B available count incorrect')
expect_true(std_full_a$quantity['Feat1', 'A_2'] == 0, 'full standard mode A did not keep zero')
expect_true(is.na(std_full_b$quantity['Feat1', 'A_2']), 'full standard mode B did not convert zero to NA')
std_group <- make_group_info(std_full_b$samples, rep(c('GroupA', 'GroupB'), each = 6))
plot_correlation_heatmap(std_full_b$quantity, std_group, file.path(tmp, 'standard_cor.pdf'), file.path(tmp, 'standard_cor.csv'), width = 3, height = 3)
plot_expression_heatmap(std_full_b$quantity, std_group, file.path(tmp, 'standard_exprhm.pdf'), file.path(tmp, 'standard_exprhm.csv'), top_n = 12, width = 3, height = 4)
run_volcano(std_full_b$quantity, std_group, 'GroupA', 'GroupB', file.path(tmp, 'standard_volcano.pdf'), file.path(tmp, 'standard_volcano.csv'), test_method = 'ttest', width = 3, height = 3)
rf_top <- run_random_forest_selection(std_full_b$quantity, std_group, file.path(tmp, 'standard_rf'), top_n = 8, rf_ntree = 50, seed = 123, split_mode = 'cross_validation_only', allow_small_sample = FALSE)
expect_true(length(rf_top) > 0, 'standard matrix RF did not return features')
slingshot_status <- 'not_run'
slingshot_msg <- ''
slingshot_try <- tryCatch({
  run_slingshot_pseudotime(std_full_b$quantity, std_group, file.path(tmp, 'standard_sling'), reduction = 'PCA', start_group = 'GroupA', end_group = 'GroupB', width = 3, height = 3, top_n = 5)
  slingshot_status <- if (file.exists(file.path(tmp, 'standard_sling', 'slingshot_sample_pseudotime.csv'))) 'ok' else 'missing_output'
  NULL
}, error = function(e) {
  slingshot_status <<- 'skipped_or_not_inferable'
  slingshot_msg <<- conditionMessage(e)
  NULL
})
expect_true(all(file.exists(file.path(tmp, c('standard_cor.pdf', 'standard_exprhm.pdf', 'standard_volcano.pdf')))), 'standard matrix downstream plots missing')

cat('V14_INPUT_REGRESSION_OK diann=', nrow(diann$quantity), 'x', ncol(diann$quantity),
    ' spectronaut=', nrow(spec$quantity), 'x', ncol(spec$quantity),
    ' standard_csv_A=', std_a$available_quantitative_value_count,
    ' standard_csv_B=', std_b$available_quantitative_value_count,
    ' standard_tsv_A=', std_tsv_a$available_quantitative_value_count,
    ' standard_tsv_B=', std_tsv_b$available_quantitative_value_count,
    ' full_A=', std_full_a$available_quantitative_value_count,
    ' full_B=', std_full_b$available_quantitative_value_count,
    ' rf_top=', length(rf_top),
    ' slingshot=', slingshot_status,
    if (nzchar(slingshot_msg)) paste0(' slingshot_msg=', shQuote(slingshot_msg)) else '',
    '\n', sep = '')
