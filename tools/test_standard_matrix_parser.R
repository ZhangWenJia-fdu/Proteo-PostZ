coalesce_null <- function(a, b) if (!is.null(a)) a else b
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
script_path <- if (length(file_arg) > 0) normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE) else normalizePath(coalesce_null(sys.frame(1)$ofile, 'tools/test_standard_matrix_parser.R'), winslash = '/', mustWork = FALSE)
package_root <- normalizePath(file.path(dirname(script_path), '..'), winslash = '/', mustWork = TRUE)
app_root <- file.path(package_root, 'app')
setwd(app_root)
source('R/analysis_core.R')

tmp <- file.path(package_root, 'outputs', 'standard_matrix_parser_test')
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
write_utf8 <- function(path, lines) writeLines(lines, path, useBytes = TRUE)
expect_error <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  if (is.null(err)) stop('Expected error matching ', pattern, ' but no error was raised.')
  if (!grepl(pattern, err, ignore.case = TRUE)) stop('Expected error matching ', pattern, ' but got: ', err)
  invisible(err)
}
expect_true <- function(x, label) if (!isTRUE(x)) stop(label)

csv <- file.path(tmp, 'standard_matrix.csv')
write_utf8(csv, c(
  paste0('\ufefffeature id,control A,treat-1,Sci.Sample'),
  'ProteinA,123.4,0,1.2e3',
  'ProteinB,85.1, ,NA',
  'ProteinC,0,41.7,NaN',
  'ProteinD,-2.5,3,4'
))
tsv <- file.path(tmp, 'standard_matrix.tsv')
write_utf8(tsv, c(
  paste('gene symbol', 'sample alpha', 'sample beta', sep = '\t'),
  paste('GeneA', '1', '0', sep = '\t'),
  paste('GeneB', '', '-3.2', sep = '\t'),
  paste('GeneC', '5e-1', 'NaN', sep = '\t')
))

value_mode <- read_standard_matrix(csv, 'zero_is_value')
missing_mode <- read_standard_matrix(csv, 'zero_is_missing')
tsv_mode <- read_standard_matrix(tsv, 'zero_is_missing')
expect_true(identical(value_mode$input_source, 'standard_matrix'), 'input_source not recorded')
expect_true(identical(value_mode$feature_column_name, 'feature id'), 'feature column name not preserved')
expect_true(identical(value_mode$samples, c('control A', 'treat-1', 'Sci.Sample')), 'sample names/order not preserved')
expect_true(identical(rownames(value_mode$raw_quant_matrix), c('ProteinA', 'ProteinB', 'ProteinC', 'ProteinD')), 'feature order not preserved')
expect_true(value_mode$raw_quant_matrix['ProteinD', 'control A'] == -2.5, 'negative values not parsed')
expect_true(value_mode$raw_quant_matrix['ProteinA', 'Sci.Sample'] == 1200, 'scientific notation not parsed')
expect_true(value_mode$raw_zero_cell_count == 2, 'zero count incorrect')
expect_true(value_mode$zero_to_missing_count == 0, 'zero_is_value converted zeros')
expect_true(missing_mode$zero_to_missing_count == 2, 'zero_is_missing conversion count incorrect')
expect_true(value_mode$available_quantitative_value_count == 9, 'available value count for zero_is_value incorrect')
expect_true(missing_mode$available_quantitative_value_count == 7, 'available value count for zero_is_missing incorrect')
expect_true(!is.na(missing_mode$raw_quant_matrix['ProteinA', 'treat-1']), 'raw matrix was overwritten in zero_is_missing mode')
expect_true(is.na(missing_mode$analysis_quant_matrix['ProteinA', 'treat-1']), 'analysis matrix did not convert zero to NA')
expect_true(value_mode$missingness_matrix['ProteinB', 'treat-1'], 'blank value missingness not detected')
expect_true(value_mode$missingness_matrix['ProteinB', 'Sci.Sample'], 'NA missingness not detected')
expect_true(value_mode$missingness_matrix['ProteinC', 'Sci.Sample'], 'NaN missingness not detected')
expect_true(tsv_mode$raw_quant_matrix['GeneC', 'sample alpha'] == 0.5, 'TSV scientific notation not parsed')
expect_true(tsv_mode$raw_quant_matrix['GeneB', 'sample beta'] == -3.2, 'TSV negative value not parsed')

write_utf8(file.path(tmp, 'one_col.csv'), c('feature', 'A', 'B'))
expect_error(read_standard_matrix(file.path(tmp, 'one_col.csv'), 'zero_is_value'), 'one column|delimiter')
write_utf8(file.path(tmp, 'duplicate_samples.csv'), c('feature,s1,s1', 'A,1,2'))
expect_error(read_standard_matrix(file.path(tmp, 'duplicate_samples.csv'), 'zero_is_value'), 'sample column names must be unique')
write_utf8(file.path(tmp, 'blank_feature.csv'), c('feature,s1', ' ,1'))
expect_error(read_standard_matrix(file.path(tmp, 'blank_feature.csv'), 'zero_is_value'), 'feature identifiers cannot be blank|entirely empty')
write_utf8(file.path(tmp, 'duplicate_feature.csv'), c('feature,s1', 'A,1', 'A,2'))
expect_error(read_standard_matrix(file.path(tmp, 'duplicate_feature.csv'), 'zero_is_value'), 'feature identifiers must be unique')
write_utf8(file.path(tmp, 'bad_text.csv'), c('feature,s1', 'A,abc'))
expect_error(read_standard_matrix(file.path(tmp, 'bad_text.csv'), 'zero_is_value'), 'non-numeric text.*s1.*abc')
write_utf8(file.path(tmp, 'pos_inf.csv'), c('feature,s1', 'A,Inf'))
expect_error(read_standard_matrix(file.path(tmp, 'pos_inf.csv'), 'zero_is_value'), 'Inf or -Inf.*s1')
write_utf8(file.path(tmp, 'neg_inf.csv'), c('feature,s1', 'A,-Inf'))
expect_error(read_standard_matrix(file.path(tmp, 'neg_inf.csv'), 'zero_is_value'), 'Inf or -Inf.*s1')
write_utf8(file.path(tmp, 'all_missing.csv'), c('feature,s1,s2', 'A,0,1', 'B,0,2'))
expect_error(read_standard_matrix(file.path(tmp, 'all_missing.csv'), 'zero_is_missing'), 'entirely missing.*s1')
write_utf8(file.path(tmp, 'wrong_delimiter.csv'), c('feature\ts1', 'A\t1'))
expect_error(read_standard_matrix(file.path(tmp, 'wrong_delimiter.csv'), 'zero_is_value'), 'extension and delimiter|delimiter matches')
expect_error(read_standard_matrix(csv), 'zero handling mode')

cat('STANDARD_MATRIX_TEST_OK', value_mode$available_quantitative_value_count, missing_mode$available_quantitative_value_count, '\n')