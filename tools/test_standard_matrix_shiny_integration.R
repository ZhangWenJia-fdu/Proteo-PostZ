coalesce_null <- function(a, b) if (!is.null(a)) a else b
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
script_path <- if (length(file_arg) > 0) normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE) else normalizePath(coalesce_null(sys.frame(1)$ofile, 'tools/test_standard_matrix_shiny_integration.R'), winslash = '/', mustWork = FALSE)
package_root <- normalizePath(file.path(dirname(script_path), '..'), winslash = '/', mustWork = TRUE)
app_root <- file.path(package_root, 'app')
setwd(app_root)
source('app.R', local = TRUE)
source(file.path(package_root, 'tools', 'test_input_paths.R'))

tmp <- file.path(package_root, 'outputs', 'standard_matrix_shiny_test')
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
write_utf8 <- function(path, lines) writeLines(lines, path, useBytes = TRUE)
expect_true <- function(x, label) if (!isTRUE(x)) stop(label)

csv <- file.path(tmp, 'ui_standard.csv')
write_utf8(csv, c(
  'Feature,S1,S2,S3,S4',
  'F1,1,0,2,3',
  'F2,4,5,0,6',
  'F3,7,,8,9',
  'F4,10,11,12,13',
  'F5,14,NA,15,16'
))
tsv <- file.path(tmp, 'ui_standard.tsv')
write_utf8(tsv, c(
  paste('Any ID', 'A-1', 'B-1', 'C-1', 'D-1', sep = '\t'),
  paste('G1', '1', '0', '2', '3', sep = '\t'),
  paste('G2', '4', '5', '0', '6', sep = '\t'),
  paste('G3', '7', '', '8', '9', sep = '\t'),
  paste('G4', '10', '11', '12', '13', sep = '\t'),
  paste('G5', '14', 'NaN', '15', '16', sep = '\t')
))
bad <- file.path(tmp, 'bad_standard.csv')
write_utf8(bad, c('Feature,S1,S2', 'F1,abc,1'))

shiny::testServer(server, {
  session$setInputs(input_family = 'standard_matrix', file_path = csv, outdir = tmp, standard_zero_mode = '')
  session$setInputs(load_data = 1)
  session$flushReact()
  expect_true(is.null(rv$data), 'standard matrix loaded without zero mode')
  expect_true(grepl('select how zero', rv$load_error), 'missing zero-mode error not shown')

  session$setInputs(standard_zero_mode = 'zero_is_value')
  session$setInputs(load_data = 2)
  session$flushReact()
  expect_true(identical(rv$data$input_source, 'standard_matrix'), 'standard CSV did not load')
  expect_true(identical(rv$data$samples, c('S1', 'S2', 'S3', 'S4')), 'sample names/order not from matrix columns')
  expect_true(rv$data$available_quantitative_value_count == 18, 'zero_is_value available count incorrect')
  expect_true(!is.na(rv$data$quantity['F1', 'S2']), 'zero_is_value did not retain zero')
  expect_true(identical(names(rv$groups), rv$data$samples), 'group control names not refreshed from standard samples')

  session$setInputs(standard_zero_mode = 'zero_is_missing')
  session$setInputs(load_data = 3)
  session$flushReact()
  expect_true(rv$data$available_quantitative_value_count == 16, 'zero_is_missing available count incorrect')
  expect_true(is.na(rv$data$quantity['F1', 'S2']), 'zero_is_missing did not remove zero from analysis matrix')
  expect_true(!is.na(rv$data$raw_quant_matrix['F1', 'S2']), 'zero_is_missing overwrote raw matrix')

  session$setInputs(file_path = tsv, standard_zero_mode = 'zero_is_value')
  session$setInputs(load_data = 4)
  session$flushReact()
  expect_true(identical(rv$data$feature_column_name, 'Any ID'), 'TSV first column name not shown/preserved')
  expect_true(identical(rv$data$samples, c('A-1', 'B-1', 'C-1', 'D-1')), 'TSV sample names/order incorrect')
  expect_true(rv$data$available_quantitative_value_count == 18, 'standard TSV zero_is_value available count incorrect')

  session$setInputs(standard_zero_mode = 'zero_is_missing')
  session$setInputs(load_data = 5)
  session$flushReact()
  expect_true(rv$data$available_quantitative_value_count == 16, 'standard TSV zero_is_missing available count incorrect')

  session$setInputs(input_family = 'dia', dia_format = 'diann')
  session$flushReact()
  expect_true(is.null(rv$data), 'switching input type did not clear loaded standard matrix')
  session$setInputs(input_family = 'standard_matrix', file_path = bad, standard_zero_mode = 'zero_is_value')
  session$setInputs(load_data = 6)
  session$flushReact()
  expect_true(is.null(rv$data), 'bad standard matrix should not load')
  expect_true(grepl('non-numeric text', rv$load_error), 'bad standard matrix error was not shown cleanly')
})

source('R/analysis_core.R')
diann <- extract_protein_data(resolve_external_test_file('PROTEOPOSTZ_DIANN_TEST_FILE'), 'DIANN', 'd', 'protein_name')
diann_group <- make_group_info(diann$samples, rep(c('G1', 'G2'), each = 4))
plot_rank_abundance(diann$quantity, diann_group, file.path(tmp, 'diann_rank.pdf'), file.path(tmp, 'diann_rank.csv'), 3, 3)
expect_true(file.exists(file.path(tmp, 'diann_rank.pdf')), 'DIA-NN rank module did not produce PDF')

spectronaut <- extract_protein_data(resolve_external_test_file('PROTEOPOSTZ_SPECTRONAUT_TEST_FILE'), 'Spectronaut', 'd', 'protein_name')
spectronaut_group <- make_group_info(spectronaut$samples, rep(c('G1', 'G2'), each = 4))
plot_rank_abundance(spectronaut$quantity, spectronaut_group, file.path(tmp, 'spectronaut_rank.pdf'), file.path(tmp, 'spectronaut_rank.csv'), 3, 3)
expect_true(file.exists(file.path(tmp, 'spectronaut_rank.pdf')), 'Spectronaut rank module did not produce PDF')

cat('STANDARD_MATRIX_SHINY_TEST_OK csv_value=18 csv_missing=16 tsv_value=18 tsv_missing=16 diann=', nrow(diann$quantity), 'x', ncol(diann$quantity), ' spectronaut=', nrow(spectronaut$quantity), 'x', ncol(spectronaut$quantity), '\n', sep = '')
