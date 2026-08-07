coalesce_null <- function(a, b) if (!is.null(a)) a else b
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(coalesce_null(sys.frame(1)$ofile, "tools/test_v201_count_modes.R"), winslash = "/", mustWork = FALSE)
}
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
app_root <- file.path(package_root, "app")
setwd(app_root)
source("app.R", local = TRUE)
source("R/analysis_core.R")
source(file.path(package_root, "tools", "test_input_paths.R"))

out <- file.path(package_root, "outputs", "v201_count_modes")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_utf8 <- function(path, lines) writeLines(lines, path, useBytes = TRUE)
expect_true <- function(x, label) if (!isTRUE(x)) stop(label)

fragpipe_path <- file.path(out, "fragpipe_test.tsv")
write_utf8(fragpipe_path, c(
  paste(
    c("Protein ID", "Gene", "Entry Name", "S1 MaxLFQ Intensity", "S2 MaxLFQ Intensity", "S1 Spectral Count", "S2 Spectral Count"),
    collapse = "\t"
  ),
  paste(c("P001", "Gene1", "Prot1", "0", "100", "3", "0"), collapse = "\t"),
  paste(c("P002", "Gene2", "Prot2", "25", "0", "0", "5"), collapse = "\t")
))
frag <- extract_fragpipe_protein_data(fragpipe_path, "accession")
expect_true(identical(frag$counts$Sample, c("S1", "S2")), "FragPipe sample order changed")
expect_true(identical(frag$counts$Identified_Protein_Count, c(1, 1)), "FragPipe identified counts are incorrect")
expect_true(identical(frag$counts$Quantified_Protein_Count, c(1, 1)), "FragPipe quantified counts are incorrect")
expect_true(is.na(frag$quantity["P001", "S1"]) && !is.na(frag$qualitative["P001", "S1"]), "FragPipe did not preserve identified-only case")
expect_true(!is.na(frag$quantity["P002", "S1"]) && is.na(frag$qualitative["P002", "S1"]), "FragPipe did not preserve quantified-only case")
frag_group <- make_group_info(frag$samples, c("G1", "G2"))
plot_sample_count_bar(
  frag$counts,
  frag_group,
  "Quantified_Protein_Count",
  file.path(out, "frag_quantified.pdf"),
  file.path(out, "frag_quantified.csv"),
  3,
  3,
  y_label = "Quantified protein count",
  plot_title = "Sample quantified protein count"
)
frag_export <- data.table::fread(file.path(out, "frag_quantified_sample_counts.csv"), data.table = FALSE)
expect_true(all(c("Identified_Protein_Count", "Quantified_Protein_Count") %in% colnames(frag_export)), "FragPipe exported sample counts did not retain both count columns")

peaks_path <- file.path(out, "peaks_test.csv")
write_utf8(peaks_path, c(
  "Protein Group,Top,Accession,Gene,Area Sample1,Area Sample2,Coverage(%) Sample1,Coverage(%) Sample2,Coverage(%)",
  "1,TRUE,P10001,GeneA,0,200,30,0,55",
  "1,TRUE,P10001_alt,GeneA_alt,500,500,80,80,90",
  "2,FALSE,P20002,GeneB,100,0,0,40,20"
))
peaks <- extract_peaks_protein_data(peaks_path, "accession")
expect_true(nrow(peaks$quantity) == 2, "PEAKS Protein Group reduction changed")
expect_true(identical(rownames(peaks$quantity), c("P10001", "P20002")), "PEAKS did not keep the first Top=TRUE entry per Protein Group")
expect_true(is.na(peaks$quantity["P10001", "Sample1"]) && !is.na(peaks$qualitative["P10001", "Sample1"]), "PEAKS identified-only case failed")
expect_true(!is.na(peaks$quantity["P20002", "Sample1"]) && is.na(peaks$qualitative["P20002", "Sample1"]), "PEAKS quantified-only case failed")
expect_true(identical(peaks$counts$Identified_Protein_Count, c(1, 1)), "PEAKS identified counts are incorrect")
expect_true(identical(peaks$counts$Quantified_Protein_Count, c(1, 1)), "PEAKS quantified counts are incorrect")

maxquant_path <- file.path(out, "maxquant_test.txt")
write_utf8(maxquant_path, c(
  paste(
    c("Protein IDs", "Fasta headers", "Gene names", "Protein names", "Reverse", "Potential contaminant", "Only identified by site", "LFQ intensity S1", "LFQ intensity S2"),
    collapse = "\t"
  ),
  paste(c("sp|P30001|P30001_HUMAN", "sp|P30001|P30001_HUMAN Protein 1", "GeneC", "Protein C", "", "", "", "100", "0"), collapse = "\t"),
  paste(c("sp|P30002|P30002_HUMAN", "sp|P30002|P30002_HUMAN Protein 2", "GeneD", "Protein D", "+", "", "", "999", "999"), collapse = "\t"),
  paste(c("sp|P30003|P30003_HUMAN", "sp|P30003|P30003_HUMAN Protein 3", "GeneE", "Protein E", "", "+", "", "999", "999"), collapse = "\t"),
  paste(c("sp|P30004|P30004_HUMAN", "sp|P30004|P30004_HUMAN Protein 4", "GeneF", "Protein F", "", "", "TRUE", "999", "999"), collapse = "\t")
))
mq <- extract_maxquant_protein_data(maxquant_path, "accession")
expect_true(nrow(mq$quantity) == 1, "MaxQuant filtering changed")
expect_true(identical(mq$counts$Identified_Protein_Count, mq$counts$Quantified_Protein_Count), "MaxQuant identified and quantified counts should match in V2.0.1")
expect_true(identical(mq$counts$Quantified_Protein_Count, c(1, 0)), "MaxQuant zero-to-NA handling changed")
expect_true(grepl("approximated by non-missing LFQ intensity", mq$count_approximation_note_en), "MaxQuant approximation note missing")

std_path <- file.path(out, "standard_test.csv")
write_utf8(std_path, c(
  "Feature,S1,S2",
  "Feat1,0,1",
  "Feat2,2,0"
))
std_a <- read_standard_matrix(std_path, "zero_is_value")
std_b <- read_standard_matrix(std_path, "zero_is_missing")
expect_true(identical(std_a$counts$Quantified_Protein_Count, c(2, 2)), "Standard matrix mode A quantified counts are incorrect")
expect_true(identical(std_b$counts$Quantified_Protein_Count, c(1, 1)), "Standard matrix mode B quantified counts are incorrect")
expect_true(all(is.na(std_a$counts$Identified_Protein_Count)) && all(is.na(std_b$counts$Identified_Protein_Count)), "Standard matrix should not fabricate identification counts")

diann <- extract_protein_data(resolve_external_test_file("PROTEOPOSTZ_DIANN_TEST_FILE"), "DIANN", "d", "protein_name")
expect_true("Quantified_Protein_Count" %in% colnames(diann$counts), "DIA-NN quantified count column missing")
expect_true(identical(diann$counts$Identified_Protein_Count, diann$counts$Quantified_Protein_Count), "DIA-NN count columns unexpectedly diverged")

spec <- extract_protein_data(resolve_external_test_file("PROTEOPOSTZ_SPECTRONAUT_TEST_FILE"), "Spectronaut", "d", "protein_name")
expect_true(all(spec$counts$Identified_Protein_Count == colSums(!is.na(spec$ibaq))), "Spectronaut identified counts no longer use PG.IBAQ")
expect_true(all(spec$counts$Quantified_Protein_Count == colSums(!is.na(spec$quantity))), "Spectronaut quantified counts no longer use PG.Quantity non-missing values")

shiny::testServer(server, {
  session$setInputs(input_family = "standard_matrix", file_path = std_path, outdir = out, standard_zero_mode = "zero_is_value")
  session$setInputs(load_data = 1)
  session$flushReact()
  expect_true(!is.null(rv$data), "Shiny failed to load standard matrix test data")
  expect_true(identical(count_metric_choices(rv$data), c("Quantified protein count" = "Quantified_Protein_Count")), "Standard matrix UI should only expose quantified count")
  expect_true(identical(active_count_column(), "Quantified_Protein_Count"), "Standard matrix default count selection is incorrect")

  session$setInputs(input_family = "dda", dda_format = "maxquant", file_path = maxquant_path)
  session$setInputs(load_data = 2)
  session$flushReact()
  expect_true(!is.null(rv$data) && identical(rv$data$input_source, "maxquant"), "Shiny failed to load MaxQuant test data")
  expect_true(identical(active_count_column(), "Identified_Protein_Count"), "Non-standard default count selection should remain identification count")
  expect_true(grepl("approximated by non-missing LFQ intensity", rv$data$count_approximation_note_en), "Shiny MaxQuant approximation note missing")
})

cat(
  "V201_COUNT_MODES_OK",
  "fragpipe=", paste(frag$counts$Identified_Protein_Count, frag$counts$Quantified_Protein_Count, collapse = "/"),
  "peaks_rows=", nrow(peaks$quantity),
  "maxquant_rows=", nrow(mq$quantity),
  "standard_mode_a=", paste(std_a$counts$Quantified_Protein_Count, collapse = ","),
  "standard_mode_b=", paste(std_b$counts$Quantified_Protein_Count, collapse = ","),
  "diann_samples=", ncol(diann$quantity),
  "spectronaut_samples=", ncol(spec$quantity),
  "\n",
  sep = ""
)
