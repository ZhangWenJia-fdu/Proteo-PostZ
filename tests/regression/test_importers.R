test_importers <- function(ctx) {
  f <- ctx$canonical
  writers <- list(
    importer_diann_basic = function(p) write_synthetic_diann(f, p),
    importer_spectronaut_basic = function(p) write_synthetic_spectronaut(f, p),
    importer_fragpipe_basic = function(p) write_synthetic_fragpipe(f, p),
    importer_peaks_db_basic = function(p) write_synthetic_peaks_db(f, p),
    importer_peaks_lfq_basic = function(p) write_synthetic_peaks_lfq(f, p),
    importer_maxquant_basic = function(p) write_synthetic_maxquant(f, p)
  )
  for (id in names(writers)) {
    record_test(id, {
      path <- file.path(ctx$fixtures, paste0(id, ".tsv"))
      writers[[id]](path)
      dat <- if (grepl("diann", id)) extract_protein_data(path, "DIANN", "d", "protein_name") else if (grepl("spectronaut", id)) extract_protein_data(path, "Spectronaut", "d", "protein_name") else extract_dda_protein_data(path, if (grepl("fragpipe", id)) "fragpipe" else if (grepl("peaks", id)) "peaks" else "maxquant", "protein_name")
      expect_equal(dim(dat$quantity), dim(f$quantity), paste(id, "dimensions differ"))
      expect_equal(dat$samples, f$samples, paste(id, "sample names differ"))
      expect_equal(length(rownames(dat$quantity)), nrow(f$quantity), paste(id, "feature identifiers missing"))
    })
  }
  for (sep in c(",", "\t")) {
    record_test(paste0("importer_standard_", if (sep == ",") "csv" else "tsv"), {
      path <- file.path(ctx$fixtures, paste0("standard.", if (sep == ",") "csv" else "tsv"))
      write_synthetic_standard(f, path, sep)
      dat <- read_standard_matrix(path, "zero_is_value")
      expect_equal(dim(dat$quantity), dim(f$quantity), "Standard dimensions differ")
      expect_equal(dat$samples, f$samples, "Standard sample names differ")
    })
  }
}
