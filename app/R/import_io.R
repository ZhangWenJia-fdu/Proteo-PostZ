# ProteoPostZ input readers, format detection, and importers

read_result_file <- function(file) {
  first_line <- readLines(file, n = 1, warn = FALSE, encoding = "UTF-8")
  delim <- if (grepl("\t", first_line)) "\t" else if (grepl(";", first_line)) ";" else ","
  dat <- data.table::fread(
    file = file, sep = delim, header = TRUE,
    na.strings = c("", "NA", "NaN", "null", "NULL"), quote = "\"",
    data.table = FALSE, check.names = FALSE, encoding = "UTF-8"
  )
  colnames(dat) <- trimws(colnames(dat))
  colnames(dat)[1] <- sub("^\ufeff", "", colnames(dat)[1])
  dat
}
has_any_col <- function(dat, candidates) {
  any(candidates %in% colnames(dat))
}
grep_cols <- function(dat, pattern, ignore.case = TRUE) {
  grep(pattern, colnames(dat), value = TRUE, ignore.case = ignore.case)
}

first_present_col <- function(dat, candidates) {
  hit <- candidates[candidates %in% colnames(dat)]
  if (length(hit) > 0) hit[[1]] else NA_character_
}

diann_annotation_columns <- function() {
  c(
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description",
    "N.Sequences", "N.Proteotypic.Sequences", "Genes.Quantity", "Genes.Normalised",
    "Protein.Q.Value", "PG.Q.Value", "Global.Q.Value", "Q.Value"
  )
}

diann_quantity_columns <- function(dat) {
  cols <- colnames(dat)
  direct <- grep("\\.(d|raw|wiff)$", cols, value = TRUE, ignore.case = TRUE)
  if (length(direct) > 0) return(direct)

  annotation_idx <- which(cols %in% diann_annotation_columns())
  if (length(annotation_idx) > 0 && max(annotation_idx) < length(cols)) {
    candidates <- cols[seq.int(max(annotation_idx) + 1, length(cols))]
  } else {
    candidates <- setdiff(cols, diann_annotation_columns())
  }
  numeric_fraction <- vapply(candidates, function(col) {
    vals <- suppressWarnings(as.numeric(take_first(dat[[col]])))
    mean(!is.na(vals))
  }, numeric(1))
  candidates[numeric_fraction > 0]
}

detect_diann_signature <- function(dat) {
  has_ids <- has_any_col(dat, c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes"))
  sample_cols <- diann_quantity_columns(dat)
  if (has_ids && length(sample_cols) > 0) {
    return(paste0("DIA-NN signature: protein annotation columns plus ", length(sample_cols), " sample quantity columns."))
  }
  NULL
}

detect_spectronaut_signature <- function(dat) {
  qcols <- grep_cols(dat, "PG\\.Quantity$")
  icols <- grep_cols(dat, "PG\\.IBAQ$")
  if (length(qcols) > 0 && length(icols) > 0) {
    return(paste0("Spectronaut signature: ", length(qcols), " PG.Quantity columns and ", length(icols), " PG.IBAQ columns."))
  }
  NULL
}

detect_maxquant_signature <- function(dat) {
  qcols <- grep_cols(dat, "^LFQ intensity ")
  if (length(qcols) > 0 && has_any_col(dat, c("Gene names", "Protein names", "Majority protein IDs", "Protein IDs"))) {
    return(paste0("MaxQuant signature: ", length(qcols), " LFQ intensity columns."))
  }
  NULL
}

detect_fragpipe_signature <- function(dat) {
  qcols <- grep_cols(dat, " MaxLFQ Intensity$")
  if (length(qcols) > 0 && has_any_col(dat, c("Protein", "Protein ID", "Gene", "Description"))) {
    return(paste0("FragPipe/MSFragger signature: ", length(qcols), " sample MaxLFQ Intensity columns in combined_protein.tsv."))
  }
  NULL
}

detect_proteome_discoverer_signature <- function(dat) {
  qcols <- grep_cols(dat, "^Abundance")
  if (length(qcols) > 0 && has_any_col(dat, c("Accession", "Gene Symbol", "Description"))) {
    return(paste0("Proteome Discoverer signature: ", length(qcols), " Abundance columns."))
  }
  NULL
}

detect_peaks_signature <- function(dat) {
  schema <- peaks_schema(dat)
  if (!is.null(schema)) return(schema$evidence)
  NULL
}

input_format_registry <- function(input_family = c("dia", "dda")) {
  input_family <- match.arg(input_family)
  if (input_family == "dia") {
    return(list(
      list(id = "diann", label = "DIA-NN", detector = detect_diann_signature),
      list(id = "spectronaut", label = "Spectronaut", detector = detect_spectronaut_signature)
    ))
  }
  list(
    list(id = "fragpipe", label = "FragPipe / MSFragger", detector = detect_fragpipe_signature),
    list(id = "peaks", label = "PEAKS", detector = detect_peaks_signature),
    list(id = "maxquant", label = "MaxQuant", detector = detect_maxquant_signature)
  )
}

detect_input_format <- function(file, input_family = c("dia", "dda")) {
  input_family <- match.arg(input_family)
  raw <- read_result_file(file)
  hits <- list()
  for (entry in input_format_registry(input_family)) {
    evidence <- entry$detector(raw)
    if (!is.null(evidence)) {
      if (identical(entry$id, "peaks")) {
        schema <- peaks_schema(raw)
        if (is.null(schema) || identical(schema$id, "ambiguous")) {
          hits[[length(hits) + 1]] <- list(id = NA_character_, label = "PEAKS ambiguous", evidence = schema$evidence %||% "Unsupported or ambiguous PEAKS protein result schema.")
        } else {
          hits[[length(hits) + 1]] <- list(id = entry$id, label = schema$label, evidence = evidence)
        }
      } else {
        hits[[length(hits) + 1]] <- list(id = entry$id, label = entry$label, evidence = evidence)
      }
    }
  }
  if (length(hits) == 0) {
    return(list(id = NA_character_, label = "No supported format detected", evidence = "No registered parser signature matched the file header.", ambiguous = FALSE))
  }
  first <- hits[[1]]
  first$ambiguous <- length(hits) > 1
  first$all_hits <- vapply(hits, function(x) x$label, character(1))
  first
}

standard_matrix_zero_modes <- c("zero_is_value", "zero_is_missing")

standard_matrix_zero_mode_label <- function(missingness_mode) {
  missingness_mode <- match.arg(missingness_mode, standard_matrix_zero_modes)
  if (missingness_mode == "zero_is_value") "0 is a valid quantitative value" else "0 represents missing or unquantified"
}

read_standard_matrix_file <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (!ext %in% c("csv", "tsv")) stop("Standard matrix input must be a .csv or .tsv file.")
  first_line <- readLines(file, n = 1, warn = FALSE, encoding = "UTF-8")
  first_line <- sub("^\ufeff", "", first_line)
  expected_sep <- if (ext == "tsv") "\t" else ","
  wrong_sep <- if (ext == "tsv") "," else "\t"
  if (!grepl(expected_sep, first_line, fixed = TRUE) && grepl(wrong_sep, first_line, fixed = TRUE)) {
    stop("File extension and delimiter do not match: .", ext, " input must use ", if (ext == "tsv") "tab" else "comma", " delimiters.")
  }
  dat <- data.table::fread(
    file = file, sep = expected_sep, header = TRUE, na.strings = character(),
    quote = "\"", data.table = FALSE, check.names = FALSE, encoding = "UTF-8",
    colClasses = list(character = 1)
  )
  colnames(dat)[1] <- sub("^\ufeff", "", colnames(dat)[1])
  if (ncol(dat) == 1) {
    stop("Standard matrix input was read as one column. Check that the file delimiter matches the .", ext, " extension.")
  }
  dat
}

parse_standard_quant_value <- function(x) {
  raw <- as.character(x)
  trimmed <- trimws(raw)
  explicit_missing <- is.na(raw) | trimmed == "" | tolower(trimmed) %in% c("na", "nan")
  numeric <- rep(NA_real_, length(trimmed))
  candidates <- !explicit_missing
  suppressWarnings(numeric[candidates] <- as.numeric(trimmed[candidates]))
  invalid_text <- candidates & is.na(numeric)
  infinite <- candidates & !is.na(numeric) & !is.finite(numeric)
  list(value = numeric, explicit_missing = explicit_missing, invalid_text = invalid_text, infinite = infinite, raw = raw)
}

summarize_bad_values <- function(values, max_values = 5) {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]
  paste(head(values, max_values), collapse = ", ")
}

read_standard_matrix <- function(file, missingness_mode = NULL) {
  load_required_packages()
  if (is.null(missingness_mode) || !nzchar(as.character(missingness_mode))) {
    stop("Standard matrix zero handling mode must be explicitly selected: zero_is_value or zero_is_missing.")
  }
  missingness_mode <- match.arg(missingness_mode, standard_matrix_zero_modes)
  raw <- read_standard_matrix_file(file)
  if (ncol(raw) < 2) stop("Standard matrix input must contain at least two columns: one feature identifier column and at least one sample column.")

  feature_colname <- colnames(raw)[1]
  sample_names <- colnames(raw)[-1]
  if (any(is.na(sample_names) | trimws(sample_names) == "")) stop("Standard matrix sample column names cannot be empty.")
  dup_samples <- unique(sample_names[duplicated(sample_names)])
  if (length(dup_samples) > 0) stop("Standard matrix sample column names must be unique. Duplicates: ", paste(dup_samples, collapse = ", "))

  features <- as.character(raw[[1]])
  feature_blank <- is.na(features) | trimws(features) == ""
  if (all(feature_blank)) stop("The first feature identifier column cannot be entirely empty.")
  if (any(feature_blank)) stop("Standard matrix feature identifiers cannot be blank. Blank rows: ", paste(which(feature_blank), collapse = ", "))
  dup_features <- unique(features[duplicated(features)])
  if (length(dup_features) > 0) stop("Standard matrix feature identifiers must be unique. Duplicates: ", paste(head(dup_features, 5), collapse = ", "), ". Please resolve duplicate features before analysis.")

  qdat <- raw[, -1, drop = FALSE]
  parsed <- lapply(qdat, parse_standard_quant_value)
  invalid_cols <- names(parsed)[vapply(parsed, function(z) any(z$invalid_text), logical(1))]
  if (length(invalid_cols) > 0) {
    details <- vapply(invalid_cols, function(nm) summarize_bad_values(parsed[[nm]]$raw[parsed[[nm]]$invalid_text]), character(1))
    stop("Standard matrix quantitative columns contain non-numeric text. ", paste(paste0(invalid_cols, ": ", details), collapse = "; "))
  }
  inf_cols <- names(parsed)[vapply(parsed, function(z) any(z$infinite), logical(1))]
  if (length(inf_cols) > 0) {
    details <- vapply(inf_cols, function(nm) summarize_bad_values(parsed[[nm]]$raw[parsed[[nm]]$infinite]), character(1))
    stop("Standard matrix quantitative columns contain Inf or -Inf, which are invalid values. ", paste(paste0(inf_cols, ": ", details), collapse = "; "))
  }

  raw_quant_matrix <- do.call(cbind, lapply(parsed, `[[`, "value"))
  raw_quant_matrix <- as.matrix(raw_quant_matrix)
  storage.mode(raw_quant_matrix) <- "numeric"
  rownames(raw_quant_matrix) <- features
  colnames(raw_quant_matrix) <- sample_names

  explicit_missing_matrix <- do.call(cbind, lapply(parsed, `[[`, "explicit_missing"))
  explicit_missing_matrix <- as.matrix(explicit_missing_matrix)
  rownames(explicit_missing_matrix) <- features
  colnames(explicit_missing_matrix) <- sample_names

  zero_matrix <- !is.na(raw_quant_matrix) & raw_quant_matrix == 0
  analysis_quant_matrix <- raw_quant_matrix
  zero_to_missing_count <- 0L
  if (missingness_mode == "zero_is_missing") {
    zero_to_missing_count <- sum(zero_matrix)
    analysis_quant_matrix[zero_matrix] <- NA_real_
  }
  missingness_matrix <- is.na(analysis_quant_matrix)
  all_missing_samples <- colnames(analysis_quant_matrix)[colSums(!is.na(analysis_quant_matrix)) == 0]
  if (length(all_missing_samples) > 0) {
    stop("Standard matrix sample columns cannot be entirely missing under the selected zero handling mode. Columns: ", paste(all_missing_samples, collapse = ", "))
  }

  counts <- make_sample_count_table(
    analysis_quant_matrix,
    analysis_quant_matrix,
    identified_available = FALSE
  )
  counts$Available_Quantitative_Value_Count <- counts$Quantified_Protein_Count
  feature_counts <- data.frame(
    Metric = "number of features with at least one available quantitative value",
    Value = sum(rowSums(!is.na(analysis_quant_matrix)) > 0),
    stringsAsFactors = FALSE
  )
  meta <- data.frame(RowID = features, AnalysisID = features, FeatureID = features, stringsAsFactors = FALSE)

  list(
    input_source = "standard_matrix",
    software = "standard_matrix",
    feature_column_name = feature_colname,
    features = features,
    samples = sample_names,
    raw = raw,
    meta = meta,
    raw_quant_matrix = raw_quant_matrix,
    analysis_quant_matrix = analysis_quant_matrix,
    quantity = analysis_quant_matrix,
    qualitative = analysis_quant_matrix,
    ibaq = NULL,
    missingness_matrix = missingness_matrix,
    explicit_missing_matrix = explicit_missing_matrix,
    missingness_mode = missingness_mode,
    missingness_mode_label = standard_matrix_zero_mode_label(missingness_mode),
    raw_zero_cell_count = sum(zero_matrix),
    zero_to_missing_count = zero_to_missing_count,
    explicit_missing_value_count = sum(explicit_missing_matrix),
    counts = counts,
    feature_counts = feature_counts,
    available_quantitative_value_count = sum(!is.na(analysis_quant_matrix))
  )
}

make_sample_count_table <- function(quantity, qualitative = NULL, identified_available = TRUE) {
  quantity <- as.matrix(quantity)
  if (is.null(qualitative)) qualitative <- quantity
  qualitative <- as.matrix(qualitative)
  if (!identical(colnames(quantity), colnames(qualitative))) {
    stop("Quantity and qualitative matrices must share the same sample columns for count summaries.")
  }
  data.frame(
    Sample = colnames(quantity),
    Identified_Protein_Count = if (isTRUE(identified_available)) colSums(!is.na(qualitative)) else rep(NA_integer_, ncol(quantity)),
    Quantified_Protein_Count = colSums(!is.na(quantity)),
    stringsAsFactors = FALSE
  )
}

take_first <- function(x) {
  x <- trimws(as.character(x))
  x <- sub(";.*$", "", x)
  x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
  x
}

clean_accession_value <- function(x) {
  x <- take_first(x)
  x <- ifelse(grepl("^(sp|tr)\\|[^|]+\\|", x), sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", x), sub("\\|.*$", "", x))
  x
}

entry_name_from_accession_field <- function(x) {
  x <- take_first(x)
  out <- ifelse(grepl("^[^|]+\\|[^|]+$", x), sub("^[^|]+\\|([^|]+)$", "\\1", x), NA_character_)
  out
}

entry_name_from_fasta_header <- function(x) {
  x <- take_first(x)
  out <- ifelse(grepl("^(sp|tr)\\|[^|]+\\|[^ ]+", x), sub("^(sp|tr)\\|[^|]+\\|([^ ]+).*$", "\\2", x), NA_character_)
  out
}

take_numeric <- function(x) {
  x <- take_first(x)
  suppressWarnings(as.numeric(x))
}

clean_sample_name <- function(x, suffix = c("diann", "d", "raw", "spectronaut")) {
  suffix <- match.arg(suffix)
  x <- gsub("^\\[\\d+\\]\\s*", "", x)
  x <- gsub("\\\\", "/", x)
  x <- basename(x)
  if (suffix == "spectronaut") {
    x <- gsub("\\.d\\.PG\\.(Quantity|IBAQ)$", "", x, ignore.case = TRUE)
    x <- gsub("\\.raw\\.PG\\.(Quantity|IBAQ)$", "", x, ignore.case = TRUE)
    x <- gsub("\\.PG\\.(Quantity|IBAQ)$", "", x, ignore.case = TRUE)
  }
  x <- gsub("\\.(d|raw|wiff)$", "", x, ignore.case = TRUE)
  x
}

id_columns <- function(software) {
  if (software == "DIANN") {
    c(protein_name = "Protein.Names", gene_name = "Genes", accession = "Protein.Group")
  } else {
    c(protein_name = "PG.ProteinNames", gene_name = "PG.Genes", accession = "PG.ProteinAccessions")
  }
}

extract_protein_data <- function(file, software = c("DIANN", "Spectronaut"), diann_type = c("d", "raw"), row_id = c("protein_name", "gene_name", "accession")) {
  load_required_packages()
  software <- match.arg(software)
  diann_type <- match.arg(diann_type)
  row_id <- match.arg(row_id)
  raw <- read_result_file(file)
  ids <- id_columns(software)
  missing_ids <- setdiff(unname(ids), colnames(raw))
  if (length(missing_ids) > 0) stop("Missing required columns: ", paste(missing_ids, collapse = ", "))

  row_values <- take_first(raw[[ids[[row_id]]]])
  accession_values <- take_first(raw[[ids[["accession"]]]])
  protein_values <- take_first(raw[[ids[["protein_name"]]]])
  gene_values <- take_first(raw[[ids[["gene_name"]]]])

  if (software == "DIANN") {
    qcols <- diann_quantity_columns(raw)
    if (length(qcols) == 0) stop("No DIA-NN sample quantity columns were found.")
    qmat <- as.data.frame(lapply(raw[, qcols, drop = FALSE], take_numeric), check.names = FALSE)
    if (all(grepl("\\.(d|raw|wiff)$", qcols, ignore.case = TRUE))) {
      colnames(qmat) <- clean_sample_name(colnames(qmat), suffix = "diann")
    } else {
      colnames(qmat) <- make.unique(qcols)
    }
    ident_mat <- qmat
    ibaq_mat <- NULL
  } else {
    qcols <- grep("PG\\.Quantity$", colnames(raw), value = TRUE)
    icols <- grep("PG\\.IBAQ$", colnames(raw), value = TRUE)
    if (length(qcols) == 0) stop("No Spectronaut PG.Quantity columns were found.")
    if (length(icols) == 0) stop("No Spectronaut PG.IBAQ columns were found.")
    qmat <- as.data.frame(lapply(raw[, qcols, drop = FALSE], take_numeric), check.names = FALSE)
    ibaq_mat <- as.data.frame(lapply(raw[, icols, drop = FALSE], take_numeric), check.names = FALSE)
    colnames(qmat) <- clean_sample_name(colnames(qmat), suffix = "spectronaut")
    colnames(ibaq_mat) <- clean_sample_name(colnames(ibaq_mat), suffix = "spectronaut")
    if (!setequal(colnames(qmat), colnames(ibaq_mat))) stop("PG.Quantity and PG.IBAQ sample names do not match.")
    if (!identical(colnames(qmat), colnames(ibaq_mat))) ibaq_mat <- ibaq_mat[, match(colnames(qmat), colnames(ibaq_mat)), drop = FALSE]
    ident_mat <- ibaq_mat
  }

  keep <- !is.na(row_values) & row_values != ""
  meta <- data.frame(RowID = row_values[keep], ProteinName = protein_values[keep], GeneName = gene_values[keep], Accession = accession_values[keep], stringsAsFactors = FALSE)
  quantity <- qmat[keep, , drop = FALSE]
  ident <- ident_mat[keep, , drop = FALSE]
  if (!is.null(ibaq_mat)) ibaq_mat <- ibaq_mat[keep, , drop = FALSE]
  analysis_ids <- make.unique(meta$RowID)
  meta$AnalysisID <- analysis_ids
  rownames(quantity) <- analysis_ids
  rownames(ident) <- analysis_ids
  if (!is.null(ibaq_mat)) rownames(ibaq_mat) <- analysis_ids

  counts <- make_sample_count_table(quantity, ident)
  list(raw = raw, meta = meta, quantity = as.matrix(quantity), qualitative = as.matrix(ident), ibaq = if (is.null(ibaq_mat)) NULL else as.matrix(ibaq_mat), counts = counts, software = software, samples = colnames(quantity))
}

clean_dda_sample_name <- function(x, software = c("fragpipe", "proteome_discoverer", "peaks", "maxquant")) {
  software <- match.arg(software)
  x <- trimws(as.character(x))
  if (software == "maxquant") {
    x <- sub("^LFQ intensity\\s+", "", x, ignore.case = TRUE)
    x <- sub("^Intensity\\s+", "", x, ignore.case = TRUE)
  } else if (software == "fragpipe") {
    x <- sub("\\s+MaxLFQ Unique Intensity$", "", x, ignore.case = TRUE)
    x <- sub("\\s+MaxLFQ Intensity$", "", x, ignore.case = TRUE)
  } else if (software == "proteome_discoverer") {
    x <- sub("^Abundance\\s*[:]?\\s*", "", x, ignore.case = TRUE)
    x <- sub("^Ratio\\s*[:]?\\s*", "", x, ignore.case = TRUE)
  } else if (software == "peaks") {
    x <- sub("^Area\\s*[:]?\\s*", "", x, ignore.case = TRUE)
    x <- sub("\\s+Area$", "", x, ignore.case = TRUE)
    x[x == "Area"] <- paste0("Sample", seq_len(sum(x == "Area")))
  }
  x <- gsub("^F[0-9]+\\s*[:]\\s*", "", x)
  x <- trimws(x)
  x[x == ""] <- paste0("Sample", seq_len(sum(x == "")))
  make.unique(x)
}

dda_id_candidates <- function(software = c("fragpipe", "proteome_discoverer", "peaks", "maxquant")) {
  software <- match.arg(software)
  if (software == "maxquant") {
    return(list(
      protein_name = c("Protein names", "Fasta headers", "Majority protein IDs", "Protein IDs"),
      gene_name = c("Gene names", "Gene Names"),
      accession = c("Majority protein IDs", "Protein IDs", "Accession")
    ))
  }
  if (software == "fragpipe") {
    return(list(
      protein_name = c("Entry Name", "Protein", "Protein ID", "Description"),
      gene_name = c("Gene", "Gene Name", "Mapped Genes"),
      accession = c("Protein ID", "Protein", "Entry Name")
    ))
  }
  if (software == "proteome_discoverer") {
    return(list(
      protein_name = c("Description", "Protein Description", "Accession"),
      gene_name = c("Gene Symbol", "Gene", "Gene ID"),
      accession = c("Accession", "Master Protein Accessions")
    ))
  }
  list(
    protein_name = c("Accession", "Description", "Protein ID"),
    gene_name = c("Gene", "Gene Name"),
    accession = c("Accession", "Protein ID", "Protein Group")
  )
}

dda_quantity_columns <- function(raw, software = c("fragpipe", "proteome_discoverer", "peaks", "maxquant")) {
  software <- match.arg(software)
  if (software == "maxquant") return(grep_cols(raw, "^LFQ intensity "))
  if (software == "fragpipe") return(grep_cols(raw, " MaxLFQ Intensity$"))
  if (software == "proteome_discoverer") return(grep_cols(raw, "^Abundance"))
  if (software == "peaks") {
    schema <- peaks_schema(raw)
    if (is.null(schema) || identical(schema$id, "ambiguous")) return(character())
    return(if (identical(schema$id, "db")) peaks_db_quantity_columns(raw) else peaks_lfq_quantity_columns(raw))
  }
  grep_cols(raw, "(^Area[: ]| Area$|^Area$)")
}

maxquant_filter_keep <- function(raw) {
  keep <- rep(TRUE, nrow(raw))
  for (col in c("Reverse", "Potential contaminant", "Only identified by site")) {
    if (col %in% colnames(raw)) {
      val <- trimws(as.character(raw[[col]]))
      keep <- keep & !(val %in% c("+", "TRUE", "True", "true", "1"))
    }
  }
  keep
}

fragpipe_spectral_count_columns <- function(raw) {
  cols <- grep_cols(raw, " Spectral Count$")
  cols[!grepl("(^Combined | Unique Spectral Count$| Total Spectral Count$)", cols, ignore.case = TRUE)]
}

extract_maxquant_protein_data <- function(file, row_id = c("protein_name", "gene_name", "accession")) {
  row_id <- match.arg(row_id)
  raw_all <- read_result_file(file)
  qcols <- dda_quantity_columns(raw_all, "maxquant")
  if (length(qcols) == 0) stop("No MaxQuant LFQ intensity sample columns were found.")

  keep_filter <- maxquant_filter_keep(raw_all)
  raw <- raw_all[keep_filter, , drop = FALSE]
  qcols <- dda_quantity_columns(raw, "maxquant")

  accession_col <- first_present_col(raw, c("Protein IDs", "Majority protein IDs"))
  gene_col <- first_present_col(raw, c("Gene names", "Gene Names"))
  fasta_col <- first_present_col(raw, c("Fasta headers", "Fasta Headers"))
  description_col <- first_present_col(raw, c("Protein names", "Protein Names"))
  if (is.na(accession_col)) stop("No MaxQuant Protein IDs or Majority protein IDs column was found.")
  if (is.na(fasta_col)) stop("No MaxQuant Fasta headers column was found for UniProt entry-name extraction.")

  accession_values <- clean_accession_value(raw[[accession_col]])
  protein_values <- entry_name_from_fasta_header(raw[[fasta_col]])
  fallback_protein <- if (!is.na(description_col)) take_first(raw[[description_col]]) else accession_values
  protein_values[is.na(protein_values) | protein_values == ""] <- fallback_protein[is.na(protein_values) | protein_values == ""]
  gene_values <- if (!is.na(gene_col)) take_first(raw[[gene_col]]) else rep(NA_character_, length(accession_values))
  row_values <- switch(row_id, protein_name = protein_values, gene_name = gene_values, accession = accession_values)

  qmat <- as.data.frame(lapply(raw[, qcols, drop = FALSE], take_numeric), check.names = FALSE)
  qmat[qmat == 0] <- NA_real_
  colnames(qmat) <- clean_dda_sample_name(qcols, "maxquant")

  keep <- !is.na(row_values) & row_values != ""
  if (!any(keep)) stop("No non-empty feature identifiers were found in the selected MaxQuant identifier field.")
  meta <- data.frame(
    RowID = row_values[keep],
    ProteinName = protein_values[keep],
    GeneName = gene_values[keep],
    Accession = accession_values[keep],
    SourceFormat = "maxquant",
    IdentifierColumn = if (row_id == "accession") accession_col else if (row_id == "gene_name") gene_col else fasta_col,
    ProteinIDs = if ("Protein IDs" %in% colnames(raw)) raw[["Protein IDs"]][keep] else NA,
    MajorityProteinIDs = if ("Majority protein IDs" %in% colnames(raw)) raw[["Majority protein IDs"]][keep] else NA,
    ProteinDescriptions = if (!is.na(description_col)) raw[[description_col]][keep] else NA,
    FastaHeaders = raw[[fasta_col]][keep],
    stringsAsFactors = FALSE
  )
  quantity <- qmat[keep, , drop = FALSE]
  analysis_ids <- make.unique(meta$RowID)
  meta$AnalysisID <- analysis_ids
  rownames(quantity) <- analysis_ids
  quantity <- as.matrix(quantity)
  storage.mode(quantity) <- "numeric"
  counts <- make_sample_count_table(quantity, quantity)
  list(
    raw = raw,
    raw_all = raw_all,
    meta = meta,
    quantity = quantity,
    qualitative = quantity,
    ibaq = NULL,
    counts = counts,
    software = "maxquant",
    input_family = "dda",
    input_source = "maxquant",
    data_level = "protein",
    format_evidence = "MaxQuant proteinGroups.txt: filtered reverse, contaminant, and only-identified-by-site rows; sample LFQ intensity columns are used for label-free quantification; zeros are treated as missing.",
    count_approximation_note_en = "MaxQuant identification count is approximated by non-missing LFQ intensity in the current proteinGroups.txt importer.",
    features = rownames(quantity),
    samples = colnames(quantity)
  )
}

extract_fragpipe_protein_data <- function(file, row_id = c("protein_name", "gene_name", "accession")) {
  row_id <- match.arg(row_id)
  raw <- read_result_file(file)
  qcols <- dda_quantity_columns(raw, "fragpipe")
  if (length(qcols) == 0) stop("No FragPipe sample MaxLFQ Intensity columns were found.")
  scols <- fragpipe_spectral_count_columns(raw)
  if (length(scols) == 0) stop("No FragPipe sample Spectral Count columns were found for identification counts.")

  ids <- dda_id_candidates("fragpipe")
  selected_id_col <- first_present_col(raw, ids[[row_id]])
  accession_col <- first_present_col(raw, ids[["accession"]])
  protein_col <- first_present_col(raw, ids[["protein_name"]])
  gene_col <- first_present_col(raw, ids[["gene_name"]])
  if (is.na(selected_id_col)) stop("No suitable ", row_id, " identifier column was found for FragPipe/MSFragger.")

  row_values <- if (row_id == "accession") clean_accession_value(raw[[selected_id_col]]) else take_first(raw[[selected_id_col]])
  accession_values <- if (!is.na(accession_col)) clean_accession_value(raw[[accession_col]]) else row_values
  protein_values <- if (!is.na(protein_col)) take_first(raw[[protein_col]]) else row_values
  gene_values <- if (!is.na(gene_col)) take_first(raw[[gene_col]]) else rep(NA_character_, length(row_values))
  if (!is.na(protein_col) && protein_col == "Entry Name") protein_values <- take_first(raw[[protein_col]])

  qmat <- as.data.frame(lapply(raw[, qcols, drop = FALSE], take_numeric), check.names = FALSE)
  qmat[qmat == 0] <- NA_real_
  colnames(qmat) <- clean_dda_sample_name(qcols, "fragpipe")

  ident_mat <- as.data.frame(lapply(raw[, scols, drop = FALSE], take_numeric), check.names = FALSE)
  ident_mat[ident_mat == 0] <- NA_real_
  colnames(ident_mat) <- clean_dda_sample_name(sub(" Spectral Count$", " MaxLFQ Intensity", scols), "fragpipe")
  common_samples <- intersect(colnames(qmat), colnames(ident_mat))
  if (length(common_samples) == 0) stop("FragPipe MaxLFQ Intensity columns and Spectral Count columns could not be matched by sample name.")
  qmat <- qmat[, common_samples, drop = FALSE]
  ident_mat <- ident_mat[, common_samples, drop = FALSE]

  keep <- !is.na(row_values) & row_values != ""
  if (!any(keep)) stop("No non-empty feature identifiers were found in the selected FragPipe identifier column: ", selected_id_col)
  meta <- data.frame(
    RowID = row_values[keep],
    ProteinName = protein_values[keep],
    GeneName = gene_values[keep],
    Accession = accession_values[keep],
    SourceFormat = "fragpipe",
    IdentifierColumn = selected_id_col,
    Protein = if ("Protein" %in% colnames(raw)) raw[["Protein"]][keep] else NA,
    IndistinguishableProteins = if ("Indistinguishable Proteins" %in% colnames(raw)) raw[["Indistinguishable Proteins"]][keep] else NA,
    stringsAsFactors = FALSE
  )
  quantity <- qmat[keep, , drop = FALSE]
  ident <- ident_mat[keep, , drop = FALSE]
  analysis_ids <- make.unique(meta$RowID)
  meta$AnalysisID <- analysis_ids
  rownames(quantity) <- analysis_ids
  rownames(ident) <- analysis_ids
  quantity <- as.matrix(quantity)
  ident <- as.matrix(ident)
  storage.mode(quantity) <- "numeric"
  storage.mode(ident) <- "numeric"
  counts <- make_sample_count_table(quantity, ident)
  list(
    raw = raw,
    meta = meta,
    quantity = quantity,
    qualitative = ident,
    ibaq = NULL,
    counts = counts,
    software = "fragpipe",
    input_family = "dda",
    input_source = "fragpipe",
    data_level = "protein",
    format_evidence = "FragPipe combined_protein.tsv: rows are protein groups; sample MaxLFQ Intensity columns are used for label-free quantification; zeros are treated as missing; sample Spectral Count > 0 defines identification.",
    features = rownames(quantity),
    samples = colnames(quantity)
  )
}

peaks_coverage_columns <- function(raw) {
  cols <- grep_cols(raw, "^Coverage\\(%\\)\\s+")
  cols[nzchar(trimws(sub("^Coverage\\(%\\)\\s+", "", cols)))]
}

peaks_db_quantity_columns <- function(raw) {
  grep_cols(raw, "^Area\\s+.+$", ignore.case = FALSE)
}

peaks_lfq_quantity_columns <- function(raw) {
  cols <- grep_cols(raw, "^.+\\s+Area$", ignore.case = FALSE)
  cols[!grepl("^Group\\s+[0-9]+\\s+Area$", cols, ignore.case = TRUE)]
}

peaks_schema <- function(raw) {
  cols <- colnames(raw)
  if (!all(c("Protein Group", "Top", "Accession") %in% cols)) return(NULL)
  db_q <- peaks_db_quantity_columns(raw)
  coverage <- peaks_coverage_columns(raw)
  lfq_q <- peaks_lfq_quantity_columns(raw)
  lfq_markers <- c("Sample Profile (Ratio)", "Group Profile (Ratio)", "Significance")
  has_lfq_marker <- any(lfq_markers %in% cols) || any(grepl("^Group\\s+[0-9]+\\s+Area$", cols, ignore.case = TRUE))
  is_db <- length(db_q) > 0 && length(coverage) > 0
  is_lfq <- length(lfq_q) > 0 && has_lfq_marker
  if (is_db) {
    mixed_note <- if (is_lfq) " LFQ-style fields are also present; DB schema takes precedence." else ""
    return(list(id = "db", label = "PEAKS DB protein result", mixed_db_lfq = is_lfq, evidence = paste0("PEAKS DB protein result: ", length(db_q), " Area <sample> columns and ", length(coverage), " sample-specific Coverage(%) columns.", mixed_note)))
  }
  if (is_lfq) return(list(id = "lfq", label = "PEAKS LFQ protein result", evidence = paste0("PEAKS LFQ protein result: ", length(lfq_q), " <sample> Area columns; LFQ profile/group fields detected.")))
  NULL
}

peaks_top_keep_index <- function(raw) {
  if (!"Protein Group" %in% colnames(raw)) return(seq_len(nrow(raw)))
  group <- as.character(raw[["Protein Group"]])
  missing_group <- is.na(group) | trimws(group) == ""
  group[missing_group] <- paste0("__missing_group_", which(missing_group))
  top <- if ("Top" %in% colnames(raw)) {
    toupper(trimws(as.character(raw[["Top"]]))) %in% c("TRUE", "T", "YES", "Y", "1")
  } else {
    rep(FALSE, nrow(raw))
  }
  idx <- seq_len(nrow(raw))
  groups_in_order <- unique(group)
  unlist(lapply(groups_in_order, function(g) {
    ii <- idx[group == g]
    top_ii <- ii[top[ii]]
    if (length(top_ii) > 0) top_ii[[1]] else ii[[1]]
  }), use.names = FALSE)
}

extract_peaks_protein_data <- function(file, row_id = c("protein_name", "gene_name", "accession")) {
  row_id <- match.arg(row_id)
  raw_all <- read_result_file(file)
  schema <- peaks_schema(raw_all)
  if (is.null(schema)) stop("Unsupported or ambiguous PEAKS protein result schema. Expected DB-style Area/Coverage columns or LFQ-style sample Area/profile fields.")
  if (identical(schema$id, "ambiguous")) stop(schema$evidence)
  qcols <- if (identical(schema$id, "db")) peaks_db_quantity_columns(raw_all) else peaks_lfq_quantity_columns(raw_all)
  coverage_cols <- if (identical(schema$id, "db")) peaks_coverage_columns(raw_all) else character()
  keep_idx <- peaks_top_keep_index(raw_all)
  raw <- raw_all[keep_idx, , drop = FALSE]

  ids <- dda_id_candidates("peaks")
  selected_id_col <- first_present_col(raw, ids[[row_id]])
  accession_col <- first_present_col(raw, ids[["accession"]])
  protein_col <- first_present_col(raw, ids[["protein_name"]])
  gene_col <- first_present_col(raw, ids[["gene_name"]])
  if (is.na(selected_id_col)) stop("No suitable ", row_id, " identifier column was found for PEAKS.")

  row_values <- if (row_id == "accession") clean_accession_value(raw[[selected_id_col]]) else take_first(raw[[selected_id_col]])
  accession_values <- if (!is.na(accession_col)) clean_accession_value(raw[[accession_col]]) else row_values
  protein_values <- if (!is.na(protein_col) && protein_col == "Accession") entry_name_from_accession_field(raw[[protein_col]]) else if (!is.na(protein_col)) take_first(raw[[protein_col]]) else row_values
  protein_values[is.na(protein_values) | protein_values == ""] <- row_values[is.na(protein_values) | protein_values == ""]
  gene_values <- if (!is.na(gene_col)) take_first(raw[[gene_col]]) else rep(NA_character_, length(row_values))

  qmat <- as.data.frame(lapply(raw[, qcols, drop = FALSE], take_numeric), check.names = FALSE)
  qmat[qmat == 0] <- NA_real_
  colnames(qmat) <- clean_dda_sample_name(qcols, "peaks")

  if (identical(schema$id, "db")) {
    coverage_mat <- as.data.frame(lapply(raw[, coverage_cols, drop = FALSE], take_numeric), check.names = FALSE)
    coverage_mat[coverage_mat == 0] <- NA_real_
    colnames(coverage_mat) <- clean_dda_sample_name(sub("^Coverage\\(%\\)\\s+", "Area ", coverage_cols), "peaks")
    common_samples <- intersect(colnames(qmat), colnames(coverage_mat))
    if (length(common_samples) == 0) stop("PEAKS DB Area sample columns and sample-specific Coverage(%) columns could not be matched by sample name.")
    qmat <- qmat[, common_samples, drop = FALSE]
    coverage_mat <- coverage_mat[, common_samples, drop = FALSE]
  } else {
    coverage_mat <- matrix(NA_real_, nrow = nrow(qmat), ncol = ncol(qmat), dimnames = dimnames(qmat))
  }

  keep <- !is.na(row_values) & row_values != ""
  if (!any(keep)) stop("No non-empty feature identifiers were found in the selected PEAKS identifier column: ", selected_id_col)
  meta <- data.frame(
    RowID = row_values[keep],
    ProteinName = protein_values[keep],
    GeneName = gene_values[keep],
    Accession = accession_values[keep],
    SourceFormat = "peaks",
    IdentifierColumn = selected_id_col,
    ProteinGroup = if ("Protein Group" %in% colnames(raw)) raw[["Protein Group"]][keep] else NA,
    Top = if ("Top" %in% colnames(raw)) raw[["Top"]][keep] else NA,
    OriginalRowIndex = keep_idx[keep],
    stringsAsFactors = FALSE
  )
  quantity <- qmat[keep, , drop = FALSE]
  ident <- coverage_mat[keep, , drop = FALSE]
  analysis_ids <- make.unique(meta$RowID)
  meta$AnalysisID <- analysis_ids
  rownames(quantity) <- analysis_ids
  rownames(ident) <- analysis_ids
  quantity <- as.matrix(quantity)
  ident <- as.matrix(ident)
  storage.mode(quantity) <- "numeric"
  storage.mode(ident) <- "numeric"
  counts <- make_sample_count_table(quantity, ident, identified_available = identical(schema$id, "db"))
  list(
    raw = raw,
    raw_all = raw_all,
    meta = meta,
    quantity = quantity,
    qualitative = ident,
    ibaq = NULL,
    counts = counts,
    software = "peaks",
    peaks_subtype = schema$id,
    input_family = "dda",
    input_source = "peaks",
    data_level = "protein",
    format_evidence = if (identical(schema$id, "db")) "PEAKS DB protein result: one row retained per Protein Group using Top == TRUE, first TRUE wins; Area zeros treated as missing; sample-specific Coverage(%) > 0 defines identification." else "PEAKS LFQ protein result: one row retained per Protein Group using Top == TRUE, first TRUE wins; sample-level <sample> Area values are used directly as LFQ abundance and zeros are treated as missing; sample-level identification evidence is unavailable.",
    count_approximation_note_en = if (identical(schema$id, "lfq")) "Sample-level identification counts are unavailable because this PEAKS LFQ protein result does not contain sample-specific identification evidence." else NULL,
    features = rownames(quantity),
    samples = colnames(quantity)
  )
}

extract_dda_protein_data <- function(file, software = c("fragpipe", "proteome_discoverer", "peaks", "maxquant"), row_id = c("protein_name", "gene_name", "accession")) {
  load_required_packages()
  software <- match.arg(software)
  row_id <- match.arg(row_id)
  if (software == "maxquant") return(extract_maxquant_protein_data(file, row_id))
  if (software == "fragpipe") return(extract_fragpipe_protein_data(file, row_id))
  if (software == "peaks") return(extract_peaks_protein_data(file, row_id))
  raw <- read_result_file(file)
  qcols <- dda_quantity_columns(raw, software)
  if (length(qcols) == 0) stop("No supported DDA quantitative columns were found for ", software, ".")

  ids <- dda_id_candidates(software)
  selected_id_col <- first_present_col(raw, ids[[row_id]])
  accession_col <- first_present_col(raw, ids[["accession"]])
  protein_col <- first_present_col(raw, ids[["protein_name"]])
  gene_col <- first_present_col(raw, ids[["gene_name"]])
  if (is.na(selected_id_col)) {
    stop("No suitable ", row_id, " identifier column was found for ", software, ".")
  }

  row_values <- take_first(raw[[selected_id_col]])
  accession_values <- if (!is.na(accession_col)) take_first(raw[[accession_col]]) else row_values
  protein_values <- if (!is.na(protein_col)) take_first(raw[[protein_col]]) else row_values
  gene_values <- if (!is.na(gene_col)) take_first(raw[[gene_col]]) else rep(NA_character_, length(row_values))
  qmat <- as.data.frame(lapply(raw[, qcols, drop = FALSE], take_numeric), check.names = FALSE)
  colnames(qmat) <- clean_dda_sample_name(qcols, software)

  keep <- !is.na(row_values) & row_values != ""
  if (!any(keep)) stop("No non-empty feature identifiers were found in the selected DDA identifier column: ", selected_id_col)
  meta <- data.frame(
    RowID = row_values[keep],
    ProteinName = protein_values[keep],
    GeneName = gene_values[keep],
    Accession = accession_values[keep],
    SourceFormat = software,
    IdentifierColumn = selected_id_col,
    stringsAsFactors = FALSE
  )
  quantity <- qmat[keep, , drop = FALSE]
  analysis_ids <- make.unique(meta$RowID)
  meta$AnalysisID <- analysis_ids
  rownames(quantity) <- analysis_ids
  quantity <- as.matrix(quantity)
  storage.mode(quantity) <- "numeric"
  counts <- make_sample_count_table(quantity, quantity)
  list(
    raw = raw,
    meta = meta,
    quantity = quantity,
    qualitative = quantity,
    ibaq = NULL,
    counts = counts,
    software = software,
    input_family = "dda",
    input_source = software,
    data_level = "protein",
    features = rownames(quantity),
    samples = colnames(quantity)
  )
}

load_input_dataset <- function(file, input_family = c("dia", "dda", "standard_matrix"), selected_format = "auto", diann_type = c("d", "raw"), row_id = c("protein_name", "gene_name", "accession"), standard_zero_mode = NULL) {
  input_family <- match.arg(input_family)
  row_id <- match.arg(row_id)
  if (input_family == "standard_matrix") {
    dat <- read_standard_matrix(file, standard_zero_mode)
    dat$input_family <- "standard_matrix"
    dat$format_label <- "Standard quantitative matrix"
    dat$format_evidence <- "User selected standard quantitative matrix input."
    return(dat)
  }
  if (identical(selected_format, "auto")) {
    detected <- detect_input_format(file, input_family)
    if (is.na(detected$id)) stop(detected$evidence)
    selected_format <- detected$id
  } else {
    detected <- detect_input_format(file, input_family)
  }
  if (input_family == "dia") {
    if (selected_format == "diann") {
      dat <- extract_protein_data(file, "DIANN", diann_type, row_id)
    } else if (selected_format == "spectronaut") {
      dat <- extract_protein_data(file, "Spectronaut", diann_type, row_id)
    } else {
      stop("Unknown DIA format: ", selected_format)
    }
  } else {
    dat <- extract_dda_protein_data(file, selected_format, row_id)
  }
  entry <- input_format_registry(input_family)
  labels <- stats::setNames(vapply(entry, `[[`, character(1), "label"), vapply(entry, `[[`, character(1), "id"))
  dat$input_family <- input_family
  dat$input_source <- selected_format
  dat$data_level <- dat$data_level %||% "protein"
  dat$format_label <- if (identical(selected_format, "peaks") && !is.null(dat$peaks_subtype)) {
    if (identical(dat$peaks_subtype, "db")) "PEAKS DB protein result" else "PEAKS LFQ protein result"
  } else labels[[selected_format]] %||% selected_format
  dat$format_evidence <- detected$evidence %||% "Format selected manually."
  dat
}

analysis_ids_to_accessions <- function(meta, ids) {
  if (is.null(meta) || !"Accession" %in% colnames(meta)) return(rep(NA_character_, length(ids)))
  key_col <- if ("AnalysisID" %in% colnames(meta)) "AnalysisID" else "RowID"
  meta$Accession[match(ids, meta[[key_col]])]
}
