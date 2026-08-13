# ProteoPostZ core functions
# Developed by Wenjia Zhang

`%||%` <- function(a, b) if (!is.null(a)) a else b

dependency_manifest <- list(
  cran_runtime = c(
    "data.table", "dplyr", "tidyr", "ggplot2", "pheatmap", "RColorBrewer",
    "VennDiagram", "UpSetR", "ggridges", "uwot", "randomForest", "glmnet",
    "Peptides", "pdftools", "tibble", "shiny", "bslib", "DT"
  ),
  bioconductor_runtime = c(
    "slingshot", "SingleCellExperiment", "S4Vectors", "limma"
  ),
  r_base_runtime = c(
    "grid", "grDevices", "graphics", "stats", "utils", "tools"
  ),
  cran_annotation_update = c(
    "httr"
  ),
  project_runtime = c(
    "app/R/analysis_core.R",
    "app/annotations/uniprot_all_celegans_6239_annotations.csv",
    "app/annotations/uniprot_reviewed_human_9606_annotations.csv",
    "app/annotations/uniprot_reviewed_mouse_10090_annotations.csv",
    "portable/R-4.5.1",
    "portable/Rlibs"
  ),
  system_runtime = c(
    "Windows shell for launcher scripts",
    "Microsoft .NET runtime for launcher executable",
    "local loopback TCP access to 127.0.0.1:3840",
    "native DLL dependencies bundled with portable R and compiled R packages"
  ),
  uncertain_system_runtime = c(
    "pdftools native PDF rendering dependencies when not fully bundled by the Windows binary package"
  )
)

required_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "pheatmap", "RColorBrewer",
  "VennDiagram", "UpSetR", "ggridges", "uwot", "randomForest", "glmnet",
  "Peptides", "grid", "pdftools", "slingshot", "SingleCellExperiment", "S4Vectors",
  "tibble"
)

load_required_packages <- function() {
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "), ". Please install them before running the app.")
  }
}

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

make_group_info <- function(samples, groups) {
  groups <- trimws(as.character(groups))
  groups[groups == ""] <- "Group1"
  data.frame(Sample = samples, Group = factor(groups, levels = unique(groups)), stringsAsFactors = FALSE)
}

write_matrix_csv <- function(mat, path, id_col = "ProteinID") {
  out <- data.frame(mat, check.names = FALSE)
  out <- cbind(stats::setNames(data.frame(rownames(mat), check.names = FALSE), id_col), out)
  data.table::fwrite(out, path)
}

plot_sample_count_bar <- function(counts, group_info, count_col, out_pdf, out_csv, width = 3.3, height = 3.3, palette = "npg", y_label = "Protein groups", plot_title = NULL) {
  if (!count_col %in% colnames(counts)) stop("Count column not found: ", count_col)
  df <- dplyr::left_join(counts, group_info, by = "Sample")
  summary <- df |>
    dplyr::group_by(Group) |>
    dplyr::summarise(
      Mean = if (all(is.na(.data[[count_col]]))) NA_real_ else mean(.data[[count_col]], na.rm = TRUE),
      SD = if (sum(!is.na(.data[[count_col]])) < 2) NA_real_ else stats::sd(.data[[count_col]], na.rm = TRUE),
      N = sum(!is.na(.data[[count_col]])),
      .groups = "drop"
    )
  summary$Count_Column <- count_col
  data.table::fwrite(df, sub("\\.csv$", "_sample_counts.csv", out_csv))
  data.table::fwrite(summary, out_csv)
  cols <- sci_palette(length(unique(df$Group)), palette)
  p <- ggplot2::ggplot(summary, ggplot2::aes(Group, Mean, fill = Group)) +
    ggplot2::geom_col(width = 0.65, color = "black", linewidth = 0.25) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.18, linewidth = 0.35) +
    ggplot2::geom_jitter(data = df, ggplot2::aes(x = Group, y = .data[[count_col]]), inherit.aes = FALSE, width = 0.09, size = 1.7) +
    ggplot2::scale_fill_manual(values = cols) + theme_sci() + ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = y_label, title = plot_title %||% y_label)
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
  p
}

plot_identification_bar <- function(counts, group_info, out_pdf, out_csv, width = 3.3, height = 3.3, palette = "npg") {
  plot_sample_count_bar(counts, group_info, "Identified_Protein_Count", out_pdf, out_csv, width, height, palette, "Protein groups")
}

plot_available_quantitative_bar <- function(counts, group_info, out_pdf, out_csv, width = 3.3, height = 3.3, palette = "npg") {
  plot_sample_count_bar(counts, group_info, "Available_Quantitative_Value_Count", out_pdf, out_csv, width, height, palette, "Number of available quantitative values")
}

sci_palette <- function(n, palette = "npg") {
  sets <- list(
    npg = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2"),
    lancet = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A"),
    jama = c("#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"),
    nejm = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91"),
    uchicago = c("#800000", "#767676", "#FFA319", "#8A9045", "#155F83", "#C16622", "#58593F")
  )
  pal <- sets[[palette]] %||% sets$npg
  rep(pal, length.out = n)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

theme_sci <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", size = 8),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      axis.text = ggplot2::element_text(size = 8, color = "black"),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
      strip.text = ggplot2::element_text(size = 8)
    )
}

identified_by_group <- function(qual_mat, group_info, min_reps = 1) {
  groups <- levels(group_info$Group)
  sets <- lapply(groups, function(g) {
    samples <- group_info$Sample[group_info$Group == g]
    samples <- intersect(samples, colnames(qual_mat))
    rownames(qual_mat)[rowSums(!is.na(qual_mat[, samples, drop = FALSE])) >= min_reps]
  })
  names(sets) <- groups
  sets
}

validate_group_sets <- function(sets, min_reps = 1, analysis = c("venn", "upset")) {
  analysis <- match.arg(analysis)
  set_sizes <- vapply(sets, length, integer(1))
  if (length(set_sizes) == 0 || all(set_sizes == 0)) {
    stop("All group-level protein sets are empty. Minimum replicates detected in group may be too high, or no proteins in the current data meet the detection rule.")
  }
  nonempty <- sets[set_sizes > 0]
  if (length(nonempty) < 2) {
    stop("At least 2 non-empty group-level protein sets are required for Venn/UpSet plots.")
  }
  if (analysis == "venn" && length(nonempty) > 4) {
    stop("Venn diagram is recommended for 2-4 non-empty groups. Please use UpSet plot for 5 or more groups.")
  }
  nonempty
}

make_set_membership <- function(sets) {
  membership <- unique(unlist(sets, use.names = FALSE))
  if (length(membership) == 0) {
    return(data.frame(ProteinID = character(), stringsAsFactors = FALSE))
  }
  set_df <- data.frame(ProteinID = membership, stringsAsFactors = FALSE)
  for (nm in names(sets)) set_df[[nm]] <- as.integer(membership %in% sets[[nm]])
  set_df
}

plot_venn_upset <- function(sets, outdir, width = 3.3, height = 3.3) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  sets <- validate_group_sets(sets, analysis = "upset")
  set_df <- make_set_membership(sets)
  data.table::fwrite(set_df, file.path(outdir, "group_detected_membership.csv"))
  if (length(sets) >= 2 && length(sets) <= 4) {
    grDevices::pdf(file.path(outdir, "venn.pdf"), width = width, height = height)
    grid::grid.draw(VennDiagram::venn.diagram(sets, filename = NULL, fill = sci_palette(length(sets)), alpha = 0.45, cex = 0.8, cat.cex = 0.8, margin = 0.08))
    grDevices::dev.off()
  }
  if (nrow(set_df) == 0 || ncol(set_df) - 1 < 2) {
    stop("UpSet plot requires a membership table with at least one protein and at least 2 non-empty group columns.")
  }
  grDevices::pdf(file.path(outdir, "upset.pdf"), width = max(width, 5), height = max(height, 4))
  UpSetR::upset(as.data.frame(set_df[, -1, drop = FALSE]), nsets = length(sets), order.by = "freq")
  grDevices::dev.off()
}

preprocess_expr <- function(mat, log2_transform = TRUE, min_valid_fraction = 0.5) {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  if (log2_transform) mat <- log2(mat + 1)
  keep_n <- max(1, ceiling(ncol(mat) * min_valid_fraction))
  mat <- mat[rowSums(!is.na(mat)) >= keep_n, , drop = FALSE]
  imp <- t(apply(mat, 1, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))
  imp[complete.cases(imp), , drop = FALSE]
}

plot_correlation_heatmap <- function(mat, group_info, out_pdf, out_csv, method = "pearson", order_mode = "original", cluster_distance = "one_minus_correlation", cluster_linkage = "complete", cluster_k = NULL, digits = 2, fontsize_number = 8, color_scheme = "blue_white_red", min_cor = -1, max_cor = 1, width = 3.3, height = 3.3) {
  used <- preprocess_expr(mat, TRUE, 0.3)
  cor_mat <- cor(used, method = method, use = "pairwise.complete.obs")
  if (!is.finite(min_cor) || !is.finite(max_cor) || min_cor >= max_cor) {
    stop("Correlation legend min must be smaller than max.")
  }
  if (ncol(cor_mat) < 2) stop("Correlation heatmap requires at least two samples.")

  block_gaps <- function(labels) {
    if (length(labels) < 2) return(NULL)
    gaps <- which(labels[-1] != labels[-length(labels)])
    if (length(gaps) == 0) NULL else gaps
  }

  plot_mat <- cor_mat
  plot_cluster_rows <- FALSE
  plot_cluster_cols <- FALSE
  plot_gaps <- NULL
  if (order_mode == "group") {
    ord <- unlist(lapply(levels(group_info$Group), function(g) {
      smp <- group_info$Sample[group_info$Group == g]
      intersect(smp, colnames(cor_mat))
    }))
    cor_mat <- cor_mat[ord, ord, drop = FALSE]
    plot_mat <- cor_mat
    plot_gaps <- block_gaps(as.character(group_info$Group[match(colnames(cor_mat), group_info$Sample)]))
  } else if (order_mode == "hclust") {
    if (is.null(cluster_k)) cluster_k <- length(levels(group_info$Group))
    cluster_k <- as.integer(cluster_k)
    if (!is.finite(cluster_k) || cluster_k < 1 || cluster_k > ncol(cor_mat)) {
      stop("Hierarchical clustering k must be between 1 and the number of samples.")
    }
    if (cluster_distance == "one_minus_correlation") {
      distance_matrix <- 1 - cor_mat
      diag(distance_matrix) <- 0
      if (any(!is.finite(distance_matrix[upper.tri(distance_matrix)]))) {
        stop("Hierarchical clustering requires finite pairwise correlations. Check for samples with insufficient valid quantitative values.")
      }
      distance_object <- as.dist(pmax(distance_matrix, 0))
    } else if (cluster_distance == "correlation_euclidean") {
      if (any(!is.finite(cor_mat))) {
        stop("Hierarchical clustering requires finite correlations. Check for samples with insufficient valid quantitative values.")
      }
      distance_object <- stats::dist(cor_mat, method = "euclidean")
    } else {
      stop("Unsupported hierarchical clustering distance.")
    }
    hc <- stats::hclust(distance_object, method = cluster_linkage)
    cor_mat <- cor_mat[hc$order, hc$order, drop = FALSE]
    plot_cluster_rows <- hc
    plot_cluster_cols <- hc
    plot_gaps <- cluster_k
  }
  sample_group <- as.character(group_info$Group[match(colnames(cor_mat), group_info$Sample)])
  data.table::fwrite(data.frame(Sample = rownames(cor_mat), Group = sample_group, cor_mat, check.names = FALSE), out_csv)
  long_csv <- sub("\\.csv$", "_long.csv", out_csv)
  long_df <- as.data.frame(cor_mat) |>
    tibble::rownames_to_column(".cor_row_sample") |>
    tidyr::pivot_longer(-.cor_row_sample, names_to = ".cor_col_sample", values_to = "Correlation") |>
    dplyr::mutate(
      Group1 = as.character(group_info$Group[match(.cor_row_sample, group_info$Sample)]),
      Group2 = as.character(group_info$Group[match(.cor_col_sample, group_info$Sample)]),
      Method = method,
      OrderedSample1 = match(.cor_row_sample, rownames(cor_mat)),
      OrderedSample2 = match(.cor_col_sample, colnames(cor_mat))
    ) |>
    dplyr::transmute(Method, Sample1 = .cor_row_sample, Group1, OrderedSample1, Sample2 = .cor_col_sample, Group2, OrderedSample2, Correlation)
  data.table::fwrite(long_df, long_csv)
  ann <- data.frame(Group = group_info$Group[match(colnames(plot_mat), group_info$Sample)])
  rownames(ann) <- colnames(plot_mat)
  cols <- if (color_scheme == "purple_white_orange") colorRampPalette(c("#5E3C99", "white", "#E66101"))(100) else colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  breaks <- seq(min_cor, max_cor, length.out = 101)
  nums <- matrix(sprintf(paste0("%.", digits, "f"), plot_mat), nrow(plot_mat), dimnames = dimnames(plot_mat))
  grDevices::pdf(out_pdf, width = width, height = height)
  pheatmap::pheatmap(plot_mat, color = cols, breaks = breaks, cluster_rows = plot_cluster_rows, cluster_cols = plot_cluster_cols, cutree_rows = if (order_mode == "hclust") plot_gaps else NA, cutree_cols = if (order_mode == "hclust") plot_gaps else NA, gaps_row = if (order_mode == "group") plot_gaps else NULL, gaps_col = if (order_mode == "group") plot_gaps else NULL, annotation_col = ann, annotation_row = ann, display_numbers = nums, fontsize_number = fontsize_number, border_color = NA)
  grDevices::dev.off()
  cor_mat
}
plot_rank_abundance <- function(mat, group_info, out_pdf, out_csv, width = 3.3, height = 3.3, palette = "npg") {
  df <- as.data.frame(mat) |>
    tibble::rownames_to_column("ProteinID") |>
    tidyr::pivot_longer(-ProteinID, names_to = "Sample", values_to = "Intensity") |>
    dplyr::filter(!is.na(Intensity), Intensity > 0) |>
    dplyr::group_by(Sample) |>
    dplyr::arrange(dplyr::desc(Intensity), .by_group = TRUE) |>
    dplyr::mutate(Rank = dplyr::row_number(), Log2Intensity = log2(Intensity + 1)) |>
    dplyr::ungroup() |>
    dplyr::left_join(group_info, by = "Sample")
  data.table::fwrite(df, out_csv)
  p <- ggplot2::ggplot(df, ggplot2::aes(Rank, Log2Intensity, color = Group, group = Sample)) + ggplot2::geom_line(alpha = 0.55, linewidth = 0.35) + theme_sci() + ggplot2::scale_color_manual(values = sci_palette(length(levels(group_info$Group)), palette)) + ggplot2::labs(y = "log2(intensity + 1)")
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
}

plot_cv_ridges <- function(mat, group_info, out_pdf, out_csv, cv_max = 60, width = 3.3, height = 3.3, palette = "npg") {
  rows <- lapply(levels(group_info$Group), function(g) {
    smp <- intersect(group_info$Sample[group_info$Group == g], colnames(mat))
    if (length(smp) < 2) return(NULL)
    vals <- mat[, smp, drop = FALSE]
    cv <- apply(vals, 1, function(x) if (sum(!is.na(x)) >= 2 && mean(x, na.rm = TRUE) > 0) 100 * sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) else NA_real_)
    data.frame(ProteinID = rownames(mat), Group = g, CV = cv)
  })
  df <- dplyr::bind_rows(rows) |> dplyr::filter(!is.na(CV), CV >= 0, CV <= cv_max)
  med <- df |> dplyr::group_by(Group) |> dplyr::summarise(MedianCV = median(CV), .groups = "drop")
  data.table::fwrite(df, out_csv)
  data.table::fwrite(med, sub("\\.csv$", "_median.csv", out_csv))
  p <- ggplot2::ggplot(df, ggplot2::aes(CV, Group, fill = Group)) + ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) + ggplot2::geom_segment(data = med, ggplot2::aes(x = MedianCV, xend = MedianCV, y = as.numeric(factor(Group)) - 0.35, yend = as.numeric(factor(Group)) + 0.35), inherit.aes = FALSE, linetype = "dashed") + ggplot2::coord_cartesian(xlim = c(0, cv_max)) + ggplot2::scale_fill_manual(values = sci_palette(length(levels(group_info$Group)), palette)) + theme_sci() + ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
}

run_pca_umap <- function(mat, group_info, outdir, prefix = "all", top_features = NULL, width = 3.3, height = 3.3, palette = "npg", n_neighbors = 10, min_dist = 0.1, seed = 123) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  used <- preprocess_expr(mat, TRUE, 0.5)
  if (!is.null(top_features)) used <- used[intersect(top_features, rownames(used)), , drop = FALSE]
  sample_mat <- t(used)
  pca <- prcomp(sample_mat, center = TRUE, scale. = TRUE)
  var <- summary(pca)$importance[2, 1:2] * 100
  pca_df <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2]) |> dplyr::left_join(group_info, by = "Sample")
  data.table::fwrite(pca_df, file.path(outdir, paste0(prefix, "_PCA_coordinates.csv")))
  cols <- sci_palette(length(levels(group_info$Group)), palette)
  p1 <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2, color = Group)) + ggplot2::geom_point(size = 2.4) + ggplot2::scale_color_manual(values = cols) + theme_sci() + ggplot2::labs(x = sprintf("PC1 (%.2f%%)", var[1]), y = sprintf("PC2 (%.2f%%)", var[2]))
  ggplot2::ggsave(file.path(outdir, paste0(prefix, "_PCA.pdf")), p1, width = width, height = height)
  set.seed(seed)
  nn <- min(n_neighbors, max(2, nrow(sample_mat) - 1))
  um <- uwot::umap(sample_mat, n_neighbors = nn, min_dist = min_dist, metric = "euclidean", verbose = FALSE)
  um_df <- data.frame(Sample = rownames(sample_mat), UMAP1 = um[,1], UMAP2 = um[,2]) |> dplyr::left_join(group_info, by = "Sample")
  data.table::fwrite(um_df, file.path(outdir, paste0(prefix, "_UMAP_coordinates.csv")))
  p2 <- ggplot2::ggplot(um_df, ggplot2::aes(UMAP1, UMAP2, color = Group)) + ggplot2::geom_point(size = 2.4) + ggplot2::scale_color_manual(values = cols) + theme_sci()
  ggplot2::ggsave(file.path(outdir, paste0(prefix, "_UMAP.pdf")), p2, width = width, height = height)
}

heatmap_scheme_colors <- function(scheme = "blue_white_red", n = 100) {
  anchors <- switch(scheme,
    blue_white_red = c("#2166AC", "white", "#B2182B"),
    red_white_blue = c("#B2182B", "white", "#2166AC"),
    purple_white_orange = c("#5E3C99", "white", "#E66101"),
    green_black_red = c("#1B9E77", "#000000", "#D73027"),
    stop("Unsupported heatmap color scheme.")
  )
  grDevices::colorRampPalette(anchors)(n)
}

validate_annotation_colors <- function(colors, labels, annotation_name) {
  if (length(labels) == 0) return(character())
  colors <- as.character(colors)
  if (is.null(names(colors))) names(colors) <- labels
  colors <- colors[labels]
  if (any(is.na(colors)) || any(!grepl("^#[0-9A-Fa-f]{6}$", colors))) {
    stop(annotation_name, " colors must be valid #RRGGBB values for every displayed label.")
  }
  colors
}

plot_expression_heatmap <- function(mat, group_info, out_pdf, out_csv, top_n = 100, row_cluster = "hclust", col_cluster = "hclust", row_kmeans_k = 4, col_kmeans_k = 4, row_hclust_k = 2, col_hclust_k = 2, color_scheme = "blue_white_red", group_annotation_colors = NULL, cluster_annotation_colors = NULL, row_labels = NULL, width = 4.5, height = 4.5) {
  used <- preprocess_expr(mat, TRUE, 0.5)
  if (nrow(used) < 2 || ncol(used) < 2) stop("Expression heatmap requires at least two proteins and two samples after filtering.")
  vars <- apply(used, 1, var, na.rm = TRUE)
  selected <- rownames(used)[order(vars, decreasing = TRUE)[seq_len(min(top_n, length(vars)))]]
  used <- used[rownames(used) %in% selected, , drop = FALSE]
  scaled <- t(scale(t(used)))
  scaled[!is.finite(scaled)] <- 0
  group_labels <- as.character(group_info$Group[match(colnames(scaled), group_info$Sample)])
  if (anyNA(group_labels)) stop("Some heatmap samples do not have group labels.")
  group_levels <- levels(group_info$Group)
  group_colors <- validate_annotation_colors(group_annotation_colors %||% stats::setNames(sci_palette(length(group_levels), "npg"), group_levels), group_levels, "Group annotation")

  safe_k <- function(k, max_k, label) {
    k <- suppressWarnings(as.integer(k))
    if (!is.finite(k) || k < 1) stop(label, " must be at least 1.")
    min(k, max_k)
  }
  cluster_gaps <- function(cluster) {
    if (length(cluster) < 2) return(NULL)
    gaps <- which(cluster[-1] != cluster[-length(cluster)])
    if (length(gaps) == 0) NULL else gaps
  }
  set.seed(123)
  row_km <- if (row_cluster == "kmeans") kmeans(scaled, centers = safe_k(row_kmeans_k, nrow(scaled), "Row K-means k"))$cluster else rep(NA_integer_, nrow(scaled))
  set.seed(123)
  col_km <- if (col_cluster == "kmeans") kmeans(t(scaled), centers = safe_k(col_kmeans_k, ncol(scaled), "Column K-means k"))$cluster else rep(NA_integer_, ncol(scaled))
  row_order <- if (row_cluster == "kmeans") rownames(scaled)[order(row_km, seq_along(row_km))] else rownames(scaled)
  col_order <- if (col_cluster == "kmeans") colnames(scaled)[order(col_km, seq_along(col_km))] else colnames(scaled)
  plot_scaled <- scaled[row_order, col_order, drop = FALSE]
  plot_row_labels <- NULL
  if (!is.null(row_labels)) {
    if (is.null(names(row_labels))) stop("Heatmap row labels must be named by the current protein identifiers.")
    plot_row_labels <- as.character(row_labels[match(rownames(plot_scaled), names(row_labels))])
    plot_row_labels[is.na(plot_row_labels)] <- ""
  }
  row_kmeans_gaps <- if (row_cluster == "kmeans") cluster_gaps(row_km[match(rownames(plot_scaled), rownames(scaled))]) else NULL
  col_kmeans_gaps <- if (col_cluster == "kmeans") cluster_gaps(col_km[match(colnames(plot_scaled), colnames(scaled))]) else NULL

  row_hc <- if (row_cluster == "hclust") stats::hclust(stats::dist(plot_scaled), method = "complete") else FALSE
  col_hc <- if (col_cluster == "hclust") stats::hclust(stats::dist(t(plot_scaled)), method = "complete") else FALSE
  row_hc_k <- if (row_cluster == "hclust") safe_k(row_hclust_k, nrow(plot_scaled), "Row hierarchical clustering k") else NA_integer_
  col_hc_k <- if (col_cluster == "hclust") safe_k(col_hclust_k, ncol(plot_scaled), "Column hierarchical clustering k") else NA_integer_
  row_hc_cluster <- if (row_cluster == "hclust") stats::cutree(row_hc, k = row_hc_k) else rep(NA_integer_, nrow(plot_scaled))
  col_hc_cluster <- if (col_cluster == "hclust") stats::cutree(col_hc, k = col_hc_k) else rep(NA_integer_, ncol(plot_scaled))
  col_cluster_values <- if (col_cluster == "hclust") col_hc_cluster else if (col_cluster == "kmeans") col_km[match(colnames(plot_scaled), colnames(scaled))] else rep(NA_integer_, ncol(plot_scaled))
  cluster_levels <- if (col_cluster == "none") character() else paste0("Cluster ", sort(unique(col_cluster_values)))
  ann <- data.frame(Group = factor(group_labels[match(colnames(plot_scaled), colnames(scaled))], levels = group_levels), stringsAsFactors = FALSE)
  if (length(cluster_levels) > 0) ann$Cluster <- factor(paste0("Cluster ", col_cluster_values), levels = cluster_levels)
  rownames(ann) <- colnames(plot_scaled)
  annotation_colors <- list(Group = group_colors)
  if (length(cluster_levels) > 0) {
    cluster_colors <- validate_annotation_colors(cluster_annotation_colors %||% stats::setNames(sci_palette(length(cluster_levels), "lancet"), cluster_levels), cluster_levels, "Cluster annotation")
    annotation_colors$Cluster <- cluster_colors
  }
  # Reserve vertical space for column annotations and plot margins, then keep each
  # visible label below the available per-row height to avoid overlap.
  row_label_fontsize <- if (is.null(plot_row_labels)) {
    10
  } else {
    available_row_height_pt <- max(1, height * 72 - 115) / nrow(plot_scaled)
    max(1.5, min(10, 0.75 * available_row_height_pt))
  }
  grDevices::pdf(out_pdf, width = width, height = height)
  ph <- pheatmap::pheatmap(plot_scaled, scale = "none", color = heatmap_scheme_colors(color_scheme), cluster_rows = row_hc, cluster_cols = col_hc, cutree_rows = row_hc_k, cutree_cols = col_hc_k, gaps_row = row_kmeans_gaps, gaps_col = col_kmeans_gaps, kmeans_k = NA, annotation_col = ann, annotation_colors = annotation_colors, labels_row = plot_row_labels, show_rownames = !is.null(plot_row_labels), fontsize_row = row_label_fontsize, border_color = NA)
  grDevices::dev.off()
  row_order <- if (row_cluster == "hclust" && !is.null(ph$tree_row)) {
    rownames(plot_scaled)[ph$tree_row$order]
  } else {
    row_order
  }
  col_order <- if (col_cluster == "hclust" && !is.null(ph$tree_col)) colnames(plot_scaled)[ph$tree_col$order] else col_order
  ordered_scaled <- scaled[row_order, col_order, drop = FALSE]
  row_hc_ordered <- if (row_cluster == "hclust") row_hc_cluster[match(row_order, rownames(plot_scaled))] else rep(NA_integer_, length(row_order))
  col_hc_ordered <- if (col_cluster == "hclust") col_hc_cluster[match(col_order, colnames(plot_scaled))] else rep(NA_integer_, length(col_order))
  data.table::fwrite(data.frame(ProteinID = rownames(ordered_scaled), KmeansCluster = row_km[match(rownames(ordered_scaled), rownames(scaled))], HierarchicalCluster = row_hc_ordered, ordered_scaled, check.names = FALSE), out_csv)
  data.table::fwrite(data.frame(RowOrder = seq_along(row_order), ProteinID = row_order, KmeansCluster = row_km[match(row_order, rownames(scaled))], HierarchicalCluster = row_hc_ordered), sub("\\.csv$", "_row_order.csv", out_csv))
  data.table::fwrite(data.frame(ColOrder = seq_along(col_order), Sample = col_order, Group = as.character(group_info$Group[match(col_order, group_info$Sample)]), KmeansCluster = col_km[match(col_order, colnames(scaled))], HierarchicalCluster = col_hc_ordered), sub("\\.csv$", "_col_order.csv", out_csv))
}
run_volcano <- function(mat, group_info, group_a, group_b, out_pdf, out_csv, log2fc_cutoff = 1, adj_p_cutoff = 0.05, raw_p_cutoff = 0.05, fc_method = c("log2_then_diff", "ratio_then_log2"), test_method = c("limma", "ttest"), sig_metric = c("adj_p", "raw_p"), width = 3.3, height = 3.3) {
  fc_method <- match.arg(fc_method)
  test_method <- match.arg(test_method)
  sig_metric <- match.arg(sig_metric)
  used_log2 <- log2(mat + 1)
  a <- intersect(group_info$Sample[group_info$Group == group_a], colnames(used_log2)); b <- intersect(group_info$Sample[group_info$Group == group_b], colnames(used_log2))
  if (length(a) < 2 || length(b) < 2) stop("Volcano requires at least two samples in each selected group.")
  mean_a_log2 <- rowMeans(used_log2[, a, drop = FALSE], na.rm = TRUE)
  mean_b_log2 <- rowMeans(used_log2[, b, drop = FALSE], na.rm = TRUE)
  mean_a_raw <- rowMeans(mat[, a, drop = FALSE], na.rm = TRUE)
  mean_b_raw <- rowMeans(mat[, b, drop = FALSE], na.rm = TRUE)
  log2fc <- if (fc_method == "log2_then_diff") {
    mean_b_log2 - mean_a_log2
  } else {
    suppressWarnings(log2(mean_b_raw / mean_a_raw))
  }
  log2fc[!is.finite(log2fc)] <- NA_real_
  res <- data.frame(
    ProteinID = rownames(used_log2),
    MeanA_Log2 = mean_a_log2,
    MeanB_Log2 = mean_b_log2,
    MeanA_Raw = mean_a_raw,
    MeanB_Raw = mean_b_raw,
    log2FC = log2fc,
    stringsAsFactors = FALSE
  )
  if (test_method == "limma") {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("Volcano plot with limma requires the limma package.")
    }
    sample_group <- factor(c(rep(group_a, length(a)), rep(group_b, length(b))), levels = c(group_a, group_b))
    design <- stats::model.matrix(~ 0 + sample_group)
    colnames(design) <- c(group_a, group_b)
    contrast <- limma::makeContrasts(contrasts = paste0("`", group_b, "`-`", group_a, "`"), levels = design)
    fit <- limma::lmFit(used_log2[, c(a, b), drop = FALSE], design)
    fit2 <- limma::contrasts.fit(fit, contrast)
    fit2 <- limma::eBayes(fit2)
    res$P.Value <- fit2$p.value[, 1]
    p_value_method <- "limma moderated p value"
  } else {
    res$P.Value <- apply(used_log2, 1, function(x) tryCatch(t.test(x[b], x[a])$p.value, error = function(e) NA_real_))
    p_value_method <- "two-sample t-test p value"
  }
  res$BH.Adjusted.P.Value <- p.adjust(res$P.Value, method = "BH")
  sig_values <- if (sig_metric == "adj_p") res$BH.Adjusted.P.Value else res$P.Value
  sig_cutoff <- if (sig_metric == "adj_p") adj_p_cutoff else raw_p_cutoff
  sig_label <- if (sig_metric == "adj_p") "BH-adjusted p value" else "p value"
  res$FCMethod <- fc_method
  res$TestMethod <- test_method
  res$PValueMethod <- p_value_method
  res$AdjustedPValueMethod <- "BH-adjusted p value"
  res$SignificanceMetric <- if (sig_metric == "adj_p") "BH-adjusted p value" else "p value"
  res$SignificanceValue <- sig_values
  res$Regulation <- dplyr::case_when(!is.na(sig_values) & sig_values < sig_cutoff & res$log2FC >= log2fc_cutoff ~ "Up", !is.na(sig_values) & sig_values < sig_cutoff & res$log2FC <= -log2fc_cutoff ~ "Down", TRUE ~ "NotSig")
  data.table::fwrite(res, out_csv)
  plot_df <- res
  plot_df$NegLog10Significance <- -log10(plot_df$SignificanceValue)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(log2FC, NegLog10Significance, color = Regulation)) + ggplot2::geom_point(size = 1.2, alpha = 0.75, na.rm = TRUE) + ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") + ggplot2::geom_hline(yintercept = -log10(sig_cutoff), linetype = "dashed") + ggplot2::scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC", NotSig = "grey75")) + theme_sci() + ggplot2::labs(x = "log2FC", y = paste0("-log10(", sig_label, ")"))
  ggplot2::ggsave(out_pdf, p, width = width, height = height)
}

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
  if (length(training$test) > 0) {
    x_train <- x[training$train, , drop = FALSE]
    y_train <- droplevels(y[training$train])
    foldid_train <- make_stratified_foldid(y_train, min(nfolds, min(table(y_train))), seed)
    fit <- glmnet::cv.glmnet(x_train, y_train, family = "multinomial", alpha = l1_alpha, type.measure = "class", nfolds = max(foldid_train), foldid = foldid_train)
    probs <- predict(fit, x[training$test, , drop = FALSE], s = lambda_selection, type = "response")
    probs <- as.data.frame(probs[, , 1, drop = FALSE][, , 1], check.names = FALSE)
    test_y <- droplevels(y[training$test])
    pred <- factor(colnames(probs)[max.col(probs, ties.method = "first")], levels = levels(test_y))
    write_classification_evaluation(test_y, pred, probs, outdir, "l1_test")
  }
  foldid <- make_cv_folds(y, seed, max_folds = min(5, nfolds))
  pred <- rep(NA_character_, length(y))
  probs_all <- matrix(NA_real_, nrow = length(y), ncol = length(levels(y)), dimnames = list(names(y), levels(y)))
  for (fold in sort(unique(foldid))) {
    train_idx <- which(foldid != fold)
    test_idx <- which(foldid == fold)
    y_train <- droplevels(y[train_idx])
    inner_folds <- min(nfolds, min(table(y_train)))
    inner_foldid <- make_stratified_foldid(y_train, inner_folds, seed + fold)
    set.seed(seed + fold)
    fit <- glmnet::cv.glmnet(x[train_idx, , drop = FALSE], y_train, family = "multinomial", alpha = l1_alpha, type.measure = "class", nfolds = inner_folds, foldid = inner_foldid)
    probs <- predict(fit, x[test_idx, , drop = FALSE], s = lambda_selection, type = "response")
    probs <- as.data.frame(probs[, , 1, drop = FALSE][, , 1], check.names = FALSE)
    pred[test_idx] <- colnames(probs)[max.col(probs, ties.method = "first")]
    probs_all[test_idx, colnames(probs)] <- as.matrix(probs)
  }
  write_classification_evaluation(y, factor(pred, levels = levels(y)), as.data.frame(probs_all, check.names = FALSE), outdir, "l1_cross_validation")
  data.table::fwrite(data.frame(Sample = rownames(x), Fold = foldid, Actual = as.character(y), Predicted = pred, probs_all, check.names = FALSE), file.path(outdir, "l1_cross_validation_predictions.csv"))
}

write_rf_stability <- function(x, y, outdir, top_n = 50, rf_ntree = 500, mtry = NA_integer_, seed = 123, repeats = 20, sample_fraction = 0.8) {
  selected <- list()
  for (i in seq_len(repeats)) {
    set.seed(seed + 1000 + i)
    idx <- unlist(lapply(levels(y), function(cls) {
      cls_idx <- which(y == cls)
      sample(cls_idx, max(2, floor(length(cls_idx) * sample_fraction)), replace = FALSE)
    }), use.names = FALSE)
    args <- list(x = x[idx, , drop = FALSE], y = droplevels(y[idx]), ntree = rf_ntree, importance = TRUE)
    if (!is.na(mtry)) args$mtry <- mtry
    fit <- do.call(randomForest::randomForest, args)
    imp <- randomForest::importance(fit)
    col <- if ("MeanDecreaseGini" %in% colnames(imp)) "MeanDecreaseGini" else tail(colnames(imp), 1)
    selected[[i]] <- rownames(imp)[order(imp[, col], decreasing = TRUE)][seq_len(min(top_n, nrow(imp)))]
  }
  freq <- sort(table(unlist(selected, use.names = FALSE)) / repeats, decreasing = TRUE)
  out <- data.frame(ProteinID = names(freq), SelectionFrequency = as.numeric(freq), Repeats = repeats, stringsAsFactors = FALSE)
  data.table::fwrite(out, file.path(outdir, "random_forest_stability_selection.csv"))
  invisible(out)
}

write_l1_stability <- function(x, y, outdir, top_n = 50, l1_alpha = 1, lambda_selection = "lambda.1se", seed = 123, repeats = 20, sample_fraction = 0.8) {
  selected <- list()
  for (i in seq_len(repeats)) {
    set.seed(seed + 2000 + i)
    idx <- unlist(lapply(levels(y), function(cls) {
      cls_idx <- which(y == cls)
      sample(cls_idx, min(length(cls_idx), max(2, floor(length(cls_idx) * sample_fraction))), replace = FALSE)
    }), use.names = FALSE)
    y_sub <- droplevels(y[idx])
    if (min(table(y_sub)) < 2) next
    nfolds <- min(5, min(table(y_sub)))
    foldid <- make_stratified_foldid(y_sub, nfolds, seed + i)
    fit <- glmnet::cv.glmnet(x[idx, , drop = FALSE], y_sub, family = "multinomial", alpha = l1_alpha, type.measure = "class", nfolds = nfolds, foldid = foldid)
    co <- coef(fit, s = lambda_selection)
    ids <- unique(unlist(lapply(co, function(m) rownames(as.matrix(m))[as.numeric(as.matrix(m)[, 1]) != 0]), use.names = FALSE))
    ids <- setdiff(ids, "(Intercept)")
    selected[[i]] <- head(ids, top_n)
  }
  freq <- sort(table(unlist(selected, use.names = FALSE)) / repeats, decreasing = TRUE)
  out <- data.frame(ProteinID = names(freq), SelectionFrequency = as.numeric(freq), Repeats = repeats, stringsAsFactors = FALSE)
  data.table::fwrite(out, file.path(outdir, "l1_stability_selection.csv"))
  invisible(out)
}

run_random_forest_selection <- function(mat, group_info, outdir, top_n = 50, rf_ntree = 500, seed = 123, split_mode = "auto", train_prop = 0.7, rf_mtry = NA, allow_small_sample = FALSE) {
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
  write_rf_stability(x[training$train, , drop = FALSE], droplevels(y[training$train]), outdir, top_n = top_n, rf_ntree = rf_ntree, mtry = mtry, seed = seed, repeats = 20)
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
    EvaluationOutputs = "confusion_matrix, class_metrics, summary_metrics, ROC/AUC for binary groups, stability_selection",
    SmallSampleExploratoryML = allow_small_sample,
    ReliabilityNote = if (allow_small_sample) "Exploratory only: no independent test set; feature selection may be unstable." else "",
    AutoNote = training$auto_note %||% ""
  ))
  top <- head(imp$ProteinID, min(top_n, nrow(imp)))
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
    if (allow_small_sample) min(5, min_class_n) else min(10, min_class_n)
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

run_l1_selection <- function(mat, group_info, outdir, top_n = 50, l1_alpha = 1, seed = 123, split_mode = "auto", train_prop = 0.7, lambda_selection = "lambda.1se", cv_folds = "Auto", allow_small_sample = FALSE) {
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
  cv_expr <- quote(glmnet::cv.glmnet(x_train, y_train, family = "multinomial", alpha = l1_alpha, type.measure = "class", nfolds = nfolds, foldid = foldid, grouped = cv_settings$grouped))
  cv <- if (allow_small_sample) {
    withCallingHandlers(eval(cv_expr), warning = function(w) {
      glmnet_warnings <<- c(glmnet_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  } else {
    eval(cv_expr)
  }
  write_l1_evaluation(x, y, training, l1_alpha, lambda_selection, nfolds, outdir, seed)
  co <- coef(cv, s = lambda_selection)
  coef_df <- dplyr::bind_rows(lapply(names(co), function(cls) {
    m <- as.matrix(co[[cls]])
    data.frame(ProteinID = rownames(m), Class = cls, Coefficient = as.numeric(m[, 1]), row.names = NULL)
  })) |>
    dplyr::filter(ProteinID != "(Intercept)", Coefficient != 0)
  data.table::fwrite(coef_df, file.path(outdir, "l1_nonzero_coefficients.csv"))
  scores <- coef_df |>
    dplyr::group_by(ProteinID) |>
    dplyr::summarise(L1Score = sum(abs(Coefficient)), NonzeroClasses = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(L1Score))
  data.table::fwrite(scores, file.path(outdir, "l1_feature_scores.csv"))
  write_l1_stability(x_train, y_train, outdir, top_n = top_n, l1_alpha = l1_alpha, lambda_selection = lambda_selection, seed = seed, repeats = 20)
  top <- head(scores$ProteinID, min(top_n, nrow(scores)))
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
    CrossValidationFolds = nfolds,
    RequestedCrossValidationFolds = cv_folds,
    GroupedCV = cv_settings$grouped,
    EvaluationOutputs = "confusion_matrix, class_metrics, summary_metrics, ROC/AUC for binary groups, stability_selection",
    SmallSampleExploratoryML = allow_small_sample,
    GlmnetWarnings = if (length(glmnet_warnings) > 0) paste(unique(glmnet_warnings), collapse = " | ") else "",
    ReliabilityNote = if (allow_small_sample) "Exploratory only: no independent test set; feature selection may be unstable." else "",
    AutoNote = training$auto_note %||% ""
  ))
  top
}

run_feature_selection <- function(mat, group_info, outdir, top_n = 50, rf_ntree = 500, l1_alpha = 1, seed = 123, split_mode = "auto", train_prop = 0.7, rf_mtry = NA, lambda_selection = "lambda.1se", cv_folds = "Auto", allow_small_sample = FALSE) {
  rf_top <- run_random_forest_selection(mat, group_info, outdir, top_n, rf_ntree, seed, split_mode, train_prop, rf_mtry, allow_small_sample)
  if (length(rf_top) < 2) stop("RF + L1 combined requires at least 2 RF-selected candidate proteins before running the L1 stage.")
  l1_top <- run_l1_selection(mat, group_info, outdir, top_n, l1_alpha, seed, split_mode, train_prop, lambda_selection, cv_folds, allow_small_sample)
  top <- unique(c(rf_top, l1_top))[seq_len(min(length(unique(c(rf_top, l1_top))), top_n))]
  data.table::fwrite(data.frame(ProteinID = top), file.path(outdir, paste0("top", length(top), "_rf_l1_union_features.csv")))
  if (length(top) > 0) write_matrix_csv(mat[top, , drop = FALSE], file.path(outdir, paste0("top", length(top), "_rf_l1_union_quantity_matrix.csv")))
  top
}
run_physicochemical <- function(sets, annotation_file, outdir, width = 3.3, height = 3.3, palette = "npg") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ann <- data.table::fread(annotation_file, data.table = FALSE)
  required <- c("Accession", "GRAVY", "TM_helices", "Subcellular_class", "MW", "pI", "Length")
  miss <- setdiff(required, colnames(ann))
  if (length(miss) > 0) stop("Annotation table missing columns: ", paste(miss, collapse = ", "))
  long <- dplyr::bind_rows(lapply(names(sets), function(g) data.frame(Group = g, Accession = sets[[g]], stringsAsFactors = FALSE))) |>
    dplyr::left_join(ann, by = "Accession")
  tm_status <- if ("AnnotationStatus" %in% colnames(long)) as.character(long$AnnotationStatus) else rep(NA_character_, nrow(long))
  if ("TM_helices" %in% colnames(long)) long$TM_helices[!is.na(tm_status) & tm_status == "not_found_in_uniprot"] <- NA
  data.table::fwrite(long, file.path(outdir, "detected_proteins_annotations_long.csv"))
  props <- c("GRAVY", "MW", "pI", "Length")
  for (prop in props) {
    df <- long |> dplyr::filter(!is.na(.data[[prop]]))
    if (nrow(df) == 0) next
    p <- ggplot2::ggplot(df, ggplot2::aes(Group, .data[[prop]], fill = Group)) + ggplot2::geom_violin(trim = TRUE, alpha = 0.65) + ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.3) + ggplot2::scale_fill_manual(values = sci_palette(length(unique(df$Group)), palette)) + theme_sci() + ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(file.path(outdir, paste0(prop, "_distribution.pdf")), p, width = width, height = height)
  }
  tm_df <- long |>
    dplyr::mutate(
      TM_helices_numeric = suppressWarnings(as.numeric(TM_helices)),
      TM_valid_annotation = !is.na(TM_helices_numeric) & (is.na(tm_status) | tm_status != "not_found_in_uniprot"),
      TM_category = dplyr::case_when(
        TM_valid_annotation & TM_helices_numeric == 0 ~ "0",
        TM_valid_annotation & TM_helices_numeric == 1 ~ "1",
        TM_valid_annotation & TM_helices_numeric == 2 ~ "2",
        TM_valid_annotation & TM_helices_numeric >= 3 ~ ">=3",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(TM_valid_annotation, !is.na(TM_category))
  if (nrow(tm_df) > 0) {
    tm_levels <- c("0", "1", "2", ">=3")
    tm_tab <- tm_df |>
      dplyr::count(Group, TM_category, name = "protein_count") |>
      tidyr::complete(Group, TM_category = tm_levels, fill = list(protein_count = 0)) |>
      dplyr::group_by(Group) |>
      dplyr::mutate(
        denominator = sum(protein_count),
        proportion = dplyr::if_else(denominator > 0, protein_count / denominator, NA_real_)
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(TM_category = factor(TM_category, levels = tm_levels)) |>
      dplyr::arrange(Group, TM_category)
    data.table::fwrite(tm_tab, file.path(outdir, "TM_helices_category_proportions.csv"))
    cols <- sci_palette(length(tm_levels), palette)
    names(cols) <- tm_levels
    p <- ggplot2::ggplot(tm_tab, ggplot2::aes(Group, proportion, fill = TM_category)) +
      ggplot2::geom_col(width = 0.75, color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      theme_sci() +
      ggplot2::labs(x = NULL, y = "Proportion", fill = "TM helices")
    ggplot2::ggsave(file.path(outdir, "TM_helices_category_proportions.pdf"), p, width = width, height = height)
  }
  if (any(!is.na(long$Subcellular_class))) {
    tab <- long |>
      dplyr::filter(!is.na(Subcellular_class), Subcellular_class != "") |>
      dplyr::count(Group, Subcellular_class) |>
      dplyr::group_by(Group) |>
      dplyr::mutate(
        denominator = sum(n),
        proportion = dplyr::if_else(denominator > 0, n / denominator, NA_real_),
        label = dplyr::if_else(proportion >= 0.06, sprintf("%.1f%%\nn=%s", proportion * 100, n), NA_character_)
      ) |>
      dplyr::ungroup()
    data.table::fwrite(tab, file.path(outdir, "Subcellular_class_counts.csv"))
    # Labels below 6% are hidden because narrow stacked segments tend to overlap.
    p <- ggplot2::ggplot(tab, ggplot2::aes(Group, proportion, fill = Subcellular_class)) +
      ggplot2::geom_col(width = 0.75, color = "white", linewidth = 0.2) +
      ggplot2::geom_text(ggplot2::aes(label = label), position = ggplot2::position_stack(vjust = 0.5), size = 2.2, na.rm = TRUE) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      theme_sci() +
      ggplot2::labs(y = "Fraction")
    ggplot2::ggsave(file.path(outdir, "Subcellular_class_distribution.pdf"), p, width = width, height = height)
  }
}

run_slingshot_pseudotime <- function(mat, group_info, outdir, reduction = c("PCA", "UMAP"), start_group = NULL, end_group = NULL, width = 3.3, height = 3.3, palette = "npg", n_neighbors = 10, min_dist = 0.1, seed = 123, top_n = 50, heatmap_width = 4.17, heatmap_height = 5.56) {
  reduction <- match.arg(reduction)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  missing <- c(
    if (!requireNamespace("uwot", quietly = TRUE)) "uwot" else character(),
    if (!requireNamespace("slingshot", quietly = TRUE)) "slingshot" else character(),
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) "SingleCellExperiment" else character(),
    if (!requireNamespace("pheatmap", quietly = TRUE)) "pheatmap" else character()
  )
  if (length(missing) > 0) stop("Slingshot pseudotime requires missing package(s): ", paste(missing, collapse = ", "), ".")
  used <- preprocess_expr(mat, TRUE, 0.5)
  if (nrow(used) < 3) stop("Slingshot requires at least three proteins after preprocessing.")
  sample_mat <- t(used)
  sample_groups <- group_info$Group[match(rownames(sample_mat), group_info$Sample)]
  if (any(is.na(sample_groups))) stop("Some samples do not have group labels.")
  if (length(unique(sample_groups)) < 2) stop("Slingshot requires at least two groups/clusters.")
  if (nrow(sample_mat) < 4) stop("Slingshot requires at least four samples for a meaningful trajectory.")
  if (!is.null(start_group) && nzchar(start_group) && !start_group %in% as.character(sample_groups)) stop("Start group is not present in the current group labels: ", start_group)
  if (!is.null(end_group) && nzchar(end_group) && end_group != "None" && !end_group %in% as.character(sample_groups)) stop("End group is not present in the current group labels: ", end_group)
  pca <- prcomp(sample_mat, center = TRUE, scale. = TRUE)
  pca_coords <- pca$x[, 1:2, drop = FALSE]
  colnames(pca_coords) <- c("PC1", "PC2")
  set.seed(seed)
  nn <- min(n_neighbors, max(2, nrow(sample_mat) - 1))
  umap_coords <- uwot::umap(pca_coords, n_neighbors = nn, min_dist = min_dist, metric = "euclidean", verbose = FALSE)
  rownames(umap_coords) <- rownames(sample_mat)
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  pca_df <- data.frame(Sample = rownames(pca_coords), Group = as.character(sample_groups), PC1 = pca_coords[, 1], PC2 = pca_coords[, 2], check.names = FALSE)
  umap_df <- data.frame(Sample = rownames(umap_coords), Group = as.character(sample_groups), UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2], check.names = FALSE)
  data.table::fwrite(pca_df, file.path(outdir, "slingshot_input_PCA_coordinates.csv"))
  data.table::fwrite(umap_df, file.path(outdir, "slingshot_input_UMAP_coordinates.csv"))
  coords <- if (reduction == "PCA") pca_coords else umap_coords
  colnames(coords) <- c("Dim1", "Dim2")
  start_group <- if (is.null(start_group) || !nzchar(start_group)) as.character(sample_groups[1]) else start_group
  end_group <- if (is.null(end_group) || !nzchar(end_group) || end_group == "None") NULL else end_group
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = used),
    colData = data.frame(Sample = rownames(sample_mat), Group = factor(sample_groups, levels = levels(factor(sample_groups))), Cluster = factor(sample_groups), row.names = rownames(sample_mat))
  )
  SingleCellExperiment::reducedDims(sce)$PCA <- pca_coords
  SingleCellExperiment::reducedDims(sce)$UMAP <- umap_coords
  sce <- tryCatch(
    slingshot::slingshot(sce, clusterLabels = "Cluster", reducedDim = reduction, start.clus = start_group, end.clus = end_group),
    error = function(e) {
      stop("Slingshot could not infer a trajectory from the current group labels and ", reduction, " coordinates. Check that stages form a clear progression, or try the other reduction. Original error: ", conditionMessage(e), call. = FALSE)
    }
  )
  pst <- slingshot::slingPseudotime(sce)
  pst_df <- data.frame(Sample = colnames(sce), Group = as.character(sce$Group), Cluster = as.character(sce$Cluster), pst, check.names = FALSE)
  data.table::fwrite(pst_df, file.path(outdir, "slingshot_pseudotime_by_sample.csv"))
  data.table::fwrite(pst_df, file.path(outdir, "slingshot_sample_pseudotime.csv"))
  data.table::fwrite(data.frame(Group = levels(factor(sample_groups))), file.path(outdir, "slingshot_group_clusters.csv"))
  capture.output(slingshot::slingLineages(sce), file = file.path(outdir, "slingshot_lineages.txt"))
  cols <- sci_palette(length(levels(factor(sample_groups))), palette)
  names(cols) <- levels(factor(sample_groups))
  plot_coords <- if (reduction == "PCA") pca_coords else umap_coords
  plot_xlab <- if (reduction == "PCA") "PC1" else "UMAP1"
  plot_ylab <- if (reduction == "PCA") "PC2" else "UMAP2"
  grDevices::pdf(file.path(outdir, "Slingshot_PCA_curve.pdf"), width = width, height = height)
  graphics::plot(plot_coords[, 1], plot_coords[, 2], col = cols[as.character(sample_groups)], pch = 16, xlab = plot_xlab, ylab = plot_ylab, main = paste("Slingshot trajectory on", reduction))
  graphics::text(plot_coords[, 1], plot_coords[, 2], labels = rownames(plot_coords), pos = 3, cex = 0.55)
  graphics::legend("topright", legend = names(cols), col = cols, pch = 16, cex = 0.75, bty = "n")
  graphics::lines(slingshot::SlingshotDataSet(sce), lwd = 2, col = "black")
  grDevices::dev.off()
  file.copy(file.path(outdir, "Slingshot_PCA_curve.pdf"), file.path(outdir, "slingshot_pseudotime_trajectory.pdf"), overwrite = TRUE)
  p_stage <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, color = Group)) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(values = cols) +
    theme_sci() +
    ggplot2::labs(title = "UMAP by stage")
  ggplot2::ggsave(file.path(outdir, "UMAP_by_stage.pdf"), p_stage, width = width, height = height)
  pt_col <- colnames(pst)[1]
  umap_pt_df <- dplyr::left_join(umap_df, pst_df[, c("Sample", pt_col), drop = FALSE], by = "Sample")
  colnames(umap_pt_df)[colnames(umap_pt_df) == pt_col] <- "Pseudotime"
  p_pt <- ggplot2::ggplot(umap_pt_df, ggplot2::aes(UMAP1, UMAP2, color = Pseudotime)) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_gradient(low = "#3B4CC0", high = "#B40426", na.value = "grey80") +
    theme_sci() +
    ggplot2::labs(title = "UMAP by pseudotime")
  ggplot2::ggsave(file.path(outdir, "UMAP_by_pseudotime.pdf"), p_pt, width = width, height = height)
  pt_vec <- umap_pt_df$Pseudotime
  names(pt_vec) <- umap_pt_df$Sample
  valid_samples <- names(pt_vec)[!is.na(pt_vec)]
  if (length(valid_samples) < 3) stop("Slingshot produced too few samples with finite pseudotime for protein association analysis.")
  expr_for_cor <- used[, valid_samples, drop = FALSE]
  pt_vec <- pt_vec[valid_samples]
  cor_res <- data.frame(ProteinID = rownames(expr_for_cor), SpearmanR = NA_real_, Pvalue = NA_real_)
  for (i in seq_len(nrow(expr_for_cor))) {
    x <- expr_for_cor[i, ]
    if (length(unique(x)) > 2 && length(unique(pt_vec)) > 2) {
      ct <- suppressWarnings(stats::cor.test(x, pt_vec, method = "spearman"))
      cor_res$SpearmanR[i] <- unname(ct$estimate)
      cor_res$Pvalue[i] <- ct$p.value
    }
  }
  cor_res$FDR <- p.adjust(cor_res$Pvalue, method = "BH")
  cor_res <- cor_res[order(cor_res$FDR, -abs(cor_res$SpearmanR)), , drop = FALSE]
  data.table::fwrite(cor_res, file.path(outdir, "pseudotime_associated_proteins_spearman.csv"))
  top_proteins <- head(cor_res$ProteinID[!is.na(cor_res$FDR)], top_n)
  if (length(top_proteins) >= 2) {
    ordered_samples <- names(sort(pt_vec))
    heat_mat <- expr_for_cor[top_proteins, ordered_samples, drop = FALSE]
    heat_z <- t(scale(t(heat_mat)))
    heat_z[!is.finite(heat_z)] <- 0
    ann_col <- data.frame(Group = umap_pt_df$Group[match(colnames(heat_z), umap_pt_df$Sample)], Pseudotime = pt_vec[colnames(heat_z)])
    rownames(ann_col) <- colnames(heat_z)
    grDevices::pdf(file.path(outdir, "Top_pseudotime_proteins_heatmap.pdf"), width = heatmap_width, height = heatmap_height)
    pheatmap::pheatmap(heat_z, color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100), breaks = seq(-2, 2, length.out = 101), cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = ann_col, show_colnames = FALSE, fontsize_row = 7, border_color = NA, main = paste0("Top ", length(top_proteins), " pseudotime proteins"))
    grDevices::dev.off()
  }
  invisible(pst_df)
}
