# ProteoPostZ Formal V2.0 Release Notes

ProteoPostZ Formal V2.0 is the next major source release after ProteoDIAPostZ Formal V1.4. The name removes `DIA` because the application now supports protein-level DIA outputs, protein-level DDA quantitative tables, and user-prepared standard quantitative matrices.

## Input Expansion

V2.0 separates format recognition from downstream analysis and converts supported inputs into one shared internal protein quantitative matrix.

DIA inputs:

- DIA-NN `pg_matrix` tables remain supported.
- DIA-NN sample column recognition now supports `.d`, `.raw`, `.wiff`, and direct sample-name columns when vendor suffixes are absent.
- Spectronaut report tables remain supported; `PG.Quantity` is used for quantitative analysis and `PG.IBAQ` is used as identification evidence for qualitative summaries.

DDA inputs:

- FragPipe/MSFragger `combined_protein.tsv`: MaxLFQ intensity columns are imported as label-free quantitative values. Protein groups are normalized to a single representative accession for the shared matrix.
- PEAKS Online protein result tables: protein groups are deduplicated by group number, preferring the first `Top = true` row, with table order used to break multiple-top cases. Area columns are imported as quantitative values; zero Area is treated as missing or unquantified. Sample-specific Coverage columns are used as identification evidence, while the unsuffixed global Coverage column is not treated as a sample.
- MaxQuant `proteinGroups.txt`: `LFQ intensity` columns are imported as label-free quantitative values. Zero LFQ values are treated as missing or unquantified. The first accession before a semicolon in `Protein IDs` is used as the default accession-level identifier.

DDA import in V2.0 is limited to protein-level quantitative result tables. Peptide-level and modification-site-level DDA workflows are not enabled in this source release.

Standard quantitative matrices:

- Users can import prepared CSV or TSV matrices with one identifier column followed by sample columns.
- Two explicit zero modes are available: `0 is a valid quantitative value` and `0 represents missing or unquantified`.
- Explicit missing markers such as `NA`, `NaN`, blank cells, and empty strings are handled as missing.

## Internal Data Model

All accepted input formats are converted into a unified protein-level matrix model with sample names, quantitative values, qualitative/identification evidence, protein metadata, and format-detection notes. Accession is the recommended and default row identifier because physicochemical annotation and mapping are accession-based, and not every accession has a unique gene-name mapping in all FASTA/database designs.

## Dimension Reduction

V2.0 adds t-SNE alongside the existing PCA and UMAP workflow. PCA, UMAP, and t-SNE outputs use the same loaded matrix and user grouping definitions.

## Machine Learning

The random forest, L1, and RF + L1 modes are retained. V2.0 expands the evaluation outputs where applicable:

- cross-validation details
- confusion matrix
- per-class metrics
- sensitivity and specificity
- ROC/AUC for binary comparisons
- stability-selection summaries
- selected-feature quantity matrices and downstream feature UMAP/heatmap compatibility

For small sample sizes, machine-learning results should be interpreted as exploratory feature-ranking evidence rather than validated predictive performance.

## Traceability And Reporting

V2.0 adds run traceability outputs:

- `analysis_manifest.json`
- `analysis_summary.html`

These files record available input file metadata, detected format, grouping, filtering, normalization, missing-value handling, method settings, application version, and output files. Existing vector PDF plots and CSV exports are retained.

## Interface Updates

The input panel was redesigned around three input categories: DIA result, DDA result, and Standard quantitative matrix. The software-format selector changes with the selected category. The input area uses a responsive two-column layout where screen width permits it and collapses cleanly on narrower layouts.

The current visible application title is **ProteoPostZ Formal V2.0**.

## Offline Runtime Scope

Normal application use remains offline. Built-in annotation CSV files are included in the source tree. Optional annotation refresh scripts are kept under `tools/` and are not part of normal offline analysis.

## Source And Portable Package Boundary

The Git source repository contains R/Shiny source code, launcher source, startup scripts, tests, documentation, and lightweight annotation CSV files.

The repository does not commit portable R, `portable/Rlibs`, outputs, logs, ZIP archives, compiled executables, or analysis-generated PDF/CSV/TSV/manifest/HTML files. The full Windows portable package is distributed separately through a GitHub Release asset.