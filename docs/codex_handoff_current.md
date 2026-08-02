# ProteoPostZ Current Handoff

Current development version: ProteoPostZ Formal V2.0.1.

## V2.0.1 Count Definitions

V2.0.1 adds `Quantified_Protein_Count` while retaining the existing `Identified_Protein_Count` field for compatibility.

- `Quantified_Protein_Count`: per-sample count of non-missing values in the internal quantitative matrix.
- `Identified_Protein_Count`: per-sample count of non-missing values in the internal qualitative or identification-evidence matrix.
- Exported sample count data keeps both `Identified_Protein_Count` and `Quantified_Protein_Count`.

## Input-Specific Semantics

- FragPipe/MSFragger: identification uses sample `Spectral Count > 0`; quantification uses non-missing `MaxLFQ Intensity` after zero-to-`NA` handling.
- PEAKS: identification uses sample-specific `Coverage(%) sample` values greater than zero; quantification uses non-missing sample `Area` after zero-to-`NA` handling. Protein Group reduction is unchanged: first `Top = TRUE` row wins, and if multiple rows are top entries the earliest source-table row is retained.
- MaxQuant: the current `proteinGroups.txt` importer has no independent sample-level identification evidence. `Identified_Protein_Count` is therefore approximated by non-missing `LFQ intensity`, and currently matches `Quantified_Protein_Count`.
- DIA-NN: counts are computed from the existing DIA-NN qualitative and quantitative matrices. In the current importer both matrices use the same sample quantity evidence.
- Spectronaut: identification counts use the existing `PG.IBAQ` evidence matrix, while quantified counts use non-missing `PG.Quantity`.
- Standard quantitative matrix: counts represent available quantitative values under the selected zero-handling mode. The app does not fabricate an independent strict identification count for standard matrices.

## UI And Export State

- The input preview count table displays both count fields.
- The sample protein count barplot module supports choosing identification count or quantified protein count.
- Standard matrix barplot selection exposes quantified count only.
- MaxQuant displays the approximation note in English where relevant.
- Count-related CSV exports retain both count fields; generated plot titles and axis labels follow the selected count type.

## Validation Completed

- `Rscript -e "parse('E:/Project/ProteomicsDIAApp_Project/repo/app/app.R')"`: passed with `APP_PARSE_OK`.
- `Rscript -e "parse('E:/Project/ProteomicsDIAApp_Project/repo/app/R/analysis_core.R')"`: passed with `CORE_PARSE_OK`.
- `Rscript E:\Project\ProteomicsDIAApp_Project\repo\tools\test_v201_count_modes.R`: passed with `V201_COUNT_MODES_OK`.
- `Rscript E:\Project\ProteomicsDIAApp_Project\repo\tools\test_v14_input_regression.R`: passed with `V14_INPUT_REGRESSION_OK`.
- Source app startup on `http://127.0.0.1:3840/`: HTTP 200 verified.

Known non-blocking local warnings during tests: Shiny package build-version warning, systemfonts/textshaping Freetype warning, and an existing bslib `layout_columns` breakpoint warning. These are environment or existing UI warnings and were not introduced by the V2.0.1 count logic.

## Not Completed In This Source Commit

- No GitHub push.
- No tag.
- No GitHub Release.
- No formal release ZIP.
- No SHA256 checksum for a formal release asset.
