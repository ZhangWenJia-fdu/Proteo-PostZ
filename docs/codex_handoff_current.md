# ProteoPostZ Current Handoff

Current development version: ProteoPostZ Formal V2.0.2.

## V2.0.2 Release Scope

V2.0.2 is a source and portable-runtime update focused on PEAKS Online protein-result input and targeted quantitative-plot usability refinements. The release package is prepared separately from the source repository and excludes local test data, generated outputs, logs, and test records.

The bilingual V2.0.2 user guides identify the two user-facing additions: schema-detected PEAKS LFQ protein-result input while retaining PEAKS DB compatibility, and refined Sample correlation / Expression / Feature-protein heatmap settings. They also provide balanced qualitative and quantitative analysis overviews, a concise statement that module settings are adjustable, and PDF/CSV export guidance.

## V2.0.1 Count Definitions

V2.0.1 adds `Quantified_Protein_Count` while retaining the existing `Identified_Protein_Count` field for compatibility.

- `Quantified_Protein_Count`: per-sample count of non-missing values in the internal quantitative matrix.
- `Identified_Protein_Count`: per-sample count of non-missing values in the internal qualitative or identification-evidence matrix.
- Exported sample count data keeps both `Identified_Protein_Count` and `Quantified_Protein_Count`.

## V2.0.2 PEAKS Input Update

- One PEAKS protein file is schema-detected as either `PEAKS DB protein result` or `PEAKS LFQ protein result`; the filename is not used for subtype detection. When DB and LFQ-style fields are mixed, DB schema takes precedence.
- DB uses `Area <sample>` for quantitative values and sample-specific `Coverage(%) <sample> > 0` for identification.
- LFQ uses `<sample> Area` for quantitative values and excludes `Group N Area`, ratio/profile, significance, and other non-sample fields. LFQ sample-level identification evidence is unavailable, so `Identified_Protein_Count` is `NA` while `Quantified_Protein_Count` counts non-missing LFQ Area values.
- Both subtypes retain one representative row per `Protein Group`: the first `Top == TRUE` row in source order, or the first group row when no Top row exists. No within-group abundance aggregation is performed.
- The detected subtype is shown in the input preview. LFQ identification counts remain unavailable (`NA`) and the UI explains that sample-level identification evidence is absent; choosing that metric for the barplot returns a clear unavailable message.

## Input-Specific Semantics

- FragPipe/MSFragger: identification uses sample `Spectral Count > 0`; quantification uses non-missing `MaxLFQ Intensity` after zero-to-`NA` handling.
- PEAKS: identification uses sample-specific `Coverage(%) sample` values greater than zero; quantification uses non-missing sample `Area` after zero-to-`NA` handling. Protein Group reduction is unchanged: first `Top = TRUE` row wins, and if multiple rows are top entries the earliest source-table row is retained.
- MaxQuant: the protein-level importer has no independent sample-level identification evidence. `Identified_Protein_Count` is therefore approximated by non-missing `LFQ intensity`, and currently matches `Quantified_Protein_Count`.
- DIA-NN: counts are computed from the existing DIA-NN qualitative and quantitative matrices. In the current importer both matrices use the same sample quantity evidence.
- Spectronaut: identification counts use the existing `PG.IBAQ` evidence matrix, while quantified counts use non-missing `PG.Quantity`.
- Standard quantitative matrix: counts represent available quantitative values under the selected zero-handling mode. The app does not fabricate an independent strict identification count for standard matrices.

## UI And Export State

- The input preview count table displays both count fields.
- The input Sample count preview retains the pre-existing DT table behavior, including its existing `Show entries` and search controls.
- The sample protein count barplot module supports choosing identification count or quantified protein count.
- Standard matrix barplot selection exposes quantified count only.
- MaxQuant displays the approximation note in English where relevant.
- Count-related CSV exports retain both count fields; generated plot titles and axis labels follow the selected count type.
- Venn and UpSet each show their group-level set-definition note directly above the minimum-replicate control.
- Quantitative plots are arranged across five fixed-viewport pages. Analysis cards retain full-screen expansion and independently scrollable controls/results.
- Sample correlation heatmap supports `Original` order (default), `By group`, and global `Hierarchical clustering`. Hierarchical clustering supports `1 - correlation` or Euclidean distance of correlation profiles, `complete` or `average` linkage, and a cut-tree `k` with visible cluster gaps. The selected Pearson or Spearman method is the correlation basis for clustering.
- Expression and Feature-protein heatmaps independently support row and column `Hierarchical`, `K-means`, and `None` choices. Row hierarchical/K-means `k` defaults to 1; column hierarchical/K-means `k` defaults to the configured Group count. K-means and cut-tree hierarchical clusters have visible gaps.
- Both heatmaps support Blue-White-Red, Red-White-Blue, Purple-White-Orange, and Green-Black-Red expression color schemes; independent Group and sample-cluster annotation palettes; and custom per-label `#RRGGBB` colors. Sample-cluster annotation applies to column clustering.
- Both heatmaps support optional row labels from Accession, Protein name, or Gene name, with `None` as the default. Label font size adapts to PDF height and displayed row count. Default heatmap PDF size is 600 x 800 pt.

## Validation Completed

- Current source version and visible app title are `ProteoPostZ Formal V2.0.2`.
- PEAKS DB, PEAKS LFQ, and mixed DB/LFQ schema precedence checks passed using local development inputs supplied outside the repository.
- PEAKS LFQ preview/export/barplot handling passed with unavailable identification counts and valid quantified counts.
- The unpacked package at `F:\ProteoPostZ_v2.0.2_windows_x86_release` contains 34,775 files and was verified at approximately 1.24 GB before user-side compression.
- The refreshed unpacked package was independently cold-started from its own portable R runtime and passed HTTP 200 / Formal V2.0.2 page verification. Its bilingual user guides were synchronized with the current source; the final package audit found no generated outputs, logs, test records, or Git metadata.
- Package launch from a stopped state returned HTTP 200 at `http://127.0.0.1:3840/` and exposed Formal V2.0.2.
- The package audit found no application `logs`, `outputs`, `Rplots.pdf`, analysis manifests/summaries, test directories, or local test records.

- From the repository root, `Rscript -e "parse(file = 'app/app.R'); parse(file = 'app/R/analysis_core.R')"` passed.
- From the repository root, `Rscript tools/test_v201_count_modes.R` passed with `V201_COUNT_MODES_OK`.
- From the repository root, `Rscript tools/test_v14_input_regression.R` passed with `V14_INPUT_REGRESSION_OK`.
- DIA-NN and Spectronaut regression inputs are supplied through `PROTEOPOSTZ_DIANN_TEST_FILE` and `PROTEOPOSTZ_SPECTRONAUT_TEST_FILE`; local file paths are not stored in the repository.
- Expression heatmap row/column combinations: 9/9 passed, with PDF and CSV outputs generated for every combination.
- Feature-protein heatmap row/column combinations: 9/9 passed, with PDF and CSV outputs generated for every combination.
- Current V2.0.2 correlation and heatmap refinements passed on locally supplied real importer inputs and a temporary multi-group simulated matrix. Coverage included all supported importer families, PEAKS DB/LFQ detection, all correlation order modes, both clustering distances/linkages, heatmap clustering modes, color schemes, custom annotations, and row-label modes. No local input paths, source data, fixtures, or derived data were added to the repository.
- Correlation long-format export now safely handles a legitimate sample label that matches its former internal `Sample1` temporary-column name.
- Source app startup on `http://127.0.0.1:3840/`: HTTP 200 verified.

Known non-blocking local warnings during tests: Shiny package build-version warning, systemfonts/textshaping Freetype warning, and an existing bslib `layout_columns` breakpoint warning. These are environment or existing UI warnings and were not introduced by the V2.0.2 PEAKS input logic.

## Publication Boundary

- The V2.0.2 source commit is recorded locally only unless a later push is explicitly requested.
- No tag or GitHub Release is created by this handoff.
- The unpacked Windows package is prepared at `F:\ProteoPostZ_v2.0.2_windows_x86_release`; the user will compress or publish it separately.
