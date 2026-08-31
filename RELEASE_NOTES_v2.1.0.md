# ProteoPostZ Formal V2.1.0 Release Notes

ProteoPostZ Formal V2.1.0 focuses on quantitative quality control, internal code organization and reproducibility, targeted machine-learning refinements, and practical improvements across several downstream analysis modules.

## Highlights

### Quantitative QC

V2.1.0 introduces a dedicated **Quantitative QC** module before the main quantitative-analysis workflow.

The module summarizes the current imported protein × sample quantitative matrix without applying an additional normalization or imputation layer. It respects the missing-value and zero-value semantics defined by the selected importer or standard-matrix mode.

Outputs include:

- sample missingness
- protein missingness distribution
- missingness heatmap
- missingness versus abundance
- sample abundance distribution
- sample-level QC summary CSV
- protein-level QC summary CSV

The QC workflow is descriptive and does not automatically remove samples or proteins.

### Internal modular refactor and reproducibility

The analysis code has been reorganized into dedicated modules to improve maintainability, regression testing, and future feature development.

V2.1.0 also strengthens reproducibility and traceability through source-level dependency records and additional runtime information in analysis manifests. The Windows portable release continues to bundle its own R runtime and required R libraries for ordinary users.

### Machine-learning workflow refinements

Random forest, L1, and RF + L1 retain the established stability-selection framework and user-adjustable parameter structure, with targeted refinements to reproducibility and model logic.

These include improved CV-fold handling, binary/multiclass L1 behavior, RF stability-score handling, and current-session ML feature-source handling for downstream Feature UMAP and heatmap.

Detailed parameter definitions, defaults, recommendations, and a concise V2.1.0-versus-V2.0.3 algorithmic note are provided in `README_ML_parameters_v2.1.0_Bilingual.txt`.

### Analysis-module usability and robustness

Several analysis modules received targeted improvements to parameter handling, visualization behavior, current-session result handling, preview/output mapping, and edge-case robustness. These changes improve consistency and practical usability without changing the overall analysis workflow.

Rank-abundance plots additionally support both `log2(abundance + 1)` and `log10(abundance + 1)` display transforms. Missing values remain excluded from ranking, and legitimate zero values remain valid when the selected input semantics define zero as quantitative.

## Supported inputs

V2.1.0 continues to support:

- DIA-NN
- Spectronaut
- FragPipe/MSFragger
- PEAKS DB/LFQ
- MaxQuant
- user-prepared standard quantitative CSV/TSV matrices

## Validation and reproducibility

V2.1.0 includes a reusable synthetic regression suite covering input parsing, zero/missing-value semantics, quantitative QC, major visualization and dimension-reduction paths, differential analysis, heatmaps, annotation, machine learning, ML downstream feature handling, and Slingshot safeguards.

The release candidate passed the standard regression suite and fresh application startup checks before release preparation.

## Compatibility

ProteoPostZ V2.1.0 remains a local Windows/R Shiny application. The portable Windows release includes the runtime and required R libraries and is intended to run without a separate system R installation.

For detailed usage, see the bundled English and Chinese guides. For machine-learning settings, see `README_ML_parameters_v2.1.0_Bilingual.txt`.
