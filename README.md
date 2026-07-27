# ProteoDIAPostZ

ProteoDIAPostZ is a Windows/R/Shiny desktop-style application for post-processing protein-level DIA proteomics results from DIA-NN, Spectronaut, and preprocessed standard quantitative matrices.

Current source version: **Formal V1.4**

Developed by Wenjia Zhang, Department of Chemistry, Fudan University, Modern Chromatography Separation and Analysis.

## Source Repository Boundary

This repository contains the application source code, tests, documentation, launcher source, and lightweight startup scripts.

The source repository intentionally does **not** include:

- portable R runtime folders such as `portable/R-4.5.1/`;
- packaged R library folders such as `portable/Rlibs/`;
- generated analysis outputs under `outputs/`;
- runtime logs under `logs/`;
- built zip archives or installer/package artifacts;
- compiled launcher executables.

The complete portable Windows software package is distributed separately as a **GitHub Release asset**. Users who only want to run the software should download the release asset rather than cloning this source repository.

## Main Features

Formal V1.4 keeps the V1.3 protein-level post-processing workflow and adds standard quantitative matrix input.

Core V1.3 features retained in V1.4:

- DIA-NN and Spectronaut protein-level result post-processing.
- Protein identification summaries, Venn/UpSet plots, physicochemical property plots, correlation heatmaps, rank-abundance plots, CV ridgelines, PCA, UMAP, volcano plots, and expression heatmaps.
- Machine-learning analyses: Random forest, L1, and RF + L1 combined.
- Feature UMAP/heatmap and Slingshot pseudotime analysis.
- Per-module controls for plot generation, dimensions, palettes, CSV export, previews, and output paths.

New in V1.4:

- Standard quantitative matrix input through the internal source type `standard_matrix`.
- CSV and TSV matrix reading with arbitrary first-column and sample-column names.
- Explicit zero-handling modes:
  - `zero_is_value`: zero is retained as a valid quantitative value;
  - `zero_is_missing`: zero is converted to missing only in the analysis matrix.
- Raw standard-matrix values are retained separately from the downstream `analysis_quant_matrix`.
- Standard matrix summaries use "available quantitative values" terminology rather than protein-identification terminology.
- DIA-NN/Spectronaut physicochemical accession mapping is more stable when Protein name or Gene name row identifiers are duplicated and internally made unique.
- The Input page uses a shorter `ProteoDIAPostZ Formal V1.4` title and a fixed developer identity area at the bottom of the left sidebar.

## Input Types

### DIA-NN

DIA-NN input behavior is preserved from V1.3. The app uses the DIA-NN quantity matrix for identification and quantitative downstream analyses.

### Spectronaut

Spectronaut input behavior is preserved from V1.3. `PG.IBAQ` is used for identification/qualitative summaries, and `PG.Quantity` is used for quantitative analyses.

### Standard Quantitative Matrix

A standard quantitative matrix is a CSV or TSV feature-by-sample matrix:

- first row: column names;
- first column: protein, gene, accession, or other feature identifiers;
- second column onward: sample quantitative values.

The quantitative region accepts numeric values and explicit missing-value entries: blank cells, whitespace-only cells, `NA`, and `NaN`. `Inf`, `-Inf`, and other nonnumeric text are rejected.

Users must explicitly choose how numeric zero is interpreted. The app does not infer this from the data distribution or zero proportion.

## Repository Layout

- `app/app.R`: Shiny UI and server wiring.
- `app/R/analysis_core.R`: input parsing and analysis functions.
- `app/annotations/`: tracked annotation tables used by built-in physicochemical annotation.
- `tools/`: regression and smoke-test scripts.
- `docs/`: handoff and V1.4 development notes.
- `ProteoDIAPostZ_v1.4_launcher.cs`: V1.4 launcher source for building the Windows launcher executable.
- `Run_ProteoDIAPostZ_v1.4.cmd`: fallback command launcher for a built portable package.

## Development Checks

Representative V1.4 checks used during this update:

```powershell
Rscript tools/test_standard_matrix_parser.R
Rscript tools/test_standard_matrix_shiny_integration.R
Rscript tools/test_standard_matrix_downstream_smoke.R
Rscript tools/test_v14_input_regression.R
```

The full portable package was also validated from `F:\ProteoDIAPostZ_v1.4_windows_x86_release` using its own bundled portable R and R libraries, with no runtime dependency on the E-drive development repository or V1.3 release folder.

## Running From Source

A source checkout requires a suitable local R installation and required R packages. For normal end users, use the complete portable Windows release asset instead.

From a configured development environment:

```powershell
Rscript run_app.R
```

Then open:

```text
http://127.0.0.1:3840/
```