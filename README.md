# ProteoPostZ

ProteoPostZ, formerly named ProteoDIAPostZ, is a Windows/R/Shiny desktop-style application for post-processing protein-level DIA and DDA proteomics results and user-prepared standard quantitative matrices.

Current source version: **ProteoPostZ Formal V2.0.1**.

## V2.0.1 Update Note

- Added `Quantified_Protein_Count` alongside the existing `Identified_Protein_Count`.
- Sample count tables retain both count fields, and the sample count barplot supports switching between identification count and quantified protein count.
- FragPipe/MSFragger and PEAKS now keep identification counts separate from quantitative non-missing counts by using their existing sample-level identification evidence.
- MaxQuant still lacks an independent sample-level identification evidence field in the current `proteinGroups.txt` importer, so identification count is currently approximated by non-missing `LFQ intensity`.

## Why The Name Changed

ProteoDIAPostZ V1.4 focused on DIA protein-level post-processing. V2.0 expands the input scope to DIA, DDA, and standard quantitative matrices, so the product name is shortened to **ProteoPostZ** while preserving the original version lineage.

## Supported Inputs

ProteoPostZ V2.0 uses a separated input recognition layer and converts all accepted formats into a shared internal protein quantitative matrix. The recommended and default row identifier is **accession**, because physicochemical mapping and annotation are accession-based and some FASTA/database choices do not provide a unique gene name for every accession.

DIA protein-level inputs:

- DIA-NN `pg_matrix` tables, including sample columns from `.d`, `.raw`, `.wiff`, and direct sample-name columns when suffixes are absent.
- Spectronaut report tables. `PG.Quantity` is used for quantitative analyses, and `PG.IBAQ` is used as identification evidence for qualitative summaries.

DDA protein-level inputs:

- FragPipe/MSFragger `combined_protein.tsv`; MaxLFQ intensity columns are used as label-free quantitative values.
- PEAKS Online protein result tables; top proteins are selected within each protein group, Area columns are used for quantitative values, and sample-specific Coverage columns are used as identification evidence.
- MaxQuant `proteinGroups.txt`; `LFQ intensity` columns are used as label-free quantitative values, and the first accession before a semicolon in `Protein IDs` is used by default.

Standard quantitative matrices:

- First column is the feature/protein identifier, remaining columns are sample quantitative values.
- Two explicit zero-value modes are available: `0 is a valid quantitative value` and `0 represents missing or unquantified`.
- Missing values represented as `NA`, `NaN`, blank cells, or empty strings are handled as missing.

After loading, users can still rename samples and assign groups manually. Identical group text means the same biological group.

## Main Analysis Modules

Qualitative and quantitative plots preserve the V1.x workflow with vector PDF output and per-plot CSV export. V2.0 also keeps feature UMAP/heatmap and Slingshot pseudotime support.

Dimension reduction:

- PCA
- UMAP
- t-SNE

Machine learning:

- Random forest
- L1-regularized model
- RF + L1 combined feature selection

V2.0 adds more model evaluation outputs where applicable, including cross-validation details, confusion matrices, class metrics, ROC/AUC for binary comparisons, sensitivity, specificity, and stability-selection summaries. For small sample sizes, these outputs should be interpreted as exploratory evidence rather than independent clinical or biological validation.

## Traceability Outputs

Each analysis module can write module outputs to the selected output directory. V2.0 adds run traceability files:

- `analysis_manifest.json`: input file, detected category/format, grouping, filtering, normalization, missing-value handling, method settings, app version, and output file list where available.
- `analysis_summary.html`: a lightweight browser-readable summary of the same analysis run.

These files are generated outputs and are intentionally excluded from the source repository.

## Repository And Release Asset Boundary

This Git repository contains source code, scripts, documentation, tests, and lightweight built-in annotation CSV files required by the application source.

The source repository does **not** commit:

- portable R runtime
- `portable/Rlibs`
- `outputs`
- `logs`
- ZIP archives
- compiled `.exe` launchers
- analysis-generated PDF, CSV, TSV, manifest, or HTML output files
- local temporary configuration files that may contain machine-specific absolute paths

The complete Windows portable package, including portable R/Rlibs and `ProteoPostZ_v2.0.exe`, is distributed separately as a GitHub Release asset.

## Run From Source

Install R and the required R packages, or use the separate Windows portable package for a self-contained runtime. From the repository root:

```bat
Rscript run_app.R
```

Then open:

```text
http://127.0.0.1:3840/
```

The current Windows source launcher is:

```text
ProteoPostZ_v2.0_launcher.cs
Run_ProteoPostZ_v2.0.cmd
```

The compiled launcher executable is not committed to this repository.

## Repository Layout

```text
app/                         Shiny application and R analysis code
app/R/analysis_core.R         Core input parsers, plotting, statistics, ML, and helper functions
app/annotations/              Lightweight offline annotation CSV files
tools/                       Regression and smoke-test scripts
docs/                        Development notes and historical archived launchers
README.md                    Current V2.0 source overview
RELEASE_NOTES_v2.0.0.md      V2.0 release notes
```

Historical V1.4 launcher files are archived under `docs/archive/v1.4/`; historical release notes retain the old ProteoDIAPostZ name where that is factually correct.

## V2.0 Test Scope

The V2.0 source checks cover:

- R syntax parsing for `app/app.R` and `app/R/analysis_core.R`
- DIA-NN and Spectronaut protein-level input regression
- FragPipe/MSFragger, PEAKS, and MaxQuant protein-level DDA input smoke tests
- standard matrix zero-as-value and zero-as-missing modes
- PCA, UMAP, and t-SNE smoke tests
- random forest, L1, and RF + L1 machine-learning smoke tests
- manifest and HTML summary generation-path checks
- Slingshot compatibility smoke tests when local dependencies are available
- static inspection of the Windows launcher and CMD startup scripts

## 中文简要说明

ProteoPostZ Formal V2.0.1 是 ProteoDIAPostZ V1.4 的后续版本。由于软件范围已经从 DIA 扩展到 DIA、DDA 以及用户整理的标准定量矩阵，软件名称去掉了 DIA。当前版本仍以蛋白水平结果为核心，推荐默认使用 accession 作为行标识，保留 PDF 矢量图和 CSV 导出，并新增 DDA 蛋白定量表输入、t-SNE、机器学习评估细节、`analysis_manifest.json` 和 `analysis_summary.html`。
