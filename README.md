# ProteoPostZ

**ProteoPostZ** is a Windows/R/Shiny application for post-processing protein-level proteomics results from **DIA**, **DDA**, and user-prepared standard quantitative matrices.

ProteoPostZ was formerly named **ProteoDIAPostZ**. The name was changed in V2.0 when the supported input scope was expanded beyond DIA workflows.

**Latest release: ProteoPostZ Formal V2.0.1**

[Download the latest Windows release](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/latest)

> For normal use, download the complete Windows portable package from the **Assets** section of the GitHub Release.
> Do not use GitHub's automatically generated **Source code (zip)** or **Source code (tar.gz)** archives as the runnable application package.

---

## Supported Inputs

ProteoPostZ converts supported protein-level results into a shared internal quantitative matrix for downstream analysis.

### DIA protein-level inputs

* **DIA-NN** `pg_matrix` tables
* **Spectronaut** protein-level report tables

For Spectronaut input, `PG.Quantity` is used for quantitative analyses, while `PG.IBAQ` can provide identification evidence for qualitative summaries.

### DDA protein-level inputs

* **FragPipe/MSFragger** `combined_protein.tsv`

  * MaxLFQ Intensity is used for label-free quantification.
* **PEAKS Online** protein result tables

  * Area is used for quantitative analysis.
  * Sample-specific Coverage is used as identification evidence.
* **MaxQuant** `proteinGroups.txt`

  * LFQ intensity is used for label-free quantification.

### Standard quantitative matrices

User-prepared `.csv` or `.tsv` quantitative matrices are also supported.

* The first column contains the feature/protein identifier.
* Remaining columns contain sample quantitative values.
* Two explicit zero-value modes are available:

  * `0 is a valid quantitative value`
  * `0 represents missing or unquantified`
* `NA`, `NaN`, blank cells, and empty strings are treated as missing values.

After data loading, users can rename samples and assign biological groups manually.

The recommended and default protein identifier is **accession**, because physicochemical-property mapping and built-in annotations are accession-based.

---

## Main Analysis Features

### Protein identification and quantification summaries

* Identified protein counts
* Quantified protein counts
* Sample protein-count barplots
* Venn diagrams
* UpSet plots

### Protein physicochemical and quantitative characteristics

* Physicochemical-property analysis
* Rank-abundance plots
* Correlation analysis
* CV distribution/ridgeline plots

### Dimension reduction

* PCA
* UMAP
* t-SNE

### Differential and expression analysis

* Volcano plots
* Expression heatmaps

### Machine learning and feature selection

* Random forest
* L1-regularized feature selection
* RF + L1 combined feature selection
* Feature-protein UMAP
* Feature-protein heatmaps

Where applicable, model evaluation outputs include:

* cross-validation details
* confusion matrices
* class metrics
* ROC/AUC for binary comparisons
* sensitivity
* specificity
* stability-selection summaries

For small sample sizes, machine-learning results should be interpreted as exploratory evidence rather than independent biological or clinical validation.

### Pseudotime analysis

* Slingshot pseudotime analysis

### Export

Analysis modules support:

* vector PDF figure export
* CSV result export
* run-level analysis traceability outputs

---

## What's New in V2.0.1

ProteoPostZ Formal V2.0.1 is a maintenance update to V2.0.0. The major supported input formats and downstream analysis framework remain unchanged.

Main changes include:

* Added `Quantified_Protein_Count` while retaining `Identified_Protein_Count`.
* Sample protein-count plots can switch between identified and quantified protein counts.
* FragPipe/MSFragger and PEAKS now distinguish sample-level identification evidence from available quantitative values.
* Expression and Feature-protein heatmaps now support independent row and column clustering.
* Row and column clustering can independently use `Hierarchical`, `K-means`, or `None`.
* Separate K-means `k` settings are available for rows and columns where applicable.

---

## Identification and Quantification Count Definitions

Where independent sample-level identification evidence is available, ProteoPostZ distinguishes **identified proteins** from **quantified proteins**.

### FragPipe/MSFragger

For `combined_protein.tsv`:

* **Identified protein count:** sample-specific `Spectral Count > 0`
* **Quantified protein count:** non-missing `MaxLFQ Intensity`

Zero MaxLFQ values are treated as missing quantitative values.

### PEAKS

For PEAKS protein-level input:

* **Identified protein count:** sample-specific `Coverage > 0`
* **Quantified protein count:** non-missing `Area`

Zero Area values are treated as missing quantitative values.

### MaxQuant

The current `proteinGroups.txt` importer does not contain an independent sample-level identification-evidence field.

Therefore:

* quantitative values are based on `LFQ intensity`;
* identification count is currently approximated using non-missing `LFQ intensity`.

Accordingly, identified and quantified counts currently coincide for MaxQuant input.

### Standard quantitative matrices

A user-prepared standard quantitative matrix does not necessarily contain independent protein-identification evidence.

ProteoPostZ therefore reports the number of available quantitative values according to the selected zero-handling mode rather than creating an artificial strict identification count.

---

## Heatmap Clustering

Expression heatmaps and Feature-protein heatmaps allow row and column clustering to be configured independently.

Available choices are:

* `Hierarchical`
* `K-means`
* `None`

When `K-means` is selected:

* the row `k` value can be configured independently;
* the column `k` value can be configured independently;
* the corresponding `k` control is displayed only when K-means is selected for that direction.

Selecting `None` preserves the corresponding input order.

---

## Windows Portable Release

The recommended way to run ProteoPostZ on Windows is to download the complete portable package from:

[GitHub Releases](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases)

For Formal V2.0.1, download:

`ProteoPostZ_v2.0.1_windows_x86_release.zip`

Do **not** download GitHub's automatically generated **Source code (zip)** or **Source code (tar.gz)** archives if you want the runnable Windows application. These archives contain the source repository only and do not include the complete portable R runtime and packaged R libraries.

### Running the portable version

After downloading:

1. Extract the ZIP archive completely.
2. Double-click `ProteoPostZ_v2.0.1.exe`.
3. If the executable launcher cannot be used, run `Run_ProteoPostZ_v2.0.1.cmd`.
4. ProteoPostZ will open locally at:

`http://127.0.0.1:3840/`

The browser interface runs only on the current computer and is not an online web service.

Normal use of the complete portable package does not require an internet connection.


---

## Traceability Outputs

ProteoPostZ can generate run-level traceability files in the selected output directory:

* `analysis_manifest.json`
* `analysis_summary.html`

`analysis_manifest.json` can record information including:

* input file
* detected input category and format
* sample grouping
* filtering
* normalization
* missing-value handling
* analysis settings
* application version
* generated output files where available

`analysis_summary.html` provides a lightweight browser-readable summary of the analysis run.

These files are generated analysis outputs and are intentionally excluded from the source repository.

---

## Run From Source

Users who wish to run ProteoPostZ directly from source can install R and the required R packages.

From the repository root:

```bat
Rscript run_app.R
```

Then open:

```text
http://127.0.0.1:3840/
```

The current Windows source launcher files are:

```text
ProteoPostZ_v2.0.1_launcher.cs
Run_ProteoPostZ_v2.0.1.cmd
```

The compiled `.exe` launcher and complete portable R environment are distributed with the GitHub Release package rather than committed to the source repository.

---

## Repository Layout

```text
app/                         Shiny application and R analysis code
app/R/analysis_core.R        Core input parsing, plotting, statistics, ML, and helper functions
app/annotations/             Built-in lightweight annotation tables
tools/                       Regression and smoke-test scripts
docs/                        Development documentation and archived files
README.md                    Current software overview
RELEASE_NOTES_v2.0.0.md      V2.0.0 release notes
RELEASE_NOTES_v2.0.1.md      V2.0.1 release notes
```

Historical documentation retains the **ProteoDIAPostZ** name where that name was used in the corresponding software version.

---

## 中文简介

**ProteoPostZ Formal V2.0.1** 是一个基于 R/Shiny 的蛋白质组学蛋白水平结果后处理软件。

软件原名 **ProteoDIAPostZ**。自 V2.0 起，输入范围由 DIA 扩展至 **DIA、DDA 以及用户整理的标准定量矩阵**，因此软件正式更名为 **ProteoPostZ**。

目前支持的主要输入包括：

* DIA-NN
* Spectronaut
* FragPipe/MSFragger
* PEAKS
* MaxQuant
* 用户整理的 CSV/TSV 标准定量矩阵

主要分析功能包括：

* 蛋白鉴定和定量数量统计
* Venn 和 UpSet 分析
* 蛋白理化性质分析
* 相关性分析
* PCA、UMAP 和 t-SNE
* 差异分析和表达 heatmap
* Random forest
* L1 特征筛选
* RF + L1 联合特征筛选
* Feature-protein UMAP 和 heatmap
* Slingshot 拟时序分析

V2.0.1 进一步区分了样本水平的**鉴定蛋白数**与**定量蛋白数**，并完善了 Expression heatmap 和 Feature-protein heatmap 的行、列独立聚类设置。
