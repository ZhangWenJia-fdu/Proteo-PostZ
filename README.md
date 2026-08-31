# ProteoPostZ

ProteoPostZ, formerly ProteoDIAPostZ, is a Windows/R Shiny application for post-processing protein-level DIA and DDA proteomics results, as well as user-prepared standard quantitative matrices.

**Latest release: ProteoPostZ Formal V2.1.0.**

## Download the Windows portable application

Download `ProteoPostZ_v2.1.0_windows_x86_release.zip` from the V2.1.0 Release Assets after the release is published.

Expected portable-package name: `ProteoPostZ_v2.1.0_windows_x86_release.zip`.

The GitHub-generated **Source code (zip)** and **Source code (tar.gz)** downloads contain source only; they are not the complete Windows portable application. Windows users should download the portable ZIP from Release Assets.

## What's New in V2.1.0

ProteoPostZ Formal V2.1.0 focuses on quantitative quality control, code maintainability and reproducibility, targeted refinements to the machine-learning workflow, and practical improvements across several downstream analysis modules.

- **Quantitative QC module.** A new quantitative quality-control workflow summarizes sample- and protein-level missingness and abundance before downstream quantitative analysis. It produces five QC visualizations and two summary tables while respecting the missing-value and zero-value semantics established by the selected importer.
- **Internal modular refactor and reproducibility.** Core analysis logic has been reorganized into dedicated R modules to improve maintainability, testing, and future development. Reproducibility metadata and dependency records have also been strengthened.
- **Machine-learning workflow refinements.** Random forest, L1, and RF + L1 stability-selection workflows received targeted methodological and reproducibility refinements while retaining the established user-facing parameter framework. Detailed parameter guidance is provided separately.
- **Analysis-module usability and robustness.** Several qualitative, quantitative, annotation, trajectory, visualization, and ML-downstream modules received practical improvements to parameter handling, current-session result handling, plotting behavior, and edge-case robustness.
- **Rank-abundance transform option.** Rank-abundance plots now support both `log2(abundance + 1)` and `log10(abundance + 1)` display transforms. Missing values remain excluded, and legitimate zero values are retained when the selected input semantics define zero as quantitative.

See [README_ML_parameters_v2.1.0_Bilingual.txt](README_ML_parameters_v2.1.0_Bilingual.txt) for the detailed bilingual machine-learning parameter guide and [RELEASE_NOTES_v2.1.0.md](RELEASE_NOTES_v2.1.0.md) for the concise change record.

## Supported inputs

ProteoPostZ converts accepted files into a shared internal protein quantitative matrix. The recommended default row identifier is **accession**, because offline physicochemical annotation is accession-based and gene names are not universally unique.

### DIA protein-level results

- **DIA-NN** `pg_matrix` tables. The importer recognizes common instrument-file suffixes and direct sample-name columns.
- **Spectronaut** report tables. `PG.Quantity` supplies quantitative values and `PG.IBAQ` supplies identification evidence for qualitative summaries.

### DDA protein-level results

- **FragPipe/MSFragger** protein-level results: `MaxLFQ intensity` columns are used as label-free quantitative values.
- **PEAKS** protein-level results: `Area` is used for quantitative analysis. Sample-specific `Coverage` is used as identification evidence.
- **MaxQuant** protein-level results: `LFQ intensity` columns are used as label-free quantitative values; the first accession in `Protein IDs` is used by default.

### Note: PEAKS DB/LFQ

- **PEAKS** protein-level results remain one upload entry: **one PEAKS protein file is uploaded at a time**. The subtype is determined from the internal column schema, never from the filename.
- **PEAKS DB:** uses `Area <sample>` as sample quantitative values; sample-specific `Coverage(%) <sample> > 0` supplies independent sample-level identification evidence.
- **PEAKS LFQ:** uses `<sample> Area` as the reported protein-level LFQ abundance; `Group N Area`, ratio/profile columns, significance, and other non-sample fields are excluded. PEAKS LFQ does not provide independent sample-specific identification evidence, so `Identified_Protein_Count` is unavailable while `Quantified_Protein_Count` is calculated from non-missing LFQ Area values.

For either PEAKS subtype, one representative row is retained per `Protein Group`: the first source-order row with `Top == TRUE`; if no such row exists, the first row in the group. Areas from different members of a group are never summed or averaged. If a file contains both DB- and LFQ-style fields, the DB schema takes precedence.

### Standard quantitative matrix

- The first column is the protein/feature identifier and remaining columns are sample quantitative values.
- Users explicitly choose whether zero is a valid quantitative value or denotes missing/unquantified data.
- `NA`, `NaN`, blank, and empty values are treated as missing.

After loading, sample names and biological groups can be edited in the interface; identical group text denotes the same group.

## Main analysis features

### Qualitative analysis

- **Sample protein counts:** identification and quantified protein counts are kept separate whenever the importer provides independent identification evidence.
- **Venn diagram and UpSet plot:** both use group-level protein sets rather than sample-level sets. The *Minimum replicates detected in group* setting determines how many detected replicates are required for a protein to enter each group set.
- **Physicochemical-property analysis:** compares group protein sets using accession-based offline annotations, including GRAVY, molecular weight, pI, length, transmembrane-helices, and subcellular-class summaries.

### Quantitative QC

The Quantitative QC module operates directly on the current imported protein × sample quantitative matrix and is intended to describe the data before downstream quantitative analysis.

It provides:

- sample missingness
- protein missingness distribution
- missingness heatmap
- missingness versus abundance
- sample abundance distribution

Two summary tables are exported:

- `quantitative_qc_sample_summary.csv`
- `quantitative_qc_protein_summary.csv`

Quantitative QC respects importer-defined missing-value semantics. It does not reinterpret zero, perform normalization, imputation, scaling replacement, batch correction, or automatically filter samples/proteins based on QC results.

### Quantitative overview and dimension reduction

- **Sample correlation heatmap:** Pearson or Spearman correlation with configurable ordering, clustering, annotation, color range, precision, and PDF dimensions.
- **Rank-abundance plot:** supports `log2(abundance + 1)` and `log10(abundance + 1)`. Missing values are excluded from ranking; legitimate zero values remain zero when zero is defined as quantitative.
- **Within-group CV ridgeline:** summarizes protein-level coefficient-of-variation distributions within assigned groups.
- **PCA:** produces scaled PC1/PC2 sample coordinates after valid-feature filtering. Constant/zero-variance features are excluded from PCA scaling when necessary while leaving the underlying quantitative matrix unchanged.
- **UMAP:** sample-level UMAP with configurable `n_neighbors` and `min_dist` and a reproducible seed.
- **t-SNE:** two-dimensional sample embeddings with automatic/manual perplexity, configurable iterations, valid-feature filtering, and reproducible seeding.

### Differential analysis

**Volcano plot** compares a selected reference group with a selected comparison group and supports limma moderated testing or two-sample t tests, selectable fold-change calculation, raw/BH-adjusted p-value thresholds, and configurable significance cutoffs. Missing values are not silently converted to zero.

### Expression heatmap

The expression heatmap selects variable quantitative proteins, applies row-wise display scaling, and provides independent row and column choices of **Hierarchical**, **K-means**, or **None**. Users can choose row labels, clustering parameters, annotation palettes, and vector-PDF dimensions.

### Machine-learning feature selection

- **Random forest** uses repeated stability selection with stratified resampling, training-set-only variance prefiltering, configurable tree number and `mtry`, and RF frequency/Gini stability summaries.
- **L1** supports LASSO/elastic-net behavior through `alpha`, `lambda.1se` or `lambda.min`, automatic/manual CV folds, repeated selection-frequency summaries, and reproducible fold handling. Binary analyses use binomial modeling; multiclass analyses use ungrouped multinomial modeling.
- **RF + L1** uses a shared stratified split plan within each stability repeat and combines RF top-feature frequencies, L1 selection frequency, and mean Gini using configurable weights.

The workflows retain training-only variance filtering and training-derived imputation to reduce information leakage. Detailed settings and recommendations are documented in [README_ML_parameters_v2.1.0_Bilingual.txt](README_ML_parameters_v2.1.0_Bilingual.txt).

### Feature-protein analysis

- **Feature-protein UMAP** can use proteins selected by standalone RF, standalone L1, RF + L1 combined selection, or the union of currently available final ML results.
- **Feature-protein heatmap** uses the same selectable ML feature sources and the expression-heatmap framework.
- ML source availability is based on valid current-session results rather than historical output files.

### Slingshot pseudotime

Slingshot pseudotime uses the current protein quantitative matrix and user-assigned groups as clusters. PCA or UMAP can be selected for trajectory representation; users select a start group and may select an optional end group. It exports pseudotime coordinates, trajectory visualization, stage/group visualization, pseudotime-colored sample visualization, and top pseudotime-associated protein results when supported by the data and runtime packages.

### Analysis controls and output

Across modules, V2.1.0 includes practical refinements to parameter controls, plotting behavior, preview handling, output mapping, and edge-case safeguards. Plot modules generate vector PDF output and can export corresponding CSV data.

## Export, traceability, reproducibility, and offline use

Each completed analysis records traceability files in the selected output directory:

- `analysis_manifest.json` records the input/detected format, grouping, preprocessing, relevant settings, app version, runtime information, and generated outputs.
- `analysis_summary.html` is a browser-readable run summary.

V2.1.0 also includes source-level dependency and reproducibility records. The Windows portable release bundles the runtime and required R libraries for ordinary local use.

Normal application use is offline. Internet access is only relevant when deliberately updating annotations or installing missing R packages.

## Release notes

### ProteoPostZ Formal

- ProteoPostZ Formal Release V2.1.0
- [ProteoPostZ Formal Release V2.0.3](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.3)
- [ProteoPostZ Formal Release V2.0.2](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.2)
- [ProteoPostZ Formal Release V2.0.1](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.1)
- [ProteoPostZ Formal Release V2.0](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.0)

### ProteoDIAPostZ Formal

- [ProteoDIAPostZ Formal Release V1.4](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v1.4.0)
- [ProteoDIAPostZ Formal Release V1.3](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v1.3.0)

> **Version lineage:** ProteoPostZ is the continuation of ProteoDIAPostZ. The product was renamed beginning with V2.0 when support was expanded from DIA-focused workflows to DIA, DDA, and user-prepared standard quantitative matrices.

## 中文简要说明

ProteoPostZ Formal V2.1.0 是基于 R/Shiny 的 Windows 蛋白质组学后处理软件，支持 DIA、DDA 和用户整理的标准定量矩阵。

V2.1.0 主要加入 Quantitative QC 模块，用于在下游分析前检查样品和蛋白层面的缺失情况与丰度分布；同时完成内部代码模块化整理和可复现性增强，对 RF、L1 和 RF + L1 稳定性筛选流程进行了小幅方法学完善，并优化了部分分析模块的参数交互、结果预览、输出映射和边界条件处理。

软件支持 DIA-NN、Spectronaut、FragPipe/MSFragger、PEAKS DB/LFQ、MaxQuant 和标准 CSV/TSV 定量矩阵，并提供定性分析、Quantitative QC、相关性热图、丰度排序、CV、PCA、UMAP、t-SNE、差异分析、表达热图、机器学习特征筛选、Feature-protein 分析、理化性质分析和 Slingshot 等功能。

机器学习模块的详细参数说明见 [README_ML_parameters_v2.1.0_Bilingual.txt](README_ML_parameters_v2.1.0_Bilingual.txt)。

软件支持矢量 PDF、CSV、`analysis_manifest.json` 和 `analysis_summary.html` 输出。Windows portable release 可在本地离线运行。
