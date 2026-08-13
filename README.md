# ProteoPostZ

ProteoPostZ, formerly ProteoDIAPostZ, is a Windows/R Shiny application for post-processing protein-level DIA and DDA proteomics results, as well as user-prepared standard quantitative matrices.

**Latest release: ProteoPostZ Formal V2.0.2.**

## Download the Windows portable application

Download [`ProteoPostZ_v2.0.2_windows_x86_release.zip`](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.2) from the V2.0.2 Release Assets. 

The GitHub-generated **Source code (zip)** and **Source code (tar.gz)** downloads contain source only; they are not the complete Windows portable application. Windows users should download the portable ZIP from Release Assets.

## What's New in V2.0.2

- **PEAKS LFQ protein-result input.** Upload one PEAKS protein file at a time. The application detects PEAKS DB and PEAKS LFQ result schemas from their columns, retains DB compatibility, and supports PEAKS LFQ protein abundance as the recommended label-free quantitative input when an LFQ result is available.
- **Sample correlation heatmap.** Sample ordering, hierarchical clustering, and layout controls have been refined.
- **Expression and Feature-protein heatmaps.** Clustering and ordering controls, heatmap palettes, annotation colors, row labels, and output layout have been expanded.

See [RELEASE_NOTES_v2.0.2.md](RELEASE_NOTES_v2.0.2.md) for the concise V2.0.2 change record.

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
- **PEAKS** protein-level results remains one upload entry: **one PEAKS protein file is uploaded at a time**. The subtype is determined from the internal column schema, never from the filename.

- **PEAKS DB protein result**
- DB uses `Area <sample>` as the sample quantitative values.
- Sample-specific `Coverage(%) <sample> > 0` supplies independent sample-level identification evidence.
- This preserves the V2.0.1 PEAKS DB workflow.

- **PEAKS LFQ protein result**
- LFQ uses `<sample> Area` as the PEAKS-reported protein-level LFQ abundance.
- `Group N Area`, ratio/profile columns, significance, and other non-sample fields are excluded from the quantitative matrix.
- ProteoPostZ does not recalculate Top3, MaxLFQ, or protein abundance.
- PEAKS LFQ exports do not provide independent sample-specific identification evidence. Therefore `Identified_Protein_Count` is shown as unavailable, while `Quantified_Protein_Count` is calculated normally from non-missing LFQ Area values.
- If LFQ analysis has already been performed in PEAKS, the LFQ protein result is the preferred PEAKS quantitative input.

For either PEAKS subtype, one representative row is retained per `Protein Group`: the first source-order row with `Top == TRUE`; if no such row exists, the first row in the group. Areas from different members of a group are never summed or averaged. If a file contains both DB- and LFQ-style fields, the DB schema takes precedence.

### Standard quantitative matrix

- The first column is the protein/feature identifier and remaining columns are sample quantitative values.
- Users explicitly choose whether zero is a valid quantitative value or denotes missing/unquantified data.
- `NA`, `NaN`, blank, and empty values are treated as missing.

After loading, sample names and biological groups can be edited in the interface; identical group text denotes the same group.

## Main analysis features

### Qualitative analysis

- **Sample protein counts:** identification and quantified protein counts are kept separate whenever the selected importer provides independent identification evidence. Count availability and meaning follow the importer schema rather than forcing one definition on every format.
- **Venn diagram:** available for 2–4 non-empty groups. **UpSet plot:** supports group-set intersections and is recommended for 5 or more groups. Both use group-level protein sets, not sample-level sets. The *Minimum replicates detected in group* setting determines how many detected replicates are required for a protein to enter each group set.
- **Physicochemical-property analysis:** compares group protein sets using accession-based offline annotations, including GRAVY, molecular weight, pI, length, transmembrane-helices, and subcellular-class summaries. Built-in annotations and a compatible user annotation table can be selected; unmatched accessions remain appropriately unavailable.

### Quantitative overview and dimension reduction

- **Sample correlation heatmap:** Pearson or Spearman correlation; sample order can remain original, follow the assigned groups with group gaps, or be rearranged by hierarchical clustering. Hierarchical ordering supports `1 - correlation` or Euclidean distance of correlation profiles, complete or average linkage, and a dendrogram-cut `k`. Group annotation, correlation color scheme, displayed precision, legend range, and vector-PDF dimensions are available.
- **Rank-abundance plot:** visualizes the per-sample distribution of protein abundance against abundance rank.
- **Within-group CV ridgeline:** calculates distributions of protein-level coefficient of variation within each assigned group, for groups with sufficient quantitative replicates; the displayed CV range is adjustable.
- **PCA:** produces scaled PC1/PC2 sample coordinates after valid-feature filtering, with selectable palette and vector-PDF/CSV export.
- **UMAP:** provides sample-level UMAP using valid-feature filtering, configurable `n_neighbors` and `min_dist`, and a fixed reproducible seed.
- **t-SNE:** provides two-dimensional sample embeddings with valid-feature filtering, configurable or automatic perplexity, configurable iteration count, and a fixed reproducible seed. Input-size safeguards are applied.

### Differential analysis

**Volcano plot** compares a selected reference group with a selected comparison group and requires at least two samples in each group. It supports limma moderated testing or two-sample t tests on `log2(x + 1)` values; log2 fold change calculated either as the difference of log2 means or as the log2 raw-mean ratio; raw-p-value or Benjamini-Hochberg adjusted-p-value thresholds; and independently adjustable log2FC and significance cutoffs. The exported table records the chosen calculation and test method, raw and BH-adjusted p values, and the comparison direction. Missing values remain unavailable to the relevant per-protein calculation rather than being silently converted to zero.

### Expression heatmap

The expression heatmap selects the most variable quantitative proteins, applies row-wise display scaling, and provides independent row and column choices of **Hierarchical**, **K-means**, or **None**. Hierarchical and K-means modes have their own row/column cluster-number controls; `None` preserves the applicable input order. Users can choose protein accession, protein name, gene name, or no row labels.

Four heatmap color schemes are provided: Blue–White–Red, Red–White–Blue, Purple–White–Orange, and Green–Black–Red. Group and sample-cluster annotations have separate palette choices and optional custom `#RRGGBB` colors. The default heatmap PDF is 600 × 800 pt, with adaptive row-label sizing and adjustable output dimensions.

### Machine-learning feature selection

- **Random forest feature selection** ranks proteins by model importance (MeanDecreaseGini when available).
- **L1 feature selection** supports LASSO or elastic-net behavior through `alpha`, with cross-validated `lambda.1se` or `lambda.min` selection.
- **RF + L1 combined feature selection** exports the selected-feature union and its corresponding quantity matrix.

The models support reproducible seeds, automatic/cross-validation-only/train-test modes, train proportion, and model-specific settings. Outputs include cross-validation predictions, confusion matrices, class metrics, and summary metrics; binary analyses additionally produce ROC/AUC when calculable. Stability-selection summaries are produced for RF and L1. Default strict mode requires adequate replication; optional small-sample mode is explicitly exploratory and should not be interpreted as independent validation.

### Feature-protein analysis

- **Feature-protein UMAP** can use RF-selected proteins, L1-selected proteins, RF + L1 proteins, or their available union, with top-feature, neighbor, and minimum-distance settings.
- **Feature-protein heatmap** uses the same selectable feature sources and the same independent row/column clustering, palettes, annotation colors, row-label choices, CSV export, and vector-PDF layout controls as the expression heatmap.

### Slingshot pseudotime

Slingshot pseudotime analysis uses the current protein quantitative matrix and user-assigned groups as clusters. PCA or UMAP can be selected for the trajectory representation; users select a start group and may select an optional end group. It exports pseudotime coordinates, trajectory plots, and top pseudotime-associated protein outputs when the required runtime packages and data conditions are available.

## Export, traceability, and offline use

Plot modules generate vector PDF output and can export the corresponding CSV data. Each completed analysis also records traceability files in the chosen output directory:

- `analysis_manifest.json` records the input/detected format, grouping, preprocessing, relevant settings, app version, and generated outputs.
- `analysis_summary.html` is a browser-readable run summary.

Normal application use is offline. Internet access is only relevant when deliberately updating annotations or installing missing R packages.

---

## Release Notes

### ProteoPostZ Formal

* [ProteoPostZ Formal Release V2.0.2](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.2)
* [ProteoPostZ Formal Release V2.0.1](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.1)
* [ProteoPostZ Formal Release V2.0](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.0)

### ProteoDIAPostZ Formal

* [ProteoDIAPostZ Formal Release V1.4](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v1.4.0)
* [ProteoDIAPostZ Formal Release V1.3](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v1.3.0)

> **Version lineage:** ProteoPostZ is the continuation of ProteoDIAPostZ. The product was renamed beginning with V2.0 when support was expanded from DIA-focused workflows to DIA, DDA, and user-prepared standard quantitative matrices.

---

## 中文简要说明

**ProteoPostZ Formal V2.0.2** 是一个基于 R/Shiny 的 Windows 蛋白质组学后处理软件，主要用于 **DIA、DDA 蛋白水平结果以及用户整理的标准定量矩阵**的后续统计分析与可视化。

ProteoPostZ 原名 **ProteoDIAPostZ**。自 V2.0 起，软件的输入范围由 DIA 扩展至 DIA、DDA 和标准定量矩阵，因此正式更名为 ProteoPostZ，并延续原有版本体系。

### 支持的输入

目前主要支持：

* **DIA-NN** 蛋白水平结果
* **Spectronaut** 蛋白水平结果
* **FragPipe/MSFragger** 蛋白水平结果
* **PEAKS DB** protein result
* **PEAKS LFQ** protein result
* **MaxQuant** `proteinGroups.txt`
* 用户整理的 CSV/TSV 标准定量矩阵

对于 PEAKS，用户每次只需上传一个 protein result 文件。软件根据文件内部 column schema 自动识别 **PEAKS DB** 或 **PEAKS LFQ**，不依赖文件名判断。

PEAKS DB 保留原有兼容模式，以 `Area <sample>` 作为样品定量值，并利用 sample-specific `Coverage(%)` 作为独立的样品水平鉴定证据。

PEAKS LFQ 使用 `<sample> Area` 作为 PEAKS 输出的蛋白水平 LFQ 定量值，并自动排除 `Group N Area`、ratio/profile 等非样品定量字段。如果用户已经在 PEAKS 中完成 LFQ 分析，建议优先使用 LFQ protein result 进行后续定量分析。

### 主要分析功能

定性分析主要包括：

* 样品蛋白鉴定数和定量数统计
* Venn diagram
* UpSet plot
* 蛋白理化性质分析

定量分析主要包括：

* Pearson / Spearman 样品相关性热图
* Rank-abundance plot
* 组内 CV ridgeline
* PCA
* UMAP
* t-SNE
* Volcano plot
* Expression heatmap
* Random forest 特征筛选
* L1 特征筛选
* RF + L1 联合特征筛选
* Feature-protein UMAP
* Feature-protein heatmap
* Slingshot pseudotime analysis

V2.0.2 进一步优化了 **Sample correlation heatmap** 的样品排序和层次聚类设置，并扩展了 **Expression heatmap** 和 **Feature-protein heatmap** 的参数、颜色及输出布局选项。

### 输出与离线使用

各主要分析模块支持：

* PDF 矢量图导出
* 对应 CSV 数据导出
* `analysis_manifest.json` 分析过程记录
* `analysis_summary.html` 分析摘要

完整 Windows portable package 已包含运行所需的 portable R、R packages 和内置 annotation tables，正常分析过程不需要互联网连接。

---

