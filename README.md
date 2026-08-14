# ProteoPostZ

ProteoPostZ, formerly ProteoDIAPostZ, is a Windows/R Shiny application for post-processing protein-level DIA and DDA proteomics results, as well as user-prepared standard quantitative matrices.

**Latest release: ProteoPostZ Formal V2.0.3.**

## Download the Windows portable application

Download [`ProteoPostZ_v2.0.3_windows_x86_release.zip`](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.3) from the V2.0.3 Release Assets after the release is published.

Expected portable-package name: `ProteoPostZ_v2.0.3_windows_x86_release.zip`.

The GitHub-generated **Source code (zip)** and **Source code (tar.gz)** downloads contain source only; they are not the complete Windows portable application. Windows users should download the portable ZIP from Release Assets.

## What's New in V2.0.3

ProteoPostZ Formal V2.0.3 is a focused refinement of the three protein-level machine-learning feature-selection modules.

- **Random forest, L1, and RF + L1.** All three modules now expose repeated stability selection, training-set proportion control, training-only variance prefiltering, configurable Top features, and model-performance/stability summaries.
- **Random forest.** RF stability ranking combines top-20 frequency, top-50 frequency, and mean Gini according to user-configurable weights; tree number and `mtry` can be adjusted.
- **L1.** L1 supports `alpha`, `lambda.1se`/`lambda.min`, automatic or manual cross-validation folds, multinomial modeling, and feature-selection frequency summaries.
- **RF + L1.** RF and L1 use the same stratified train/test split within each stability repeat. The combined score uses RF top-20 frequency, RF top-50 frequency, L1 selection frequency, and mean Gini with configurable component weights.

See [README_ML_parameters_v2.0.3_Bilingual.txt](README_ML_parameters_v2.0.3_Bilingual.txt) for the detailed bilingual parameter guide and [RELEASE_NOTES_v2.0.3.md](RELEASE_NOTES_v2.0.3.md) for the concise change record.

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

- **Random forest feature selection** uses repeated stability selection. Each repeat uses a stratified training/test split, applies a training-set-only variance prefilter, fits the configured number of trees and `mtry`, and records top-20 frequency, top-50 frequency, mean Gini, feature rankings, and model-performance summaries.
- **L1 feature selection** supports LASSO or elastic-net behavior through `alpha`, cross-validated `lambda.1se` or `lambda.min` selection, automatic/manual CV folds, multinomial modeling, selection frequency, and model-performance summaries.
- **RF + L1 combined feature selection** uses one shared stratified split per stability repeat, combines RF top-20 frequency, RF top-50 frequency, L1 selection frequency, and mean Gini with configurable weights, and exports the selected-feature union and corresponding quantity matrix.

The models support reproducible seeds, cross-validation-only or stratified train/test modes, training proportion, repeated stability selection, training-only variance prefiltering, final Top features, and model-specific settings. Outputs include cross-validation predictions, confusion matrices, class metrics, summary metrics, stability-selection summaries, and resample-performance summaries; binary analyses additionally produce ROC/AUC when calculable. Default strict mode requires adequate replication; optional small-sample mode is explicitly exploratory and should not be interpreted as independent validation.

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

## Release notes

### ProteoPostZ Formal

- [ProteoPostZ Formal Release V2.0.3](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.3)
- [ProteoPostZ Formal Release V2.0.2](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.2)
- [ProteoPostZ Formal Release V2.0.1](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.1)
- [ProteoPostZ Formal Release V2.0](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v2.0.0)

### ProteoDIAPostZ Formal

- [ProteoDIAPostZ Formal Release V1.4](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v1.4.0)
- [ProteoDIAPostZ Formal Release V1.3](https://github.com/ZhangWenJia-fdu/Proteo-PostZ/releases/tag/v1.3.0)

> **Version lineage:** ProteoPostZ is the continuation of ProteoDIAPostZ. The product was renamed beginning with V2.0 when support was expanded from DIA-focused workflows to DIA, DDA, and user-prepared standard quantitative matrices.

## 中文简要说明

ProteoPostZ Formal V2.0.3 是基于 R/Shiny 的 Windows 蛋白质组学后处理软件，支持 DIA、DDA 和用户整理的标准定量矩阵。

V2.0.3 主要优化三个机器学习模块：Random forest、L1 和 RF + L1 联合特征筛选。三个模块支持稳定性重复筛选、训练集比例、训练集方差预筛选、Top features 和模型表现汇总；RF + L1 在每次稳定性重复中使用同一套分层训练集/测试集拆分；L1 使用多分类 multinomial 建模并汇总特征筛选频率。

软件支持 DIA-NN、Spectronaut、FragPipe/MSFragger、PEAKS DB/LFQ、MaxQuant 和标准 CSV/TSV 定量矩阵，并提供相关性热图、表达热图、Feature-protein 热图、火山图、PCA、UMAP、t-SNE、ML-based Feature-protein 分析和 Slingshot 等功能。

机器学习模块的详细参数说明见 [README_ML_parameters_v2.0.3_Bilingual.txt](README_ML_parameters_v2.0.3_Bilingual.txt)。

软件支持 PDF 矢量图、CSV 数据、`analysis_manifest.json` 和 `analysis_summary.html` 导出，正常使用不需要联网。
