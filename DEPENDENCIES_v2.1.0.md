# ProteoPostZ V2.1.0 source dependencies

This file records the environment used for the V2.1.0 source validation. The
portable Windows application continues to use its bundled R and Rlibs; this
document and `renv.lock` are for source-code reproducibility only.

## Runtime

- ProteoPostZ: Formal V2.1.0
- R: 4.5.1 (2025-06-13 ucrt)
- Platform: x86_64-w64-mingw32
- Portable runtime: bundled `portable/R-4.5.1` and `portable/Rlibs`
- renv: not activated; no `.Rprofile` or automatic restore is required

## R packages

| Package | Version | Source |
|---|---:|---|
| data.table | 1.18.2.1 | CRAN |
| dplyr | 1.1.4 | CRAN |
| tidyr | 1.3.1 | CRAN |
| ggplot2 | 4.0.1 | CRAN |
| pheatmap | 1.0.13 | CRAN |
| RColorBrewer | 1.1.3 | CRAN |
| VennDiagram | 1.8.2 | CRAN |
| UpSetR | 1.4.0 | CRAN |
| ggridges | 0.5.7 | CRAN |
| uwot | 0.2.4 | CRAN |
| randomForest | 4.7.1.2 | CRAN |
| glmnet | 4.1.10 | CRAN |
| Peptides | 2.4.6 | CRAN |
| pdftools | 3.9.0 | CRAN |
| tibble | 3.3.0 | CRAN |
| shiny | 1.14.0 | CRAN |
| bslib | 0.9.0 | CRAN |
| DT | 0.34.0 | CRAN |
| slingshot | 2.16.0 | Bioconductor |
| SingleCellExperiment | 1.30.1 | Bioconductor |
| S4Vectors | 0.46.0 | Bioconductor |
| limma | 3.64.3 | Bioconductor |
| grid | 4.5.1 | base/recommended R |
| grDevices | 4.5.1 | base/recommended R |
| graphics | 4.5.1 | base/recommended R |
| stats | 4.5.1 | base/recommended R |
| utils | 4.5.1 | base/recommended R |
| tools | 4.5.1 | base/recommended R |

The versions above were read from installed package metadata with
`packageVersion()` and from `R.version`; optional packages are represented as
unavailable in runtime manifests if they are absent.
