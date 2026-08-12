ProteoPostZ Formal V2.0.2
Developed by Wenjia Zhang

1. How to Start
1. Recommended: double-click ProteoPostZ_v2.0.2.exe.
2. Alternative: run Run_ProteoPostZ_v2.0.2.cmd.
3. Local browser address: http://127.0.0.1:3840/

The exe is a launcher for the local R/Shiny application. It is not an online website. The bundled portable R runtime and R libraries are used locally, and the browser interface is available only on the current computer.

2. Supported Inputs
- DIA-NN and Spectronaut protein-level DIA results.
- FragPipe/MSFragger, PEAKS, and MaxQuant protein-level DDA results.
- User-prepared standard quantitative matrices.

PEAKS accepts one protein file and detects its subtype from the internal schema. DB uses Area <sample> and sample-specific Coverage(%) columns. LFQ uses <sample> Area and excludes Group N Area and ratio/profile fields. If DB and LFQ fields are mixed, DB takes precedence. LFQ sample-level identification counts are unavailable; quantified counts remain available.

3. Main Features
Sample counts, Venn/UpSet, physicochemical properties, correlation heatmaps, rank-abundance, CV, PCA, UMAP, t-SNE, volcano, expression heatmaps, RF, L1, RF + L1, feature UMAP/heatmap, and Slingshot pseudotime.

PDF vector output and per-module CSV export are supported.

4. Offline Use
Normal use does not require internet access. Internet is needed only for explicit annotation updates or package installation.

5. Package Boundary
This package excludes local test files, test records, generated outputs, logs, and development history.
