ProteoPostZ Formal V2.0.2
Developed by Wenjia Zhang

1. How to Start
1. Recommended: double-click ProteoPostZ_v2.0.2.exe.
2. Alternative: run Run_ProteoPostZ_v2.0.2.cmd.
3. Local browser address: http://127.0.0.1:3840/

The exe is a launcher for the local R/Shiny application. It is not an online website. The bundled portable R runtime and R libraries are used locally, and the browser interface is available only on the current computer.

2. What's New in V2.0.2
- PEAKS LFQ protein-result input is now supported.
- Upload one PEAKS protein file; the app automatically detects PEAKS DB protein result or PEAKS LFQ protein result from the internal column schema, without relying on the filename.
- PEAKS DB remains supported for compatibility. If DB and LFQ-style fields are mixed in one file, the file is read as PEAKS DB.
- Sample correlation heatmap, Expression heatmap, and Feature-protein heatmap analysis settings have also been refined, including ordering or clustering choices, annotation colors, row labels, and improved heatmap output layout.

3. Supported Inputs
- DIA-NN and Spectronaut protein-level DIA results.
- FragPipe/MSFragger, PEAKS, and MaxQuant protein-level DDA results.
- User-prepared standard quantitative matrices.

PEAKS accepts one protein file and detects its subtype from the internal schema. DB uses Area <sample> and sample-specific Coverage(%) columns. LFQ uses <sample> Area and excludes Group N Area and ratio/profile fields. If DB and LFQ fields are mixed, DB takes precedence. LFQ sample-level identification counts are unavailable; quantified counts remain available.

4. Main Features
Qualitative analysis includes sample protein counts, Venn diagrams, UpSet plots, and physicochemical-property distributions. Group-level protein or feature sets are controlled by the selected minimum replicate requirement.

Quantitative analysis includes sample correlation heatmaps, rank-abundance plots, within-group CV ridgelines, PCA, UMAP, t-SNE, volcano plots, expression heatmaps, Random forest feature selection, L1 feature selection, RF + L1 combined feature selection, feature-protein UMAP, feature-protein heatmaps, and Slingshot pseudotime. The app also provides sample grouping, input preview, count summaries, and output-file browsing.

Analysis settings are adjustable within each module, including applicable thresholds, clustering, dimensions, statistical or machine-learning settings, palettes, labels, and PDF dimensions. Vector PDF figures and corresponding CSV data can be exported per module. Analysis cards can be expanded to full screen; quantitative pages use fixed viewports with scrollable module controls and results.

5. Offline Use
Normal use does not require internet access. Internet is needed only for explicit annotation updates or package installation.
