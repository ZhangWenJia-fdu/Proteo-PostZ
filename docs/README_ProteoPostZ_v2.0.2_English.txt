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
Qualitative analysis includes sample protein counts, Venn/UpSet group-level sets, and physicochemical properties. Quantitative analysis includes correlation heatmaps, rank-abundance, within-group CV, PCA, UMAP, t-SNE, volcano plots, expression heatmaps, Random forest selection, L1 selection, RF + L1 combined selection, feature-protein UMAP/heatmap, and Slingshot pseudotime.

Correlation heatmap: choose Pearson or Spearman, then retain Original input order, show By group, or use global hierarchical clustering. Hierarchical clustering provides 1 - correlation or Euclidean distance of correlation profiles, complete or average linkage, and a dendrogram cut k with cluster gaps.

Expression heatmap and Feature-protein heatmap: choose Hierarchical, K-means, or None independently for protein rows and sample columns. Row k defaults to 1; sample-column k defaults to the configured Group count. Both heatmaps provide four expression color schemes (Blue-White-Red, Red-White-Blue, Purple-White-Orange, Green-Black-Red), separate Group and sample-cluster annotation palettes, optional custom #RRGGBB annotation colors, and optional row labels from Accession, Protein name, or Gene name. The default heatmap PDF size is 600 x 800 pt.

Vector PDF output and per-module CSV export are supported. Analysis cards can be expanded to full screen; quantitative pages use fixed viewports with scrollable module controls and results.

4. Offline Use
Normal use does not require internet access. Internet is needed only for explicit annotation updates or package installation.

5. Package Boundary
This package excludes local test files, test records, generated outputs, logs, and development history.
