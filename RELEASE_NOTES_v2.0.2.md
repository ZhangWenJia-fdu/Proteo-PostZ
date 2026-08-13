# ProteoPostZ Formal V2.0.2 Release Notes

ProteoPostZ Formal V2.0.2 updates PEAKS Online protein-result input while preserving the V2.0.1 analysis framework and other importers.

## PEAKS input

- Upload one PEAKS protein result file. The app detects DB versus LFQ from the internal column schema, not from the filename.
- DB schema uses `Area <sample>` for quantitative values and sample-specific `Coverage(%) <sample> > 0` for identification evidence.
- LFQ schema uses `<sample> Area` for quantitative values and excludes `Group N Area`, profile/ratio fields, significance, and other non-sample fields.
- If DB and LFQ-style fields are mixed, DB schema takes precedence.
- Both subtypes retain the first `Top == TRUE` row per `Protein Group`; if no Top row exists, the first source row in the group is retained. No within-group Area aggregation is performed.

## Count semantics

- PEAKS DB: `Identified_Protein_Count` uses sample-specific Coverage > 0; `Quantified_Protein_Count` uses non-missing Area.
- PEAKS LFQ: `Identified_Protein_Count` is unavailable (`NA`); `Quantified_Protein_Count` uses non-missing LFQ Area.
- Existing FragPipe/MSFragger, MaxQuant, DIA-NN, Spectronaut, and standard-matrix definitions remain unchanged.

## Runtime

- The Windows package includes the portable R runtime, portable R libraries, app source, annotations, launcher, scripts, and bilingual instructions.
- Normal operation is offline and opens the local Shiny interface at `http://127.0.0.1:3840/`.
- Local test files, test records, logs, generated outputs, and source-repository history are not included in the package.

## Quantitative plot refinements

- Venn and UpSet display their group-level set-definition note within each module.
- Quantitative plots are organized into five fixed-viewport pages; each analysis card supports independent scrolling and full-screen expansion.
- Sample correlation heatmap has three ordering modes: Original input order, By group, and global hierarchical clustering. Hierarchical clustering can use `1 - correlation` or Euclidean distance of correlation profiles, complete or average linkage, and a user-selected dendrogram cut `k` with cluster gaps. Pearson or Spearman is selected before the clustering distance is calculated.
- Expression heatmap and Feature-protein heatmap independently support hierarchical, K-means, or no row/column clustering. Row clustering k defaults to 1; sample/column clustering k defaults to the configured Group count. Cluster boundaries are displayed for both cut-tree hierarchical and K-means modes.
- Both heatmaps provide four expression color schemes (Blue-White-Red, Red-White-Blue, Purple-White-Orange, Green-Black-Red), independent Group and sample-cluster annotation palettes, optional per-label custom `#RRGGBB` colors, and optional row labels from Accession, Protein name, or Gene name.
- Heatmap PDFs default to 600 x 800 pt. Row-label font size adapts to output height and the number of displayed proteins.
