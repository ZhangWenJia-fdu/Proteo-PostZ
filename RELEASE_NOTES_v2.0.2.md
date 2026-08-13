# ProteoPostZ Formal V2.0.2 Release Notes

ProteoPostZ Formal V2.0.2 adds PEAKS LFQ protein-result support and refines selected quantitative-heatmap controls while preserving the V2.0.1 analysis framework and other importers.

## PEAKS LFQ protein-result input

- PEAKS still accepts **one protein file at a time**. The importer recognizes PEAKS DB and PEAKS LFQ result types from their internal column schema, not from the filename.
- **PEAKS DB protein result:** `Area <sample>` supplies quantitative values; sample-specific `Coverage(%) <sample> > 0` remains the independent identification evidence.
- **PEAKS LFQ protein result:** `<sample> Area` supplies the protein-level LFQ abundance. `Group N Area`, ratio/profile columns, significance, and other non-sample fields are excluded.
- If DB- and LFQ-style fields occur in the same PEAKS file, the DB schema takes precedence.
- Each `Protein Group` retains one representative protein: the first source-order `Top == TRUE` row, or the first source row when no Top row is present. Protein-group Areas are not summed, averaged, or recalculated.
- PEAKS LFQ does not invent a sample-level `Identified_Protein_Count`, because the exported LFQ schema has no independent sample-specific identification evidence. It reports `Identified_Protein_Count` as unavailable and still generates `Quantified_Protein_Count` from non-missing LFQ Area values.
- PEAKS DB count behavior remains compatible with V2.0.1. When available, an LFQ protein result is the recommended PEAKS label-free quantitative input; ProteoPostZ uses the abundance already reported by PEAKS and does not recompute Top3 or MaxLFQ.

## Sample correlation heatmap

- Added clearer sample-order choices: original input order, assigned-group order with group gaps, or global hierarchical clustering.
- Hierarchical clustering can use `1 - correlation` or Euclidean distance of correlation profiles, complete or average linkage, and a user-defined dendrogram cut `k` with cluster gaps.
- Pearson or Spearman correlation is selected before the correlation matrix and its clustering distance are calculated.
- The module retains adjustable correlation-number precision, heatmap color scheme, legend range, vector-PDF dimensions, and CSV export.

## Expression heatmap and Feature-protein heatmap

- Both heatmaps now offer independent row and sample-column clustering modes: Hierarchical, K-means, or None, with corresponding cluster-number controls and visible cluster boundaries.
- Added four heatmap color schemes: Blue–White–Red, Red–White–Blue, Purple–White–Orange, and Green–Black–Red.
- Added separate Group and sample-cluster annotation palettes, plus optional dynamic custom `#RRGGBB` colors.
- Added row-label choices: accession, protein name, gene name, or no labels.
- Improved output layout with a default 600 × 800 pt heatmap PDF and adaptive row-label sizing. Output dimensions and CSV export remain configurable.
- Feature-protein heatmap continues to support RF, L1, RF + L1, or available-union feature sources.

## Interface and compatibility

- Quantitative analyses are arranged across fixed-viewport pages; individual analysis cards can scroll independently and support full-screen expansion.
- Venn and UpSet modules each display the same group-level set-definition note near their replicate setting.
- Existing FragPipe/MSFragger, MaxQuant, DIA-NN, Spectronaut, and standard-matrix importer/count semantics are unchanged.
