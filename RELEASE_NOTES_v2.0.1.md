# ProteoPostZ Formal V2.0.1 Release Notes

ProteoPostZ Formal V2.0.1 is a patch release of Formal V2.0. It keeps the existing input formats and downstream analysis framework while clarifying the distinction between identified and quantified proteins.

## Changes

- Added `Quantified_Protein_Count` while retaining `Identified_Protein_Count` in previews, sample-count summaries, and CSV exports.
- Sample protein count barplots can switch between identified protein count and quantified protein count.
- Expression and Feature-protein heatmaps now provide independent row and column choices for Hierarchical clustering, K-means, or None. Row and column K-means `k` values are shown separately when applicable.
- FragPipe/MSFragger identification counts use sample `Spectral Count > 0`; quantified counts use non-missing `MaxLFQ Intensity` after zero values are treated as missing.
- PEAKS identification counts use sample-specific `Coverage > 0`; quantified counts use non-missing `Area` after zero values are treated as missing. The existing Protein Group deduplication rule is unchanged.
- MaxQuant has no independent sample-level identification evidence in the current `proteinGroups.txt` importer, so its identification count is approximated by non-missing `LFQ intensity` and currently matches the quantified count.
- Standard quantitative matrices report available quantitative values under the selected zero-value mode and do not fabricate independent strict identification evidence.

The existing DIA-NN and Spectronaut parsing and downstream analysis behavior is retained. Formal V2.0.0 release notes remain unchanged.
