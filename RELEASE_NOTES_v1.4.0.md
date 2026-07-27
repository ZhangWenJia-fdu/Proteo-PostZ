# ProteoDIAPostZ Formal V1.4 Release Notes

This file documents the source-level changes for ProteoDIAPostZ Formal V1.4. The complete portable Windows package is distributed separately as a GitHub Release asset and is not committed to the source repository.

## Scope

V1.4 is built from the V1.3 application line. It preserves DIA-NN and Spectronaut input semantics and adds support for standard quantitative matrix input.

## New Features

- Added `standard_matrix` as a third input source.
- Added CSV/TSV standard quantitative matrix parsing.
- Added explicit zero-handling modes:
  - `zero_is_value`: numeric zero is retained as a valid quantitative value.
  - `zero_is_missing`: numeric zero is converted to `NA` in `analysis_quant_matrix` while the raw matrix is retained.
- Added validation for duplicate/blank feature identifiers, duplicate/blank sample names, illegal numeric text, `Inf`/`-Inf`, delimiter mismatch, and all-missing sample columns.
- Added standard-matrix summaries based on the number of available quantitative values.
- Connected standard matrices to generic downstream quantitative workflows that operate on feature-by-sample matrices.

## Compatibility Preserved

- DIA-NN input parsing remains unchanged from V1.3.
- Spectronaut input parsing remains unchanged from V1.3.
- Spectronaut continues to use `PG.IBAQ` for identification/qualitative summaries and `PG.Quantity` for quantitative analyses.
- Existing plot parameters and downstream statistical defaults were not intentionally changed for DIA-NN or Spectronaut workflows.

## Fixes

- Improved physicochemical accession mapping for DIA-NN and Spectronaut when selected row identifiers such as Protein name or Gene name contain duplicates and are internally made unique.
- Standard matrix feature identifiers are matched to annotation `Accession` as-is; accession-based annotation may not match if users provide gene symbols, protein names, or other identifiers.

## Interface Changes

- Shortened the visible title to `ProteoDIAPostZ Formal V1.4`.
- Added `Modern Chromatography Separation and Analysis` to the developer identity block.
- Fixed the Input-page left sidebar so upper input controls scroll independently while Output directory, Load data, error text, and developer identity stay in a fixed bottom section.
- Adjusted the developer identity block background to match the sidebar instead of appearing as a white card.

## Source/Package Boundary

The GitHub source repository includes source code, tests, documentation, annotation tables, launcher source, and lightweight scripts. It does not include portable R, packaged R libraries, generated outputs, runtime logs, zip packages, or compiled exe files.

The complete portable Windows software package should be uploaded as a GitHub Release asset.

## Validation Summary

Representative checks passed during the V1.4 update:

```text
STANDARD_MATRIX_TEST_OK 9 7
STANDARD_MATRIX_SHINY_TEST_OK csv_value=18 csv_missing=16 tsv_value=18 tsv_missing=16 diann=4862x8 spectronaut=3734x8
STANDARD_MATRIX_DOWNSTREAM_SMOKE_OK mode_a=22 mode_b=20 zero_mode_b_is_na=TRUE
V14_INPUT_REGRESSION_OK diann=4862x8 spectronaut=3734x8 standard_csv_A=9 standard_csv_B=7 standard_tsv_A=4 standard_tsv_B=3 full_A=286 full_B=284 rf_top=8 slingshot=ok
E_REPO_PARSE_OK
F_V14_PARSE_OK
DEPENDENCY_CHECK_OK
```