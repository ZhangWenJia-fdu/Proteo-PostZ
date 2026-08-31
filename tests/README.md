# Regression tests

Run the reusable, version-independent synthetic regression suite from the repository root:

```text
Rscript tests/run_all_regression.R
```

The suite creates canonical synthetic fixtures in `tempdir()` and does not require files from `F:\test`, network access, `renv::restore()`, or a packaged runtime. It covers supported importer shapes, Standard matrix zero semantics, quantitative QC, preprocessing and dimension reduction, core analysis outputs, offline annotation, ML paths, and trajectory safeguards. Optional analyses are reported as `SKIP` when their package is unavailable; required failures return a non-zero exit status.

Generated PDFs and CSVs are temporary by default. Set `PROTEOPOSTZ_KEEP_TEST_ARTIFACTS=1` to retain them for inspection. The suite checks output existence, schemas, dimensions, finite values, and deterministic results; it is not a replacement for manual Shiny visual checks or real-data validation. Existing legacy regression scripts remain separate and may use their configured real test inputs.

The fixture and test names are deliberately version-agnostic so this suite can be reused by later releases without renaming its infrastructure.

For daily development, use the standard suite above. Before a release, enable the slower high-fidelity ML validation from PowerShell:

```powershell
$env:PROTEOPOSTZ_HIGH_FIDELITY="1"
Rscript tests/run_all_regression.R
Remove-Item Env:PROTEOPOSTZ_HIGH_FIDELITY
```

This mode adds deterministic 320-protein/24-sample binary and 320-protein/30-sample multiclass ML checks, including full ranking invariance, score formulas, fold fields, seed sensitivity, and repeated-run reproducibility. Its seed `2026` is validation-only; production GUI and function defaults remain `123`.
