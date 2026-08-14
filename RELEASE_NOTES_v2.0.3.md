# ProteoPostZ Formal V2.0.3 Release Notes

ProteoPostZ Formal V2.0.3 is a focused refinement of the three protein-level machine-learning feature-selection modules. Existing V2.0.2 input formats and quantitative-analysis features remain supported.

## Machine-learning workflow refinements

- Random forest, L1, and RF + L1 each support repeated stability selection, stratified evaluation, training-set proportion control, training-set-only variance prefiltering, configurable final Top features, and model-performance summaries.
- Each stability repeat calculates preprocessing from its training set: top-variance features are selected from training data, training medians are used for imputation, and those medians are applied to both training and test data.
- Stability-selection CSV files report feature-selection frequencies or component scores, while resample-performance CSV files summarize model behavior across repeats.

## Random forest

- Configurable stability repeats, training proportion, top-variance prefilter, tree number, and automatic or manual `mtry`.
- Stability ranking uses RF top-20 frequency, RF top-50 frequency, and mean Gini with configurable component weights.
- Mean Gini is accumulated across the full number of stability repeats, with zero contribution for a feature absent from a repeat.
- Final output is controlled by Top RF selected features.

## L1

- Configurable stability repeats, training proportion, top-variance prefilter, `alpha`, `lambda.1se` or `lambda.min`, and automatic or manual cross-validation folds.
- Multiclass feature selection uses multinomial modeling with ungrouped coefficients.
- Stability output reports L1 feature-selection frequency and resample performance.
- Final output is controlled by Top L1 features.

## RF + L1 combined feature selection

- RF and L1 use the same stratified train/test split within each stability repeat.
- The combined stability score consists of RF top-20 frequency, RF top-50 frequency, L1 selection frequency, and mean Gini.
- The four component weights are user-configurable and normally sum to 1.
- Final output is controlled by Top combined features and includes the selected-feature union with its corresponding quantitative matrix.

## Documentation

`README_ML_parameters_v2.0.3_Bilingual.txt` provides detailed Chinese and English explanations of the RF, L1, and RF + L1 parameters and practical starting settings.

## Compatibility

- PEAKS DB/LFQ support from V2.0.2 remains available.
- DIA-NN, Spectronaut, FragPipe/MSFragger, MaxQuant, and standard quantitative-matrix inputs remain available.
- Correlation heatmap, expression heatmap, and Feature-protein heatmap improvements from V2.0.2 remain available.
