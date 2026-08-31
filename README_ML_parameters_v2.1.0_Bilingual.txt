ProteoPostZ V2.1.0
Machine-learning parameter guide / 机器学习参数说明

This guide covers the three protein-level feature-selection modules: Random forest, L1, and RF + L1 combined feature selection, together with the ML-selected feature sources used by Feature-protein UMAP and Feature-protein heatmap.
本说明适用于三个蛋白质水平特征筛选模块：单独随机森林、单独 L1、RF + L1 联合特征筛选，以及 Feature-protein UMAP 和 Feature-protein heatmap 使用的机器学习特征来源。

GENERAL SETTINGS / 通用设置

Random seed / 随机种子
Sets the reproducible random-number stream for data splitting and model fitting. With the same input data and the same settings, the same seed is intended to reproduce the same resampling sequence and result. Changing the seed produces another valid resampling sequence and may change feature-selection frequencies or rankings.
控制数据拆分和模型拟合所使用的可复现随机数序列。在输入数据和其他参数完全相同的情况下，相同随机种子用于复现相同的重复拆分和结果。更换随机种子会产生另一套有效的随机拆分序列，特征筛选频率或排序可能随之变化。

Default / 默认值: 123
Recommendation / 建议：比较不同参数设置时建议固定随机种子。只有在需要评估随机重采样敏感性时才更换种子。

Train/test split mode / 训练集/测试集拆分模式
Cross-validation only uses cross-validation for model assessment. Train/test split creates stratified training and test subsets for performance assessment and stability-selection repeats.
Cross-validation only 使用交叉验证进行模型评估。Train/test split 按组别进行分层训练集/测试集拆分，并用于模型表现评估和稳定性重复筛选。

Training set proportion / 训练集比例
The fraction of samples assigned to the training set in each stratified split. The remaining samples form the test set. In RF + L1, RF and L1 use the same outer split in each stability repeat.
每次分层拆分中分配给训练集的样本比例，剩余样本作为测试集。RF + L1 模块中，每次稳定性重复的 RF 与 L1 共用同一套外层拆分。

Default / 默认值: 0.7
Recommendation / 建议：0.7 可作为一般默认值。需要更大的训练集时可以使用 0.8，但相应测试集会更小。比较不同特征筛选参数时建议保持该值一致。

Allow small-sample exploratory ML / 允许小样本探索性机器学习
Strict mode applies minimum group-size checks before running the model. Enable the exploratory small-sample option only when the dataset does not satisfy the strict requirements and the analysis is intended for exploratory use.
严格模式会在模型运行前检查每组最小样本量。只有当数据无法满足严格要求且分析目的明确为探索性时，才建议启用小样本模式。

Recommendation / 建议：常规分析尽量使用严格模式。小样本模式不应被解释为与具有充分生物学重复的数据提供同等程度的独立验证。

Stability selection repeats / 稳定性筛选次数
Repeats the complete stability-selection procedure. Each repeat creates a reproducible stratified split, performs training-only preprocessing, fits the corresponding model, and records selected features and performance.
重复执行完整的稳定性筛选流程。每次重复都会产生可复现的分层拆分，只在训练集上完成相应预处理，拟合模型，并记录特征筛选和模型表现。

Default / 默认值: 50
Recommendation / 建议：50 次适合作为一般分析起点。如果特征排序不稳定且计算时间允许，可适当增加重复次数。

Stability prefilter: top variance features / 稳定性预筛选：Top variance features
Within each stability repeat, feature variance is calculated using the training data only. Only the specified number of highest-variance features is passed to the machine-learning model. This is a computational/noise-reduction prefilter and is not the final Top features output.
每次稳定性重复中，仅使用该次训练数据计算特征方差，并将方差最高的指定数量特征送入机器学习模型。该参数用于计算量和噪声控制，不等于最终输出的 Top features 数量。

Default / 默认值: 200
Recommendation / 建议：200 适合作为常规蛋白水平数据的起始设置。希望更多特征参与竞争且样本量和计算时间允许时可增加；数据很小或计算受限时可谨慎降低。

Top features / Top features 数
Controls the final number of ranked selected features displayed/exported by the standalone RF or L1 module after the relevant stability ranking has been calculated.
控制单独 RF 或 L1 模块完成稳定性排序后最终展示和导出的 Top 特征数量。

Typical starting value / 常用起始值: 50
Recommendation / 建议：根据后续解释目的选择适当数量。调整最终 Top-N 不应被理解为改变底层稳定性重复筛选过程本身。

Random forest trees / Random forest tree 数
The number of trees fitted in each random forest. Larger values generally make importance estimates less noisy but increase runtime.
每次随机森林拟合使用的树数量。增加树数通常可使重要性估计更平稳，但会增加运行时间。

Default / 默认值: 500
Recommendation / 建议：500 适合常规分析。只有在需要更稳定的重要性估计且计算时间允许时才建议增加。

Random forest mtry / Random forest mtry
The number of candidate features randomly considered at each tree split. Auto lets the software choose the data-dependent value used by the current RF workflow.
每棵树每次分裂时随机考虑的候选特征数量。Auto 由当前 RF 流程根据数据选择相应数值。

Default / 默认值: Auto
Recommendation / 建议：除非需要复现明确指定的比较或专门研究 mtry 的影响，否则建议使用 Auto。

PDF width / PDF width (pt) and PDF height / PDF height (pt)
Control exported vector-PDF dimensions in points. These settings affect figure layout only and do not change the machine-learning model.
控制导出矢量 PDF 的宽度和高度，单位为 pt。这些参数只影响图形布局，不改变机器学习模型。

RANDOM FOREST ONLY / 单独 RANDOM FOREST 模块

Stability weight: RF top 20 frequency / RF top 20 频率权重
Controls the contribution of the frequency with which a feature appears among the top 20 RF features across stability repeats.
控制某特征在稳定性重复中进入 RF 前 20 名的频率对最终 RF 稳定性评分的贡献。

Stability weight: RF top 50 frequency / RF top 50 频率权重
Controls the contribution of the frequency with which a feature appears among the top 50 RF features across stability repeats.
控制某特征在稳定性重复中进入 RF 前 50 名的频率对最终 RF 稳定性评分的贡献。

Stability weight: mean Gini / 平均 Gini 权重
Controls the contribution of mean RF decrease-in-Gini importance across stability-selection repeats. If a feature is absent from the relevant per-repeat RF subset, that repeat contributes zero before the mean is summarized.
控制稳定性重复中 RF Mean Decrease in Gini 重要性对最终评分的贡献。某特征若在某次重复的相关 RF 子集中不存在，则该次贡献按 0 处理后再汇总平均值。

The three RF components determine the RF stability score. Their weights are normally set to sum to 1. RF top-20 and top-50 frequencies are calculated independently of the final requested Top features value.
三个 RF 组成部分共同决定 RF 稳定性评分，三个权重通常应加和为 1。RF top-20 和 top-50 出现频率的计算独立于最终请求输出的 Top features 数量。

Recommendation / 建议：如果没有预先设定的分析理由，一般建议保留默认权重。调整权重时应保持权重和具有清晰解释，并尽量一次只改变一个组成部分。

L1 ONLY / 单独 L1 模块

L1 alpha / L1 alpha
Controls the balance between LASSO and elastic net. `alpha = 1` gives pure LASSO; values between 0 and 1 add a ridge component.
控制 LASSO 与 elastic net 的比例。`alpha = 1` 表示纯 LASSO；0 到 1 之间的值会加入 ridge 成分。

Default / 默认值: 1
Recommendation / 建议：标准稀疏 L1 筛选建议使用 1。如果预期多个相关蛋白共同携带信号，并希望模型不完全依赖纯稀疏选择，可考虑小于 1 的值。

Lambda selection / Lambda selection
`lambda.1se` selects the more regularized model within one standard error of the minimum cross-validation error and generally gives a simpler model. `lambda.min` selects the lambda with minimum cross-validation error and often retains more features.
`lambda.1se` 在最小交叉验证误差一个标准误范围内选择正则化更强、通常更简单的模型。`lambda.min` 选择交叉验证误差最小的 lambda，通常会保留更多特征。

Default / 默认值: lambda.1se
Recommendation / 建议：建议将 `lambda.1se` 作为较保守的默认选择。需要主动保留更广的特征集合，且更大的模型仍具有合理解释时，可使用 `lambda.min`。

Cross-validation folds / 交叉验证折数
Controls the number of folds used for lambda selection. Auto selects a valid fold count based on the minimum class size available in the current training/resampling subset. In V2.1.0, Auto uses up to 5 folds; a manually requested value is reduced when necessary to remain valid for the current class sizes.
控制选择 lambda 时使用的交叉验证折数。Auto 会根据当前训练/重复子集中最小类别样本数选择有效折数。V2.1.0 中 Auto 最多使用 5 折；手动指定的折数在当前类别样本量不足时会相应降低，以保持有效。

Default / 默认值: Auto
Recommendation / 建议：常规分析建议使用 Auto。需要复现特定分析且各类别样本量允许时，可明确设为 5 等固定值。对于较小训练集，不建议使用过大的折数。

Model family / 模型 family
For two-class analyses, V2.1.0 uses binomial modeling. For analyses with three or more classes, it uses multinomial modeling with `type.multinomial = "ungrouped"`.
两分类分析使用 binomial 模型；三分类及以上分析使用 multinomial 模型，并设置 `type.multinomial = "ungrouped"`。

This choice is automatic and is not exposed as a routine user parameter.
该设置由软件根据类别数自动确定，不作为常规用户参数开放。

L1 selection frequency / L1 筛选频率
The stability summary records how often each feature receives a non-zero selected coefficient across repeated resampling. The final standalone L1 output is ranked from this stability evidence and limited by Top features.
稳定性汇总记录每个特征在重复重采样中获得非零入选系数的频率。单独 L1 的最终输出依据该稳定性证据排序，并由 Top features 限制数量。

RF + L1 COMBINED / RF + L1 联合模块

The combined score uses four stability components:
1. RF top-20 frequency
2. RF top-50 frequency
3. L1 selection frequency
4. mean Gini

联合评分使用四项稳定性指标：
1. RF 前 20 频率
2. RF 前 50 频率
3. L1 筛选频率
4. 平均 Gini

Default weights / 默认权重:
RF top-20 frequency: 0.35
RF top-50 frequency: 0.25
L1 selection frequency: 0.30
Mean Gini: 0.10

The four weights normally sum to 1. Larger weight gives the corresponding component greater influence on the final combined ranking.
四个权重通常加和为 1。某项权重越大，该指标对最终联合排序的影响越大。

Top combined features / 联合 Top features 数
Controls the final number of features returned after the combined RF + L1 stability score has been calculated and ranked.
控制 RF + L1 联合稳定性评分完成并排序后最终返回的特征数量。

Typical starting value / 常用起始值: 50

Shared outer split / 共享外层拆分
Within each stability repeat, RF and L1 use the same stratified training/test split.
每次稳定性重复中，RF 和 L1 使用同一套分层训练集/测试集拆分。

Training-only preprocessing / 仅训练集预处理
Variance filtering and imputation parameters used by the stability workflow are derived from the training data within each repeat. Test data do not determine the training medians or variance-ranking step.
稳定性流程中的方差筛选和缺失值填补参数均在每次重复的训练数据内确定。测试数据不参与训练中位数或方差排序的计算。

MEAN GINI SCALING / 平均 GINI 标准化

For stability scoring, mean Gini importance is summarized over the full stability-eligible protein universe. If a protein is not present in the relevant per-repeat RF subset, its Gini contribution for that repeat is zero. The resulting mean Gini values are then scaled across the full eligible universe before entering the stability score.
用于稳定性评分时，平均 Gini 重要性在完整的 stability-eligible 蛋白集合上进行汇总。某蛋白若未进入某次重复的相关 RF 子集，则该次 Gini 贡献为 0。随后在完整 eligible universe 上对平均 Gini 进行标准化，再进入稳定性评分。

This keeps the Gini scale independent of the final requested Top-N output size.
这样可以使 Gini 标准化尺度不依赖最终请求输出的 Top-N 数量。

ML-SELECTED FEATURE SOURCES / 机器学习下游特征来源

Feature-protein UMAP and Feature-protein heatmap can use:
- Random forest selected proteins
- L1 selected proteins
- RF + L1 selected proteins
- Union of currently available ML selected proteins

Feature-protein UMAP 和 Feature-protein heatmap 可以使用：
- Random forest selected proteins
- L1 selected proteins
- RF + L1 selected proteins
- 当前 session 中可用机器学习结果的 union

Source availability / 来源可用性
A source is considered available only when the corresponding standalone/final ML analysis has been completed in the current application session.
只有在当前应用 session 中已经完成相应的独立/最终 ML 分析时，该来源才视为可用。

Running RF + L1 does not make its internal RF and L1 components count as standalone RF or standalone L1 analyses.
运行 RF + L1 不会把其内部 RF 和 L1 组成部分视为已经独立运行了 RF 或 L1。

Union / Union
Union is the set union of final selected-protein lists currently available in the session. If only one source is available, union is identical to that source.
Union 是当前 session 中真正可用的最终 selected-protein 列表的集合并集。如果当前只有一个来源可用，则 union 与该来源完全一致。

When more than one source is available, the current implementation preserves a deterministic source order (`rf` → `l1` → `rfl1`), removes duplicates, and then applies the requested downstream Top-N limit. Scores from different model types are not directly compared.
当多个来源同时可用时，当前实现按固定来源顺序（`rf` → `l1` → `rfl1`）保持各来源内部排序，合并后去重，再应用下游请求的 Top-N 限制。不同模型类型的分数不会被直接进行跨模型数值比较。

Feature-protein UMAP / Feature-protein UMAP
Feature-protein UMAP performs a sample-level UMAP using proteins selected from the chosen ML source. It supports Top-N selection, a minimum valid fraction, `n_neighbors`, and `min_dist`.
Feature-protein UMAP 使用所选 ML 来源的蛋白构建样本水平 UMAP。支持 Top-N、minimum valid fraction、`n_neighbors` 和 `min_dist` 设置。

Recommendation / 建议：比较两个分析时，应保持相同的 feature source 和相同的 UMAP 参数。如果当前只运行了 RF + L1，则 `rfl1` 与 `union` 应返回相同的 selected feature 列表。

Feature-protein heatmap / Feature-protein heatmap
Feature-protein heatmap uses proteins from the selected ML source and otherwise follows the expression-heatmap framework for clustering, scaling, labeling, palettes, and PDF output.
Feature-protein heatmap 使用所选 ML 来源的蛋白，其他聚类、缩放、标签、配色和 PDF 输出逻辑沿用表达热图框架。

PRACTICAL RECOMMENDATIONS / 实用建议

A practical starting configuration for ordinary datasets is:
- Random seed: 123
- Training set proportion: 0.7
- Stability repeats: 50
- Stability prefilter: Top 200 variance features
- RF trees: 500
- RF mtry: Auto
- L1 alpha: 1
- Lambda: lambda.1se
- Cross-validation folds: Auto
- Final Top features: 50
- RF + L1 weights: 0.35 / 0.25 / 0.30 / 0.10

常规数据可从以下设置开始：
- Random seed：123
- Training set proportion：0.7
- Stability repeats：50
- Stability prefilter：Top 200 variance features
- RF trees：500
- RF mtry：Auto
- L1 alpha：1
- Lambda：lambda.1se
- Cross-validation folds：Auto
- Final Top features：50
- RF + L1 权重：0.35 / 0.25 / 0.30 / 0.10

When evaluating the effect of a parameter, keep the random seed and other settings fixed and change one parameter at a time. For deliberate 80/20 analyses, set Training set proportion to 0.8. For deliberate fold-matched comparisons, set the same explicit fold value and ensure current class sizes support it.
评估某个参数影响时，建议保持随机种子和其他参数固定，并尽量一次只修改一个参数。需要明确进行 80/20 分析时，将 Training set proportion 设置为 0.8。需要固定折数进行比较时，在被比较分析中使用相同的显式 folds，并确保当前类别样本量能够支持该设置。

Interpret ML-selected proteins as model-supported candidates rather than automatic biological conclusions. Stability-selection frequency and predictive performance should be considered together with biological context and independent validation.
机器学习筛选出的蛋白应理解为模型支持的候选特征，而不是自动等同于生物学结论。建议结合稳定性筛选频率、模型表现、生物学背景以及独立验证进行解释。

V2.1.0 ALGORITHMIC NOTES COMPARED WITH V2.0.3
V2.1.0 与 V2.0.3 的算法逻辑差异说明

The overall RF/L1/RF + L1 stability-selection framework remains the same: reproducible stratified resampling, training-only variance filtering, training-derived imputation, and a shared RF/L1 split plan in the combined workflow are retained.
RF/L1/RF + L1 的整体稳定性筛选框架保持不变：仍采用可复现的分层重采样、仅训练集方差筛选、由训练集确定的缺失值填补，以及联合模块中的 RF/L1 共享拆分。

Compared with V2.0.3, V2.1.0 makes the following targeted methodological refinements:
- Binary L1 analyses use the binomial family.
- Multiclass L1 analyses explicitly use multinomial modeling with `type.multinomial = "ungrouped"`.
- L1 cross-validation fold handling is unified across standalone and stability-selection paths; Auto uses up to 5 folds according to the current minimum class size, while explicit requests are constrained to valid class sizes.
- RF top-20 and top-50 stability frequencies are calculated independently of the final requested Top-N output.
- Mean Gini scaling is calculated over the full stability-eligible protein universe, with zero contribution in repeats where a protein is absent from the relevant RF subset.
- Downstream ML feature-source availability is based on current-session final analysis results rather than historical output files.

与 V2.0.3 相比，V2.1.0 进行了以下针对性方法学调整：
- 两分类 L1 使用 binomial family。
- 多分类 L1 明确使用 multinomial，并设置 `type.multinomial = "ungrouped"`。
- 单独 L1 和稳定性筛选路径采用统一的交叉验证折数处理；Auto 根据当前最小类别样本量使用最多 5 折，手动指定值也会被限制在当前样本量可支持的有效范围内。
- RF top-20 和 top-50 稳定性频率独立于最终请求输出的 Top-N 数量计算。
- Mean Gini 在完整 stability-eligible 蛋白集合上进行标准化；某蛋白未进入某次重复的相关 RF 子集时，该次贡献按 0 计算。
- 下游 ML feature source 的可用性依据当前 session 的最终分析结果，而不再依据历史输出文件。
