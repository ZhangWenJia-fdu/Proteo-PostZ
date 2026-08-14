ProteoPostZ V2.0.3
Machine-learning parameter guide / 机器学习参数说明

This guide covers the three protein-level feature-selection modules: Random forest, L1, and RF + L1 combined.
本说明适用于三个蛋白质水平特征筛选模块：单独随机森林、单独 L1，以及 RF + L1 联合特征筛选。

GENERAL SETTINGS / 通用设置

Random seed / 随机种子
Sets the reproducible random-number stream for splitting and model fitting. The same seed with the same input and settings is intended to reproduce the same result; changing it produces another valid random resampling sequence and may change selected features.
用于控制拆分和模型拟合的随机数序列。在输入文件和其他参数相同的情况下，相同种子用于复现相同结果；更换种子会产生另一套有效的随机重复拆分，筛选结果可能改变。

Train/test split mode / 训练集/测试集拆分模式
Cross-validation only uses cross-validation for model assessment. Train/test split creates stratified training and test subsets for performance evaluation and stability repeats.
Cross-validation only 只使用交叉验证评估模型。Train/test split 会按组分层划分训练集和测试集，并用于模型表现评估及稳定性重复筛选。

Training set proportion / 训练集比例
The fraction of samples used for training in each stratified split. The remaining samples are the test set. The default is 0.7; entering 0.8 means 80% training and 20% testing. In RF + L1, RF and L1 use the same split in each repeat.
每次分层拆分中用于训练的样本比例，剩余样本作为测试集。默认值为 0.7；输入 0.8 即表示 80% 训练、20% 测试。RF + L1 模块中，每次重复的 RF 和 L1 使用同一套拆分。

Allow small-sample exploratory ML / 允许小样本探索性机器学习
The strict mode applies minimum group-size checks. Enable this option only when the dataset is too small for the strict requirements; results are exploratory and less reliable.
严格模式会检查每组最小样本数。只有在数据量不足以满足严格要求时才建议勾选；此时结果仅用于探索，可靠性较低。

Stability selection repeats / 稳定性筛选次数
Repeats the complete selection procedure. Each repeat makes a new stratified split, fits the model on the training data, and records selected features and performance. A larger value generally gives a more stable frequency estimate but takes longer.
重复执行完整的特征筛选过程。每次重复都会重新分层拆分，在训练集上拟合模型，并记录特征是否入选及模型表现。次数越多，筛选频率通常越稳定，但运行时间也越长。

Stability prefilter: top variance features / 稳定性预筛选：Top variance features
Before each repeat, features are ranked by variance using the training set only. Only the specified number is passed to the machine-learning model. This reduces noise and computation; it is not the final Top features output.
每次重复前，仅使用该次训练集计算特征方差并排序，只有指定数量的特征进入机器学习模型。该参数用于减少噪声和计算量，不等于最终输出的 Top features 数量。

Top features / Top features 数
The final number of selected features exported and displayed by the module. It is applied after the stability score or model score ranking.
最终展示和导出的特征数量，在稳定性评分或模型评分排序之后应用。

Random forest tree / Random forest tree 数
The number of trees in the random forest. More trees usually make the importance estimate less noisy but increase runtime. 500 is a reasonable general-purpose value.
随机森林包含的树数量。树越多通常重要性估计越平滑，但运行时间越长。500 是适合一般分析的设置。

Random forest mtry / Random forest mtry
The number of candidate features randomly considered at each tree split. Auto lets the software choose a data-dependent value. A smaller value increases tree diversity; a larger value may improve individual-tree strength but can make trees more similar. Auto is recommended unless a specific comparison is needed.
每棵树每次分裂时随机抽取并考虑的候选特征数量。Auto 会根据数据自动选择。较小的 mtry 会增加树之间的差异；较大的 mtry 可能增强单棵树，但也会使树更相似。除非需要进行特定比较，建议使用 Auto。

PDF width / PDF width (pt) and PDF height / PDF height (pt)
Control the exported vector-PDF dimensions in points. They change the figure layout, not the machine-learning result.
控制导出矢量 PDF 的宽度和高度，单位为 pt。它们只影响图形版式，不改变机器学习结果。

RANDOM FOREST ONLY / 单独 RANDOM FOREST 模块

Stability weight: RF top 20 frequency / RF top 20 频率权重
Weight of how often a feature appears among the top 20 RF features across repeats.
特征在每次重复的 RF 前 20 名中出现频率的权重。

Stability weight: RF top 50 frequency / RF top 50 频率权重
Weight of how often a feature appears among the top 50 RF features across repeats.
特征在每次重复的 RF 前 50 名中出现频率的权重。

Stability weight: mean Gini / 平均 Gini 权重
Weight of the feature's mean decrease in Gini importance across all stability repeats. The mean is calculated over the total number of repeats, so a feature absent from a repeat contributes zero for that repeat.
特征在所有稳定性重复中的平均 Gini 重要性权重。平均值按总重复次数计算；某特征在某次重复未入选时，该次对其贡献按 0 计。

The three RF weights determine the RF stability score. They are normally set to sum to 1. The final RF feature list is controlled by Top features.
三个 RF 权重共同决定 RF 稳定性评分，通常应使三者之和为 1。最终输出数量由 Top features 控制。

L1 ONLY / 单独 L1 模块

L1 alpha / L1 alpha
Controls the balance between LASSO and elastic net. 1 means pure LASSO; values between 0 and 1 add an elastic-net ridge component.
控制 LASSO 和 elastic net 的比例。1 表示纯 LASSO；0 到 1 之间的值表示加入 elastic-net 的 ridge 成分。

Lambda selection / Lambda selection
lambda.1se chooses the more regularized, simpler model within one standard error of the best cross-validated result. lambda.min chooses the lambda with the minimum cross-validated error and generally retains more features.
lambda.1se 在最佳交叉验证结果一个标准误范围内选择正则化更强、通常更简单的模型。lambda.min 选择交叉验证误差最小的 lambda，通常会保留更多特征。

Cross-validation folds / 交叉验证折数
Controls the number of folds used to select lambda. Auto lets the software choose a valid value; a manual value such as 5 is useful when reproducing a specific analysis. More folds can use data more efficiently but take longer.
控制选择 lambda 时的交叉验证折数。Auto 会自动选择有效值；如果要复现特定分析，可以手动输入 5 等固定值。折数增加通常能更充分利用数据，但运行时间更长。

L1 selection records how often each feature is selected across repeats. The final L1 feature list is controlled by Top features.
L1 模块记录每个特征在多次重复中被选中的频率。最终输出数量由 Top features 控制。

RF + L1 COMBINED / RF + L1 联合模块

The RF + L1 score combines four normalized stability measures: RF top 20 frequency, RF top 50 frequency, L1 selection frequency, and mean Gini. The four weights should normally sum to 1. Larger weight means that component has greater influence on the final ranking.
RF + L1 评分由四项标准化稳定性指标组成：RF 前 20 频率、RF 前 50 频率、L1 筛选频率和平均 Gini。四个权重通常应加和为 1。某项权重越大，该项对最终排序的影响越大。

RF top 20 frequency weight / RF 前 20 频率权重
Controls the contribution of RF top-20 occurrence frequency.
控制 RF 前 20 出现频率的贡献。

RF top 50 frequency weight / RF 前 50 频率权重
Controls the contribution of RF top-50 occurrence frequency.
控制 RF 前 50 出现频率的贡献。

L1 selection frequency weight / L1 筛选频率权重
Controls the contribution of the proportion of repeats in which L1 selects the feature.
控制 L1 在多少比例的重复中选中特征这一指标的贡献。

Mean Gini weight / 平均 Gini 权重
Controls the contribution of mean RF Gini importance across all repeats.
控制所有重复中平均 RF Gini 重要性的贡献。

The combined module uses one shared split plan per stability repeat, then fits RF and L1 on the corresponding training set. Its final feature list is the union-ranking result from the combined score and is limited by Top combined features.
联合模块在每次稳定性重复中使用一套共享拆分，然后分别在对应训练集上拟合 RF 和 L1。最终特征列表按联合评分排序，并由 Top combined features 限制数量。

PRACTICAL RECOMMENDATIONS / 实用建议

For ordinary datasets, start with 50 stability repeats, training proportion 0.7, variance prefilter 200, RF trees 500, RF mtry Auto, L1 alpha 1, lambda.1se, and Top features 50. For a deliberate 80/20 analysis, set Training set proportion to 0.8 in the selected module. Keep the random seed fixed when comparing parameter changes, and change only one parameter at a time when investigating its effect.
一般数据可从以下设置开始：稳定性筛选 50 次、训练集比例 0.7、方差预筛选 200、RF 树数 500、RF mtry 使用 Auto、L1 alpha 为 1、lambda.1se、Top features 为 50。如果需要明确进行 80/20 分析，在相应模块将 Training set proportion 改为 0.8。比较参数影响时应固定随机种子，并尽量一次只改变一个参数。
