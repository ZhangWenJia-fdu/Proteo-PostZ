ProteoPostZ Formal V2.0.2
Developed by Wenjia Zhang

一、启动方式
1. 推荐：双击 ProteoPostZ_v2.0.2.exe。
2. 备用：运行 Run_ProteoPostZ_v2.0.2.cmd。
3. 本地浏览器地址：http://127.0.0.1:3840/

exe 是启动本地 R/Shiny 应用的启动器，不是在线网站。软件使用随包提供的 portable R 和 R 包，浏览器界面只在当前电脑本地可访问。

二、V2.0.2 新增内容
- 新增对 PEAKS LFQ protein result 的输入支持。
- PEAKS 用户只需上传一个 protein 文件；软件根据内部列 schema 自动识别 PEAKS DB protein result 或 PEAKS LFQ protein result，不依赖文件名。
- 继续兼容 PEAKS DB；如果同一文件混合 DB 与 LFQ 风格字段，则按 PEAKS DB 读取。
- 同时优化了 Sample correlation heatmap、Expression heatmap 和 Feature-protein heatmap 的分析设置，包括排序或聚类选择、注释颜色、行标签及热图输出布局。

三、支持的输入
- DIA-NN 和 Spectronaut 蛋白水平 DIA 结果。
- FragPipe/MSFragger、PEAKS 和 MaxQuant 蛋白水平 DDA 结果。
- 用户整理的标准定量矩阵。

PEAKS 用户只需上传一个 protein 文件。软件根据内部列结构自动识别 DB 或 LFQ，不根据文件名判断。DB 使用 Area <sample> 和 sample-specific Coverage(%)；LFQ 使用 <sample> Area，并排除 Group N Area、Ratio/Profile 等非样品定量字段。如果 DB 与 LFQ 字段混合，优先按 DB 读取。LFQ 没有样品级独立鉴定证据，因此 Identified_Protein_Count 不可用；Quantified_Protein_Count 仍然正常生成。

四、主要功能
定性分析包括样品蛋白计数、Venn 图、UpSet 图和理化性质分布。组水平蛋白或特征集合可通过 Minimum replicates detected in group 设置进行控制。

定量分析包括样品相关性热图、丰度排序图、组内 CV ridgeline、PCA、UMAP、t-SNE、火山图、表达热图、随机森林特征选择、L1 特征选择、RF + L1 联合特征选择、特征蛋白 UMAP、特征蛋白热图和 Slingshot 伪时间分析。软件还提供样品分组、输入预览、样品计数汇总和输出文件浏览。

各分析模块均可根据分析目的调整适用的阈值、聚类、维度、统计或机器学习设置、配色、标签和 PDF 尺寸。软件支持各模块的 PDF 矢量图和对应 CSV 数据导出。分析卡片可全屏放大；定量分析页面采用固定视口，模块内的控制区和结果区可独立滚动。

五、离线使用
正常使用不需要网络。只有在用户明确更新注释或安装缺失 R 包时才需要网络。
