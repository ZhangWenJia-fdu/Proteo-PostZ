ProteoPostZ Formal V2.0.2
Developed by Wenjia Zhang

一、启动方式
1. 推荐：双击 ProteoPostZ_v2.0.2.exe。
2. 备用：运行 Run_ProteoPostZ_v2.0.2.cmd。
3. 本地浏览器地址：http://127.0.0.1:3840/

exe 是启动本地 R/Shiny 应用的启动器，不是在线网站。软件使用随包提供的 portable R 和 R 包，浏览器界面只在当前电脑本地可访问。

二、支持的输入
- DIA-NN 和 Spectronaut 蛋白水平 DIA 结果。
- FragPipe/MSFragger、PEAKS 和 MaxQuant 蛋白水平 DDA 结果。
- 用户整理的标准定量矩阵。

PEAKS 用户只需上传一个 protein 文件。软件根据内部列结构自动识别 DB 或 LFQ，不根据文件名判断。DB 使用 Area <sample> 和 sample-specific Coverage(%)；LFQ 使用 <sample> Area，并排除 Group N Area、Ratio/Profile 等非样品定量字段。如果 DB 与 LFQ 字段混合，优先按 DB 读取。LFQ 没有样品级独立鉴定证据，因此 Identified_Protein_Count 不可用；Quantified_Protein_Count 仍然正常生成。

三、主要功能
样品计数、Venn/UpSet、理化性质、相关性热图、丰度排序、CV、PCA、UMAP、t-SNE、火山图、表达热图、RF、L1、RF + L1、特征蛋白 UMAP/热图和 Slingshot 伪时间分析。

软件支持 PDF 矢量图输出和各模块 CSV 导出。

四、离线使用
正常使用不需要网络。只有在用户明确更新注释或安装缺失 R 包时才需要网络。

五、软件包边界
本软件包不包含本地测试文件、测试记录、生成的输出、日志和开发历史。
