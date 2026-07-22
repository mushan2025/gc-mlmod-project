# F1 正式执行计划（SCTransform 主线，审核稿）

更新日期：2026-07-22
当前状态：F0准入、F1必需依赖、总入口dry-run和脚本静态检查均已通过；**尚未执行真实数据，正式设备与资源pilot仍待批准**。

## 1. F1 要回答什么

F1只完成四件与后续研究直接相关的工作：

1. 从F0确认过的40个公开raw count矩阵建立可追溯的Seurat对象。
2. 按冻结规则去除明显低质量细胞和scDblFinder主判双细胞，并评估retained-cell ambient RNA风险。
3. 用SCTransform v2、PCA、Harmony和Leiden建立全细胞图谱，完成主要谱系及上皮细胞注释。
4. 在上皮细胞内结合inferCNV、CopyKAT、marker和组织来源，形成F2使用的恶性上皮主对象。

F1不计算MLMOD或UCell分数，不使用预后信息，也不根据后续结果改变QC、聚类或恶性判定。

## 2. 正式启动条件

总入口在真正运行前必须同时确认：

- `data/metadata/processed_input_manifest.tsv`、`sample_info.tsv`和`data_audit.tsv`来自正式F0执行；
- `results/F0_audit/F0_gate_checklist.tsv`没有blocking FAIL；
- `F0_execution_report.md`中的`F0_scRNA_F1_gate`为`PASS`或`PASS_WITH_NOTED_LIMITATIONS`；
- `sample_info.tsv`与`data_audit.tsv`只纳入`include_in_f1=true`且`audit_decision=enter_full_F1_independent_reQC`的样本；
- `environment/F1/required_packages.tsv`中的本阶段必需包均可用；
- 用户已批准本次运行设备和范围。

F0已正式执行并得到`PASS_WITH_NOTED_LIMITATIONS`，不存在blocking FAIL。2026-07-22已在固定R 4.4.3/Bioconductor 3.20环境补齐全部F1必需依赖：`scDblFinder 1.20.2`、`celda 1.22.1`、`DoubletFinder 2.0.6`、`glmGamPoi 1.18.0`和`leidenbase 0.1.36`；逐包`loadNamespace()`、关键函数检查、总入口只读模式和静态契约均通过。`DoubletFinder`固定记录Git commit `1B244D8F0D54B4B1CB4365639931BBB16F01E1CD`。

本次环境补齐采用以下等价安装路线；版本以实际加载结果为准，不以安装命令返回代替验证：

```r
BiocManager::install(c("scDblFinder", "celda", "glmGamPoi"), ask = FALSE, update = FALSE)
install.packages(c("leidenbase", "fields"), type = "binary")
# DoubletFinder按上方记录的commit浅克隆后使用R CMD INSTALL安装。
```

启动前实测可用内存仅3.89 GB，低于F1.1预计8-20 GB及后续阶段需求，因此当前只判定“软件环境就绪”，不判定“本机可立即正式运行”。正式执行前需释放内存并重测，再对最大样本做资源pilot；正式运行时仍将实际版本和sessionInfo写入`F1_parameters_and_versions.tsv`。

## 3. 六个连续步骤

### F1.1 导入与初始对象

脚本：`scripts/F1_single_cell/F1_01_import.R`

- 按manifest逐样本读取gene × cell的`.csv.gz`矩阵，不自动转置。
- 实测行列数、方向、gene和barcode唯一性必须与F0一致；不一致即停止。
- 细胞名改为`sample_id__original_barcode`，同时保留原barcode、GEO号、patient和组织分组。
- 完整26,571行feature和全部公开细胞的RNA raw counts保存到`01_all_cells_raw_or_initial.rds`。

### F1.2 固定QC、doublet与ambient RNA

脚本：`scripts/F1_single_cell/F1_02_qc_doublet_ambient.R`

- 每个样本先保留至少在3个细胞中检出的feature，再只计算一套`nCount_RNA`、`nFeature_RNA`、`mt_percent`和`HB_percent`。
- 唯一主规则：`nFeature_RNA >= 500`、`nFeature_RNA < 6000`、`nCount_RNA > 1000`、`mt_percent <= 20`、`HB_percent < 5`。
- scDblFinder按sample/capture运行并作为主删除依据；DoubletFinder按同一输入运行，只作敏感性标记，不机械取并集删除。
- 删除scDblFinder阳性后按样本运行DecontX。corrected counts逐样本独立保存，不自动合并为主对象assay；raw counts仍是主矩阵，污染分数不用于删细胞。
- DecontX只运行一次，后续复用其score和corrected矩阵，在四处检查：F1保留细胞基线、主要谱系/上皮候选、CNV reference与上皮观察细胞、F2评分/关键DE及F4通讯对象。
- 每处只看四类信息：对象内contamination中位数/P90、异谱系marker泄漏变化、本谱系marker保留、raw/corrected目标结果是否实质改变。污染分数无跨样本统一硬阈值，也不单独触发删除或校正。
- 输出逐细胞决定、逐样本过滤影响、doublet与ambient摘要，并建立唯一跨阶段汇总表`ambient_decision_summary.tsv`；保存`02_all_cells_qc_filtered.rds`。

### F1.3 SCTransform、Harmony与全细胞聚类

脚本：`scripts/F1_single_cell/F1_03_sct_harmony_cluster.R`

- 在相同QC后细胞上按`sample_id`建立SCTransform v2模型。
- 固定`vars.to.regress=NULL`，不回归`mt_percent`、`nCount_RNA`或细胞周期；HVG数为3,000。
- 在SCT assay上运行50个PC，在PCA上仅按`sample_id`运行Harmony。
- 主dims为1:30；Leiden预生成0.2、0.4、0.6、0.8、1.0供marker审核，执行前默认0.6。
- 保存未整合SCT-PCA UMAP和Harmony UMAP。RNA assay另建LogNormalize `data`层，仅用于marker和表达图。
- 输出`03a_all_cells_sct_harmony_clustered.rds`、参数表和整合诊断。

### F1.4 主要细胞类型注释

脚本：`scripts/F1_single_cell/F1_04_annotation.R`

- 用RNA LogNormalize `data`层计算cluster marker，并结合冻结marker panel制作DotPlot/FeaturePlot。
- 第一次运行生成`F1_cluster_annotation_template.tsv`后停止；研究者审核并形成`F1_cluster_annotation_approved.tsv`。
- 再次运行时校验每个cluster均有唯一批准标签，才写入major/minor/state/confidence并保存`03_all_cells_integrated_annotated.rds`。
- 对marker冲突、稀有、边界cluster和上皮候选复核raw/corrected谱系证据；主要谱系标签一致率低于95%只触发人工复核，身份改变或主要marker依据消失才视为实质变化。
- 这是必要的生物学判断点，不把自动打分当作最终细胞身份。

### F1.5 上皮提取与二次聚类

脚本：`scripts/F1_single_cell/F1_05_epithelial_recluster.R`

- 只提取高/中置信Epithelial，排除已批准的mixed/doublet或明显污染群。
- 从RNA raw counts重建上皮对象，再次按样本运行相同SCTransform v2 + Harmony流程。
- 生成上皮marker和污染复核模板，保存`04_epithelial_reclustered.rds`。
- F1.6启动前需要批准`F1_epithelial_cluster_review_approved.tsv`，决定哪些上皮cluster可进入CNV判定。

### F1.6 恶性上皮判定

脚本：`scripts/F1_single_cell/F1_06_malignancy_inference.R`

- inferCNV和CopyKAT主分析均使用RNA raw counts，按样本运行；SCT residual和Harmony坐标不作输入。DecontX corrected counts不作默认输入，只允许按下述条件用于inferCNV敏感性。
- inferCNV优先使用同一样本高置信T/NK，必要时加入B/Plasma；不足时才登记使用pooled同谱系reference。运行前检查reference和上皮观察细胞ambient风险，少数污染reference优先删除或替换，无合格reference则该样本记为不可评估。
- 仅当污染可信且可能改变CNV结论时运行一次corrected inferCNV敏感性，reference与观察细胞必须同时使用同一corrected矩阵；禁止只校正reference。大片段支持等级或06a纳入改变时暂停裁决，否则raw结果保持主结果。CopyKAT不默认接受corrected输入。
- inferCNV辅助cell burden定义为每个细胞相对于reference逐基因中心值的平均绝对偏差，并记录reference细胞burden的P95；该数值只帮助定位，不自动判恶性。最终inferCNV支持等级由sample × epithelial cluster热图中的连续大片段模式人工审核。
- CopyKAT的aneuploid支持恶性，diploid不能排除近二倍体肿瘤，uncalled按未知。
- 两种方法都来自RNA表达，只是互补方法稳健性，不是独立DNA证据。
- 先生成`F1_malignancy_cluster_review_template.tsv`供人工审阅；批准后按主线规则形成`05`、`06a`和`06b`对象。

## 4. 表达矩阵用途

| 数据层 | 允许的主要用途 | 禁止的替代用途 |
|---|---|---|
| RNA raw counts | QC、doublet、DecontX输入、inferCNV、CopyKAT、pseudobulk、F2 UCell | 不直接用于PCA/Harmony主聚类 |
| SCT `scale.data` | HVG、PCA、Harmony、UMAP、Leiden | 不用于UCell、inferCNV、CopyKAT或pseudobulk |
| RNA LogNormalize `data` | marker检验和表达展示 | 不作为第二条主聚类线 |
| DecontX corrected assay | ambient风险描述和获批敏感性 | 不替换raw counts，不用于DESeq2主分析 |

## 5. 运行与恢复

仅查看计划和缺失条件，不执行：

```powershell
& 'D:\software\R\R443\R-4.4.3\bin\x64\Rscript.exe' scripts\F1_single_cell\run_F1.R
```

F0、依赖、环境和审核文件均就绪后，才允许显式执行：

```powershell
& 'D:\software\R\R443\R-4.4.3\bin\x64\Rscript.exe' scripts\F1_single_cell\run_F1.R --execute
```

各阶段分别保存RDS和TSV，失败后从最近完成阶段恢复，不需要重跑已经确认的上游步骤。总入口不会自动安装软件包，也不会自动把计算切到服务器。

## 6. 资源与停止条件

SCTransform全细胞步骤和inferCNV/CopyKAT最可能成为内存瓶颈。正式执行前先以最大样本或批准的代表样本做资源pilot；出现明显换页、峰值内存超过设备安全范围或运行失败时，只切换相应阶段到台式机/服务器。

只有以下情况暂停请用户决定：F0输入事实不一致、固定QC造成无法解释的样本极端损失、主doublet算法失败、SCTransform/Harmony明显破坏主要谱系结构、主要谱系无法可靠注释、CNV reference不可用、恶性证据大范围冲突，或可信ambient污染使目标科学结论在raw/corrected间发生实质变化。普通格式和可恢复运行问题修复后继续，不新增多层gate。
