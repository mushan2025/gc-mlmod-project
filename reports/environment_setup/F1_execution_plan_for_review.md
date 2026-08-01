# F1 正式执行计划（SCTransform 主线，审核稿）

更新日期：2026-08-01
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

### F1.2 固定QC与doublet

脚本：`scripts/F1_single_cell/F1_02_qc_doublet.R`

- 每个样本先保留至少在3个细胞中检出的feature，再只计算一套`nCount_RNA`、`nFeature_RNA`、`mt_percent`和`HB_percent`。
- 唯一主规则：`nFeature_RNA >= 500`、`nFeature_RNA < 6000`、`nCount_RNA > 1000`、`mt_percent <= 20`、`HB_percent < 5`。
- scDblFinder按sample/capture运行并作为主删除依据；DoubletFinder按同一输入运行，只作敏感性标记，不机械取并集删除。
- 输出逐细胞决定、逐样本过滤影响和doublet摘要，保存`02_all_cells_qc_filtered.rds`。本节不运行DecontX，因为此时尚无经过marker审核的可靠粗谱系标签。
- 同时在逐细胞审计表中输出预注册的`fixed_qc_pass_no_ncount_sensitivity`：只去掉`nCount_RNA > 1000`，其余四条QC规则完全不变。该列不改变主对象，也不复用主mask上的doublet结果冒充平行分析。
- 正式敏感性臂由该mask重跑其余相同步骤和参数，比较主要谱系及注释、恶性上皮群与06a纳入，以及进入F2后的MLMOD-high/low成员和主要方向。主分析始终使用含`nCount_RNA > 1000`的mask，不按哪套结果更好看选择；仅细胞数和少量边界细胞变化不算核心结论改变，结果并入现有QC/F1/F2报告，不增设gate。

### F1.3 SCTransform、Harmony与全细胞聚类

脚本：`scripts/F1_single_cell/F1_03_sct_harmony_cluster.R`

- 以相同QC后RNA raw counts为起点，按`sample_id`建立SCTransform v2模型，形成供大谱系注释的第一轮粗聚类。
- 固定`vars.to.regress=NULL`，不回归`mt_percent`、`nCount_RNA`或细胞周期；HVG数为3,000。
- 在SCT assay上运行50个PC，在PCA上仅按`sample_id`运行Harmony。
- 主dims为1:30；Leiden预生成0.2、0.4、0.6、0.8、1.0供marker审核，执行前默认0.6。
- 保存未整合SCT-PCA UMAP和Harmony UMAP。RNA assay另建LogNormalize `data`层，仅用于marker和表达图。
- 输出`03a_all_cells_sct_harmony_clustered.rds`、参数表和整合诊断。

### F1.4 主要细胞类型注释与DecontX评估

脚本：`scripts/F1_single_cell/F1_04_annotation.R`

- 用RNA LogNormalize `data`层计算cluster marker，并结合冻结marker panel制作DotPlot/FeaturePlot。
- 第一次运行生成`F1_cluster_annotation_template.tsv`后停止；研究者审核并形成`F1_cluster_annotation_approved.tsv`。
- 再次运行时先校验每个cluster均有唯一批准标签并写入major/minor/state/confidence；这些粗谱系标签在DecontX前冻结，不由校正结果反向迭代。
- 随后在每个sample/capture内，对固定QC且scDblFinder为singlet的raw integer counts运行一次DecontX，以批准的`cell_type_major`作为`z`。不使用精细上皮亚型、恶性标签、MLMOD或预后信息。
- 只有样本内至少存在两个可靠粗谱系时才把DecontX记为可评价；否则记为`not_evaluable_fewer_than_two_reliable_lineages`。由于无raw droplets，`background=NULL`，不能把结果解释为空液滴支持的污染真值。
- corrected counts逐样本独立保存，contamination score写入metadata；raw counts仍为主矩阵，污染分数不用于删细胞，也不自动改写注释。
- DecontX只运行一次。四个检查点均在F1.4产出之后：本次运行先建立注释后全细胞ambient基线；随后复核上皮候选和F1.5二次聚类；F1.6前复核CNV reference与观察细胞；F2评分/关键DE及F4通讯前复核相应对象。每处只看对象内contamination中位数/P90、异谱系marker泄漏、本谱系marker保留及raw/corrected目标结果是否实质改变。
- 输出`ambient_rna_cell_estimates.tsv`、`ambient_rna_summary_by_sample.tsv`和逐样本corrected矩阵，并保存`03_all_cells_integrated_annotated.rds`。这是必要的生物学判断点，不把自动打分当作最终细胞身份。

### F1.5 上皮提取与二次聚类

脚本：`scripts/F1_single_cell/F1_05_epithelial_recluster.R`

- 只提取高/中置信Epithelial，排除已批准的mixed/doublet或明显污染群。
- 从RNA raw counts重建上皮对象，再次按样本运行相同SCTransform v2 + Harmony流程。
- 生成上皮marker和污染复核模板，保存`04_epithelial_reclustered.rds`。
- F1.6启动前需要批准`F1_epithelial_cluster_review_approved.tsv`，决定哪些上皮cluster可进入CNV判定。

### F1.6 恶性上皮判定

脚本：`scripts/F1_single_cell/F1_06_malignancy_inference.R`

- inferCNV和CopyKAT主分析均使用RNA raw counts，按样本运行；SCT residual和Harmony坐标不作输入。DecontX corrected counts不作默认输入，只允许按下述条件用于inferCNV敏感性。
- inferCNV reference按固定顺序选择，且不混合不同层级：先用当前样本高/中置信T/NK，必要时在同层加入B/Plasma；不足50个时，改用`patient_id`完全一致的配对`Normal_Gastric`样本免疫细胞；仍不足时，使用其他患者`Normal_Gastric`样本的免疫细胞。配对只依据F0核对后的精确`patient_id`和`group_analysis`，不根据样本名猜测。reference最多500个；候选超过500个时在来源样本间轮流均衡抽取，三层均不足50个时该样本inferCNV记为`not_evaluable`，继续其他样本及CopyKAT。
- 运行前分别检查reference和上皮观察细胞的ambient风险，少数污染reference优先删除或替换。使用外部正常胃reference时保留每个reference细胞的来源样本、患者、谱系和置信度；这属于同一数据集内的条件性基线，不消除跨样本技术差异。
- inferCNV显式使用`analysis_mode="subclusters"`、`tumor_subcluster_partition_method="leiden"`、`k_nn=20`、`leiden_resolution="auto"`、`leiden_method="PCA"`、`leiden_function="CPM"`、`inspect_subclusters=TRUE`和`cluster_by_groups=TRUE`。每个sample的全部候选上皮只设一个`observations`组；原上皮cluster仅保留为cell-level注释与subcluster组成背景。快速分组检查可用`HMM=FALSE`；正式运行固定`HMM=TRUE`、`HMM_type="i6"`、`HMM_report_by="subcluster"`、`BayesMaxPNormal=0.5`和`reassignCNVs=TRUE`。
- inferCNV subcluster只是帮助从热图中定位具有相似CNV模式的细胞，不称为真实肿瘤亚克隆，也不预设目标数量。少数小组记为不可单独评价；只有广泛过度拆分且组间没有可区分大片段模式时才比较一个相邻更低分辨率，不进行参数扫描或按恶性细胞产量选择。正式恶性证据以热图及同谱系Normal_Gastric比较为主，HMM和burden只作辅助。
- CopyKAT每次仅输入当前样本的全部QC后singlet，不混入其他样本，也不只截取候选上皮。A主臂不传`norm.cell.names`；若预测diploid群中候选上皮超过50%，保留原始call但标记`not_evaluable_baseline_suspect`。B敏感性臂只使用同样本高/中置信T/NK（必要时加B/Plasma，至少50个）作known normal。C敏感性臂只在Normal_Gastric中用同样本免疫细胞加一半正常胃上皮作known normal、评价另一半，两折交换并固定seed=42；不跨样本借reference。
- 仅当审核表确认污染可信且可能改变CNV结论时运行一次corrected inferCNV敏感性；脚本强制reference与观察细胞同时来自各自保存的DecontX corrected矩阵，禁止只校正一方。大片段支持等级、最终恶性标签或06a纳入任一改变时写出比较表并暂停裁决，否则raw结果保持主结果。CopyKAT不默认接受corrected输入。
- inferCNV辅助cell burden定义为每个细胞相对于reference逐基因中心值的平均绝对偏差，并记录reference细胞burden的P95；该数值只帮助定位，不自动判恶性。最终inferCNV支持等级由sample × inferCNV subcluster热图中的连续大片段模式人工审核，原`epithelial_cluster_id`以组成比例形式保留为注释背景。
- 每个样本保留完整inferCNV输出目录、最终RDS、写出的热图表达矩阵、输入annotation、gene order、最终热图细胞顺序、细胞—subcluster对应表和作图来源manifest；最终RDS中的`expr.data`作为后续发表级热图重绘的数值来源，不依赖默认图片反推数据。服务器结果用`package_F1_06_results.sh`压成`tar.zst`并生成大写SHA256后再下载。
- CopyKAT的aneuploid支持恶性，diploid不能排除近二倍体肿瘤，uncalled按未知。
- 两种方法都来自RNA表达，只是互补方法稳健性，不是独立DNA证据。
- 结论边界固定为：CNV定义可能低估近二倍体恶性细胞，包括潜在的基因组稳定型胃癌细胞；后续MLMOD结论适用于本流程可由CNV证据识别的恶性上皮，不能排除未纳入的近二倍体恶性亚群。该限制不新增恶性类别，也不改变06a/06b组成。
- 先生成`F1_malignancy_cluster_review_template.tsv`供人工审阅。`epithelial_subtype`使用marker panel中的固定名称：`Pit_mucous_epithelial`、`Mucous_neck_epithelial`、`Chief_epithelial`、`Parietal_epithelial`、`Enteroendocrine_epithelial`或`Intestinal_like_epithelial`；暂时不能归类时填`epithelial_uncertain`。模板列出父级上皮cluster marker、候选sample × inferCNV subcluster与同谱系`Normal_Gastric`的marker检出比例、正常参照细胞数及样本数、组织/样本构成和DecontX摘要；不能映射或没有同谱系正常参照时明确记为不可评估，不自动按阴性处理。
- 谱系marker只说明细胞像哪类胃上皮，不能直接证明非恶性。`normal_program_support`要求整体接近同谱系`Normal_Gastric`；`tumor_program_support`要求相对同谱系正常细胞出现成套异常变化，不采用单个“癌marker”硬判。增殖、高MT、缺氧、应激、凋亡、intestinal-like状态、ambient或肿瘤组织来源均不能单独形成肿瘤样支持；也不使用MLMOD或预后信息。
- 两个程序字段允许同时为TRUE或同时为FALSE；冲突或均不明确时保持`epithelial_uncertain`，不为追求纳入量强行裁决。
- 标签实现与主线一致：`malignant_probable_infercnv`和`malignant_probable_copykat`均要求无强正常样程序；`non_malignant_epithelial`要求无肿瘤样程序。`malignant_probable_copykat`只保留在`05`和探索性输出；`06a`只纳入至少有weak inferCNV大片段支持的`malignant_high_confidence`与`malignant_probable_infercnv`，`06b`只纳入前者。

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

只有以下情况暂停请用户决定：F0输入事实不一致、固定QC造成无法解释的样本极端损失、主doublet算法失败、SCTransform/Harmony明显破坏主要谱系结构、主要谱系无法可靠注释、CNV reference在多数或关键样本中不可用、恶性证据大范围冲突，或可信ambient污染使目标科学结论在raw/corrected间发生实质变化。单一样本reference不足只记为不可评估，不阻断其他样本。普通格式和可恢复运行问题修复后继续，不新增多层gate。
