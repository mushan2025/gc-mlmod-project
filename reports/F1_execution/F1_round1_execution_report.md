# F1 第一轮执行报告：F1.1-F1.5 已完成

更新日期：2026-07-30

## 1. 执行范围与当前停止位置

本轮已完成：

1. F1.1：40个公开样本矩阵导入、方向和整数原始计数检查。
2. F1.2：冻结QC、scDblFinder主双细胞删除和DoubletFinder敏感性分析。
3. F1.3：逐样本SCTransform v2、SCT PCA、Harmony、UMAP和Leiden聚类。
4. F1.4：24个cluster的研究者注释批准、按样本运行DecontX、ambient RNA影响复核和结果回收。
5. F1.5：上皮候选提取、逐样本SCTransform、Harmony、Leiden二次聚类、输出验证和ambient复核。

19个上皮cluster已于2026-07-30完成研究者审核。当前允许进入F1.6，但尚未运行inferCNV/CopyKAT，也未生成恶性上皮对象。

## 2. 数据导入与QC结果

- 公开输入：40个样本、158,641个细胞，归档矩阵共有26,571个feature。
- 矩阵方向、非负整数值、样本数和细胞数检查均通过。
- 固定QC排除5,321个细胞，153,320个细胞通过固定QC。
- scDblFinder判定并排除10,670个双细胞。
- 最终保留142,650个细胞，占公开细胞的89.92%。
- 去掉本项目经验性`nCount_RNA > 1000`条件时，会额外保留2,871个边界细胞；该结果只记录为预注册敏感性范围，不改变主对象。

固定QC严格使用已批准规则：每个样本先保留至少在3个细胞中检出的基因，再要求`nFeature_RNA >= 500`、`nFeature_RNA < 6000`、`mt_percent <= 20`、`HB_percent < 5`和`nCount_RNA > 1000`。

## 3. 双细胞结果

scDblFinder是主方法，DoubletFinder只用于方法敏感性检查。40个样本中39个完成DoubletFinder；sample37在固定QC后仅172个细胞，按预登记规则记为`not_evaluable_small_sample`。

在可评价细胞中，scDblFinder与DoubletFinder的Jaccard一致度为0.130。最终保留细胞中有8,209个被DoubletFinder单独标记；这些细胞没有自动删除，而是在注释时作为风险信息复核。cluster 22、15、8和4的DoubletFinder敏感性阳性比例相对较高，但其谱系身份结合marker、跨样本覆盖和内部一致性后仍可解释，因此没有追加删除。

## 4. SCTransform、Harmony与聚类

- 主归一化：逐样本SCTransform v2，`vars.to.regress = NULL`。
- 候选高变基因：3,000个；实际进入PCA的共同已缩放基因：2,141个。
- PCA：50维；Harmony仅按`sample_id`校正；主聚类使用前30个Harmony维度。
- 主Leiden分辨率：0.6，共24个cluster。
- 分辨率0.2、0.4、0.8和1.0分别得到15、20、27和29个cluster，主分辨率处于可解释的中间尺度。
- 03a对象回读检查全部通过：142,650个细胞顺序不变，RNA raw counts指纹一致，40个SCT模型、PCA、Harmony和两个UMAP均完整。
- 未使用MLMOD、预后或生存信息参与聚类和注释。

PCA与`nCount_RNA`、`nFeature_RNA`、`mt_percent`的最大绝对相关分别约为0.476、0.377和0.386。这说明部分主成分仍含细胞复杂度或线粒体状态信息，但marker和主要谱系结构清楚，不支持反向修改冻结QC或回归这些变量。

## 5. 研究者批准的24簇注释

24个cluster均有可用marker，并于2026-07-30由研究者批准。142,650个细胞被归入8个粗谱系：

- T/NK：50,064个。
- B/Plasma：29,259个。
- Epithelial：26,828个。
- Myeloid：12,385个。
- Endothelial/Pericyte：11,192个。
- Fibroblast/CAF：8,946个。
- Mast：2,703个。
- Mesothelial：1,273个。

上皮候选为c3、c9、c16、c19和c21。边界簇按“谱系与状态分开”处理：c18保留B/Plasma母谱系并记录stress-high状态；c22保留Mesothelial身份；c24标为Cycling T；c3保留为混合胃型/肠型上皮；c19标为chief-like epithelial并记录metallothionein-high状态。恶性身份留待F1.6的CNV证据判定。

## 6. DecontX分组修正

第一次DecontX技术上完成了40个样本，但将8个粗谱系直接作为模型内的`z`分组。该做法把内部差异很大的B细胞和浆细胞等合并，导致c12常规B细胞的污染中位数约为0.986，明显不符合marker和表达结构。

随后在sample14、sample15和sample34进行小规模对照。把`z`改为24个研究者审核后的Seurat cluster后，c12污染中位数分别由0.9864、0.9794和0.9883降至0.0250、0.0389和0.0408。由此确认问题来自模型分组过粗，而不是c12本身几乎全部由ambient RNA构成。

最终规则为：

- 每个样本独立运行DecontX。
- 输入为QC和scDblFinder过滤后保留细胞的RNA raw integer counts。
- `z`使用研究者审核后的24个Seurat cluster，帮助模型区分同一样本内不同表达群体。
- 8个粗谱系只用于汇总和解释，不再作为模型分组。
- 第一次粗谱系结果降级为无效诊断，完整归档但不参与后续决策。

## 7. 最终ambient RNA结果

最终DecontX成功完成40/40个样本，142,650个细胞均获得有限污染估计；40个corrected矩阵均通过非负值、细胞顺序、总量一致性和原始对象保护检查。

全体细胞污染比例：

- 均值：0.0760。
- 中位数：0.00984。
- P90：0.205。
- P95：0.505。
- 污染估计不低于0.10的细胞：15.0%。
- 污染估计不低于0.20的细胞：10.2%。

不同谱系的中位数/P90分别为：

- T/NK：0.0105 / 0.0621。
- B/Plasma：0.0300 / 0.5868。
- Epithelial：0.0045 / 0.3253。
- Myeloid：0.0076 / 0.2190。
- Endothelial/Pericyte：0.0081 / 0.1378。
- Fibroblast/CAF：0.0081 / 0.1739。
- Mast：0.0038 / 0.0336。
- Mesothelial：0.0071 / 0.0620。

污染估计具有明显的样本和cluster异质性。需要后续重点留意sample35、sample28、sample22、sample36和sample9；其中sample22只有385个保留细胞，小簇结果尤其不稳定。上皮cluster的中位污染估计均较低，但部分样本或细胞存在高污染长尾，因此不能只看全局中位数。

## 8. raw与corrected结果是否发生实质变化

使用冻结的粗谱系marker panel和每个cluster的前30个非伪影marker，对raw与DecontX corrected结果进行了对照：

- 24/24个cluster的粗谱系panel最高得分方向在校正前后不变。
- 22/24个cluster的预期谱系panel在raw和corrected中均为最高。c15和c22存在预先已知的谱系panel重叠，但校正前后均不变，不能归因于ambient RNA。
- 所有cluster的前列marker仍可检出。
- 前列marker计数保留率最低的是c18，为0.891。
- 五个上皮cluster中最低的是c16，为0.985，其余约为0.989至0.999。
- c12的谱系panel计数保留率为0.668，但前列marker计数保留率为0.996，说明校正主要移除了浆细胞或非本簇来源信号，并未破坏其常规B细胞身份。

因此，ambient校正没有改变24簇粗注释，也没有改变c3、c9、c16、c19和c21作为上皮候选的判断。局部高污染样本仍需在后续关键结论处做raw/corrected针对性敏感性检查，但不支持全局替换主矩阵或删除细胞。

## 9. F1.4决策

F1.4结论为`continue_with_raw_main_and_targeted_corrected_sensitivity`：

- RNA raw counts继续作为主矩阵。
- DecontX corrected矩阵仅用于后续明确检查点的条件敏感性分析。
- 不根据DecontX污染比例新增细胞删除阈值。
- corrected值含浮点数，不能作为DESeq2主分析的原始整数counts。
- 允许进入F1.5上皮细胞重聚类。

## 10. F1.5输入修正

F1.4共识别26,828个高/中置信上皮候选。F1.5第一次正式运行发现sample37只有7个上皮候选；按SCTransform固定的`min_cells=5`计算，该样本只有283个SCT可用基因，使所有样本共同PCA特征不足500，流程按预设检查停止。

研究者于2026-07-30批准：

- sample37的7个细胞不进入F1.5独立样本SCT模型及后续恶性判定主对象。
- 这7个细胞仍保留在F1.4全细胞对象，并逐细胞记录排除原因。
- 规则按数据写为：样本需有至少500个SCT可用基因，即每个基因至少在5个该样本上皮细胞中检出。
- 其余39个样本全部保留；其中SCT可用基因最少的样本仍有7,821个，不受影响。

该处理解决的是极小样本无法支持独立SCT模型的问题，不代表7个细胞被重新判为非上皮或低质量。

## 11. F1.5二次聚类结果

- 正式输入：26,821个上皮细胞、39个样本。
- 归一化：逐样本SCTransform v2，`vars.to.regress=NULL`，`min_cells=5`。
- Harmony：仅按`sample_id`校正；使用前30个Harmony维度。
- 主Leiden分辨率0.6，得到19个上皮cluster。
- 分辨率0.2、0.4、0.8和1.0分别得到10、15、20和24个cluster；0.6与0.8只相差1簇，主尺度没有处于明显跳变点。
- 共同PCA特征1,666个；39个SCT模型完整。
- 19个cluster均覆盖至少24个样本；其中16个覆盖37至39个样本。没有单一样本垄断的新cluster。
- 输出对象保留RNA raw counts、SCT/PCA/Harmony/UMAP、所有分辨率字段和原全细胞来源cluster。
- 11项输出验证全部通过，未发现MLMOD、预后或生存字段参与F1.5。

## 12. F1.5 ambient与注释复核

19个新cluster的上皮panel最高方向在raw与corrected中均保持一致。多数cluster前列marker计数保留率为0.925至0.995；需要重点复核3簇：

- c1：污染中位数0.004、P90为0.682，说明是少数细胞形成高污染长尾；前列marker保留率0.969且上皮/亚型方向不变，因此建议保留并标记定向敏感性，而不是整簇删除。
- c13：污染中位数0.277、P90为0.909、前列marker保留率0.560；PTPRC、TYROBP、LST1、CD37等免疫marker成套表达，DoubletFinder敏感性阳性比例10.6%，为19簇最高。建议作为免疫混合/双细胞疑似簇排除出F1.6。
- c16：污染中位数0.490、P90为0.959、前列marker保留率0.463；MZB1、JCHAIN、IGHA1/2、PRDM1等浆细胞marker成套表达，而广义上皮marker很弱。建议作为浆细胞混合或ambient主导簇排除出F1.6。

其余17簇均可解释为上皮亚型或上皮状态，包括pit/surface mucous、mucous neck、chief、parietal、enteroendocrine、intestinal-like，以及增殖、金属硫蛋白高、热休克、炎症和HLA-II状态。状态标签不替代上皮亚型，更不替代F1.6的CNV恶性判定。

## 13. F1.5批准结论

研究者于2026-07-30批准19簇决策草案并冻结正式`F1_epithelial_cluster_review_approved.tsv`：17簇进入F1.6，c13和c16排除。c13、c16仍保留在F1.5完整对象中，但不进入inferCNV/CopyKAT观察细胞及后续恶性主对象。F1.5检查点通过，允许开始F1.6输入和资源预检。

## 14. 可复现性与资源记录

- 服务器：Ubuntu 22.04，48 vCPU，96 GiB内存。
- R 4.4.3，Seurat 5.5.0。
- F1.3主要并行设置：future 12 workers、scDblFinder 8 workers、BLAS/OpenMP 4线程。
- 最终cluster分组DecontX约12分钟；ambient影响复核约95秒。
- DecontX函数按单样本运行时主要使用单核；本轮在已完成31/40时没有中断重排。后续按样本独立的inferCNV和CopyKAT任务将预先拆分并行，充分使用服务器CPU。
- F1.5使用24个外层future进程、每进程1个BLAS/OpenMP线程；最终成功运行约19分钟，实测内存约30至32 GiB。
- 运行中曾出现SeuratObject由Matrix 1.7-2构建、当前环境为Matrix 1.7-5的ABI提示；最终对象、raw counts、降维和校正矩阵验证均通过。F1.5启动前已从本地上传SeuratObject 5.4.0源码并针对Matrix 1.7-5重新编译，包加载和最小Seurat对象smoke test通过，ABI提示不再出现。

关键结果、日志、作图数据和小型表已同步回本地。最终注释对象为18,492,611,572 bytes，SHA256为`AB99F8E681C60D2AE9006E6A51E8889BE339E8C29D2A1B79387FE44DCA3227B9`；40个corrected矩阵合计2,826,560,016 bytes。主对象和40个矩阵共41项均已按服务器清单完成大小与SHA256逐文件校验，41/41通过。

F1.5上皮二次聚类对象为3,658,906,739 bytes，SHA256为`AB503E7058679E750E8B559DBFF1222DF63095930E8D7341500340D373C1F5DE`。对象、结果表、图、脚本和日志共40项已按F1.5服务器清单完成本地逐文件校验，40/40通过；服务器端传输压缩包和补字段前临时备份已在本地验证通过后删除。

## 15. F1.6 inferCNV/CopyKAT第一轮计算

F1.6对17个获准上皮cluster中的25,636个细胞、39个样本完成了逐样本inferCNV和CopyKAT计算。39/39个样本检查点完整，inferCNV、CopyKAT和子群成员表的细胞ID一一闭合。正式资源设置为2个样本并行、每样本10个inferCNV线程或10个CopyKAT核；4样本并行试运行因物理内存耗尽而停止，未用于正式结果。

sample34的CopyKAT算法已完成10,634个输入细胞的prediction，但包内最终CNA表写出发生维度错误。恢复分支仅在prediction细胞ID唯一且与输入集合完全相等时复用官方prediction；候选上皮1,467个细胞均获得调用。该样本的最终CopyKAT CNA表和默认热图记为不可评估。

首轮结果暴露出两个需要在恶性判定前解决的问题：

- inferCNV `leiden_resolution="auto"`产生1,183个子群；29.4%不超过5个细胞，且127个单细胞子群来自总细胞数大于20的父群，提示实际分辨率过细。
- Normal Gastric中CopyKAT明确调用细胞的`aneuploid`比例为72.2%，高于Primary Tumor的54.3%。因此CopyKAT二分类在本数据中不能单独证明恶性，只能作为联合证据。

已修正联合判定代码：高置信恶性也必须来自肿瘤或腹膜转移样本，并满足肿瘤程序支持且正常程序不支持；另设硬检查禁止Normal Gastric进入`06a/06b`。

当前状态为`COMPUTATION_COMPLETE_REVIEW_PENDING`。现有1,183行模板尚未批准，建议先做代表性样本的Leiden分辨率诊断，再冻结正式审核粒度。详见`reports/F1_execution/F1_06_computational_review.md`。
