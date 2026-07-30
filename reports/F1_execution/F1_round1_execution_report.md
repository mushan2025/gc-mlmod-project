# F1 第一轮执行报告：F1.1-F1.4 人工注释审核点

## 1. 执行范围与当前停止位置

本轮已完成：

1. F1.1：40个公开样本矩阵导入、方向和整数原始计数检查。
2. F1.2：冻结QC、scDblFinder主双细胞删除和DoubletFinder敏感性分析。
3. F1.3：逐样本SCTransform v2、SCT PCA、Harmony、UMAP和Leiden聚类。
4. F1.4第一部分：全细胞cluster marker、粗谱系marker图和人工注释模板。

当前按计划停在人工粗谱系注释审核点。尚未运行DecontX、F1.5上皮重聚类或F1.6 inferCNV/CopyKAT，也未生成恶性上皮对象。

## 2. 数据导入与QC结果

- 公开输入：40个样本，158,641个细胞，归档矩阵共有26,571个feature。
- 矩阵方向、非负整数值、样本数和细胞数检查均通过。
- 固定QC排除5,321个细胞，153,320个细胞通过固定QC。
- scDblFinder判定并排除10,670个双细胞。
- 最终保留142,650个细胞，占公开细胞的89.92%。
- 去掉经验性`nCount_RNA > 1000`条件时，会比主规则额外保留4,871个边界细胞；该结果只记录为预注册敏感性范围，不改变主对象。

固定QC严格使用已批准规则：每个样本先保留至少在3个细胞中检出的基因，再要求`nFeature_RNA >= 500`、`nFeature_RNA < 6000`、`mt_percent <= 20`、`HB_percent < 5`和`nCount_RNA > 1000`。

## 3. 双细胞结果

scDblFinder是主方法，DoubletFinder只用于方法敏感性检查。40个样本中39个完成DoubletFinder；sample37在固定QC后仅172个细胞，按预登记规则记为`not_evaluable_small_sample`。

在可评价细胞中，scDblFinder与DoubletFinder的Jaccard一致度为0.130。最终保留细胞中有5,209个被DoubletFinder单独标记；这些细胞没有被自动删除，而是在注释时作为风险信息复核。cluster 22、15、6和5的DoubletFinder敏感性阳性比例相对较高，需结合marker混合情况重点查看。

## 4. SCTransform、Harmony与聚类

- 主归一化：逐样本SCTransform v2，`vars.to.regress = NULL`。
- 候选高变基因：3,000个；实际进入PCA的共同已缩放基因：2,141个。
- PCA：50维；Harmony仅按`sample_id`校正；主聚类使用前30个Harmony维度。
- 主Leiden分辨率：0.6，共24个cluster。
- 分辨率0.2、0.4、0.8和1.0分别得到15、20、27和29个cluster，主分辨率处于可解释的中间尺度。
- 03a对象回读检查全部通过：142,650个细胞顺序不变，RNA raw counts指纹一致，40个SCT模型、PCA/Harmony/两个UMAP均完整且数值有效。
- 未在对象中发现MLMOD、预后或生存信息参与聚类。

PCA与`nCount_RNA`、`nFeature_RNA`、`mt_percent`的最大绝对相关分别约为0.476、0.377和0.386。这说明部分主成分仍包含细胞复杂度或线粒体状态信息，但目前marker和主要谱系结构清楚，不足以支持反向修改已冻结QC或回归这些变量。

## 5. F1.4首轮注释材料

- 共得到23,574条阳性marker记录。
- 24个cluster均有可用marker，每个cluster有93-2,155条阳性marker。
- 已生成24行人工注释模板及粗谱系marker DotPlot。
- 当前只形成待审核判断，未写入批准标签。

初步可辨认的主要大类包括T/NK细胞、B/浆细胞、髓系细胞、肥大细胞、成纤维细胞、内皮细胞、周细胞/平滑肌细胞和多组上皮细胞。cluster 18主要表现为热休克/即时早期应激程序，cluster 22呈间皮样或特殊基质样特征且双细胞敏感性风险较高，cluster 24为增殖细胞；这三类需要先判断其母谱系，不能只按状态marker命名。

## 6. 方法边界与下一步

公开输入是作者公开的called/retained-cell基因计数矩阵，不是raw droplet矩阵，因此不能重做真实barcode calling、EmptyDrops或以空液滴为背景的SoupX。DecontX将在粗谱系标签经研究者批准后，按样本对保留细胞的raw integer counts运行；其结果先用于ambient风险评估，不默认替换raw主矩阵。

下一步只需要审核24个cluster的粗谱系标签，重点复核cluster 18、22、24以及DoubletFinder敏感性比例较高的cluster。标签批准后再继续DecontX和F1.5；F1.6 CNV分析仍需结合实际上皮细胞与参考细胞规模做单独资源确认。

## 7. 可复现性记录

- 服务器：Ubuntu 22.04.5 LTS，24个物理核/48个逻辑线程，实测总内存90.79 GiB。
- R 4.4.3，Seurat 5.5.0。
- 主要并行设置：future 12 workers、scDblFinder 8 workers、BLAS/OpenMP 4线程；后续CNV预留16线程。
- F1.3实测进程峰值RSS约60.23 GiB，未使用swap。
- 最终成功运行耗时：F1.1约7分钟，补齐DoubletFinder的F1.2重跑约43分钟，F1.3约29分钟，F1.4 marker与模板生成约9分钟。
- 运行代码对应Git提交：`7f71797ff718ea413865d7d90b3f75520a99c9bd`。

关键RDS、结果表、图、作图数据、参数版本和日志已经回收到本地：

- `02_all_cells_qc_filtered.rds`：SHA256 `9D33A1676176F43FEE8248E5EAA5B66E2785169CFE901C9CA25543FE0F1C892C`。
- `03a_all_cells_sct_harmony_clustered.rds`：SHA256 `2D2BB752AACA39267E1DDA91527142A9DB7374EB016E3C8A68B080CB82A8D383`。
- 结果、日志和脚本便携包：SHA256 `49BE2AD30CBF2BEDFA082FD4332F5A386BFE027F3C3A5BB6A7AFD00D08668E2F`。
- 清单中119项本地文件的大小和SHA256全部通过；唯一未下载的是可由本地40个公开矩阵重建的`01_all_cells_raw_or_initial.rds`，其服务器指纹已保存在清单中。

因此，后续可以直接从本地`02`对象重做QC后流程，或从`03a`对象重新注释、重绘聚类图和marker图，无需依赖当前服务器继续存在。
