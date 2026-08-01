# F1.6 inferCNV/CopyKAT 计算结果审核

生成日期：2026-07-30

当前状态：`HISTORICAL_RUN_ARCHIVED_NEW_METHOD_RERUN_PENDING`

> 2026-08-01方法修订说明：本报告以下数字完整保留为第一次计算的历史结果，不改写为新结果。该次inferCNV把观察细胞按原上皮cluster拆成631个父组且使用`HMM=FALSE`；CopyKAT相当于后续定义的B臂。经审核后，正式方案已改为每个样本一个`observations`组、正式i6 HMM、热图主判和CopyKAT A/B/C三臂，并禁止CopyKAT单独阳性细胞进入`06a`。因此旧checkpoint、1,183个subcluster及旧审核模板不得用于最终恶性判定，需按修订脚本重跑；旧文件仅作方法问题追溯。

## 1. 实际完成范围

- 输入为F1.5批准的17个上皮cluster，共25,636个候选上皮细胞、39个样本。
- inferCNV使用RNA raw counts，`cutoff=0.1`、`denoise=TRUE`、`HMM=FALSE`、`analysis_mode="subclusters"`。
- inferCNV内部子聚类使用Leiden，`k_nn=20`、`leiden_method="PCA"`、`leiden_function="CPM"`、`leiden_resolution="auto"`。
- CopyKAT仅使用当前样本的QC后singlet细胞，并优先指定当前样本高/中置信T/NK细胞为known normal。
- 正式运行采用2个样本并行；每个inferCNV任务10线程，每个CopyKAT任务10核。4个样本并行的试运行曾耗尽物理内存，因此未用于正式主运行。

39/39个样本均形成有效检查点，失败文件为0。inferCNV、CopyKAT和inferCNV子群成员表均含25,636个唯一细胞，三者细胞ID集合完全一致。39个inferCNV参考集合均达到最低数量要求，参考细胞数为123至500。

## 2. CopyKAT结果

25,636个候选上皮细胞中：

- `aneuploid`：12,605。
- `diploid`：8,299。
- `uncalled`：4,732。

sample34的CopyKAT已完成细胞预测，但包内在写出最终CNA结果时发生维度不一致错误。恢复分支只复用了CopyKAT已经保存的官方prediction；10,634个输入细胞与10,634行prediction逐一匹配，候选上皮中的1,467个细胞均得到调用。sample34的最终`CNA_results`和默认CopyKAT热图不可评估，该限制已单独记录。

## 3. 正常胃阴性对照揭示的限制

正常胃组织中的8,451个候选上皮细胞来自10个样本。CopyKAT结果为：

- `aneuploid`：4,785。
- `diploid`：1,845。
- `uncalled`：1,821。
- 在有明确调用的细胞中，`aneuploid`比例为72.2%。

作为对照，Primary Tumor中有明确调用的细胞被判为`aneuploid`的比例为54.3%。运行日志确认CopyKAT确实识别并使用了当前样本known normal，因此该结果不能简单归因于参考细胞未传入。

结论是：本数据中的CopyKAT二分类缺乏足够的正常上皮特异性，不能把`aneuploid`直接等同于恶性，也不能把inferCNV与CopyKAT的一致称为独立DNA证据。CopyKAT最多作为辅助证据，必须同时满足肿瘤来源、肿瘤程序支持、正常程序不支持及可解释的大片段CNV模式。

代码已修正一个与此相关的逻辑漏洞：`malignant_high_confidence`现在也必须来自Primary Tumor或Peritoneal Metastasis，并满足肿瘤程序支持且正常程序不支持；另加硬检查，禁止任何Normal Gastric细胞进入`06a/06b`。

## 4. inferCNV结果与子聚类粒度

- 所有25,636个候选细胞均有inferCNV结果。
- 16,525个细胞高于同一样本免疫参考的CNV burden P95。
- Normal Gastric、Primary Tumor和Peritoneal Metastasis中高于参考P95的比例分别为43.8%、75.3%和54.1%。

这说明inferCNV对肿瘤上皮有一定区分信息，但免疫参考不能完全消除正常胃上皮与免疫细胞之间的谱系表达差异。最终恶性判定仍需结合正常胃同谱系上皮对照和大片段、同群一致的CNV模式，不能只用单一burden阈值。

自动Leiden产生1,183个inferCNV子群，对应631个`sample × epithelial cluster`组合。子群细胞数中位数为16：

- 348个子群不超过5个细胞，占29.4%。
- 148个为单细胞子群，占12.5%。
- 127个单细胞子群来自总细胞数大于20的父群。
- 108个单细胞子群来自总细胞数大于50的父群。

因此，过细拆分不能只用“原父群本来很小”解释。当前1,183行模板不适合直接逐行人工审核，也与inferCNV官方关于需要按实际数据调整subclustering分辨率的提醒一致。

## 5. Ambient RNA检查点

sample35的ambient污染估计最高，但它属于Normal Gastric，不会进入恶性主对象。sample22的候选上皮污染存在高尾部，后续若其恶性判定依赖边界CNV模式，应运行raw与DecontX-corrected inferCNV定向敏感性；不进行全体样本统一校正。

## 6. 本地保存与校验

- 小型审核包包含结果表、日志、检查点、inferCNV图、CopyKAT JPEG、prediction、聚类对象、实际脚本及阶段报告，共749个文件；压缩包和包内文件均已完成SHA256校验，0个缺失、0个不一致。
- 39个inferCNV最终对象、39个CopyKAT结果对象和sample34原始CNA恢复文件共79项，解压后合计13,989,564,741 bytes；已保存到`C:\Users\14799\gc-mlmod-server-artifacts\F1_06_20260730\replot_objects`，逐文件SHA256校验为0个错误。
- 服务器临时`/swapfile_f1_16g`从未实际使用，计算与下载校验完成后已安全移除。
- 服务器上的完整原始输出和冗余中间矩阵暂未删除。

## 7. 当前建议

不批准现有1,183行恶性判定模板。修订后的正式重跑先使用每样本单一`observations`组和`leiden_resolution="auto"`，不设目标subcluster数；只有新热图仍出现广泛过度拆分且组间没有可区分大片段模式时，才比较一个相邻更低分辨率，不做参数网格，也不按恶性细胞产量选参数。

CopyKAT改为A/B/C三臂方法诊断。A臂只有在自估diploid组成和CNA输出可检查时才提供辅助支持；B/C用于解释基线敏感性。CopyKAT二分类不得单独决定恶性，也不根据其与inferCNV的一致性升级为独立验证。正式`06a`必须具有inferCNV热图大片段支持。
