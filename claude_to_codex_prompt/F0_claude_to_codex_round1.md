# Codex修改指令：主线TXT第一轮第三方审查修订

## 背景

我请了一个独立于我们讨论的第三方（Cowork中的Claude）审查了我们讨论好的主线TXT（`胃癌MLMOD亚群主线研究方案.txt`）。以下是审查后我同意的修改点。请逐条修改主线TXT，修改完成后告诉我每条改了什么、改在哪一行附近。

注意：第4-9条涉及F1和F2，目前这些部分还是大纲阶段，后续我们还会详细讨论。但这几条是方向性的重要修正，请先写入，后续细化时再调整措辞。

---

## 修改清单

### 修改1：F0.5 增加UMI count分布特征检查（F0改动）

当前F0.5只要求检查`integer_value_rate`来判断矩阵是否为count-like。但"非负整数"不等于"UMI count"——有些处理流程会输出整数但实际是normalized后round的结果。

请在F0.5中，在integer_value_rate相关段落之后，增加以下要求：

> 除integer_value_rate外，正式审计还应检查数值分布是否符合UMI count的典型模式：
> （1）per-cell total count（nCount）分布：应呈现连续分布且有合理范围（通常几百到几万），不应出现所有细胞total count完全相同或高度集中在少数round number的情况（这提示normalized后取整）。
> （2）per-gene mean expression分布：应呈现zero-inflated右偏分布，大量基因均值接近零，少数基因均值较高。
> （3）若发现所有细胞total count完全一致、或数值分布呈现明显的normalization artifact，则不能当作raw UMI count使用，即使integer_value_rate >= 0.99。
> （4）审计结论中应记录per-cell nCount的median、min、max、Q1、Q3，以及是否发现normalization artifact。

### 修改2：F0.7 增加precheck与正式audit一致性比对要求（F0改动）

当前F0.7第143行写了"若正式审计结果与当前预审计不一致，必须暂停"，但没有说明Biomni应如何比较。

请在F0.7的这句话附近，将其扩展为：

> Biomni生成data_audit.tsv后，必须与研究者提供的gse183904_csv_structure_precheck.tsv逐样本比对以下关键字段：matrix_rows_genes、matrix_cols_cells、mt_gene_count、hb_gene_count、gene_order_sha256、suspected_matrix_type和audit_decision_precheck。若任一样本的上述字段出现不一致，Biomni必须暂停，在analysis_log.md中记录差异内容和可能原因，并等待研究者确认后才可继续。不一致不一定意味着错误（例如正式审计使用了更严格的方法），但必须有明确解释。

### 修改3：区分全细胞注释和上皮二次聚类的marker panel使用范围（F0+F1改动）

当前cell_type_marker_panel.tsv包含22个细胞类型，其中有5种上皮亚型（Pit_mucous_epithelial、Mucous_neck_epithelial、Chief_epithelial、Intestinal_like_epithelial、Enteroendocrine_epithelial）。这些亚型在全细胞聚类阶段很难区分，只有在上皮二次聚类后才有意义。

请做以下修改：

（a）在F0.6（marker panel相关段落），增加一段说明：

> cell_type_marker_panel.tsv中的细胞类型分为两个使用层级：
> 第一层级（全细胞注释使用）：Epithelial、T_cell、Treg、NK_cell、B_cell、Plasma_cell、Myeloid、Monocyte、Macrophage、Dendritic_cell、Mast_cell、Fibroblast、CAF、Endothelial、Pericyte、Smooth_muscle_cell、Cycling。全细胞注释时，上皮细胞统一标注为Epithelial，不在此阶段区分上皮亚型。
> 第二层级（上皮二次聚类后使用）：Pit_mucous_epithelial、Mucous_neck_epithelial、Chief_epithelial、Intestinal_like_epithelial、Enteroendocrine_epithelial。这些亚型仅在F1.5上皮细胞二次聚类后的精细注释中使用。
> Biomni在F1.4全细胞注释时不应尝试区分上皮亚型。

（b）在F1.4（细胞类型注释段落），增加一句：

> 全细胞注释阶段只使用cell_type_marker_panel.tsv中的第一层级细胞类型；上皮细胞统一标为Epithelial，不区分上皮亚型。

（c）在F1.5（上皮二次聚类段落），增加一句：

> 上皮二次聚类后，使用cell_type_marker_panel.tsv中的第二层级上皮亚型marker进行精细注释，结合cluster marker genes和可视化验证。

### 修改4：明确F2.4 Bulk预后验证属于主线而非补充（F2改动）

F2.4已经在主线TXT中，这一点保持不变。但请在F2.4开头增加一句逻辑说明：

> 通路打分只反映生物学活性，不足以证明该亚群具有临床研究价值。需在Bulk数据中验证高MLMOD特征基因集能否预测患者预后——预后显著才能确认亚群值得进入后续深度分析。

### 修改5：F2.3高MLMOD亚群定义需结合分布形态判断（F2改动）

当前F2.3写"取MLMOD_Score上四分位数（Q3）定义高活性细胞，或K-means聚类（K=2/3）"。这缺少对分布形态的判断。

请将F2.3相关段落修改为：

> 打分分布：在恶性上皮细胞中可视化MLMOD_Score分布（小提琴图/密度图）。
> 阈值设定：先检查MLMOD_Score分布形态。若分布呈现明显bimodal或multimodal结构（可用mixture model或Hartigan's dip test辅助判断），则按分布谷值或mixture model交叉点定义高活性群体；若分布为连续单峰，不应硬切，建议取MLMOD_Score上四分位数（Q3）定义高活性细胞作为探索性分析，同时以连续变量形式做相关性和回归分析作为主分析。K-means聚类（K=2/3）可作为辅助参考，但不能作为唯一分组依据。

### 修改6：F1.2 明确decontX的"评估"与"校正"边界（F1改动）

当前F1.2写"主流程使用decontX评估/校正"，没有明确到底用不用corrected counts。

请将F1.2中ambient RNA相关段落修改为：

> ambient RNA：主流程使用decontX评估环境RNA污染水平并记录每个样本的estimated contamination rate。是否使用decontX corrected counts需根据污染估计水平决定：若整体污染水平较低（例如estimated contamination rate中位数<5%），可在主流程中使用原始counts，保存decontX评估结果作为记录；若污染水平较高，则使用decontX corrected counts，并保存校正前后对比图，同时在decision_evidence_log.tsv中记录此决策。SoupX仅在存在可用raw droplet/background信息时作为补充。GSE183904公开文件是filtered count-like矩阵，缺少empty droplets背景；因此不能把SoupX作为默认强制步骤。

### 修改7：F1.5 上皮二次聚类增加细胞周期回归对比（F1改动）

当前F1.5没有讨论细胞周期回归问题。肿瘤上皮中cycling cells可能形成独立cluster，影响亚群定义。

请在F1.5中，在聚类参数讨论附近增加：

> 上皮二次聚类时需对比两种方案：（1）不回归细胞周期；（2）回归S.Score和G2M.Score。若两种方案的聚类结构和恶性判定结果基本一致，则主流程不回归（保留增殖信号作为真实生物学特征）；若回归后cluster结构发生明显变化（例如cycling cluster消失并重分配），则需在结果中展示两种方案对比，并在decision_evidence_log.tsv中记录最终选择及理由。

### 修改8：F1.6 inferCNV阈值方法一致性检查（F1改动）

当前F1.6写阈值为"max(reference cells的99%分位数, reference median + 3*MAD)"。

请在该阈值描述之后增加：

> 若99%分位数与median + 3*MAD给出的阈值差异超过2倍，提示reference cell分布可能存在异常（例如reference中混入了可疑上皮细胞、某个谱系reference细胞数量过少导致分布不稳、或存在批次效应）。此时必须先检查reference cell质量，排除异常后重新计算阈值，并在decision_evidence_log.tsv中记录。

### 修改9：F2.2 综合打分标准化方法修改（F2改动）

当前F2.2写"对5种算法结果标准化（0-1缩放）后取均值"。min-max scaling对outlier非常敏感。

请将该句修改为：

> 综合打分：对5种算法结果进行rank-based normalization（即对每种算法的分数在所有细胞中取秩次后除以细胞总数，映射到0-1区间），再取均值作为最终MLMOD_Score。不使用min-max scaling，因为其对outlier敏感，单个极端值可能压缩其余细胞的得分范围。

---

## 修改后请确认

1. 每条修改的具体位置（哪个小节、哪一行附近）。
2. 是否有任何修改与现有内容冲突，如果有请指出。
3. 修改后把文件同步到WPS路径。
