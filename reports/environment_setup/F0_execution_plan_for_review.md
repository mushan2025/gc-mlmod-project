# F0 正式执行计划（审核版）

状态：执行计划和脚本已经写好，等待 Claude Code 审核及用户批准；尚未正式执行完整 F0。

## 1. 目标与边界

F0 的任务是在进入 F1 前建立可信的数据与可复现性基线，具体包括：

- 核验 GSE183904 压缩包、40 个样本矩阵、GEO 样本映射关系和文件校验和；
- 对每个矩阵中的全部数值进行流式结构与数值审计；
- 在当前公开细胞上重新计算已经批准的固定 QC 规则，仅用于确认规则可计算并估计其影响；
- 记录从组织获取到公开数据导出的已知和未知处理历史；
- 生成项目数据清单、各 F 节数据就绪状态、方法前提决策和 F0 gate 报告。

F0 **不会**删除细胞、创建过滤后的 Seurat 对象、运行双细胞或环境 RNA 算法、启动 F1，也不会作出生物学结论。只有 F0 gate 经独立审核并由用户批准后，F1 才会真正开始排除细胞。

## 2. 技术判断

### 2.1 前提验证

- 本地压缩包预期包含 40 个“基因 × 细胞”的 `csv.gz` 矩阵；合计 26,571 个存档 feature 行和 158,641 个公开细胞。
- 当前公开输入不包括 FASTQ、raw/empty droplets、Cell Ranger cell-calling 记录、作者排除的 barcode 清单或完整 DoubletFinder 参数。
- QC 定义已在主线方案中冻结，因此 F0 必须按固定规则确定性执行，不能根据 F0 结果临时挑选阈值。

### 2.2 定量依据

对 sample1 完整文件的 pilot 检查得到：

- 2,685 个公开细胞，26,571 个存档 feature 行；
- 列名全部符合本数据的 10x 16 nt barcode 加数字后缀格式，行名中 barcode-like 数为 0，方向状态为 `pass_gene_by_cell`；
- 每个样本应用 `min.cells=3` 后有 19,294 个工作 feature；
- 2,684 个细胞通过原研究报告的 `nFeature_RNA >= 500`、`nFeature_RNA < 6000`和过滤`mt_percent > 20`规则；
- 按顺序统计时，本项目经验性`nCount_RNA>1000`额外排除53个细胞，经验性`HB_percent<5`再额外排除0个细胞；
- 最终有 2,631 个细胞通过全部固定规则；
- 实际参与 HB 百分比计算的精确 globin 基因交集为 `HBA1,HBA2,HBB,HBD`；
- `HBEGF,HBP1,HBS1L` 虽以 `HB` 开头，但不是 globin 基因，已明确排除。

只读验证脚本还使用合成细胞检查每个不等式边界，并使用合成矩阵检查 globin 基因的精确匹配。它还把一个合成的 40 样本审计结果送入 Step3 的处理史构建和 Step4 的十项 gate，用于防止字段改名后下游决策在无提示的情况下失效。

方向专项反例使用人工构造的 `cell × gene` 矩阵：列名为基因、行名为 10x barcode。Step2 将其判定为 `fail_not_gene_by_cell`，证明方向检查不是根据扩展名或预期值直接写死。F0 环境锁验证同时确认当前 Python 3.10.11、NumPy 2.2.6 与两份依赖规格一致。

### 2.3 结论

F0 适合在默认笔记本上通过 Windows 原生 Python 流式执行，不需要 R、WSL、GPU 或临时服务器。完整执行仍须等待 Claude Code 审核和用户明确批准。

正式运行前必须再次测量可用内存和 D 盘剩余空间。当前资源记录只能支持执行环境规划，不能代替正式运行前的实时检查。

## 3. 输入文件

完整清单见 `reports/environment_setup/F0_input_file_checklist.tsv`。关键输入包括：

- `data/public_downloads/GSE183904_RAW.tar`
- `data/public_downloads/GEO_metadata/GSE183904_series_matrix.txt.gz`
- `docs/source_verification/GSE183904_processing_history_source_audit.tsv`
- `data/metadata/` 下已有的 manifest 文件
- `results/F0_audit/` 下已登记的结构预检查结果
- `environment/` 下的执行环境基线文件
- F0 独立依赖规格 `environment/F0/requirements.txt` 和 `environment/F0/environment.yml`
- 只读参考文件 `data/metadata/cell_type_marker_panel.tsv`

本次 F0 执行不包含任何联网下载。

## 4. 脚本及执行顺序

总入口脚本把以下四个阶段作为一个完整的 F0 审批单元顺序执行：

1. `F0_step1_structure_and_extract.py`
   - 检查所有必需路径和压缩包 SHA256；
   - 确认压缩包中恰好有 40 个文件名不重复的 `csv.gz` 成员；
   - 将压缩矩阵解压到 `data/processed_input/GSE183904/`；
   - 仅将候选文件登记为 `expected_gene_by_cell_matrix_pending_validation`；文件扩展名不能证明矩阵方向；
   - 写出处理后输入清单，所有 SHA256 统一使用大写。

2. `F0_step2_sample_info_and_audit.py`
   - 分别依据文件名和 GEO title 建立样本及分组信息，并核对两种证据是否一致；
   - 流式读取 40 个矩阵中的每个数值，不建立占用大量内存的合并矩阵；
   - 根据左上角、列名 10x barcode、行名基因格式及行名是否反而像 barcode，正式验证 `gene × cell` 方向；
   - 检查维度、整数性、缺失/非法/负值、重复基因或 barcode、基因顺序、稀疏性和标准化伪影哨兵；
   - 仅为判断输入矩阵形态，保留基于全部原始行计算的 `raw_full_nCount` 摘要；
   - 对每个样本建立唯一的 `min.cells=3` QC 工作 feature 空间；
   - 在该工作空间内计算唯一一套细胞 QC 指标；
   - 统计固定规则的通过/不通过细胞数，并执行冻结的 sample1 回归检查。

3. `F0_step3_inventory_and_markers.py`
   - 建立数据集、文件、metadata 字段和外部资源清单；
   - 在 `F0_file_manifest.tsv` 中记录所有正式 F0 脚本、只读验证脚本及两个 F0 环境锁文件的大写 SHA256；
   - 审计 marker panel 的结构，但不修改只读 panel；
   - 合并原始来源记录的处理史与 F0 实测结果；
   - 明确区分作者报告、F0 重算、跨来源推断和仍未知的处理环节。

4. `F0_step4_decisions_and_gate.py`
   - 写出数据就绪状态、方法前提、决策证据和样本排除表；
   - 评估十项 F0 gate checklist；
   - 写出全局数据摸底报告和执行报告；
   - 完成后停止，不自动启动 F1。

所有脚本在未提供 `--execute` 时默认只做不写文件的 dry run。总入口脚本一旦遇到异常或阻断条件就立即停止，不会跳过错误继续执行。

## 5. F0 对 F1.2 冻结 QC 的影响测量

QC方法学的唯一权威位置是主线方案F1.2。F0不定义新的QC方法、不删除细胞、不建立F1过滤对象；这里只规定F0如何在当前公开矩阵上复算同一规则并把逐样本影响写入`data_audit.tsv`，用于血缘和实现核对。

### 5.1 不可修改的原始输入

存档矩阵中的 26,571 行全部原样保留。`raw_full_nCount_*` 只用于判断 count 矩阵形态以及是否可能经历过标准化，绝不参与细胞删除决策。

### 5.2 用于影响测量的工作 feature 空间

每个样本独立保留至少在 3 个细胞中检出的 feature，然后仅在这一工作空间内计算：

- `nCount_RNA`：该细胞的总 UMI/count 数；
- `nFeature_RNA`：该细胞中检出的基因数；
- `mt_percent`：线粒体基因 count 占总 count 的百分比；
- `HB_percent`：冻结 globin 基因 count 占总 count 的百分比。

F0不另建第二套细胞QC指标空间；正式解释以F1.2为准。

### 5.3 与 F1.2 的一致性

F0脚本读取并复算F1.2冻结的单一mask：

```text
500 <= nFeature_RNA < 6000
nCount_RNA > 1000
mt_percent <= 20
HB_percent < 5
```

其中feature/nFeature和线粒体边界来自原研究报告；`nCount_RNA>1000`和`HB_percent<5`是本项目预注册的经验性操作阈值，不归因于原作者，也不声称为文献或本数据集最优。脚本只记录每条规则各自不通过的细胞数和全部规则合并后的唯一不通过细胞数，避免一个细胞违反多条规则时被重复计数。

### 5.4 globin 定义核对

F0使用F1.2冻结的精确globin panel，不在本计划重复定义。每个样本记录预期panel、存档矩阵中实际存在的基因、经`min.cells=3`后实际参与QC的基因、缺失基因，以及因不是panel成员而被排除的`^HB`伪匹配基因；禁止使用宽泛前缀匹配计算`HB_percent`。

### 5.5 下游边界

每个样本的 `min.cells=3` 仅用于建立 F0/F1 的 QC 工作空间，不是差异表达、pseudobulk 或评分分析的永久基因过滤条件。这些下游方法必须采用各自获批的表达覆盖规则。

## 6. 处理历史的解释原则

以下情况可以同时成立，并不互相矛盾：

- 原始来源报告使用 Cell Ranger v3.0/hg38 和 Seurat QC 阈值；
- 公开矩阵是非负整数形式的已调用/已保留细胞 count 矩阵；
- 公开矩阵有 158,641 个细胞，而论文最终 40 个组织样本对象有 152,423 个细胞；
- 作者 QC、DoubletFinder 和公开导出的准确先后顺序及其具体 barcode 影响仍无法完全还原。

F0 不得把 6,218 个细胞的差值归因于某一个特定步骤。重算得到的通过细胞数只说明获批规则应用到“当前公开矩阵”时会产生什么结果，不能重建公开前已经被排除的细胞。

在 `F0_author_processing_audit.tsv` 中：

- `author_reported_status` 只记录经审阅来源明确报告的内容；
- `record_status` 单独记录 F0 重算或跨来源核对结果；
- F0 自动生成的记录在作者状态字段中明确写为 `not_applicable_*`，防止把本项目观察误写成作者声明。

## 7. 运行命令

### 7.1 F0 独立环境规格

F0 依赖与 F1-F8 分开管理。以下是两条**可替代**路线，不应重复混装：

- pip 路线：Python 3.10.11 + `environment/F0/requirements.txt`；
- Conda 路线：`environment/F0/environment.yml`。

当前审核和只读验证使用已经存在的 Windows Python 3.10.11 / NumPy 2.2.6，
不需要为了本次审核重新创建环境。正式执行前只需确认实际解释器版本与锁文件一致。

### 7.2 只读检查输入和输出契约

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/run_F0_full_audit.py --project-root .
```

该命令不写正式输出，只检查必需输入是否存在并列出计划生成的文件。

### 7.3 只读代码验证，包括方向测试和真实 sample1 回归

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/validate_F0_readonly.py --project-root .
```

该命令只在系统临时目录中写测试文件，不生成正式 F0 产物。

该验证会确认 Step1 方向状态仍为 pending、gene × cell 合成矩阵可通过、
cell × gene 合成矩阵被拒绝，并核对 F0 环境锁、QC 边界和 sample1 冻结数字。

### 7.4 正式执行

仅在 Claude Code 审核通过并得到用户批准后允许运行：

```powershell
& 'C:\Users\14799\AppData\Local\Programs\Python\Python310\python.exe' scripts/F0_setup/run_F0_full_audit.py --project-root . --execute
```

## 8. 阻断条件

出现以下任一情况时，F0 必须暂停：

- 压缩包缺失、无法读取、校验和不一致、成员文件名重复，或 `csv.gz` 成员数不是 40；
- 压缩包文件名、GSM accession、`sample_id` 与 GEO title 的映射不一致；
- 任何样本的分组仍为 `Unclear`；
- 矩阵方向、维度、基因顺序或预登记结构字段出现无法解释且未经批准的差异；
- Step1 提前把预期方向写成已验证，或 Step2 无法确认列为 10x barcode、行为基因；
- 存在缺失值、非数值、非整数、负值、重复基因、重复 barcode，或数值分布无法评估；
- 标准化伪影哨兵被触发；
- `min.cells=3` 工作 feature 空间或任何固定 QC 指标无法计算；
- 冻结 globin 定义发生变化，或完全无法识别 panel 基因；
- sample1 未复现 19,294 个工作 feature、2,684 个原研究规则通过细胞、53 个新增 `nCount` 不通过细胞、0 个新增 HB 不通过细胞和 2,631 个最终通过细胞；
- 必需的处理历史阶段或下游约束缺失；
- 任一正式 F0 脚本未进入 `F0_file_manifest.tsv`，或缺少合法的大写 SHA256；
- 任一 F0 环境锁文件缺失、版本与实际运行时不一致，或未记录合法的大写 SHA256；
- 任一必需 F0 契约表或报告没有生成。

已确认无法获得的上游信息，只有在明确标注并映射到保守的 F1 处理措施时才不阻断执行；对应 checklist 项记为 `PASS_WITH_NOTED_ISSUES`，总体报告据此记为 `PASS_WITH_NOTED_LIMITATIONS`，不能无提示地当作完全通过。

## 9. 传递给 F1 的约束

如果 F0 获批：

- F1 保留完整原始矩阵，并按主线方案F1.2独立应用冻结的单工作空间固定QC规则；
- scDblFinder 作为逐样本双细胞识别的主方法，DoubletFinder 只作敏感性分析；
- DecontX 在固定 QC 和主双细胞排除后，对 filtered raw integer counts 运行；
- ambient score 不直接删除细胞，raw counts 始终是主分析矩阵；
- 在没有 raw droplets 时，真实 barcode knee/cell calling、emptyDrops、SoupX 和 CellBender 均记为 `not_evaluable_input_limited`；
- Normal_Peritoneum 只作参考；PM 仅 3 个样本，因此样本层结果只能作方向性观察。

## 10. 预期输出

正式清单见 `reports/environment_setup/F0_expected_output_manifest.tsv`。核心输出包括：

- `data/metadata/sample_info.tsv`
- `data/metadata/data_audit.tsv`
- `data/metadata/processed_input_manifest.tsv`
- `data/metadata/F0_author_processing_audit.tsv`
- `data/metadata/F0_data_readiness_by_F_section.tsv`
- `data/metadata/F0_method_prior_decision.tsv`
- `data/metadata/decision_evidence_log.tsv`
- `results/F0_audit/F0_gate_checklist.tsv`
- `results/F0_audit/F0_global_data_reconnaissance_report.md`
- `results/F0_audit/F0_execution_report.md`
- `logs/F0_setup/analysis_log.md`

这些文件分别回答“输入是什么、数据经历过什么、哪些信息仍未知、后续方法能否使用、F0 是否允许进入 F1”等问题。任何关键结论都必须能够追溯到相应表格、脚本和日志。

## 11. 请 Claude Code 重点审核

请核对以下内容：

- F0 只测量并报告F1.2规则的影响，没有真正执行F1的细胞删除；
- 不等式边界与F1.2完全一致，包括保留`nFeature_RNA=500`和`mt_percent=20`的细胞；
- 所有细胞 QC 指标都只在每个样本应用 `min.cells=3` 后计算；
- 原始行保持不变，且没有与 QC 工作空间混淆；
- globin 采用精确匹配且全过程可审计；
- Step1 没有根据扩展名提前断言方向，Step2 能接受 `gene × cell` 并拒绝 `cell × gene`；
- sample1 回归是强制 gate，而不是仅写在说明文字中的参考结果；
- 作者报告的处理与 F0 实测影响没有混在一起；
- 十项 gate、输出契约、校验和与暂停条件均已闭合；
- F0 的 pip/Conda 依赖规格独立于后续 F 节，并已纳入 manifest；
- 笔记本资源估算仍然可信；
- 没有任何脚本自动启动 F1，也没有使用 MLMOD 或预后信息参与 F0 判断。
