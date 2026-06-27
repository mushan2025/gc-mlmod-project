# Codex修改指令：主线TXT第二轮第三方审查修订

## 背景

第三方审查者在第一轮修订落地后重新通读了全文，发现6个仍需修正的问题。请逐条修改主线TXT（`胃癌MLMOD亚群主线研究方案.txt`），修改完成后告诉我每条改了什么、改在哪一行附近。

---

## 修改清单

### 修改1：F0.7 gate检查增加normalization artifact阻断条件

F0.5新增了normalization_artifact_flag字段和"若发现normalization artifact则不能当作raw UMI count使用"的规则，但F0.7硬性通过条件中没有检查这个字段。这是一个逻辑漏洞：Biomni可能标记了normalization artifact但gate不拦截。

请在F0.7硬性通过条件中，第133行（data_audit.tsv相关检查）附近，增加一条：

> data_audit.tsv中所有include_in_f1=true的样本，normalization_artifact_flag必须为false；若任一样本normalization_artifact_flag=true且audit_decision仍为enter_full_F1，gate判定为FAIL，必须暂停讨论。

### 修改2：F0.7 gate检查明确跨文件关联方式

当前第133行写"所有include_in_f1=true的样本必须audit_decision = enter_full_F1"，但include_in_f1在sample_info.tsv中，audit_decision在data_audit.tsv中。Biomni需要跨文件join，当前是隐含的。

请在第133行的这条检查之前或之后，增加一句说明：

> 涉及sample_info.tsv与data_audit.tsv交叉检查的gate项，通过geo_accession或sample_id将两个文件关联。Biomni必须确认两个文件的样本行数一致（均为40行）且geo_accession能一一匹配，再进行字段交叉校验。

### 修改3：明确F0作为整体执行单元，内部子节不单独审批

当前第45行写"Biomni不得连续自动执行多个F小节"，同时第38行提到"execution_report_F0.1.md、execution_report_F0.2.md"，暗示F0内部子节也需要单独报告和审批。这会导致歧义。

请做以下修改：

（a）将第45行"Biomni不得连续自动执行多个F小节"修改为：

> Biomni不得连续自动执行多个F大节（例如F0完成后不得自动进入F1）。F0内部的子节（F0.1–F0.7）作为一个整体执行单元，不需要逐个子节单独审批；但Biomni仍应在F0执行过程中按子节顺序工作，遇到需要研究者批准的事项（如下载文件、分组冲突）时暂停等待。

（b）将第38行"各小节执行报告，例如execution_report_F0.1.md、execution_report_F0.2.md"修改为：

> F0整体执行报告F0_execution_report.md（位于results/F0_audit/）。不要求F0内部每个子节单独生成execution_report；但如果某个子节出现异常、降级或需要特殊说明，应在F0_execution_report.md对应段落中记录。

（c）相应地，第45行后半段关于"每个小节开始前…每个小节执行后…"的表述，请调整为以F大节为粒度：

> 每个F大节开始前，Biomni必须先提交本节执行计划，列出拟读取的输入文件、拟生成的输出文件、关键方法、关键参数、可能风险和需要研究者批准的事项；研究者批准后才可执行。每个F大节执行后，Biomni必须导出本节执行报告、自动生成的execution trace文件夹、实际执行notebook（如.ipynb）、运行日志、关键输出表、关键图、参数文件、版本文件和文件manifest/checksum。

其余关于execution trace、executed_source vs derived_summary_script、分享链接、外部审计的规则保持不变，只把粒度从"F小节"统一为"F大节"。

### 修改4：统一SHA256大小写

当前第47行GSE183904_RAW.tar的SHA256写为大写hex（BA089D1DC186...），而precheck TSV中gene_order_sha256使用小写hex。

请在第47行SHA256值之后增加一句：

> SHA256比较应不区分大小写（case-insensitive）；本方案中SHA256统一以大写记录，但Biomni校验时应以不区分大小写的方式比对。

### 修改5：指定metadata_issue多值分隔符

当前F0.4第87行列举了metadata_issue的多种可能值，但一个样本可能同时存在多个issue，未说明分隔方式。

请在第87行末尾增加：

> 若一个样本存在多个metadata问题，以英文分号分隔，例如patient_id_missing;pairing_unknown;geo_typo_peritonium。

### 修改6：指定抽样审计最低标准

当前F0.5第99行允许降级为sampled审计，但没有规定最低抽样比例，Biomni可能只抽极少量数据就声称完成。

请在第99行"抽样结果不得伪装为full-stream审计"之前，增加最低标准：

> 降级抽样审计的最低要求：每个csv.gz文件至少抽样检查10%的数据行（基因行）或2000行（取较大者）；抽样行应均匀分布在文件头部、中部和尾部，不得只检查前N行。每个文件的per-cell nCount统计（min、Q1、median、Q3、max）仍必须基于全部细胞列计算，不可抽样。

---

## 修改后请确认

1. 每条修改的具体位置（哪个小节、哪一行附近）。
2. 是否有任何修改与现有内容冲突，如果有请指出。
3. 修改后把文件同步到WPS路径。
