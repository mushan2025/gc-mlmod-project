# 本地优先与临时服务器执行工作流

本项目不再使用网页生信agent。主线执行采用：

**Codex执行输出 → Claude Code独立审核 → 用户最终批准**

默认在本地执行，设备优先级固定为：`laptop_thinkbook16p`（默认）→ `desktop_i5_14600kf`（备用）→ 临时Linux服务器。服务器不是固定环境，只在两条本地路径均不适合、备用设备当前不可用且等待不合理，或用户基于完整评估明确批准时临时租用，并在结果完整同步到本地后删除。

## Step 1：代码、输入和资源评估

每个F节或gate执行前，Codex先在本地仓库写出并更新：

- 真实运行脚本、输入文件清单和预期输出。
- 参数、seed、版本记录方式以及失败/暂停条件。
- `reports/environment_setup/F{N}_execution_resource_assessment.tsv`：逐analysis_module先评估默认笔记本，再评估备用台式机，最后判断是否需要服务器。
- 仅在能提供额外决策价值且可安全运行时制定本地pilot；若静态资源估算已明确超出本地70%安全线，可直接提出服务器配置。纯技术smoke test与资源pilot分开记录。
- 若本地不可行，给出临时服务器最低配置、推荐配置、预计运行时间、临时磁盘和数据传输量。

资源评估至少记录：F节、gate、分析模块、输入规模、candidate_machine、machine_priority_rank、os_or_runtime（Windows原生/WSL/容器/Linux）、当前可用RAM/磁盘、CPU线程、GPU及显存、预计峰值RAM、预计临时磁盘、GPU是否被工具支持、pilot结果、预计运行时间、local_feasible、fallback_trigger、阻断原因、服务器最低/推荐CPU-RAM-磁盘-GPU、选定环境和用户审批状态。WSL必须按宿主设备记录，不得作为第三台独立设备或额外资源来源。

本地可行的判断不能只看程序能否启动。计划峰值RAM和临时磁盘必须保留安全余量；默认不应占用超过执行时可用RAM或剩余磁盘的70%。pilot峰值RAM必须优先使用OS级peak RSS监测，例如WSL/Linux的`/usr/bin/time -v`或Windows等价系统监控，并同时记录临时磁盘峰值；R的`gc()`、Python对象大小或任务管理器瞬时截图只能作为辅助信息，不能替代峰值测量。若预计超过安全线、pilot发生OOM/反复失败、依赖与本地系统不兼容，或预计运行时间明显不合理，必须暂停并重新评估。运行时间较长本身不自动否定本地执行；是否改用服务器由用户根据评估决定。

Claude Code审核代码、输入和资源计划；用户批准后进入Step 2。

## Step 2：本地基线环境

首次执行前，Codex在实际参与分析的本地设备上完成；默认先配置笔记本，只有计划使用备用台式机时才要求其达到ready状态：

- 系统、CPU、物理/逻辑核心、RAM、GPU/显存和磁盘检查。
- R、Python、Bash/WSL2或容器运行时检查；Linux-only依赖优先在经验证的WSL环境中运行。
- 项目目录、读写权限和临时目录检查。
- `renv.lock`、Python锁文件、系统依赖和外部工具版本核对。
- 小规模smoke test和随机种子可复现性检查。

Windows原生能够稳定支持且结果可审计的模块可以原生运行；需要Linux依赖、编译链或服务器迁移一致性的模块优先使用WSL2 Ubuntu 22.04或同一Linux容器。WSL使用宿主机资源，必须实测其可用RAM、CPU限制、虚拟磁盘剩余空间和跨文件系统I/O，不得按宿主机标称资源直接判定可行。若使用Windows原生R/Python，必须记录OS、编译器、BLAS、系统库和包构建差异；关键结果若迁移到Linux运行，应做方向、维度和关键统计量一致性检查。

基线信息写入 `environment/execution_environment_inventory.tsv`、`environment/environment_lock_manifest.tsv`、`data/metadata/software_versions.tsv` 和 `environment/random_seed_registry.tsv`。

## Step 3：执行路径

### Step 3A：本地执行

若资源评估以及适用时的pilot通过，Codex优先在默认笔记本执行当前gate；笔记本不适合时，须记录fallback_trigger并完成备用台式机评估后方可切换。获批本地执行后保存：

- 真实执行脚本和命令日志。
- 参数、seed和软件/数据库版本。
- 关键中间对象、结果、图和source data。
- 输入/输出manifest和SHA256。
- gate报告及资源使用实测值。

不得因另一台本地设备更方便而无记录切换设备。跨本地设备续跑时，必须先核对锁文件、输入SHA256和上游对象SHA256。

### Step 3B：临时服务器执行

只有默认笔记本和备用台式机均为 `local_feasible = false`、备用设备当前不可用且等待不合理，或用户基于两台设备资源评估明确批准使用服务器时，才能租用服务器。每次服务器均按全新环境处理：

1. 依据本gate批准的最低/推荐规格创建服务器。
2. 从版本化bootstrap脚本、锁文件或固定容器镜像重新配置环境。
3. 上传当前gate所需最小输入集，并核对本地/服务器SHA256。
4. 输出服务器系统、软件、包、数据库和容器镜像版本，与本地基线逐项比较。
5. 环境差异未解释或关键版本不一致时，不得运行正式分析。
6. 执行获批代码并记录真实命令、运行时间、峰值RAM、CPU/GPU使用和临时磁盘。

不得在服务器上临时升级包、替换数据库、修改seed或改变算法后继续运行；任何变更必须回写锁文件/参数文件并重新审批。

## Step 4：结果同步、审核和服务器删除

每个gate完成后，Codex必须把以下内容保存或同步到本地：

- 真实执行脚本、命令日志和错误日志。
- 关键中间对象、结果文件和图表source data。
- 参数文件、`data/metadata/software_versions.tsv`、`environment/random_seed_registry.tsv`和环境差异报告。
- 输入/输出manifest、文件大小和SHA256。
- gate报告、结果解读和实测资源使用记录。

服务器结果下载后必须进行本地文件数量、大小和SHA256校验，并生成 `reports/environment_setup/F{N}_{gate_id}_server_to_local_sync_manifest.tsv`。只要存在missing、size_mismatch或sha256_mismatch，就禁止删除服务器。只有同步清单完整、校验通过且没有仅存在于服务器内存或磁盘的关键产物时，才允许删除服务器。服务器删除不免除后续Claude Code审核和用户批准；审核发现问题时，应从本地保存的输入、脚本、锁文件和日志重建新服务器，不得假装旧服务器仍可恢复。

Claude Code审核执行结果；用户批准后才允许进入下一gate或下一F节。

## 不允许

- 未做资源评估就默认租服务器，或为节省服务器费用而强行在本地运行会OOM的任务。
- 不经用户批准连续执行多个gate或F节。
- 用事后整理脚本替代真实执行脚本和日志。
- 只记录包名而不记录版本、来源、系统依赖和seed。
- 在不同设备使用不同包/数据库版本而不报告。
- 服务器删除前未完成本地同步和checksum验证。
- 用agent判断替代数据、代码或文献证据。
- 因结果不显著而反复调参直到显著。
- 将大型数据和大型二进制对象直接提交到git。
