#!/usr/bin/env bash
set -Eeuo pipefail

# F1服务器正式入口。默认只跑F1.1-F1.3，随后由研究者审核聚类和marker。

if [[ $# -lt 1 || $# -gt 4 ]]; then
  echo "用法: $0 PROJECT_ROOT [FROM] [TO] [--approve-cnv-execution]" >&2
  exit 2
fi

PROJECT_ROOT="$(readlink -f "$1")"
FROM_STAGE="${2:-F1.1}"
TO_STAGE="${3:-F1.3}"
CNV_APPROVAL="${4:-}"
R_SCRIPT="${F1_RSCRIPT:-/opt/R/4.4.3/bin/Rscript}"

if [[ ! -x "${R_SCRIPT}" ]]; then
  echo "找不到可执行Rscript：${R_SCRIPT}" >&2
  exit 3
fi

# 96 GiB / 48 vCPU服务器的保守并行配置：可并行处加速，同时避免嵌套线程超过48。
export F1_FUTURE_WORKERS="${F1_FUTURE_WORKERS:-12}"
export F1_FUTURE_GLOBALS_MAX_GB="${F1_FUTURE_GLOBALS_MAX_GB:-72}"
export F1_SCDBLFINDER_WORKERS="${F1_SCDBLFINDER_WORKERS:-8}"
export F1_INFERCNV_THREADS="${F1_INFERCNV_THREADS:-16}"
export F1_COPYKAT_CORES="${F1_COPYKAT_CORES:-16}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-4}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-4}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-4}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-4}"

ulimit -n 65535
mkdir -p "${PROJECT_ROOT}/logs/F1_single_cell"
timestamp="$(date +%Y%m%d_%H%M%S)"
log_file="${PROJECT_ROOT}/logs/F1_single_cell/server_${FROM_STAGE}_to_${TO_STAGE}_${timestamp}.log"

args=(
  "${PROJECT_ROOT}/scripts/F1_single_cell/run_F1.R"
  "--execute"
  "--project-root=${PROJECT_ROOT}"
  "--from=${FROM_STAGE}"
  "--to=${TO_STAGE}"
)
if [[ "${CNV_APPROVAL}" == "--approve-cnv-execution" ]]; then
  args+=("${CNV_APPROVAL}")
elif [[ -n "${CNV_APPROVAL}" ]]; then
  echo "第四个参数只能是--approve-cnv-execution。" >&2
  exit 4
fi

{
  echo "started_at=$(date --iso-8601=seconds)"
  echo "hostname=$(hostname)"
  echo "from=${FROM_STAGE}"
  echo "to=${TO_STAGE}"
  echo "future_workers=${F1_FUTURE_WORKERS}"
  echo "scDblFinder_workers=${F1_SCDBLFINDER_WORKERS}"
  echo "inferCNV_threads=${F1_INFERCNV_THREADS}"
  echo "CopyKAT_cores=${F1_COPYKAT_CORES}"
  echo "BLAS_threads=${OPENBLAS_NUM_THREADS}"
  "${R_SCRIPT}" "${args[@]}"
  echo "finished_at=$(date --iso-8601=seconds)"
} 2>&1 | tee "${log_file}"
