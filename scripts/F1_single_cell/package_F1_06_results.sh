#!/usr/bin/env bash
set -Eeuo pipefail

# 将F1.6正式结果压成一个zstd归档，便于从临时服务器下载并校验。
# 用法：bash scripts/F1_single_cell/package_F1_06_results.sh PROJECT_ROOT [ARCHIVE_DIR]

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "用法: $0 PROJECT_ROOT [ARCHIVE_DIR]" >&2
  exit 2
fi

PROJECT_ROOT="$(readlink -f "$1")"
ARCHIVE_DIR="${2:-${PROJECT_ROOT}/server_transfer_packages}"
timestamp="$(date +%Y%m%d_%H%M%S)"
archive_name="F1_06_results_${timestamp}.tar.zst"
archive_path="${ARCHIVE_DIR}/${archive_name}"
checksum_path="${archive_path}.SHA256.txt"
listing_path="${archive_path}.contents.txt"

command -v tar >/dev/null 2>&1 || { echo "缺少tar。" >&2; exit 3; }
command -v zstd >/dev/null 2>&1 || { echo "缺少zstd。" >&2; exit 3; }
command -v sha256sum >/dev/null 2>&1 || { echo "缺少sha256sum。" >&2; exit 3; }

mkdir -p "${ARCHIVE_DIR}"
cd "${PROJECT_ROOT}"

# results/F1_malignancy内含完整inferCNV目录、HMM文件、CopyKAT三臂结果、审核表和作图来源。
items=(
  "results/F1_malignancy"
  "scripts/F1_single_cell/F1_06_malignancy_inference.R"
  "scripts/F1_single_cell/F1_config.R"
  "scripts/F1_single_cell/F1_utils.R"
  "scripts/F1_single_cell/validate_F1_static.R"
  "data/metadata/pipeline_parameters.yaml"
  "reports/environment_setup/F1_execution_plan_for_review.md"
  "胃癌MLMOD亚群主线研究方案.txt"
)
for object in \
  "objects/F1_single_cell/05_malignant_epithelial.rds" \
  "objects/F1_single_cell/06a_malignant_epithelial_main.rds" \
  "objects/F1_single_cell/06b_malignant_epithelial_high_confidence_only.rds"
do
  [[ -f "${object}" ]] && items+=("${object}")
done

for item in "${items[@]}"; do
  [[ -e "${item}" ]] || { echo "缺少必要打包项：${item}" >&2; exit 4; }
done

# -T0让zstd使用可用CPU；压缩级别10兼顾下载体积和服务器时间。
tar -I "zstd -T0 -10" -cf "${archive_path}" "${items[@]}"
tar --zstd -tf "${archive_path}" > "${listing_path}"
sha256sum "${archive_path}" | awk '{print toupper($1) "  " $2}' > "${checksum_path}"

echo "archive=${archive_path}"
echo "sha256=${checksum_path}"
echo "contents=${listing_path}"
du -h "${archive_path}"
