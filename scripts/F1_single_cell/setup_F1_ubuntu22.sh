#!/usr/bin/env bash
set -Eeuo pipefail

# Ubuntu 22.04上的F1环境安装脚本。
# R 4.4.3源码包及其.sha256文件由本地上传，避免服务器重复下载项目输入。

if [[ $# -ne 2 ]]; then
  echo "用法: $0 /path/to/R-4.4.3.tar.gz /path/to/project" >&2
  exit 2
fi

R_TARBALL="$(readlink -f "$1")"
PROJECT_ROOT="$(readlink -f "$2")"
R_VERSION="4.4.3"
R_PREFIX="/opt/R/${R_VERSION}"
R_SOURCE_DIR="/opt/src/R-${R_VERSION}"
BUILD_THREADS="${F1_BUILD_THREADS:-24}"

if [[ ! -f "${R_TARBALL}" || ! -f "${R_TARBALL}.sha256" ]]; then
  echo "缺少R源码包或配套.sha256文件。" >&2
  exit 3
fi
if [[ ! -f "${PROJECT_ROOT}/scripts/F1_single_cell/setup_F1_packages.R" ]]; then
  echo "项目目录中缺少F1包安装脚本。" >&2
  exit 4
fi

expected_sha="$(awk '{print toupper($1)}' "${R_TARBALL}.sha256")"
observed_sha="$(sha256sum "${R_TARBALL}" | awk '{print toupper($1)}')"
if [[ "${expected_sha}" != "${observed_sha}" ]]; then
  echo "R源码包SHA256不一致。" >&2
  exit 5
fi

export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
  build-essential gfortran cmake pkg-config texinfo ca-certificates curl git \
  pigz jags \
  libreadline-dev libx11-dev libxt-dev libcairo2-dev libpng-dev libjpeg-dev \
  libtiff5-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev \
  libfribidi-dev libicu-dev libpcre2-dev zlib1g-dev libbz2-dev liblzma-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev libopenblas-dev liblapack-dev \
  libgit2-dev libglpk-dev libgsl-dev libfftw3-dev libhdf5-dev \
  libudunits2-dev libgeos-dev libproj-dev libgdal-dev libmagick++-dev \
  libprotobuf-dev protobuf-compiler libjq-dev libzmq3-dev libnlopt-dev \
  libgmp3-dev libmpfr-dev libtbb-dev libsqlite3-dev

if [[ ! -x "${R_PREFIX}/bin/Rscript" ]]; then
  mkdir -p /opt/src
  if [[ ! -d "${R_SOURCE_DIR}" ]]; then
    tar -xzf "${R_TARBALL}" -C /opt/src
  fi
  cd "${R_SOURCE_DIR}"
  ./configure \
    --prefix="${R_PREFIX}" \
    --enable-R-shlib \
    --with-blas \
    --with-lapack \
    --with-x=no
  make -j"${BUILD_THREADS}"
  make install
fi

ln -sfn "${R_PREFIX}/bin/R" /usr/local/bin/R
ln -sfn "${R_PREFIX}/bin/Rscript" /usr/local/bin/Rscript

export F1_INSTALL_NCPUS="${F1_INSTALL_NCPUS:-12}"
"${R_PREFIX}/bin/Rscript" \
  "${PROJECT_ROOT}/scripts/F1_single_cell/setup_F1_packages.R" \
  "--project-root=${PROJECT_ROOT}"

echo "F1 Ubuntu环境安装完成。"
