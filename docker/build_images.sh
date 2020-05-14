#! /usr/bin/env bash
#
# build_images.sh
# Copyright (C) 2019 Maikel Nadolski <nadolski@math.fu-berlin.de>
#
# Distributed under terms of the MIT license.
#

# This function builds a clang based image which containes a compiled AMReX version.
# AMReX is built from the recent development branch on git.
#
# \param[in] COMPILER_ID  the clang compiler id, one of clang[5-8]
# \param[in] AMREX_SPACEDIM 2 or 3
build_clang_image() {
  COMPILER_ID="$1"
  AMREX_SPACEDIM="$2"
  case "${COMPILER_ID}" in
  clang5)
    CLANG_URL="https://releases.llvm.org/5.0.2/clang+llvm-5.0.2-x86_64-linux-gnu-ubuntu-16.04.tar.xz"
    CLANG_VERSION=5.0
    ;;
  clang6)
    CLANG_URL="https://releases.llvm.org/6.0.1/clang+llvm-6.0.1-x86_64-linux-gnu-ubuntu-16.04.tar.xz"
    CLANG_VERSION=6.0
    ;;
  clang7)
    CLANG_URL="https://releases.llvm.org/7.0.1/clang+llvm-7.0.1-x86_64-linux-gnu-ubuntu-16.04.tar.xz"
    CLANG_VERSION=7.0
    ;;
  clang8)
    CLANG_URL="https://releases.llvm.org/8.0.0/clang+llvm-8.0.0-x86_64-linux-gnu-ubuntu-16.04.tar.xz"
    CLANG_VERSION="8"
    ;;
  esac
  docker build \
    -t "git.imp.fu-berlin.de:5000/ag-klein/finitevolumesolver/amrex:${AMREX_SPACEDIM}d_${COMPILER_ID}" \
    --build-arg CLANG_URL="${CLANG_URL}" \
    --build-arg CLANG_VERSION="${CLANG_VERSION}" \
    --build-arg AMREX_SPACEDIM="${AMREX_SPACEDIM}" \
    --build-arg CACHEBUST="$(date +%s)" \
    -f amrex-clang_base_image .
}

build_gcc_image() {
  COMPILER_ID="$1"
  AMREX_SPACEDIM="$2"
  GCC_VERSION=$(echo "${COMPILER_ID}" | sed 's/[^[0-9]*//g')
  case "${GCC_VERSION}" in
    7)
      UBUNTU_VERSION="18.04"
      ;;
    8)
      UBUNTU_VERSION="18.10"
      ;;
    9)
      UBUNTU_VERSION="19.10"
      ;;
  esac
  docker build \
    -t "git.imp.fu-berlin.de:5000/ag-klein/finitevolumesolver/amrex:${AMREX_SPACEDIM}d_${COMPILER_ID}" \
    --build-arg UBUNTU_VERSION="${UBUNTU_VERSION}" \
    --build-arg GCC_VERSION="${GCC_VERSION}" \
    --build-arg AMREX_SPACEDIM="${AMREX_SPACEDIM}" \
    --build-arg CACHEBUST="$(date +%s)" \
    -f amrex-gcc_base_image .
}

push_docker_image() {
  COMPILER_ID="$1"
  AMREX_SPACEDIM="$2"
  docker push "git.imp.fu-berlin.de:5000/ag-klein/finitevolumesolver/amrex:${AMREX_SPACEDIM}d_${COMPILER_ID}"
}

CLANG_VERSIONS=("clang8")
GCC_VERSIONS=("gcc9")
AMREX_SPACEDIMS=("2" "3")

main() {
  COMPILER_ID=$1
  CI_JOB_TOKEN=$2
  COMPILER="$(echo "${COMPILER_ID}" | sed s/[0-9]//)"
  if [[ "${COMPILER}" == "gcc" ]]; then
    for AMREX_SPACEDIM in ${AMREX_SPACEDIMS[@]}; do
      printf "[Build] Docker Image for ${COMPILER_ID} with Dimension ${AMREX_SPACEDIM}\n" && \
      build_gcc_image "${COMPILER_ID}" "${AMREX_SPACEDIM}" "${CI_JOB_TOKEN}"              && \
      printf "[Push] Docker Image for ${COMPILER_ID} with Dimension ${AMREX_SPACEDIM}\n"  && \
      push_docker_image "${COMPILER_ID}" "${AMREX_SPACEDIM}"
      ec=$?
      if [[ "${ec}" != 0 ]]; then
        echo "[ERROR] Building or Pushing the docker image failed!"
        exit "${ec}"
      fi
    done
  elif [[ "${COMPILER}" == "clang" ]]; then
    for AMREX_SPACEDIM in ${AMREX_SPACEDIMS[@]}; do
      printf "[Build] Docker Image for ${COMPILER_ID} with Dimension ${AMREX_SPACEDIM}\n" && \
      build_clang_image "${COMPILER_ID}" "${AMREX_SPACEDIM}" "${CI_JOB_TOKEN}             && \
      printf "[Push] Docker Image for ${COMPILER_ID} with Dimension ${AMREX_SPACEDIM}\n"  && \
      push_docker_image "${COMPILER_ID}" "${AMREX_SPACEDIM}"
      ec=$?
      if [[ "${ec}" != 0 ]]; then
        echo "[ERROR] Building or Pushing the docker image failed!"
        exit "${ec}"
      fi
    done
  else 
    printf "Unknown COMPILER_ID '${COMPILER_ID}'.\n"
    exit 1;
  fi
}

if [[ "" != $1 ]]; then
  main "$1"
else
  >&2 echo "Missing an argument to '$0'"
  >&2 echo "Usage: $0 <COMPILER_ID>"
fi
