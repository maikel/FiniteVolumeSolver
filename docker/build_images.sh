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
    ;;
  clang6)
    CLANG_URL="https://releases.llvm.org/6.0.1/clang+llvm-6.0.1-x86_64-linux-gnu-ubuntu-16.04.tar.xz"
    ;;
  clang7)
    CLANG_URL="https://releases.llvm.org/7.0.1/clang+llvm-7.0.1-x86_64-linux-gnu-ubuntu-16.04.tar.xz"
    ;;
  clang8)
    CLANG_URL="https://releases.llvm.org/8.0.0/clang+llvm-8.0.0-x86_64-linux-gnu-ubuntu-16.04.tar.xz"
    ;;
  esac
  docker build \
    -t "git.imp.fu-berlin.de:5000/ag-klein/finitevolumesolver/amrex:${AMREX_SPACEDIM}d_${COMPILER_ID}" \
    --build-arg CLANG_URL="${CLANG_URL}" \
    --build-arg AMREX_SPACEDIM="${AMREX_SPACEDIM}" \
    -f amrex-clang_base_image .
}

build_gcc_image() {
  COMPILER_ID="$1"
  AMREX_SPACEDIM="$2"
  GCC_VERSION=$(echo "${COMPILER_ID}" | sed 's/[^[0-9]*//g')
  docker build \
    -t "git.imp.fu-berlin.de:5000/ag-klein/finitevolumesolver/amrex:${AMREX_SPACEDIM}d_${COMPILER_ID}" \
    --build-arg GCC_VERSION="${GCC_VERSION}" \
    --build-arg AMREX_SPACEDIM="${AMREX_SPACEDIM}" \
    -f amrex-gcc_base_image .
}

push_docker_image() {
  COMPILER_ID="$1"
  AMREX_SPACEDIM="$2"
  docker push "git.imp.fu-berlin.de:5000/ag-klein/finitevolumesolver/amrex:${AMREX_SPACEDIM}d_${COMPILER_ID}"
}

CLANG_VERSIONS=("clang5" "clang6" "clang7" "clang8")
GCC_VERSIONS=("gcc5" "gcc6" "gcc7" "gcc8")
AMREX_SPACEDIMS=("2" "3")

main() {
  for COMPILER_ID in ${CLANG_VERSIONS[@]}; do
    for AMREX_SPACEDIM in ${AMREX_SPACEDIMS[@]}; do
      build_clang_image "${COMPILER_ID}" "${AMREX_SPACEDIM}"
      push_docker_image "${COMPILER_ID}" "${AMREX_SPACEDIM}"
    done
  done

  for COMPILER_ID in ${GCC_VERSIONS[@]}; do
    for AMREX_SPACEDIM in ${AMREX_SPACEDIMS[@]}; do
      build_gcc_image "${COMPILER_ID}" "${AMREX_SPACEDIM}"
      push_docker_image "${COMPILER_ID}" "${AMREX_SPACEDIM}"
    done
  done
}

main
