#!/bin/bash -e
kokkos_path=$WORK/kokkos
[ $# -ne 2 ] && echo "Usage: $0 <opt=[Debug|Release]> <cuda=[ON|OFF]>" && exit 1
opt=$1
cuda=$2
[[ $cuda != "ON" && $cuda != "OFF" ]] && echo "cuda option $cuda is not valid" && exit 1
[[ $opt != "Debug" && $opt != "Release" ]] && echo "opt option $opt is not valid" && exit 1
set -x

cudaopts=""
if [ $cuda == "ON" ]; then 
  cudaopts="-DCMAKE_CXX_COMPILER=$kokkos_path/bin/nvcc_wrapper -DKOKKOS_ENABLE_CUDA=ON -DKOKKOS_GPU_ARCH=Pascal60"
fi
echo $cudaopts

export CC=mpicc
export CXX=mpicxx
export NVCC_WRAPPER_DEFAULT_COMPILER=$CXX

set -x
cmake ../EnGPar \
    -DCMAKE_BUILD_TYPE=$opt \
    -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DMPI_CXX_COMPILER=$CXX \
    -DMPI_C_COMPILER=$CC \
    $cudaopts \
    -DENABLE_ZOLTAN=OFF \
    -DENABLE_KOKKOS=ON \
    -DKOKKOS_HOST_ARCH=Power8 \
    -DKOKKOS_PATH=$kokkos_path \
    -DSCOREC_PREFIX=$WORK/core-build-ray-gnu/install \
    -DIS_TESTING=ON \
    -DMESHES=$HOME/develop/pumi-meshes \
    -DGRAPHS=$HOME/develop/EnGPar/EnGPar-graphs 

