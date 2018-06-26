#!/bin/bash -x

#load system modules
#source /etc/profile
export SPACK_ROOT=/lore/cwsmith/software/spack
export PATH=$SPACK_ROOT/bin:$PATH
source $SPACK_ROOT/share/spack/setup-env.sh
source /etc/profile.d/modules.sh
source /etc/profile
export MODULEPATH=$MODULEPATH:/opt/scorec/modules
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/lore/cwsmith/develop/trilinos/install-kokkos

module load pumi

cd /lore/diamog/cdash/repos/EnGPar

#update this repo
git pull

#update the mesh and graph repos
git submodule update

#cdash output root
cd /lore/diamog/cdash

#remove old compilation
rm -rf build/

#run nightly test script
ctest -VV -D Nightly -S /lore/diamog/cdash/repos/EnGPar/cdash/nightly.cmake

#EnGPar repository built by nightly.cmake
buildDir=/lore/diamog/cdash/build/master/
[ ! -e ${buildDir} ] && exit 0

cd $buildDir
make doc
if [ -d "$PWD/doc/html" ]; then
    docsdir=/net/web/engpar/
    rm -r $docsdir/*
    cp -r doc/html/* $docsdir
fi
