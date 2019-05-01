#!/bin/bash -x

#load system modules
source /etc/profile.d/modules.sh
source /etc/profile

module load gcc cmake mpich
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

#build documentation
make doc
if [ -d "$PWD/doc/html" ]; then
    docsdir=/net/web/engpar/
    rm -r $docsdir/*
    cp -r doc/html/* $docsdir
fi
