#!/bin/bash -x

source /usr/local/etc/bash_profile
module load cmake/latest
module load pumi

#cdash output root
cd /lore/diamog/cdash

#remove old compilation
rm -rf build/

#run nightly test script
ctest -VV -D Nightly -S /lore/diamog/cdash/repos/EnGPar/cdash/nightly.cmake &> /lore/diamog/cdash/repos/EnGPar/cdash/cmake_log.txt

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
