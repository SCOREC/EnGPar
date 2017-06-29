#!/bin/bash -x

module load cmake/latest
module load pumi

#cdash output root
cd /lore/diamog/cdash

#remove old compilation
rm -rf build/

#run nightly test script
ctest -VV -D Nightly -S /lore/diamog/cdash/repos/EnGPar/cdash/nightly.cmake &> /lore/diamog/cdash/repos/EnGPar/cdash/cmake_log.txt

