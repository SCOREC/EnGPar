cmake .. \
    -DCMAKE_C_COMPILER="cc" \
    -DCMAKE_C_FLAGS="-g" \
    -DCMAKE_CXX_COMPILER="CC" \
    -DCMAKE_CXX_FLAGS="-g -std=c++11 -Wl,--no-as-needed -ldl -pthread" \
    -DENABLE_ZOLTAN=OFF \
    -DSCOREC_PREFIX=$PUMI_INSTALL_DIR \
    -DENABLE_KOKKOS=OFF \
    -DKOKKOS_PREFIX=/path/to/kokkos_install \
    -DIS_TESTING=ON \
    -DMESHES=$PWD/../pumi-meshes \
    -DGRAPHS=$PWD/../EnGPar-graphs \

