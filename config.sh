cmake .. \
    -DCMAKE_C_COMPILER="mpicc" \
    -DCMAKE_C_FLAGS="-g" \
    -DCMAKE_CXX_COMPILER="mpicxx" \
    -DCMAKE_CXX_FLAGS="-g -std=c++11 -Wl,--no-as-needed -ldl -pthread" \
    -DENABLE_ZOLTAN=OFF \
    -DSCOREC_PREFIX=~/Documents/scorec/core/install \
    -DENABLE_KOKKOS=ON \
    -DKOKKOS_PREFIX=~/Documents/scorec/kokkos_install \
    -DIS_TESTING=OFF \
    -DMESHES=/path/to/meshes/ \
    -DGRAPHS=/path/to/graphs/

