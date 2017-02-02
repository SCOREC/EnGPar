cmake .. \
    -DCMAKE_C_COMPILER="mpicc" \
    -DCMAKE_C_FLAGS="-g" \
    -DCMAKE_CXX_COMPILER="mpicxx" \
    -DCMAKE_CXX_FLAGS="-g -std=c++11 -Wl,--no-as-needed -ldl -pthread" \
    -DENABLE_ZOLTAN=ON \
    -DSCOREC_PREFIX=/path/to/core/install \
    -DENABLE_KOKKOS=ON \
    -DKOKKOS_PREFIX=/path/to/kokkos_install \
    -DIS_TESTING=ON \
    -DMESHES=/path/to/meshes/ \
    -DGRAPHS=/path/to/graphs/

