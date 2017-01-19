
cmake .. \
    -DCMAKE_C_COMPILER="mpicc" \
    -DCMAKE_C_FLAGS="-g" \
    -DCMAKE_CXX_COMPILER="mpicxx" \
    -DCMAKE_CXX_FLAGS="-g -std=c++11" \
    -DENABLE_ZOLTAN=ON \
    -DSCOREC_PREFIX=/lore/diamog/core/installJessie \
    -DIS_TESTING=ON \
    -DMESHES=/users/diamog/meshes/ \
    -DGRAPHS=/users/diamog/graphs/

