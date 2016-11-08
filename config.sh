
cmake .. \
      -DCMAKE_C_COMPILER="mpicc" \
      -DCMAKE_C_FLAGS="-g" \
      -DCMAKE_CXX_COMPILER="mpicxx" \
      -DCMAKE_CXX_FLAGS="-g" \
      -DSCOREC_PREFIX=~/Documents/scorec/core/install
