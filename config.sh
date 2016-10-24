
cmake .. \
      -DCMAKE_C_COMPILER="mpicc" \
      -DCMAKE_C_FLAGS="-g -Wall" \
      -DCMAKE_CXX_COMPILER="mpicxx" \
      -DCMAKE_CXX_FLAGS="-g -Wall" \
      -DSCOREC_PREFIX=~/Documents/scorec/core/install
