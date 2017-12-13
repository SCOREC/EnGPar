#ifndef ENGPAR_SPLIT
#define ENGPAR_SPLIT

#include <ngraph.h>
#include <mpi.h>
namespace engpar {

  enum SPLIT_METHOD {
    GLOBAL_PARMETIS,
    LOCAL_PARMETIS
  };

  agi::Migration* split(agi::Ngraph* g, int split_factor, SPLIT_METHOD method);


  void expandParts(agi::Ngraph* g, MPI_Comm newComm);
}


#endif
