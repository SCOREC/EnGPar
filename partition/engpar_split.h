#ifndef ENGPAR_SPLIT
#define ENGPAR_SPLIT

#include <ngraph.h>
#include <mpi.h>
#include "engpar_input.h"
namespace engpar {

  enum SPLIT_METHOD {
    GLOBAL_PARMETIS,
    LOCAL_PARMETIS
  };

  void split(Input*, SPLIT_METHOD method);


  void expandParts(agi::Ngraph* g, MPI_Comm newComm);
}


#endif
