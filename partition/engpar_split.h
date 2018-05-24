#ifndef ENGPAR_SPLIT
#define ENGPAR_SPLIT

#include <ngraph.h>
#include <mpi.h>
#include "engpar_input.h"
#include "engpar_split_input.h"
namespace engpar {

  void expandParts(agi::Ngraph* g, MPI_Comm newComm);
}


#endif
