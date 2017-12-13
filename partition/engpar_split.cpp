#include "engpar_split.h"
#include <engpar_parmetis.h>
#include <PCU.h>
namespace engpar {

  agi::Migration* split(agi::Ngraph* g, int split_factor,SPLIT_METHOD method) {
    g->setOriginalOwners();
    if (method ==GLOBAL_PARMETIS) {
      return EnGPar_ParMETIS(g,PCU_Comm_Peers()*split_factor);
    }
    else if (method == LOCAL_PARMETIS) {
      //Set communicator to self
      return EnGPar_ParMETIS(g,split_factor);
    }
    return NULL;
  }

  void expandParts(agi::Ngraph* g, MPI_Comm newComm) {
    int s = PCU_Comm_Peers();
    int newSelf;
    MPI_Comm_rank(newComm,&newSelf);
    int* newRanks = new int[s];
    MPI_Allgather(&newSelf,1,MPI_INT,newRanks,1,MPI_INT,PCU_Get_Comm());
    g->changeOwners(newRanks);
    delete [] newRanks;
    PCU_Switch_Comm(newComm);
  }
}
