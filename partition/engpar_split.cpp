#include "engpar_split.h"
#include <engpar_parmetis.h>
#include <PCU.h>
#include <engpar_split_input.h>
namespace engpar {
  
  void split(Input* input, SPLIT_METHOD method) {
    SplitInput* inp = dynamic_cast<SplitInput*>(input);
    if (!inp){
      if (!PCU_Comm_Self()) {
        fprintf(stderr,"[ERROR] Incorrect input provided\n");
      }
      assert(false);
    }
    agi::Migration* plan = NULL;
    if (inp->isOriginal) {
      if (method ==GLOBAL_PARMETIS) {
        plan = EnGPar_ParMETIS(inp,PCU_Comm_Peers()*inp->split_factor);
      }
      else if (method == LOCAL_PARMETIS) {
        //Set communicator to self
        //return EnGPar_ParMETIS(input->g,inp->split_factor);
      }
      expandParts(inp->g,inp->largeComm);
    }
    else {
      plan = new agi::Migration(inp->g);
    }
    PCU_Switch_Comm(inp->largeComm);
    input->g->setOriginalOwners();
    inp->g->migrate(plan);
    delete input;
  }

  void expandParts(agi::Ngraph* g, MPI_Comm newComm) {
    int s = PCU_Comm_Peers();
    int newSelf;
    MPI_Comm_rank(newComm,&newSelf);
    int* newRanks = new int[s];
    MPI_Allgather(&newSelf,1,MPI_INT,newRanks,1,MPI_INT,PCU_Get_Comm());
    g->changeOwners(newRanks);
    delete [] newRanks;
  }
}
