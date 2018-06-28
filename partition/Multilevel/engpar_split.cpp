#include <engpar_parmetis.h>
#include <PCU.h>
#include <engpar_metrics.h>
#include <engpar_support.h>
#include <agiMigration.h>
#include "engpar_split_input.h"
#include "../engpar_types.h"
namespace engpar {

  void expandParts(agi::Ngraph* g, MPI_Comm newComm) {
    int s = PCU_Comm_Peers();
    int newSelf;
    MPI_Comm_rank(newComm,&newSelf);
    int* newRanks = new int[s];
    MPI_Allgather(&newSelf,1,MPI_INT,newRanks,1,MPI_INT,PCU_Get_Comm());
    g->changeOwners(newRanks);
    delete [] newRanks;
  }

  void split(Input* input, SPLIT_METHOD method) {
    SplitInput* inp = dynamic_cast<SplitInput*>(input);
    if (!inp){
      if (!PCU_Comm_Self()) {
        EnGPar_Error_Message("Incorrect input provided\n");
      }
      assert(false);
    }
    agi::Migration* plan = NULL;
    inp->self = PCU_Comm_Self();
    if (inp->isOriginal) {
      if (method ==GLOBAL_PARMETIS) {
        plan = EnGPar_ParMETIS(inp,inp->total_parts,false);
      }
      else if (method == LOCAL_PARMETIS) {
        assert(inp->split_factor!=-1);
        assert(inp->other_ranks);
        //Set communicator to self
        PCU_Switch_Comm(MPI_COMM_SELF);
        plan = EnGPar_ParMETIS(inp,inp->split_factor,true);
        agi::Migration::iterator itr;
        for (itr = plan->begin();itr!=plan->end();itr++) {
          plan->insert(std::make_pair(*itr,inp->other_ranks[plan->get(*itr)]));
        }
        PCU_Switch_Comm(inp->smallComm);
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
}
