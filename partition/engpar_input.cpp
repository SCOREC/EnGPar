#include "engpar_input.h"
#include <engpar_split_input.h>
#include <engpar_diffusive_input.h>
#include <PCU.h>
namespace engpar {
  
  Input* createDiffusiveInput(agi::Ngraph* g,double f) {
    DiffusiveInput* input = new DiffusiveInput(g);
    input->step_factor = f;
    return input;
  }

  Input* createSplitInput(agi::Ngraph* g, MPI_Comm smallComm, MPI_Comm largeComm,
                          bool isPartOfSmall, int split_factor,double tolerance,
                          agi::etype adj_type, agi::part_t* others) {
    SplitInput* input = new SplitInput(g);
    input->smallComm=smallComm;
    input->largeComm=largeComm;
    input->isOriginal = isPartOfSmall;
    input->split_factor = split_factor;
    int largeSize;
    MPI_Comm_size(largeComm,&largeSize);
    input->total_parts = largeSize;
    input->tolerance = tolerance;
    input->edge_type = adj_type;
    input->other_ranks = others;
    PCU_Switch_Comm(smallComm);
    return input;
  }
  Input* createSplitInput(agi::Ngraph* g, MPI_Comm smallComm, MPI_Comm largeComm,
                          bool isPartOfSmall, double tolerance,
                          agi::etype adj_type, agi::part_t* others) {
    SplitInput* input = new SplitInput(g);
    input->smallComm=smallComm;
    input->largeComm=largeComm;
    input->isOriginal = isPartOfSmall;
    int smallSize;
    int largeSize;
    MPI_Comm_size(smallComm,&smallSize);
    MPI_Comm_size(largeComm,&largeSize);
    input->split_factor = -1;
    if (largeSize%smallSize==0)
      input->split_factor = largeSize/smallSize;
    input->total_parts = largeSize;
    input->tolerance = tolerance;
    input->edge_type = adj_type;
    input->other_ranks = others;
    PCU_Switch_Comm(smallComm);
    return input;
  }

}
