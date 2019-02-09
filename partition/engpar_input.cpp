#include "engpar_input.h"
#include <engpar_split_input.h>
#include <engpar_diffusive_input.h>
#include <engpar_weight_input.h>
#include <engpar_coloring_input.h>
#include <PCU.h>
namespace engpar {
  
  DiffusiveInput* createDiffusiveInput(agi::Ngraph* g,double f) {
    DiffusiveInput* input = new DiffusiveInput(g);
    input->step_factor = f;
    return input;
  }

  SplitInput* createLocalSplitInput(agi::Ngraph* g, MPI_Comm smallComm, MPI_Comm largeComm,
                               bool isPartOfSmall, int split_factor,double tolerance,
                               agi::part_t* others,agi::etype adj_type) {
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
  SplitInput* createGlobalSplitInput(agi::Ngraph* g, MPI_Comm smallComm, MPI_Comm largeComm,
                                bool isPartOfSmall, double tolerance,
                                agi::etype adj_type) {
    SplitInput* input = new SplitInput(g);
    input->smallComm=smallComm;
    input->largeComm=largeComm;
    input->isOriginal = isPartOfSmall;
    input->split_factor=-1;
    int largeSize;
    MPI_Comm_size(largeComm,&largeSize);
    input->total_parts = largeSize;
    input->tolerance = tolerance;
    input->edge_type = adj_type;
    input->other_ranks = NULL;
    PCU_Switch_Comm(smallComm);
    return input;
  }

  WeightInput* createWeightInput(agi::Ngraph* g, double t, double sf, agi::etype et) {
    return new WeightInput(g,t,sf,et);
  }

  ColoringInput* createColoringInput(agi::Ngraph* g, agi::lid_t primaryType, bool vtxColoring) {
    return new ColoringInput(g,primaryType,vtxColoring);
  }
}
