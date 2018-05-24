#ifndef __ENGPAR_INPUT_H__
#define __ENGPAR_INPUT_H__

#include <mpi.h>
#include <ngraph.h>
namespace engpar {


  class Input {
  public:
  Input(agi::Ngraph* g_) : g(g_) {}
    virtual ~Input() {};

    virtual void addPriority(int, double) {}
    /** \brief The graph being balanced */
    agi::Ngraph* g;
  };
  class DiffusiveInput;
  class SplitInput;

  DiffusiveInput* createDiffusiveInput(agi::Ngraph* g, double step_factor);
  SplitInput* createLocalSplitInput(agi::Ngraph* g, MPI_Comm smallComm,
                               MPI_Comm largeComm, bool isPartOfSmall,
                               int split_factor,double tolerance,
                               agi::part_t* others, agi::etype adj_type = 0);
  SplitInput* createGlobalSplitInput(agi::Ngraph* g, MPI_Comm smallComm, MPI_Comm largeComm,
                                bool isPartOfSmall, double tolerance, agi::etype adj_type = 0);
}

#endif
