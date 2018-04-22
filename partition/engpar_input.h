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


  Input* createDiffusiveInput(agi::Ngraph* g, double step_factor);
  Input* createLocalSplitInput(agi::Ngraph* g, MPI_Comm smallComm,
                               MPI_Comm largeComm, bool isPartOfSmall,
                               int split_factor,double tolerance,
                               agi::part_t* others, agi::etype adj_type = 0);
  Input* createGlobalSplitInput(agi::Ngraph* g, MPI_Comm smallComm, MPI_Comm largeComm,
                                bool isPartOfSmall, double tolerance, agi::etype adj_type = 0);
}

#endif
