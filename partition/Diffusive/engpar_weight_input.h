#ifndef __ENGPAR_WEIGHT_INPUT_H__
#define __ENGPAR_WEIGHT_INPUT_H__

#include <ngraph.h>
#include <engpar_input.h>
namespace engpar {
  
  class WeightInput : public Input {
  public:
    WeightInput(agi::Ngraph*,double tolerance,double sf = 0.1, agi::etype et=0);

    /** \brief The maximum iterations for all load balancing */
    int maxIterations;

    /** \brief The target weight imbalance of the graph vertices */
    double tol;
    
    /** \brief The percent of difference in weight to send in each iteration 
     *
     * defaults to .1 
     */
    double step_factor;

    /** \brief The edge type used throughout the weight diffuser.
     *
     * defaults to 0
     */
    int primary_edge_type;

  };
}

#endif

