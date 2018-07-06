#ifndef ENGPAR_H__
#define ENGPAR_H__

#include <ngraph.h>
#include <engpar_diffusive_input.h>
#include <engpar_weight_input.h>
#include <engpar_split_input.h>
#include "engpar_types.h"
namespace engpar {

  /** \name Partition Routines */
  ///@{

  void split(Input*, SPLIT_METHOD method);

  ///@}
  
  /** \name Load Balancers */  
  /** \brief Performs diffusive load balancing based on the input parameters provided
   * \param input A list of input parameters provided by the user
   * \param verbosity The level of output.
   */
  void balance(Input* input, int verbosity=0);

  /** \brief Performs diffusive load balancing where vertices send weights to neighboring parts
   * \param input A collection of input parameters provided by the user.
   * \param verbosity The level of output.
   *
   * This balancer is designed with the assumption that all edges are across part 
   * boundaries and there are no interpart connections via the edge type being used to balance
   */
  void balanceWeights(WeightInput* input, int verbosity=0);

  ///@}

  /** \brief Prints an evaluation of the partition of the Ngraph
   * \param g The graph
   */
  void evaluatePartition(agi::Ngraph* g, const char* prefix = "");
}


#endif
