#ifndef ENGPAR_H__
#define ENGPAR_H__

#include <ngraph.h>
#include "Diffusive/engpar_diffusive_input.h"
#include "Diffusive/engpar_weight_input.h"
#include "Multilevel/engpar_split_input.h"
#include "engpar_types.h"
namespace engpar {

  /** \name Partition Routines */
  ///@{
  enum SPLIT_METHOD {
    GLOBAL_PARMETIS,
    LOCAL_PARMETIS
  };

  void split(Input*, SPLIT_METHOD method);

  ///@}
  
  /** \name Load Balancers */
  ///@{
  /** \brief Performs diffusive load balancing that targets the vertices of the graph
   * \param g The graph to be balanced.
   * \param tolerance The goal vertex imbalance to reach.
   * \param stepFactor Control of how much weight is migrated in each iteration.
   * \param verbosity The level of output.
   *
   * Good for initial testing. Use balance for getting results.
   */
  void balanceVertices(agi::Ngraph*& g, double tolerance, double stepFactor=0.1, int verbosity=0);
  
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

  
  /** \name Metrics */
  ///@{
  /** \brief Calculates the imbalance across all processes of the provided weight
   * \param w The weight of this process
   * \return The imbalance
   */
  double EnGPar_Get_Imbalance(wgt_t w);
  /** \brief Prints an evaluation of the partition of the Ngraph
   * \param g The graph
   */
  void evaluatePartition(agi::Ngraph* g, const char* prefix = "");

  wgt_t getWeight(agi::Ngraph*,int,bool countGhosts=false);
  
  void printImbalances(agi::Ngraph* g);

  ///@}


  
}


#endif
