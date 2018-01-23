#ifndef ENGPAR_H
#define ENGPAR_H
#include <agiBalancer.h>
#include <ngraph.h>
#include <vector>
#include "engpar_diffusive_input.h"

/** \file engpar.h
    \brief Diffusive Load Balancing APIs
*/
namespace engpar {
  // \cond
  typedef agi::wgt_t wgt_t;
  typedef agi::part_t part_t;
  typedef std::vector<agi::GraphEdge*> Queue;

  wgt_t getWeight(agi::Ngraph*,int);
  // \endcond
  /** \name Load Balancers */
  ///@{
  /** \brief Makes a diffusive load balancer that targets the vertices of the graph
   * \param g The graph to be balanced
   * \param stepFactor Control of how much weight is migrated in each iteration
   * \param verbosity The level of output.
   * \return The balancer, use balancer->run to run the operations
   */
  agi::Balancer* makeVtxBalancer(agi::Ngraph*& g, double stepFactor=0.1,
                                        int verbosity=0);
  /** \brief Makes a diffusive load balancer based on the input parameters provided
   * \param input A list of input parameters provided by the user
   * \param verbosity The level of output.
   * \return The balancer, use balancer->run to run the operations
   */
  agi::Balancer* makeBalancer(Input* input,int verbosity=0);

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
  void evaluatePartition(agi::Ngraph* g);

  void printImbalances(agi::Ngraph* g);
  ///@}
}
#endif
