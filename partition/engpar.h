#ifndef ENGPAR_H__
#define ENGPAR_H__

#include <ngraph.h>
#include "Diffusive/engpar_diffusive_input.h"
#include "Multilevel/engpar_split_input.h"
#include "engpar_types.h"
#include <agiBalancer.h>
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
  /** \brief Makes a diffusive load balancer that targets the vertices of the graph
   * \param g The graph to be balanced
   * \param stepFactor Control of how much weight is migrated in each iteration
   * \param verbosity The level of output.
   * \return The balancer, use balancer->run to run the operations
   *
   * Good for initial testing. Use makeBalancer for getting results.
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
  void evaluatePartition(agi::Ngraph* g, std::string prefix = "");

  wgt_t getWeight(agi::Ngraph*,int,bool countGhosts=false);
  
  void printImbalances(agi::Ngraph* g);

  ///@}


  
}


#endif
