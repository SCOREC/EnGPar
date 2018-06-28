#ifndef ENGPAR_METRICS_H__
#define ENGPAR_METRICS_H__

#include <ngraph.h>
#include "../engpar_types.h"
namespace engpar {
  
  /** \name Metrics */
  ///@{
  /** \brief Calculates the imbalance across all processes of the provided weight
   * \param w The weight of this process
   * \return The imbalance
   */
  double EnGPar_Get_Imbalance(wgt_t w);

  wgt_t getWeight(agi::Ngraph*,int,bool countGhosts=false);
  
  void printImbalances(agi::Ngraph* g);
  ///@}  
}


#endif
