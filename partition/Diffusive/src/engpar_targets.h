#ifndef ENGPAR_TARGETS_H
#define ENGPAR_TARGETS_H

#include <engpar_metrics.h>
#include "engpar_container.h"
#include "engpar_weights.h"
#include "engpar_diffusive_input.h"
namespace engpar {

  class Targets : public Container<wgt_t> {
  public:
    Targets(bool isHyper, double step_factor, Sides* s, Weights* targetW, int sideTol,
            Weights** completedWs,std::vector<wgt_t>& completedTolerances) {
      Sides::iterator itr;
      //For each side
      for (itr = s->begin();itr!=s->end();itr++) {
        int neighbor = itr->first;
        bool canSend=!isHyper || s->get(neighbor)<=sideTol || PCU_Comm_Peers() == 2;
        //Check if we can send to this neighbor based on previously balanced edges/vtxs
        for (unsigned int i=0;i<completedTolerances.size();i++) {
          if (completedWs[i]->get(neighbor)>=completedTolerances[i])
            canSend=false;
        }
        if (!canSend)
          continue;

        //See if we need to send weight for the target type and add if we do
        wgt_t myW = targetW->myWeight();
        wgt_t neighborW = targetW->get(neighbor);
        if (myW>neighborW) {
          wgt_t diff = myW-neighborW;
          wgt_t sideFraction = itr->second;
          sideFraction /= s->total();
          wgt_t scaledW = diff * sideFraction * step_factor;
          set(neighbor,scaledW);
        }
      }
    }

    Targets(bool isHyper, double, Sides* s, agi::WeightPartitionMap* wp_map,
            int sideTol, Weights** completedWs,std::vector<wgt_t>& completedTolerances) {
      Sides::iterator itr;
      for (itr = s->begin();itr!=s->end();itr++) {
        int neighbor = itr->first;
        bool canSend=!isHyper || s->get(neighbor)<=sideTol;
        for (unsigned int i=0;i<completedTolerances.size();i++) {
          if (completedWs[i]->get(neighbor)>=completedTolerances[i])
            canSend=false;
        }
        if (!canSend)
          continue;
        agi::WeightPartitionMap::iterator itr;
        for (itr = wp_map->begin(); itr != wp_map->end(); itr++) {
          std::unordered_map<agi::gid_t,agi::wgt_t>::iterator inner_itr;
          for (inner_itr = itr->second.begin(); inner_itr != itr->second.end(); inner_itr++) {
            part_t neighbor = inner_itr->first;
            wgt_t diff = inner_itr->second;
            wgt_t sideFraction = s->get(neighbor);
            sideFraction /= s->total();
            wgt_t scaledW = diff * sideFraction;
            set(neighbor,scaledW);
          }
        }
      }
    }
  };
  Targets* makeTargets(bool isHyper, double step_factor, Sides* s, Weights* tW,int sT,
                       Weights** cWs,std::vector<wgt_t>& cTs);

  Targets* makeTargets(DiffusiveInput* in, Sides* s, Weights* tW,int sT,
                       Weights** cWs,std::vector<wgt_t>& cTs);

  Targets* makePartWeightTargets(DiffusiveInput* in, Sides* s, agi::WeightPartitionMap* wp_map,
                                 int sT, Weights** cWs, std::vector<wgt_t>& cTs);
}

#endif
