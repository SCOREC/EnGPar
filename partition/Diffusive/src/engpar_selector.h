#ifndef ENGPAR_SELECTOR_H
#define ENGPAR_SELECTOR_H

#include <engpar_metrics.h>
#include "engpar_sides.h"
#include "engpar_targets.h"
#include "engpar_queue.h"
#include "engpar_support.h"
#include <unordered_set>
namespace engpar {

  typedef std::unordered_set<agi::GraphVertex*> Cavity;
  typedef std::vector<part_t> Peers;

  typedef std::unordered_map<part_t,wgt_t> Sending;
  typedef std::pair<agi::GraphVertex*,agi::GraphVertex*> Connection;
  //Array of weights
  typedef double* Ws;
  //Map from part id to weights
  typedef std::unordered_map<part_t, Ws> Midd;
  struct Migr;
  struct CompareMigr;
  typedef std::set<Migr,CompareMigr> MigrComm;
  typedef std::unordered_set<agi::GraphEdge*> EdgeSet;

  class Selector {
  public:
    Selector(DiffusiveInput* in_, Queue* queue,
             std::vector<int>* cd, std::vector<double>* cw);
    wgt_t select(Targets* targets,agi::Migration* plan,
                 wgt_t planW, unsigned int cavSize,int);
    wgt_t kkSelect(LIDs order, Targets* targets,agi::Migration* plan,
                 wgt_t planW, unsigned int cavSize,int);
    void getCavitiesAndPeers(agi::etype t,
        LIDs plan, LIDs vtxOwner,
        CSR& cavities, CSR& peers, CSR& eoc);
    void selectDisconnected(agi::Migration* plan, int target_dimension);
    Midd* trim(Targets* targets, agi::Migration* plan);
    void cancel(agi::Migration*& plan,Midd* capacity);

    void updatePartWeight(agi::Migration* plan, int target_dimension, agi::WeightPartitionMap* wp_map);
  protected:
    //Create the initial cavities
    void constructCavities();
    //Remove vertices that are planned to be sent in cavities
    void updateCavities(agi::Migration* plan);
    

    typedef std::map<int,EdgeSet> PeerEdgeSet;
    //Trim functions
    void insertInteriorEdges(agi::GraphVertex*,agi::part_t, EdgeSet&,int);
    void calculatePlanWeights(agi::Migration* plan, std::unordered_map<int,double>& vtx_weight,
                              PeerEdgeSet* peerEdges, std::unordered_set<part_t>& neighbors);
    void sendPlanWeight(std::unordered_map<int,double>& vtx_weight, PeerEdgeSet* peerEdges,
                        std::unordered_set<part_t>& neighbors);
    void receiveIncomingWeight(MigrComm& incoming);
    bool determineAvailability(Ws& avail);
    void acceptWeight(MigrComm& incoming, bool& isAvail, Ws& avail, Midd& accept);
    void sendAcceptedWeights(Midd& accept);
    void gatherCapacities(Midd* capacity);


    void tempInsertInteriorEdges(agi::GraphVertex*,agi::part_t,
                                           EdgeSet&,int,const EdgeSet&);
    double weight(const EdgeSet&);
    void combineSets(EdgeSet&,const EdgeSet&);

    void calculatePlanWeight(agi::Migration* plan, int target_dimension,
                             std::unordered_map<part_t,wgt_t>& weight);

    
    DiffusiveInput* in;
    agi::Ngraph* g;
    Sending sending;
    Queue* q;
    std::unordered_map<agi::GraphEdge*, Cavity> cavities;
    std::vector<int>* completed_dimensions;
    std::vector<double>* completed_weights;
  };

  Selector* makeSelector(DiffusiveInput* in,Queue* q,
                         std::vector<int>* cd,
                         std::vector<double>* cw );
}

#endif
