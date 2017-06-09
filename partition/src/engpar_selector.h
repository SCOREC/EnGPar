#ifndef ENGPAR_SELECTOR_H
#define ENGPAR_SELECTOR_H

#include "../engpar.h"
#include "engpar_sides.h"
#include "engpar_targets.h"
#include "engpar_queue.h"
#include <unordered_set>

namespace engpar {

  typedef std::vector<agi::GraphVertex*> Cavity;
  typedef std::unordered_set<part_t> Peers;
  typedef std::unordered_map<part_t,wgt_t> Sending;
  typedef double Ws[5];
  typedef std::map<int, Ws> Midd;

  class Selector {
  public:
    Selector(agi::Ngraph* graph, Queue* queue,
             std::vector<int>* cd, std::vector<double>* cw) : g(graph), q(queue),
             completed_dimensions(cd), completed_weights(cw) {}

    virtual wgt_t select(Targets* targets,agi::Migration* plan,
                         wgt_t planW, unsigned int cavSize,int);
    virtual Midd* trim(Targets* targets, agi::Migration* plan);
    virtual void cancel(agi::Migration*& plan,Midd* capacity);

  protected:
    typedef std::unordered_set<agi::GraphEdge*> EdgeSet;
    typedef std::map<int,EdgeSet> PeerEdgeSet;
    void insertInteriorEdges(agi::GraphVertex*,agi::part_t,
				       EdgeSet&,int);
    void tempInsertInteriorEdges(agi::GraphVertex*,agi::part_t,
					   EdgeSet&,int,const EdgeSet&);
    double weight(const EdgeSet&);
    void combineSets(EdgeSet&,const EdgeSet&);
				       
    agi::Ngraph* g;
    Sending sending;
    Queue* q;
    std::vector<int>* completed_dimensions;
    std::vector<double>* completed_weights;
  };
}

#endif
