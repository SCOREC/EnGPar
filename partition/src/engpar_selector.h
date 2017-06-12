#ifndef ENGPAR_SELECTOR_H
#define ENGPAR_SELECTOR_H

#include "../engpar.h"
#include "engpar_sides.h"
#include "engpar_targets.h"
#include "engpar_queue.h"
#include <unordered_set>
#include "engpar_input.h"
namespace engpar {

  typedef std::vector<agi::GraphVertex*> Cavity;
  typedef std::unordered_set<part_t> Peers;
  typedef std::unordered_map<part_t,wgt_t> Sending;
  typedef double Ws[5];
  typedef std::map<int, Ws> Midd;

  class Selector {
  public:
    Selector(Input* in_, Queue* queue,
             std::vector<int>* cd, std::vector<double>* cw) :
      in(in_), g(in_->g),
      q(queue),
      completed_dimensions(cd), completed_weights(cw) {}

    wgt_t select(Targets* targets,agi::Migration* plan,
                         wgt_t planW, unsigned int cavSize,int);
    Midd* trim(Targets* targets, agi::Migration* plan);
    void cancel(agi::Migration*& plan,Midd* capacity);

  protected:
    typedef std::unordered_set<agi::GraphEdge*> EdgeSet;
    typedef std::map<int,EdgeSet> PeerEdgeSet;
    void insertInteriorEdges(agi::GraphVertex*,agi::part_t,
				       EdgeSet&,int);
    void tempInsertInteriorEdges(agi::GraphVertex*,agi::part_t,
					   EdgeSet&,int,const EdgeSet&);
    double weight(const EdgeSet&);
    void combineSets(EdgeSet&,const EdgeSet&);

    Input* in;
    agi::Ngraph* g;
    Sending sending;
    Queue* q;
    std::vector<int>* completed_dimensions;
    std::vector<double>* completed_weights;
  };

  Selector* makeSelector(Input* in,Queue* q,
                         std::vector<int>* cd,
                         std::vector<double>* cw );
}

#endif
