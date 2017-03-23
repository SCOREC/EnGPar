#include "EdgeIterator.h"

#ifndef __HYPEREDGEITERATOR_H__
#define __HYPEREDGEITERATOR_H__

namespace agi {

  //This edge iterator is for iterating over all hyperedges in the graph
  class HyperEdgeIterator : public EdgeIterator {
    friend class Ngraph;
  public:
    ~HyperEdgeIterator() {}
  protected:
  HyperEdgeIterator(etype t, int nt, lid_t deg)
    : EdgeIterator(t,nt,NULL,deg){}
    bool isHyper() {return true;}
  };
}

#endif
