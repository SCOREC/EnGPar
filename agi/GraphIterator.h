#ifndef __GRAPH_ITERATOR__H__
#define __GRAPH_ITERATOR__H__
#include "agi.h"

namespace agi {

class Ngraph;
class EdgeIterator;
class PinIterator;
class GraphEdge;

 class GraphIterator {
  friend class Ngraph;
 private:
  bool isH;
  int count;
  EdgeIterator* eitr;
  PinIterator* pitr;
  GraphEdge* edge;
  GraphIterator(EdgeIterator* itr, bool isHyperGraph) : isH(isHyperGraph),
    count(0), eitr(itr),pitr(NULL),edge(NULL) {}
  void setEdge(GraphEdge* e, PinIterator* p) {
    edge = e;
    count = 0;
    pitr = p;
  }
};
  
}


#endif
