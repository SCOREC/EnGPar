#ifndef __ENGPAR_QUEUE_H__
#define __ENGPAR_QUEUE_H__
#include <ngraph.h>
#include <engpar_metrics.h>
#include "../engpar_diffusive_input.h"

namespace engpar {
  class ShrinkingArray {
  public:
    ShrinkingArray(int n) {
      s = 0;
      current = new agi::GraphEdge*[n];
      next = new agi::GraphEdge*[n];
      next_size = 0;
    }
    ~ShrinkingArray() {
      if (current)
        delete [] current;
      if (next)
        delete [] next;
    }
    agi::lid_t size() {return s;}
    void addElement(agi::GraphEdge* e) {
      next[next_size++] = e;
    }
    void startIteration() {
      agi::GraphEdge** tmp = current;
      current = next;
      next = tmp;
      s = next_size;
      next_size=0;
    }
    typedef agi::lid_t iterator;
    iterator begin() {return 0;}
    iterator end() {return s;}
    agi::GraphEdge* operator[](iterator pos) {return current[pos];}
    agi::GraphEdge* get(iterator pos) {return current[pos];}
  private:
    agi::lid_t s;
    agi::GraphEdge** current;
    agi::lid_t next_size;
    agi::GraphEdge** next;
  };
  typedef ShrinkingArray Queue;
  Queue* createIterationQueue(agi::Ngraph*);
  Queue* createDistanceQueue(DiffusiveInput*);
}

#endif
