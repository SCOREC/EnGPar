#include "engpar_queue.h"
#include <PCU.h>

namespace engpar {
  Queue* createIterationQueue(agi::Ngraph* g) {
    Queue* q = new Queue;
    agi::GraphEdge* edge;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((edge = g->iterate(eitr))) {
      agi::PinIterator* pitr = g->pins(edge);
      agi::lid_t degree = g->degree(edge);
      for (agi::lid_t i=0;i<degree;i++) {
        agi::GraphVertex* pin = g->iterate(pitr);
        if (g->owner(pin)!=PCU_Comm_Self()) {
          q->push_back(edge);
          break;
        }
      }
      g->destroy(pitr);
    }
    g->destroy(eitr);
    return q;
  }
  /*
  void runBFS(agi::PNgraph* pg, agi::etype t,agi::lid_t* seed,int& numSeeds) {
    int* visited = new int[pg->num_local_edges[t]];
    for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
      visited[i] = -1;
    }
    for (agi::lid_t i=0;i<numSeeds;i++)
      visited[seed[i]] = 0;
    
    int num_updates = 0;
    int level=0;
    do {
      num_updates=0;
      for (agi::lid_t i=0;i<pg->num_local_verts;i++) {
        int min_dist=-1;
        for (agi::lid_t j=degree_list[t][i];j<degree_list[t][i+1];j++)  {
          agi::lid_t edge = edge_list[j];
          min_dist = (visited[edge]!=-1&&visited[edge]<min_dist)
            ? visited[edge] : min_dist;
        }
        if (min_dist==level) {
          for (agi::lid_t j=degree_list[t][i];j<degree_list[t][i+1];j++)  {
            agi::lid_t edge = edge_list[j];
            visited[edge] = (visited[edge]==-1) ? min_dist+1 : visited[edge];
            num_updates++;
          }
        }
      }
      level++;
    }
    while (num_updates);
    delete [] visited;
  }
  */
  
  Queue* createDistanceQueue(agi::Ngraph* g) {  
    int size=0;
    agi::PNgraph* pg = g->publicize();
    agi::etype t = 0;
    agi::lid_t* first_bfs = new agi::lid_t[pg->num_local_edges[t]];


    if (pg->isHyperGraph) {
      for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
        bool isShared=false;
        for (agi::lid_t j=pg->pin_degree_list[t][i];
             j<pg->pin_degree_list[t][i+1];++j) {
          int vert = pg->pin_list[t][j];
          isShared = isShared||pg->owners[vert]!=PCU_Comm_Self();
        }
        if (isShared)
          first_bfs[size++] = i;
      }
    }
    else {
      for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
        bool isShared=false;
        int v = pg->edge_list[t][i];
        isShared = isShared||pg->owners[v]!=PCU_Comm_Self();
        if (isShared)
          first_bfs[size++] = i;
      }
    }

    //runBFS(pf,t,first_bfs,size);

    delete [] first_bfs;
    Queue* q = new Queue;
    return q;
  }
}
