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
  
  void runBFS(agi::PNgraph* pg, agi::etype t,agi::lid_t* seed,
              agi::lid_t& numSeeds, int*& visited) {
    visited = new int[pg->num_local_edges[t]];
    for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
      visited[i] = -1;
    }
    for (agi::lid_t i=0;i<numSeeds;i++) 
      visited[seed[i]] = 0;      
    int num_updates;
    int level=0;
    //Implemented for HG
    do {
      num_updates=0;
      for (agi::lid_t i=0;i<pg->num_local_verts;i++) {
        int min_dist=-1;
        for (agi::lid_t j=pg->degree_list[t][i];j<pg->degree_list[t][i+1];j++){
          agi::lid_t edge = pg->edge_list[t][j];
          min_dist = (visited[edge]!=-1&&(min_dist==-1||visited[edge]<min_dist))
            ? visited[edge] : min_dist;
        }
        if (min_dist==level) {
          for (agi::lid_t j=pg->degree_list[t][i];j<pg->degree_list[t][i+1];j++){
            agi::lid_t edge = pg->edge_list[t][j];
            if (visited[edge]==-1) {
              seed[numSeeds++] = edge;
              num_updates++;
              visited[edge] = min_dist+1;
            }
          }
        }
      }
      level++;
    }
    while (num_updates);
    //delete [] visited;
  }
  
  Queue* createDistanceQueue(agi::Ngraph* g) {  
    agi::lid_t size=0;
    agi::PNgraph* pg = g->publicize();
    agi::etype t = 0;
    agi::lid_t* first_bfs = new agi::lid_t[pg->num_local_edges[t]];

    if (PCU_Comm_Peers()==1)
      first_bfs[size++]=2;
    else if (pg->isHyperGraph) {
      for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
        bool isShared=false;
        for (agi::lid_t j=pg->pin_degree_list[t][i];
             j<pg->pin_degree_list[t][i+1];++j) {
          agi::lid_t vert = pg->pin_list[t][j];
          isShared = isShared||vert>=pg->num_local_verts;
        }
        if (isShared) {
          first_bfs[size++] = i;
        }
      }
    }
    else {
      for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
        agi::lid_t v = pg->edge_list[t][i];
        if (v>=pg->num_local_verts)
          first_bfs[size++] = i;
      }
    }
    int* visited = NULL;
    runBFS(pg,t,first_bfs,size,visited);
    agi::GraphTag* tag = g->createIntTag(0);
    agi::GraphEdge* edge;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((edge=g->iterate(eitr))) {
      agi::lid_t lid = g->localID(edge);
      g->setIntTag(tag,edge,visited[lid]);
    }
    g->destroy(eitr);
    if (g->hasCoords() &&PCU_Comm_Peers()==1) {
      std::string filename = "first_bfs";
      agi::writeVTK(g,filename.c_str(),tag,0);
    }
    delete [] visited;
    
    int* visited2 = NULL;
    agi::lid_t size2=0;
    agi::lid_t* second_bfs = new agi::lid_t[pg->num_local_edges[t]];
    int max_depth = visited[first_bfs[size-1]];
    agi::lid_t ind = size-1;
    while(visited[first_bfs[ind]]==max_depth&&ind>0) {
      second_bfs[size2++] = first_bfs[ind];
      ind--;
    }
    runBFS(pg,t,second_bfs,size2,visited2);
    agi::GraphTag* tag2 = g->createIntTag(0);
    eitr = g->begin(0);
    while ((edge=g->iterate(eitr))) {
      agi::lid_t lid = g->localID(edge);
      g->setIntTag(tag2,edge,visited2[lid]);
    }
    g->destroy(eitr);
    if (g->hasCoords()&&PCU_Comm_Peers()==1) {
      std::string filename = "second_bfs";
      agi::writeVTK(g,filename.c_str(),tag2,0);
    }    
    delete [] visited2;
    delete [] first_bfs;
    Queue* q = new Queue;
    for (int i=size2;i>0;i--) {
      agi::lid_t lid = second_bfs[i-1];
      agi::GraphEdge* edge = pg->getEdge(lid,t);
      assert(g->localID(edge)==lid);
      agi::Peers res;
      g->getResidence(edge,res);
      if (res.size()>1)
        q->push_back(edge);
    }
    delete [] second_bfs;
    return q;
  }
}
