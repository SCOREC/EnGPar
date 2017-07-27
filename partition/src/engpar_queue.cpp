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

  /* Sort algorithm adapted from:
   * http://www.algolist.net/Algorithms/Sorting/Quicksort
   */
  void bfsSort(agi::lid_t* arr,int left, int right,int* values) {
    int i = left, j = right;
    int tmp;
    int pivot = values[arr[(left + right) / 2]];
    while (i <= j) {
      while (values[arr[i]] < pivot)
        i++;
      while (values[arr[j]] > pivot)
        j--;
      if (i <= j) {
        tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
        i++;
        j--;
      }
    }
    /* recursion */
    if (left < j)
      bfsSort(arr, left, j,values);
    if (i < right)
      bfsSort(arr, i, right,values);
  }

  agi::lid_t runBFSDisjoint(agi::PNgraph* pg, agi::etype t,agi::lid_t* seed,
                    agi::lid_t& numSeeds, agi::lid_t start_seed, int* visited,
                    agi::lid_t depth=0) {
    int num_osets =numSeeds-start_seed;
    int num_sets = num_osets;
    int* parents = new int[num_sets];
    int* set_size = new int[num_sets];
    for (int i=0;i<num_sets;i++) {
      parents[i] = i;
      set_size[i] = 1;
    }
    int* labels = new int[pg->num_local_edges[t]];
    int l=0;
    for (agi::lid_t i=0;i<start_seed;i++)
      labels[seed[i]]=-1;
    for (agi::lid_t i=start_seed;i<numSeeds;i++) {
      visited[seed[i]] = depth;
      labels[seed[i]] = l++;
    }
    assert(l==num_sets);
    int num_updates;
    int level=depth;
    //Implemented for HG
    do {
      num_updates=0;
      for (agi::lid_t i=0;i<pg->num_local_verts;i++) {
        int min_dist=-1;
        int label = -1;
        for (agi::lid_t j=pg->degree_list[t][i];j<pg->degree_list[t][i+1];j++){
          agi::lid_t edge = pg->edge_list[t][j];
          if (visited[edge]!=-1&&(min_dist==-1||visited[edge]<min_dist)){
            min_dist=visited[edge];
            label = labels[edge];
            if (label ==-1)
              break;
            while(label!=parents[label])
              label = parents[label];
          }
        }
        if (label==-1)
          continue;
        if (min_dist==level) {
          for (agi::lid_t j=pg->degree_list[t][i];j<pg->degree_list[t][i+1];j++){
            agi::lid_t edge = pg->edge_list[t][j];
            if (visited[edge]==-1) {
              seed[numSeeds++] = edge;
              num_updates++;
              visited[edge] = min_dist+1;
              labels[edge] = label;
              set_size[label]++;
            }
            else {
              int l2 = labels[edge];
              while(l2!=parents[l2])
                l2=parents[l2];
              if (label!=l2) {
                num_sets--;
                set_size[label]+=set_size[l2];
                set_size[l2]=0;
                parents[l2]=label;
              }
            }
          }
        }
      }
      level++;
    }
    while (num_updates);
    int addition=0;
    if (num_sets>1) {
      for (int i=0;i<num_sets;i++) {
        int l=0;
        for (int j=1;j<num_osets;j++) {
          if (set_size[j]>set_size[l]) 
            l=j;
        }
        if (addition>0) {
          for (int j=0;j<pg->num_local_edges[t];j++) {
            int l2 = labels[j];
            while (l2!=parents[l2])
              l2 = parents[l2];
            if (l2==l)
              visited[j]+=addition;
          }
        }
        addition+=set_size[i];
        set_size[i]=0;
      }
    }
    //Sort the recent BFS in accordance to the offsets
    bfsSort(seed,start_seed,numSeeds-1,visited);
    return level+addition;
  }
  agi::lid_t runBFS(agi::PNgraph* pg, agi::etype t,agi::lid_t* seed,
                    agi::lid_t& numSeeds, agi::lid_t start_seed, int* visited) {
    for (agi::lid_t i=start_seed;i<numSeeds;i++) 
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
    return level;
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
    int* visited = new int[pg->num_local_edges[t]];
    for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
      visited[i] = -1;
    }
    runBFS(pg,t,first_bfs,size,0,visited);
    agi::GraphTag* tag = g->createIntTag(0);
    agi::GraphEdge* edge;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((edge=g->iterate(eitr))) {
      agi::lid_t lid = g->localID(edge);
      g->setIntTag(tag,edge,visited[lid]);
    }
    g->destroy(eitr);
    if (g->hasCoords()) {
      std::string filename = "first_bfs";
      agi::writeVTK(g,filename.c_str(),tag,0);
    }
    
    int* visited2 = new int[pg->num_local_edges[t]];
    agi::lid_t size2=0;
    agi::lid_t* second_bfs = new agi::lid_t[pg->num_local_edges[t]];
    for (agi::lid_t i=0;i<pg->num_local_edges[t];i++)
      visited2[i]=-1;
    agi::lid_t ind = size-1;
    agi::lid_t depth = 0;
    while (size2<size) {
      while(visited2[first_bfs[ind]]!=-1&&ind>=0)
        ind--;
      if (ind<0)
        break;
      int max_depth = visited[first_bfs[ind]];
      agi::lid_t start_seed = size2;
      while(visited[first_bfs[ind]]==max_depth&&ind>0) {
        if (visited2[first_bfs[ind]]==-1)
          second_bfs[size2++] = first_bfs[ind];
        ind--;
      }
      depth = runBFSDisjoint(pg,t,second_bfs,size2,start_seed,visited2,depth);
    }
    agi::GraphTag* tag2 = g->createIntTag(0);
    eitr = g->begin(0);
    while ((edge=g->iterate(eitr))) {
      agi::lid_t lid = g->localID(edge);
      g->setIntTag(tag2,edge,visited2[lid]);
    }
    g->destroy(eitr);
    if (g->hasCoords()) {
      std::string filename = "second_bfs";
      agi::writeVTK(g,filename.c_str(),tag2,0);
    }
    delete [] visited;
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
