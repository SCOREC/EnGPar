#include "engpar_queue.h"
#include <PCU.h>
#include <algorithm>

namespace engpar {
  Queue* createIterationQueue(agi::Ngraph* g) {
    int n=0;
    agi::GraphEdge* edge;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((edge = g->iterate(eitr))) {
      agi::PinIterator* pitr = g->pins(edge);
      agi::lid_t degree = g->degree(edge);
      for (agi::lid_t i=0;i<degree;i++) {
        agi::GraphVertex* pin = g->iterate(pitr);
        if (g->owner(pin)!=PCU_Comm_Self()) {
          n++;
          break;
        }
      }
      g->destroy(pitr);
    }
    g->destroy(eitr);
    Queue* q = new Queue(n);
    eitr = g->begin(0);
    while ((edge = g->iterate(eitr))) {
      agi::PinIterator* pitr = g->pins(edge);
      agi::lid_t degree = g->degree(edge);
      for (agi::lid_t i=0;i<degree;i++) {
        agi::GraphVertex* pin = g->iterate(pitr);
        if (g->owner(pin)!=PCU_Comm_Self()) {
          q->addElement(edge);
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

  class Inputs {
  public:
    Inputs() {
      seeds=NULL;
      numSeeds=0;
      visited=NULL;
      num_osets=0;
      num_sets=0;
      parents=NULL;
      set_size=NULL;
      labels=NULL;
    }
    ~Inputs() {
      if (seeds)
        delete [] seeds;
      if (visited)
        delete [] visited;
      if (parents) 
        delete [] parents;
      if (set_size)
        delete [] set_size;
      if (labels)
        delete [] labels;
    }
    agi::lid_t* seeds;
    agi::lid_t numSeeds;
    int* visited;

    int num_osets;//= numSeeds-start_seed;
    int num_sets; //= num_osets;
    int* parents; //= new int[num_sets];
    int* set_size;//= new int[num_sets];
    int* labels;  //= new int[pg->num_local_edges[t]];
  };
  //Visit function takes in the input, the source edge and destination edge
  typedef bool (*visitFn)(Inputs*,agi::lid_t,agi::lid_t);

  //Visit function for first traversal
  bool depth_visit(Inputs* i,agi::lid_t source,agi::lid_t dest) {
    int min_dist = i->visited[source];
    if (i->visited[dest]==-1) {
      i->seeds[i->numSeeds++] = dest;
      i->visited[dest] = min_dist+1;
      return true;
    }
    return false;
  }
  
  //Visit function for second traversal (uses disjoint sets)
  bool distance_visit(Inputs* i,agi::lid_t source, agi::lid_t dest) {
    int label = i->labels[source];
    int min_dist = i->visited[source];
    if (label==-1)
      return false;
    while(label!=i->parents[label])
      label = i->parents[label];
    if (i->visited[dest]==-1) {
      i->seeds[i->numSeeds++] = dest;
      i->visited[dest] = min_dist+1;
      i->labels[dest] = label;
      i->set_size[label]++;
      return true;
    }
    int l2 = i->labels[dest];
    while(l2!=i->parents[l2]) 
      l2=i->parents[l2];
    if (label!=l2) {
      i->num_sets--;
      i->set_size[label]+=i->set_size[l2];
      i->set_size[l2]=0;
      i->parents[l2]=label;
    }
    return false;
  }

  /*
    Push based BFS takes in :
      graph
      edge type
      first seed index into visited array
      starting depth
      the visit function
      inputs object
  */
  int bfs_push(agi::Ngraph* g, agi::etype t, agi::lid_t start_seed,
               int start_depth,visitFn visit, Inputs* in) {
    agi::GraphTag* processed = g->createIntTag(-1);
    agi::GraphVertex* v;
    agi::VertexIterator* vitr = g->begin();
    while ((v = g->iterate(vitr))) {
      g->setIntTag(processed,v,0);
    }
    for (agi::lid_t i=start_seed;i<in->numSeeds;i++) 
      in->visited[in->seeds[i]] = start_depth;
    agi::PNgraph* pg = g->publicize();
    agi::lid_t index = start_seed;
    while (index!=in->numSeeds) {
      agi::GraphEdge* source = pg->getEdge(in->seeds[index++],t);
      agi::PinIterator* pitr = g->pins(source);
      agi::GraphVertex* bridge;
      while ((bridge = g->iterate(pitr))) {
        if (g->owner(bridge)!=PCU_Comm_Self())
          continue;
        int process = g->getIntTag(processed,bridge);
        if (process)
          continue;
        agi::EdgeIterator* eitr = g->edges(bridge);
        agi::GraphEdge* dest;
        while ((dest = g->iterate(eitr))) {
          if (dest==source)
            continue;
          visit(in,g->localID(source),g->localID(dest));
        }
        g->destroy(eitr);
        g->setIntTag(processed,bridge,1);
      }
      g->destroy(pitr);
    }
    g->destroyTag(processed);
    return in->visited[in->seeds[in->numSeeds-1]];
  }
  
  /*
    Pull based BFS takes in :
      graph
      edge type
      first seed index into visited array
      starting depth
      the visit function
      inputs object
  */
  int bfs_pull(agi::Ngraph* g, agi::etype t,agi::lid_t start_seed,
               int start_depth, visitFn visit, Inputs* in) {
    agi::PNgraph* pg = g->publicize();

    int level=start_depth;
    for (agi::lid_t i=start_seed;i<in->numSeeds;i++) 
      in->visited[in->seeds[i]] = level;

    agi::lid_t* current_vtxs = new agi::lid_t[pg->num_local_verts];
    agi::lid_t num_cv = pg->num_local_verts;
    for (agi::lid_t i=0;i<pg->num_local_verts;i++)
      current_vtxs[i] = i;
    agi::lid_t* next_vtxs = new agi::lid_t[pg->num_local_verts];
    agi::lid_t num_nv = 0;
    int num_updates;
    //Implemented for HG
    do {
      num_updates=0;
      for (agi::lid_t i=0;i<num_cv;i++) {
        int source = -1;
        agi::lid_t v = current_vtxs[i];
        for (agi::lid_t j=pg->degree_list[t][v];j<pg->degree_list[t][v+1];j++){
          agi::lid_t edge = pg->edge_list[t][j];
          if (in->visited[edge] != -1 &&
              (source == -1 || in->visited[edge] < in->visited[source]))
            source = edge;
        }
        if (source!=-1&&in->visited[source]==level) {
          for (agi::lid_t j=pg->degree_list[t][v];j<pg->degree_list[t][v+1];j++){
            agi::lid_t edge = pg->edge_list[t][j];
            num_updates+=visit(in,source,edge);
          }
        }
        else {
          next_vtxs[num_nv++] = v;
        }
      }
      agi::lid_t* tmp = current_vtxs;
      current_vtxs = next_vtxs;
      next_vtxs = tmp;
      num_cv = num_nv;
      num_nv =0;
      level++;
    }
    while (num_updates);
    delete [] current_vtxs;
    delete [] next_vtxs;
    return level;
  }

  Queue* createDistanceQueue(DiffusiveInput* input) {
    agi::Ngraph* g = input->g;
    if (g->numLocalEdges()==0) {
      return new Queue(0);
    }
    agi::PNgraph* pg = g->publicize();
    agi::etype t = 0;
    //Setup Inputs for first BFS traversal
    Inputs* in1 = new Inputs;
    in1->seeds = new agi::lid_t[pg->num_local_edges[t]];
    //If run in serial use edge 2 for the seed
    if (g->numGhostVtxs()==0)
       in1->seeds[in1->numSeeds++]=2;
    //Otherwise use all edges that are shared across part boundaries
    else if (pg->isHyperGraph) {
      for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
        bool isShared=false;
        for (agi::lid_t j=pg->pin_degree_list[t][i];
             j<pg->pin_degree_list[t][i+1];++j) {
          agi::lid_t vert = pg->pin_list[t][j];
          isShared = isShared||vert>=pg->num_local_verts;
        }
        if (isShared) {
          in1->seeds[in1->numSeeds++] = i;
        }
      }
    }
    else {
      for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
        agi::lid_t v = pg->edge_list[t][i];
        if (v>=pg->num_local_verts)
          in1->seeds[in1->numSeeds++] = i;
      }
    }
    Queue* q = new Queue(in1->numSeeds);
    if (in1->numSeeds==0) {
      delete in1;
      return q;
    }
    
    in1->visited = new int[pg->num_local_edges[t]];
    for (agi::lid_t i=0;i<pg->num_local_edges[t];i++) {
      in1->visited[i] = -1;
    }
    //Run first BFS using depth visit operation
    if (input->bfsPush)
      bfs_push(g,t,0,0,depth_visit,in1);
    else
      bfs_pull(g,t,0,0,depth_visit,in1);
    //Setup inputs to second BFS traversal
    Inputs* in2 = new Inputs;
    in2->visited = new int[pg->num_local_edges[t]];
    in2->seeds = new agi::lid_t[pg->num_local_edges[t]];
    for (agi::lid_t i = 0; i < pg->num_local_edges[t]; i++)
      in2->visited[i] = -1;
    in2->labels = new int[pg->num_local_edges[t]];
    for (agi::lid_t i=0;i<pg->num_local_edges[t];i++)
      in2->labels[i]=-1;

    //Seconds BFS traversal for distance computation
    // Automatically detects disconnected components using disjoint sets
    agi::lid_t ind = in1->numSeeds-1;
    agi::lid_t depth = 0;
    while (in2->numSeeds < in1->numSeeds) {
      //Find the next deepest layer of edges
      while(ind >= 0&&in2->visited[in1->seeds[ind]] != -1)
        --ind;
      if (ind < 0)
        break;
      int max_depth = in1->visited[in1->seeds[ind]];
      agi::lid_t start_seed = in2->numSeeds;
      while(ind >= 0 && in1->visited[in1->seeds[ind]] == max_depth) {
        if (in2->visited[in1->seeds[ind]] == -1)
          in2->seeds[in2->numSeeds++] = in1->seeds[ind];
        --ind;
      }
      //Setup the disjoint set structure
      in2->num_osets = in2->numSeeds-start_seed;
      in2->num_sets = in2->num_osets;
      if (in2->num_sets==0)
        continue;
      in2->parents = new int[in2->num_sets];
      in2->set_size = new int[in2->num_sets];
      for (int i=0;i<in2->num_sets;i++) {
        in2->parents[i] = i;
        in2->set_size[i] = 1;
      }

      for (agi::lid_t i=start_seed;i<in2->numSeeds;i++) {
        in2->labels[in2->seeds[i]] = i-start_seed;
      }
      //Run the bfs using the distance visit operation
      if (input->bfsPush)
        depth = bfs_push(g,t,start_seed,depth,distance_visit,in2);
      else
        depth = bfs_pull(g,t,start_seed,depth,distance_visit,in2);
      //Correct the distance when there are multiple disjoint sets
      int addition=0;
      if (in2->num_sets>1) {
        for (int i=0;i<in2->num_sets;i++) {
          int l=0;
          for (int j=1;j<in2->num_osets;j++) {
            if (in2->set_size[j]>in2->set_size[l]) 
              l=j;
          }
          if (addition>0) {
            for (agi::lid_t j=start_seed;j<in2->numSeeds;j++) {
              int l2 = in2->labels[in2->seeds[j]];
              if (l2<0)
                continue;
              while (l2!=in2->parents[l2])
                l2 = in2->parents[l2];
              if (l2==l)
                in2->visited[in2->seeds[j]]+=addition;
            }
          }
          addition+=in2->set_size[l];
          in2->set_size[l]=0;
        }
        depth+=addition;
      }
      //Sort the recent BFS in accordance to the offsets
      bfsSort(in2->seeds,start_seed,in2->numSeeds-1,in2->visited);
      
      //Cleanup disjoint sets
      delete [] in2->parents;
      in2->parents = NULL;
      delete [] in2->set_size;
      in2->set_size = NULL;
    }
    delete in1;
    //Setup the Queue from the second bfs
    for (int i=in2->numSeeds;i>0;i--) {
      agi::lid_t lid = in2->seeds[i-1];
      agi::GraphEdge* edge = pg->getEdge(lid,t);
      agi::Peers res;
      g->getResidence(edge,res);
      if (res.size()>1)
        q->addElement(edge);
    }
    delete in2;
    return q;
  }

  LIDs getCavityOrder(agi::Ngraph* g, const agi::lid_t t, Queue* qu) {
    assert( t == 0 ); // the distance queue is hardcoded to use the 0th edge type
    agi::PNgraph* pg = g->publicize();
    const int numEdges = pg->num_local_edges[t];
    agi::lid_t* order_h = new agi::lid_t[numEdges];
    std::fill(order_h,order_h+numEdges,-1);
    int pos = 0;
    qu->startIteration();
    Queue::iterator itr;
    for (itr = qu->begin();itr!=qu->end();itr++) {
      agi::GraphEdge* e = qu->get(itr);
      const agi::lid_t eid = g->localID(e);
      assert(eid >=0 && eid < numEdges);
      order_h[eid] = pos++;
    }
    LIDs order("cavOrder", numEdges);
    hostToDevice(order,order_h);
    delete [] order_h;
    return order;
  }
}
