#include "engpar_selector.h"
#include <engpar_metrics.h>
#include <engpar_kokkosColoring.h>
#include <engpar_support.h>
#include <set>
#include <PCU.h>
#include <agiMigration.h>
namespace engpar {

  void getCavity(agi::Ngraph* g, agi::GraphEdge* edge, agi::Migration* plan,
                 Cavity& cav, Peers& peers) {
    agi::PinIterator* pitr = g->pins(edge);
    agi::lid_t deg = g->degree(edge);
    agi::GraphVertex* vtx;
    typedef std::map<agi::part_t, std::unordered_set<agi::GraphEdge*> > Peer_Map;
    Peer_Map peerMap;
    std::set<agi::GraphEdge*> edgesOfCavity;
    for (agi::lid_t i =0;i<deg;i++) {
      vtx = g->iterate(pitr);
      if (g->owner(vtx)==PCU_Comm_Self()) {
        if(!plan->has(vtx)) {
          cav.push_back(vtx);
          agi::GraphEdge* e;
          agi::EdgeIterator* eitr = g->edges(vtx);
          while ((e=g->iterate(eitr)))
            edgesOfCavity.insert(e);
          g->destroy(eitr);
        }
      }
    }
    g->destroy(pitr);
    std::set<agi::GraphEdge*>::iterator sitr;
    for (sitr = edgesOfCavity.begin(); sitr != edgesOfCavity.end(); sitr++) {
      agi::GraphVertex* v;
      pitr = g->pins(*sitr);
      while ((v = g->iterate(pitr))) {
        if (g->owner(v)!=PCU_Comm_Self())
          peerMap[g->owner(v)].insert(*sitr);
      }
      g->destroy(pitr);
    }
    //create a sorted list (in an unordered set?) of peers
    while (peers.size()!=peerMap.size()) {
      unsigned int max =0;
      Peer_Map::iterator itr;
      for (itr = peerMap.begin(); itr != peerMap.end(); itr++)
        if (itr->second.size() > max)
          max = itr->second.size();
      if (max<2&&peers.size()>0)
        break;
      for (itr = peerMap.begin(); itr != peerMap.end(); itr++)
        if (itr->second.size() == max) {
          peers.insert(itr->first);
          itr->second.clear();
        }
    }
  }

  wgt_t addCavity(agi::Ngraph* g, Cavity& cav,
                  part_t peer, agi::Migration* plan,
                  int target_dimension) {
    Cavity::iterator itr;
    wgt_t w=0.0;
    std::set<agi::GraphEdge*> target_edges;
    for (itr = cav.begin();itr!=cav.end();itr++) {
      plan->insert(std::make_pair(*itr,peer));
      if (target_dimension==-1)
        w+= g->weight(*itr);
      else {
        agi::EdgeIterator* eitr = g->edges(*itr,target_dimension);
        agi::GraphEdge* e;
        while ((e=g->iterate(eitr)))
          target_edges.insert(e);
        g->destroy(eitr);
      }
    }
    std::set<agi::GraphEdge*>::iterator sitr;
    for (sitr = target_edges.begin();sitr!=target_edges.end();++sitr) {
      agi::GraphEdge* e = *sitr;
      if (!g->isResidentOn(e,peer))
        w+=g->weight(e);
    }
    return w;
  }

  typedef std::pair<agi::GraphVertex*,agi::GraphVertex*> Connection;
  typedef std::map<Connection,int> Connection_Count;
  void getCavityConnections(agi::Ngraph* g,agi::etype con, int minCon, Cavity cav,
                            Connection_Count& conns) {
    for (unsigned int i=0;i<cav.size();i++) {
      agi::GraphVertex* u = cav[i];
      agi::GraphIterator* gitr = g->adjacent(u,con);
      agi::GraphVertex* v;
      std::map<agi::GraphVertex*, int> occ;
      while ((v = g->iterate(gitr))) {
        if (v==u) continue;
        occ[v]++;
      }
      g->destroy(gitr);
      std::map<agi::GraphVertex*, int>::iterator itr;
      for (itr=occ.begin();itr!=occ.end();itr++) {
        if (itr->second>=minCon) {
          agi::GraphVertex* v = itr->first;
          if (u<v)
            conns[std::make_pair(u,v)]++;
          else
            conns[std::make_pair(v,u)]++;
        }
      }
    }
  }

  //Returns true if the cavity isn't fully connected to the part
  //  In the example of a mesh, full connectivity is having a face connection
  bool isPartiallyConnected(agi::Ngraph* g,agi::etype con, int minCon,
                            agi::Migration* plan, Cavity cav) {
    Connection_Count conns;
    getCavityConnections(g,con,minCon,cav,conns);
    Connection_Count::iterator itr;
    for (itr=conns.begin();itr!=conns.end();itr++) {
      if (itr->second==2)
        continue;
      if (!plan->has(itr->first.first) && !plan->has(itr->first.second))
        return false;
    }
    return true;
  }

  double edgeCutGrowth(agi::Ngraph* g, Cavity cav, part_t peer) {
    int cut_pins = 0;
    int uncut_pins = 0;
    for (size_t i = 0; i < cav.size(); ++i) {
      agi::EdgeIterator* eitr = g->edges(cav[i]);
      agi::GraphEdge* e;
      while ((e = g->iterate(eitr))) {
        if (g->isResidentOn(e,peer))
          ++cut_pins;
        else
          ++uncut_pins;
      }
      g->destroy(eitr);
    }
    return uncut_pins * 1.0 / cut_pins;
  }

  /**
   * \brief a cavity is defined as the vertices adjacent to a given hyperedge
   *        that are owned and not in the plan
   * \param plan (In) size=numVerts, plan(i) = -1   : not in plan
   *                                         >= 0   : otherwise
   */
  void getCavities(int numEdges, int numVerts, int numGhostVerts,
      LIDs pin_degree, LIDs pins, LIDs plan, LIDs isVtxOwned,
      CSR& cavs) {
    //create the plan mask: 1 if in plan, 0 otherwise
    LIDs isVtxInPlan("isVtxInPlan", numVerts+numGhostVerts); //make this random access
    Kokkos::parallel_for(numVerts, KOKKOS_LAMBDA(const int v) {
      isVtxInPlan(v) = ( plan(v) != -1 );
    });

    //count the size of each cavity
    //loop over the edges - the last entry is to support the next operation where
    // we reuse this array to store the CSR offsets
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int degree = pin_degree(e+1)-pin_degree(e);
      const int firstPin = pin_degree(e);
      int size = 0;
      for(int pinIdx = 0; pinIdx < degree; pinIdx++) { //TODO nested parallel reduce?
        const int v = pins(firstPin+pinIdx);
        size += ( isVtxOwned(v) && !isVtxInPlan(v) ); //random access
      }
      cavs.off(e) = size;
    });
    //construct the cavity csr
    // -create the offsets with an exclusive scan
    degreeToOffset(cavs);
    allocateItems(cavs);
    // -fill the csr list view
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int degree = pin_degree(e+1)-pin_degree(e);
      const int firstPin = pin_degree(e);
      const int cavIdx = cavs.off(e);
      for(int pinIdx = 0; pinIdx < degree; pinIdx++) {
        const int v = pins(firstPin+pinIdx);
        if ( isVtxOwned(v) && !isVtxInPlan(v) ) { //random access
          cavs.items(cavIdx+pinIdx) = v;
        }
      }
    });
  }

  /**
   * \brief for each cavity get the list of processes that share
   *        edges adjacent to the cavity vertices
   * \remark the list of processes will contain duplicates
   */
  void getPeers(int numEdges, int numVerts, int numGhostVerts,
      LIDs edge_degree, LIDs edges,
      LIDs pin_degree, LIDs pins,
      LIDs isVtxOwned, LIDs vtxOwner,
      CSR& cavs, CSR& peers) {
    CSR eoc("edgesOfCavity",numEdges);
    //count the number of edges adjacent to cavity vertices
    // there will be duplicate edges; that should not change
    // the list of peers
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int numAdjVtx = cavs.off(e+1)-cavs.off(e);
      const int firstVtx = cavs.off(e);
      for(int vtxIdx = 0; vtxIdx < numAdjVtx; vtxIdx++) { //parallel reduce on the count?
        const int v = cavs.items(firstVtx+vtxIdx);
        eoc.off(e) += edge_degree(v+1)-edge_degree(v);
      }
    });
    degreeToOffset(eoc);
    allocateItems(eoc);
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int numAdjVtx = cavs.off(e+1)-cavs.off(e);
      const int firstVtx = cavs.off(e);
      int itemIdx = eoc.off(e);
      for(int vtxIdx = 0; vtxIdx < numAdjVtx; vtxIdx++) {
        const int v = cavs.items(firstVtx+vtxIdx);
        const int numAdjEdges = edge_degree(v+1)-edge_degree(v);
        const int firstEdge = edge_degree(v);
        for(int edgeIdx = 0; edgeIdx < numAdjEdges; edgeIdx++) {
          eoc.items(itemIdx) = edges(firstEdge+edgeIdx);
          itemIdx++;
        }
      }
    });
    //count the number of non-owned (ghost) vertices adj to each cavity edge
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int numCavEdges = eoc.off(e+1)-eoc.off(e);
      const int firstEdge = eoc.off(e);
      for(int edgeIdx = 0; edgeIdx < numCavEdges; edgeIdx++) { //parallel reduce?
        const int adjEdge = eoc.items(firstEdge+edgeIdx);
        const int pinDegree = pin_degree(adjEdge+1)-pin_degree(adjEdge);
        const int firstPin = pins(adjEdge);
        for(int pinIdx = 0; pinIdx < pinDegree; pinIdx++) {
          const int pin = pins(firstPin+pinIdx);
          peers.off(e) += !isVtxOwned(pin); //TODO make isVtxOwned random access
        }
      }
    });
    degreeToOffset(peers);
    allocateItems(peers);
    //insert the owners of the ghost vertices
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int numCavEdges = eoc.off(e+1)-eoc.off(e);
      const int firstEdge = eoc.off(e);
      int itemIdx = peers.off(e);
      for(int edgeIdx = 0; edgeIdx < numCavEdges; edgeIdx++) {
        const int adjEdge = eoc.items(firstEdge+edgeIdx);
        const int pinDegree = pin_degree(adjEdge+1)-pin_degree(adjEdge);
        const int firstPin = pins(adjEdge);
        for(int pinIdx = 0; pinIdx < pinDegree; pinIdx++) {
          const int pin = pins(firstPin+pinIdx);
          const int owner = vtxOwner(pin); //TODO make vtxOwner random access
          if(!isVtxOwned(pin)) { //TODO make isVtxOwned random access
            peers.items(itemIdx) = owner;
            itemIdx++;
          }
        }
      }
    });
  }

  void Selector::getCavitiesAndPeers(agi::etype t, LIDs plan,
      CSR& cavs, CSR& peers) {
    agi::PNgraph* pg = g->publicize();
    const int N = pg->num_local_edges[t];
    const int M = pg->num_local_verts;
    const int G = pg->num_ghost_verts;
    LIDs edge_degree("edge_degree", M+1); //convert to CSR
    LIDs edges("edges", pg->degree_list[t][M]); //convert to CSR
    LIDs pin_degree("pin_degree", N+1); //convert to CSR
    LIDs pins("pins", pg->pin_degree_list[t][N]); //convert to CSR
    hostToDevice(edge_degree, pg->degree_list[t]);
    hostToDevice(edges, pg->edge_list[t]);
    hostToDevice(pin_degree, pg->pin_degree_list[t]);
    hostToDevice(pins, pg->pin_list[t]);
    //create the ownership mask - this can be moved outside the select function
    LIDs isVtxOwned("isVtxOwned", M+G); //make this random access
    Kokkos::parallel_for(M, KOKKOS_LAMBDA(const int v) {
      isVtxOwned(v) = 1;
    });
    getCavities(N,M,G,pin_degree,pins,plan,isVtxOwned,cavs);
    //create the owners mask - this can be moved outside the select function
    LIDs vtxOwner("vtxOwner", M+G); //make this random access
    LIDs ghostOwners("ghostOwners", G);
    hostToDevice(ghostOwners,pg->owners);
    ENGPAR_LID_T self = PCU_Comm_Self();
    LIDs processId("processId", 1);
    hostToDevice(processId,&self);
    hostToDevice(pins, pg->pin_list[t]);
    Kokkos::parallel_for(M, KOKKOS_LAMBDA(const int v) {
      vtxOwner(v) = processId(0);
    });
    Kokkos::parallel_for(G, KOKKOS_LAMBDA(const int v) {
      vtxOwner(v+M) = ghostOwners(v);
    });
    getPeers(N,M,G,edge_degree,edges,pin_degree,pins,
        isVtxOwned,vtxOwner,cavs,peers);
  }

  /**
   * \brief the entry for a cavity is set to 1 if its
   *        color is 'color', 0 otherwise
   */
  LIDs buildColorMask(LIDs colors, int color) {
    const int numEdges = colors.dimension_0();
    LIDs colorMask("colorMask", numEdges);
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      colorMask(e) = ( colors(e) == color );
    });
    return colorMask;
  }

  /**
   * \brief the entry for a cavity is set to 1 if it is smaller
   *        than maxSize, 0 otherwise
   */
  LIDs buildSizeMask(CSR cavs, int maxSize) {
    const int numEdges = cavs.n;
    LIDs sizeMask("sizeMask", numEdges);
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int size = cavs.off(e+1)-cavs.off(e);
      sizeMask(e) = ( size < maxSize );
    });
    return sizeMask;
  }

  wgt_t Selector::kkSelect(Targets* targets, agi::Migration* plan,
                         wgt_t planW, unsigned int cavSize,int target_dimension) {
#ifdef KOKKOS_ENABLED
    agi::etype t = target_dimension;
    agi::PNgraph* pg = g->publicize();
    const int N = pg->num_local_edges[t];
    const int M = pg->num_local_verts;
    const int G = pg->num_ghost_verts;
    engpar::ColoringInput* in = engpar::createBdryColoringInput(g, t);
    agi::lid_t numColors;
    LIDs colors = engpar::EnGPar_KokkosColoring(in, numColors);
    LIDs plan_view("plan_view", M+G);
    //loop over colors
    for(agi::lid_t c=1; c<=numColors; c++) {
      if (planW > targets->total()) break;
      //build all the cavities and peers
      CSR cavs("cavities",N);
      CSR peers("cavityPeers",N);
      getCavitiesAndPeers(t,plan_view,cavs,peers);
      //parallel for each cavity: compute cavity color mask
      LIDs colorMask = buildColorMask(colors,c);
      //parallel for each cavity: compute cavity size mask
      LIDs sizeMask = buildSizeMask(cavs,cavSize);
      //parallel for each cavity: compute edgecutgrowth mask; size = sum_edges(peers(edge))
      //parallel for each cavity: compute neighbor target mask
      //parallel for each cavity: compute neighbor sending mask
      //parallel sort of selected cavities based on distance -> selection mask
      //create numEdges sized array of ints to store plan
      //parallel for each cavity: use logical and of masks to compute peer for each cavity entity
    }
    //build Migration object from plan array
    return planW;
#else
    fprintf(stderr,"ERROR: kokkos cavity selection disabled, recompile with ENABLE_KOKKOS\n");
    exit(1);
#endif
  }

  wgt_t Selector::select(Targets* targets, agi::Migration* plan,
                         wgt_t planW, unsigned int cavSize,int target_dimension) {
    q->startIteration();
    Queue::iterator itr;
    for (itr = q->begin();itr!=q->end();itr++) {
      if (planW > targets->total()) break;
      //Create Cavity and peers
      Cavity cav;
      Peers peers;
      getCavity(g,q->get(itr),plan,cav,peers);
      bool sent = false;
      if (cav.size() < cavSize) { //If the cavity is small enough
        Peers::iterator pitr;
        for (pitr = peers.begin(); pitr != peers.end(); ++pitr) {
          part_t peer = *pitr;
          if (targets->has(peer) && // Targeting this neighbor
              sending[*pitr]<targets->get(peer) && //Havent sent too much weight to this peer
              (in->limitEdgeCutGrowth <= 0 ||
               edgeCutGrowth(g, cav, peer) < in->limitEdgeCutGrowth)) {
                //add cavity to plan
                wgt_t w = addCavity(g,cav,peer,plan,target_dimension);
                planW+=w;
                sending[peer]+=w;
                sent=true;
                break;
          }
        }
      }
      if (!sent) {
        q->addElement(q->get(itr));
      }
    }
    return planW;
  }

  void Selector::selectDisconnected(agi::Migration* plan, int target_dimension) {
    q->startIteration();
    Queue::iterator itr;
    for (itr = q->begin();itr!=q->end();itr++) {
      //Create Cavity and peers
      Cavity cav;
      Peers peers;
      getCavity(g,q->get(itr),plan,cav,peers);
      if (in->minConnectivity > 1 &&
          isPartiallyConnected(g,in->connectivityType,in->minConnectivity,plan,cav)) {
        addCavity(g,cav,*peers.begin(),plan,target_dimension);
      }
    }
  }
  void Selector::insertInteriorEdges(agi::GraphVertex* vtx, agi::part_t dest,
                                     EdgeSet& edges, int dim) {
    
    agi::EdgeIterator* eitr = g->edges(vtx,dim);
    agi::GraphEdge* e;
    while ((e = g->iterate(eitr))) {
      if (!g->isResidentOn(e,dest))
        edges.insert(e);
    }
    g->destroy(eitr);
  }
  void Selector::tempInsertInteriorEdges(agi::GraphVertex* vtx, agi::part_t dest,
                                         EdgeSet& tmpEdges, int dim,
                                         const EdgeSet& edges) {
    agi::EdgeIterator* eitr = g->edges(vtx,dim);
    agi::GraphEdge* e;
    while ((e = g->iterate(eitr))) {
      if (edges.find(e)!=edges.end())
        continue;
      if (!g->isResidentOn(e,dest))
        tmpEdges.insert(e);
    }
    g->destroy(eitr);
  }

  double Selector::weight(const EdgeSet& edges) {
    EdgeSet::const_iterator itr;
    double w=0;
    for (itr=edges.begin();itr!=edges.end();itr++) 
      w+=g->weight(*itr);
    
    return w;
  }
  struct Migr {
    Migr(int i) :
      id(i), edge_types(0) {}
    void addWeight(double w) {ws[edge_types++]=w;}
    int id;
    double ws[5];
    int edge_types;
  };
  struct CompareMigr {
    //sort by part first completed edge type weights
    bool operator()(const Migr& a, const Migr& b) const {
      if( a.ws[0] < b.ws[0] )
        return true;
      if (a.ws[0]==b.ws[0])
        return a.id<b.id;
      return false;
    }
  };

  void Selector::calculatePlanWeights(agi::Migration* plan, std::unordered_map<int,double>& vtx_weight, PeerEdgeSet* peerEdges, std::unordered_set<part_t>& neighbors) {
    agi::Migration::iterator itr;
    for(itr = plan->begin();itr!=plan->end();itr++) {
      agi::GraphVertex* vtx = *itr;
      const int dest = plan->get(vtx);
      neighbors.insert(dest);
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        if (completed_dimensions->at(i)!=-1)
          insertInteriorEdges(vtx, dest, peerEdges[i][dest],
                              completed_dimensions->at(i));
        else
          vtx_weight[dest]+=g->weight(vtx);
      }
    }

  }

  void Selector::sendPlanWeight(std::unordered_map<int,double>& vtx_weight, PeerEdgeSet* peerEdges, std::unordered_set<part_t>& neighbors) {
    std::unordered_set<part_t>::iterator sitr;
    for (sitr = neighbors.begin();sitr!=neighbors.end();sitr++) {
      const int dest = *sitr;
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        double w;
        if (completed_dimensions->at(i)!=-1)
          w = weight(peerEdges[i][dest]);
        else
          w = vtx_weight[dest];
        PCU_COMM_PACK(dest, w);
      }
    }
    delete [] peerEdges;
  }

  void Selector::receiveIncomingWeight(MigrComm& incoming) {
    double w;
    while (PCU_Comm_Receive()) {
      Migr migr(PCU_Comm_Sender());
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        PCU_COMM_UNPACK(w);
        migr.addWeight(w);
      }
      incoming.insert(migr);        
    }

  }
  bool Selector::determineAvailability(Ws& avail) {
    bool isAvail = true;
    for (unsigned int i=0;i<completed_dimensions->size();i++) {
      double totW = getWeight(g,completed_dimensions->at(i),in->countGhosts);
      avail[i] = completed_weights->at(i) - totW;
      if (avail[i]<0)
        isAvail=false;
    }
    return isAvail;
  }

  void Selector::acceptWeight(MigrComm& incoming, bool& isAvail, Ws& avail, Midd& accept) {
    MigrComm::iterator in;
    for (in=incoming.begin();in!=incoming.end();in++) {
      const int nbor = (*in).id;
      if(isAvail) {
        bool hasSpace= true;
        for (unsigned int i=0;i<completed_dimensions->size();i++)
          if ((*in).ws[i] > avail[i])
            hasSpace=false;
        accept[nbor] = new double[completed_dimensions->size()];
        if( hasSpace ) {
          for (unsigned int i=0;i<completed_dimensions->size();i++) 
            isAvail = (avail[i]-=accept[nbor][i] = (*in).ws[i])>0;
        } else {
          for (unsigned int i=0;i<completed_dimensions->size();i++)
            isAvail = (avail[i] -=accept[nbor][i] = avail[i])>0;
        }
      } else {
        accept[nbor] = new double[completed_dimensions->size()];        
        for (unsigned int i=0;i<completed_dimensions->size();i++)
          accept[nbor][i] = 0;
      }
    }
  }

  void Selector::sendAcceptedWeights(Midd& accept) {
    Midd::iterator acc; 
    for (acc=accept.begin();acc!=accept.end();acc++) {
      for (unsigned int i=0;i<completed_dimensions->size();i++) 
        PCU_COMM_PACK(acc->first, acc->second[i]);
      delete [] acc->second;
    }
  }

  void Selector::gatherCapacities(Midd* capacity) {
    double w;
    while (PCU_Comm_Receive()) {
      int nbor = PCU_Comm_Sender();
      capacity->insert(std::make_pair(nbor,
                                      new double[completed_dimensions->size()]));
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        PCU_COMM_UNPACK(w);
        (*capacity)[nbor][i] = w;
      }
    }
  }
  //return  map<int neighbor, pair<double vtxW, double edgeW> > where
  //  neighbor is a neighbor's part id
  //  vtxW is the vtx weight capacity of neighbor
  //  edgeW is the edge weight capacity of neighbor
  Midd* Selector::trim(Targets*, agi::Migration* plan) {
    //compute the weight of the vertices and edges being sent to each peer
    std::unordered_set<part_t> neighbors;
    PeerEdgeSet* peerEdges = new PeerEdgeSet[completed_dimensions->size()];
    std::unordered_map<int,double> vtx_weight;
    calculatePlanWeights(plan,vtx_weight,peerEdges,neighbors);
    
    //send vtx and edge weight
    PCU_Comm_Begin();
    sendPlanWeight(vtx_weight,peerEdges,neighbors);
    PCU_Comm_Send();
    MigrComm incoming;
    receiveIncomingWeight(incoming);


    Ws avail = new double[completed_dimensions->size()];
    bool isAvail=determineAvailability(avail);

    Midd accept;
    acceptWeight(incoming,isAvail,avail,accept);

    PCU_Barrier();
    delete [] avail;

    PCU_Comm_Begin();
    sendAcceptedWeights(accept);
    PCU_Comm_Send();
    Midd* capacity = new Midd;
    gatherCapacities(capacity);

    return capacity;
  }

  void Selector::combineSets(EdgeSet& edges,const EdgeSet& tempEdges) {
    EdgeSet::const_iterator itr;
    for (itr=tempEdges.begin();itr!=tempEdges.end();itr++) 
      edges.insert(*itr);
    
  }
  

  void Selector::cancel(agi::Migration*& plan, Midd* capacity) {
    typedef std::pair<agi::GraphVertex*, int> PlanPair;
    std::vector<PlanPair > keep; //temporary plan container
    keep.reserve(plan->size());

    //Plan is a vector so this iterates over the plan in the order they were
    //    selected in.
    PeerEdgeSet* peerEdges = new PeerEdgeSet[completed_dimensions->size()];
    std::unordered_map<int,double> vtx_weight;
    agi::Migration::iterator pitr;
    for(pitr = plan->begin();pitr!=plan->end();pitr++) {
      agi::GraphVertex* v = *pitr;
      int dest = plan->get(v);
      bool isSpace=true;
      EdgeSet* tmpEdges = new EdgeSet[completed_dimensions->size()];
      
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        if (completed_dimensions->at(i)!=-1) {
          tempInsertInteriorEdges(v, dest, tmpEdges[i],
                                  completed_dimensions->at(i),
                                  peerEdges[i][dest]);
          if(weight(peerEdges[i][dest])+weight(tmpEdges[i])>(*capacity)[dest][i])
            isSpace=false;
        }
        else {
          if (vtx_weight[dest]+ g->weight(v) > (*capacity)[dest][i])
            isSpace = false;
        }
      }
      if (isSpace) {
        keep.push_back(PlanPair(v,dest));
        for (unsigned int i=0;i<completed_dimensions->size();i++) {
          if (completed_dimensions->at(i)!=-1)
            combineSets(peerEdges[i][dest],tmpEdges[i]);
          else
            vtx_weight[dest]+=g->weight(v);
        }
      }
      delete [] tmpEdges;
    }
    Midd::iterator acc; 
    for (acc=capacity->begin();acc!=capacity->end();acc++) {
      delete [] acc->second;
    }
    delete capacity;
    delete [] peerEdges;
    plan->clear();
    for(size_t i=0; i < keep.size(); i++)
      plan->insert(keep[i]);
  }

  void Selector::calculatePlanWeight(agi::Migration* plan, int target_dimension,
                                     std::unordered_map<part_t,wgt_t>& weights) {
    PeerEdgeSet peerEdges;
    agi::Migration::iterator itr;
    for(itr = plan->begin();itr!=plan->end();itr++) {
      agi::GraphVertex* vtx = *itr;
      const int dest = plan->get(vtx);
      if (target_dimension!=-1)
        insertInteriorEdges(vtx, dest, peerEdges[dest],
                            target_dimension);
      else
        weights[dest]+=g->weight(vtx);
    }
    if (target_dimension!=-1) {
      PeerEdgeSet::iterator pitr;
      for (pitr = peerEdges.begin(); pitr != peerEdges.end(); pitr++) {
        weights[pitr->first] = weight(pitr->second);
      }
    }
  }

  void Selector::updatePartWeight(agi::Migration* plan, int target_dimension, agi::WeightPartitionMap* wp_map) {
    std::unordered_map<part_t,wgt_t> weights;
    calculatePlanWeight(plan,target_dimension,weights);

    agi::WeightPartitionMap::iterator itr;
    for (itr = wp_map->begin(); itr != wp_map->end(); itr++) {
      std::unordered_map<part_t,wgt_t>::iterator inner_itr;
      for (inner_itr = weights.begin(); inner_itr != weights.end(); inner_itr++) {
        itr->second[inner_itr->first] -= inner_itr->second;
      }
    }
  }

  int degreeFunc(agi::Ngraph* g, agi::GraphVertex* v) {
    return g->degree(v);
  }

  Selector::Selector(DiffusiveInput* in_, Queue* queue,
           std::vector<int>* cd, std::vector<double>* cw) :
    in(in_), g(in_->g),
    q(queue),
    completed_dimensions(cd), completed_weights(cw) {
  }

  Selector* makeSelector(DiffusiveInput* in,Queue* q,
                         std::vector<int>* cd,
                         std::vector<double>* cw ) {
    return new Selector(in,q,cd,cw);
  }
}
