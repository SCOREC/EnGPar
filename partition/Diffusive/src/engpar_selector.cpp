#include "engpar_selector.h"
#include <engpar_metrics.h>
#include <engpar_kokkosColoring.h>
#include <engpar_support.h>
#include <set>
#include <PCU.h>
#include <agiMigration.h>
#include <agi_typeconvert.h>
#include <engpar_reduce.h>

#include <torch/script.h>
#include <ATen/ATen.h>

#define DEBUG_RANK 0
#define DEBUG_EDGE 50
#define DEBUG_KK 0

/** \brief define an upper limit on the number of remote 
 *          processes a hyperedge can exist on
 */
#define MAX_PEERS 80


//Adapted from https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
namespace std {
  template <>
  struct hash<engpar::Connection> {
    size_t operator()(const engpar::Connection& k) const {
      return ((hash<agi::GraphVertex*>()(k.first)
	       ^ (hash<agi::GraphVertex*>()(k.second) << 1)) >> 1);
    }
  };
}

namespace engpar {
  void getCavity(agi::Ngraph* g, agi::GraphEdge* edge, Cavity& cav) {
    //Grab all vertices connected to the edge that are on part
    agi::PinIterator* pitr = g->pins(edge);
    agi::GraphVertex* vtx;
    while ((vtx = g->iterate(pitr)))
      if (g->owner(vtx)==PCU_Comm_Self())
	cav.insert(vtx);
    g->destroy(pitr);
  }
  
  void Selector::constructCavities() {
    q->startIteration();
    int cav_sizes[100];
    for (int i =0; i < 100; ++i)
      cav_sizes[i] = 0;
    Queue::iterator itr;
    for (itr = q->begin();itr!=q->end();itr++) {
      Cavity cav;
      getCavity(g, q->get(itr), cav);
      cavities.insert(std::make_pair(q->get(itr),cav));
      q->addElement(q->get(itr));
      cav_sizes[cav.size()]++;
    }
  }
  
  void Selector::updateCavities(agi::Migration* plan) {
    std::unordered_map<agi::GraphEdge*, Cavity>::iterator c_itr;
    for (c_itr = cavities.begin(); c_itr != cavities.end(); ++c_itr) {
      agi::Migration::iterator p_itr;
      for (p_itr = plan->begin(); p_itr != plan->end(); ++p_itr) {
	c_itr->second.erase(*p_itr);
      }
    }
  }
  
  void getCavityPeers(agi::Ngraph* g, agi::Migration* plan,
		      Cavity& cav, Peers& peers) {

    //Grab all vertices connected to the edge that are on part
    typedef std::unordered_map<agi::part_t, std::unordered_set<agi::GraphEdge*> > Peer_Map;
    Peer_Map peerMap;
    EdgeSet edgesOfCavity;
    Cavity::iterator c_itr;
    //Add all of the edges of the vertex to a set
    for (c_itr = cav.begin(); c_itr != cav.end(); ++c_itr) {
      agi::GraphEdge* e;
      agi::EdgeIterator* eitr = g->edges(*c_itr);
      while ((e=g->iterate(eitr)))
	edgesOfCavity.insert(e);
      g->destroy(eitr);
    }

    //Locate all of the cut edges around the cavity
    EdgeSet::iterator sitr;
    agi::GraphVertex* v;
    agi::PinIterator* pitr;
    for (sitr = edgesOfCavity.begin(); sitr != edgesOfCavity.end(); sitr++) {
      pitr = g->pins(*sitr);
      while ((v = g->iterate(pitr))) {
        if (g->owner(v)!=PCU_Comm_Self())
          peerMap[g->owner(v)].insert(*sitr);
	else if (plan->has(v))
	  peerMap[plan->get(v)].insert(*sitr);
      }
      g->destroy(pitr);
    }

    //Reverse the Peer Map
    std::map<agi::lid_t, agi::part_t> edges_to_part;
    Peer_Map::iterator itr;
    for (itr = peerMap.begin(); itr != peerMap.end(); itr++)
      edges_to_part.insert(std::make_pair(itr->second.size(),itr->first));
    
    //Order the peers by largest surface area with cavity
    //  Iterate map backwards
    std::map<agi::lid_t, agi::part_t>::iterator etp_itr;
    for (etp_itr = edges_to_part.end(); etp_itr != edges_to_part.begin();) {
      --etp_itr;
      peers.push_back(etp_itr->second);
    }
  }

  wgt_t addCavity(agi::Ngraph* g, Cavity& cav,
                  part_t peer, agi::Migration* plan,
                  int target_dimension) {
    Cavity::iterator itr;
    wgt_t w=0.0;
    EdgeSet target_edges;
    for (itr = cav.begin();itr!=cav.end();itr++) {
      if (plan->insert(std::make_pair(*itr,peer))) {
	if (target_dimension==-1)
	  w+= g->weight(*itr);
	else {
	  agi::EdgeIterator* eitr = g->edges(*itr, target_dimension);
	  agi::GraphEdge* e;
	  while ((e = g->iterate(eitr)))
	    target_edges.insert(e);
	  g->destroy(eitr);
	}
      }
    }
    EdgeSet::iterator sitr;
    for (sitr = target_edges.begin();sitr!=target_edges.end();++sitr) {
      agi::GraphEdge* e = *sitr;
      if (!g->isResidentOn(e,peer))
        w+=g->weight(e);
    }
    return w;
  }
  typedef std::unordered_map<Connection,int> Connection_Count;
  void getCavityConnections(agi::Ngraph* g,agi::etype con, int minCon, Cavity& cav,
                            Connection_Count& conns) {
    Cavity::iterator itr;
    for (itr = cav.begin(); itr != cav.end(); ++itr) {
      agi::GraphVertex* u = *itr;
      agi::GraphIterator* gitr = g->adjacent(u,con);
      agi::GraphVertex* v;
      std::unordered_map<agi::GraphVertex*, int> occ;
      while ((v = g->iterate(gitr))) {
        if (v==u) continue;
        occ[v]++;
      }
      g->destroy(gitr);
      std::unordered_map<agi::GraphVertex*, int>::iterator itr;
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
                            agi::Migration* plan, Cavity& cav) {
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

  double edgeCutGrowth(agi::Ngraph* g, Cavity& cav, part_t peer) {
    int cut_pins = 0;
    int uncut_pins = 0;
    Cavity::iterator itr;
    for (itr = cav.begin(); itr != cav.end(); ++itr) {
      agi::EdgeIterator* eitr = g->edges(*itr);
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
    auto cavs_off = cavs.off;
    Kokkos::parallel_for("asdf", numEdges, KOKKOS_LAMBDA(const int e) {
      const int degree = pin_degree(e+1)-pin_degree(e);
      const int firstPin = pin_degree(e);
      int size = 0;
      for(int pinIdx = 0; pinIdx < degree; pinIdx++) { //TODO nested parallel reduce?
        const int v = pins(firstPin+pinIdx);
        size += ( isVtxOwned(v) && !isVtxInPlan(v) ); //random access
      }
      cavs_off(e) = size;
#if DEBUG_KK==1
      if( e == DEBUG_EDGE )
        printf("e cav.degree %4d %3d\n", e, size);
#endif
    });
    //construct the cavity csr
    // -create the offsets with an exclusive scan
    degreeToOffset(cavs);
    allocateItems(cavs);
    // -fill the csr list view
    auto cavs_items = cavs.items;
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int degree = pin_degree(e+1)-pin_degree(e);
      const int firstPin = pin_degree(e);
      const int cavIdx = cavs_off(e);
      for(int pinIdx = 0; pinIdx < degree; pinIdx++) {
        const int v = pins(firstPin+pinIdx);
        if ( isVtxOwned(v) && !isVtxInPlan(v) ) { //random access
          cavs_items(cavIdx+pinIdx) = v;
        }
      }
    });
  }

  //need to do something to inline these and make them device callable
  typedef struct PeerList {
    int n;
    int e[MAX_PEERS];
  } peerList;

  KOKKOS_INLINE_FUNCTION
  peerList makePeerList() {
    peerList p;
    p.n = 0;
    for(int i=0; i<MAX_PEERS; i++) {
      p.e[i] = -1;
    }
    return p;
  }

  KOKKOS_INLINE_FUNCTION
  int find(peerList& a, int v) {
    assert(v != -1);
    int has = 0;
    for(int i=0; i<MAX_PEERS; i++) {
      has += (v == a.e[i]);
    }
    return has;
  }

  KOKKOS_INLINE_FUNCTION
  void unionPeers(peerList& a, peerList& b) {
    for(int i=0; i<b.n; i++) {
      const int v = b.e[i];
      const int has = find(a,v); 
      if( has == 0 ) {
        a.e[a.n] = v;
        a.n++;
        assert(a.n < MAX_PEERS);
      }
    }
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
      CSR& cavs, CSR& peers, CSR& eoc) {
    const int self = PCU_Comm_Self();
    //count the number of edges adjacent to cavity vertices
    // there will be duplicate edges; that should not change
    // the list of peers
    auto cavs_off = cavs.off;
    auto cavs_items = cavs.items;
    auto eoc_off = eoc.off;
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int numAdjVtx = cavs_off(e+1)-cavs_off(e);
      const int firstVtx = cavs_off(e);
      for(int vtxIdx = 0; vtxIdx < numAdjVtx; vtxIdx++) { //parallel reduce on the count?
        const int v = cavs_items(firstVtx+vtxIdx);
        eoc_off(e) += edge_degree(v+1)-edge_degree(v);
      }
    });
    degreeToOffset(eoc);
    allocateItems(eoc);
    auto eoc_items = eoc.items;
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int numAdjVtx = cavs_off(e+1)-cavs_off(e);
      const int firstVtx = cavs_off(e);
      int itemIdx = eoc_off(e);
      for(int vtxIdx = 0; vtxIdx < numAdjVtx; vtxIdx++) {
        const int v = cavs_items(firstVtx+vtxIdx);
        const int numAdjEdges = edge_degree(v+1)-edge_degree(v);
        const int firstEdge = edge_degree(v);
        for(int edgeIdx = 0; edgeIdx < numAdjEdges; edgeIdx++) {
          eoc_items(itemIdx) = edges(firstEdge+edgeIdx);
          itemIdx++;
        }
      }
    });
    //count the number of non-owned (ghost) vertices adj to each cavity edge
    auto peers_off = peers.off;
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      peerList cavPeers = makePeerList();
      const int numCavEdges = eoc_off(e+1)-eoc_off(e);
      const int firstEdge = eoc_off(e);
      for(int edgeIdx = 0; edgeIdx < numCavEdges; edgeIdx++) { //parallel reduce?
        const int adjEdge = eoc_items(firstEdge+edgeIdx);
        const int pinDegree = pin_degree(adjEdge+1)-pin_degree(adjEdge);
        const int firstPin = pin_degree(adjEdge);
        peerList vtxPeers = makePeerList();
        for(int pinIdx = 0; pinIdx < pinDegree; pinIdx++) {
          const int pin = pins(firstPin+pinIdx);
          assert(pin>=0 && pin<numVerts+numGhostVerts);
          const int owner = vtxOwner(pin); //TODO make vtxOwner random access
          if(!isVtxOwned(pin)) { //TODO make isVtxOwned random access
            vtxPeers.e[vtxPeers.n] = owner;
            vtxPeers.n++;
            assert(vtxPeers.n<MAX_PEERS);
          }
        }
        unionPeers(cavPeers,vtxPeers); 
      }
      peers_off(e) = cavPeers.n;
    });
    degreeToOffset(peers);
#if DEBUG_KK==1
    if(self == DEBUG_RANK) {
      Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
        if( e == DEBUG_EDGE) {
          printf("e cavDegree eocDegree peerDegree %5d %5d %5d %5d\n",
            e,
            cavs_off(e+1)-cavs_off(e),
            eoc_off(e+1)-eoc_off(e),
            peers_off(e+1)-peers_off(e));
        }
      });
    }
#endif
    allocateItems(peers);
    //insert the owners of the ghost vertices
    auto peers_items = peers.items;
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      peerList cavPeers = makePeerList();
      const int numCavEdges = eoc_off(e+1)-eoc_off(e);
      const int firstEdge = eoc_off(e);
      for(int edgeIdx = 0; edgeIdx < numCavEdges; edgeIdx++) {
        const int adjEdge = eoc_items(firstEdge+edgeIdx);
        const int pinDegree = pin_degree(adjEdge+1)-pin_degree(adjEdge);
        const int firstPin = pin_degree(adjEdge);
        peerList vtxPeers = makePeerList();
        for(int pinIdx = 0; pinIdx < pinDegree; pinIdx++) {
          const int pin = pins(firstPin+pinIdx);
          const int owner = vtxOwner(pin); //TODO make vtxOwner random access
          if(!isVtxOwned(pin)) { //TODO make isVtxOwned random access
            vtxPeers.e[vtxPeers.n] = owner;
            vtxPeers.n++;
            assert(vtxPeers.n<MAX_PEERS);
          }
        }
        unionPeers(cavPeers,vtxPeers); 
      }
      const int itemIdx = peers_off(e);
      assert( cavPeers.n == peers_off(e+1)-peers_off(e) );
      for(int i=0; i<cavPeers.n; i++) {
        peers_items(itemIdx+i) = cavPeers.e[i];
      }
    });
  }

  void Selector::getCavitiesAndPeers(agi::etype t,
      LIDs plan, LIDs vtxOwner,
      CSR& cavs, CSR& peers, CSR& eoc) {
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
    getPeers(N,M,G,edge_degree,edges,pin_degree,pins,
        isVtxOwned,vtxOwner,cavs,peers,eoc);
  }

  /**
   * \brief the entry for a cavity is set to 1 if its
   *        color is 'color', 0 otherwise
   */
  LIDs buildColorMask(LIDs colors, int color) {
    const int numEdges = colors.extent(0);
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
  LIDs buildSizeMask(CSR cavs, const int maxSize) {
    const int numEdges = cavs.n;
    LIDs sizeMask("sizeMask", numEdges);
    auto cavs_off = cavs.off;
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int size = cavs_off(e+1)-cavs_off(e);
      sizeMask(e) = ( size < maxSize );
    });
    return sizeMask;
  }

  /**
   * \brief the entry for a cavity is > 0 if it is shared with
   *        tgtPeer, 0 otherwise
   */
  LIDs buildTargetMask(CSR peers, int tgtPeer) {
    const int numEdges = peers.n;
    LIDs targetMask("targetMask", numEdges);
    auto peers_off = peers.off;
    auto peers_items = peers.items;
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      const int size = peers_off(e+1)-peers_off(e);
      const int firstPeer = peers_off(e);
      for( int i=0; i<size; i++ ) {
        targetMask(e) |= ( peers_items(firstPeer+i) == tgtPeer );
      }
    });
    return targetMask;
  }

  /**
   * \brief the entry for a cavity is = 1 if the edge cut growth
   *        is less than maxGrowth, 0 otherwise
   */
  LIDs buildEdgeCutMask(CSR peers, int tgtPeer, float maxGrowth) {
    const int numEdges = peers.n;
    LIDs cutPins("edgeCutMask", numEdges);
    LIDs edgeCutMask("edgeCutMask", numEdges);
    if( maxGrowth <= 0 ) {
      // edgeCutGrowth disabled, set all cavities to 'true'
      Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
        edgeCutMask(e) = 1;
      });
      return edgeCutMask;
    } else {
      auto peers_off = peers.off;
      auto peers_items = peers.items;
      Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
        const int size = peers_off(e+1)-peers_off(e);
        const int firstPeer = peers_off(e);
        for( int i=0; i<size; i++ ) {
          cutPins(e) += ( peers_items(firstPeer+i) == tgtPeer );
        }
        const int uncutPins = size - cutPins(e);
        const float edgeCutGrowth = static_cast<float>(uncutPins)/cutPins(e);
        edgeCutMask(e) = ( edgeCutGrowth < maxGrowth );
      });
      return edgeCutMask;
    }
  }

  /**
   * \brief the entry for a cavity is set to 1 if it will 
   *        be migrated, 0 otherwise
   */
  LIDs buildMigrationMask(CSR cavs, LIDs colorMask, LIDs sizeMask,
      LIDs targetMask, LIDs edgeCutMask) {
    const int self = PCU_Comm_Self();
    const int numEdges = cavs.n;
    LIDs migrMask("migrationMask", numEdges);
    Kokkos::parallel_for(numEdges, KOKKOS_LAMBDA(const int e) {
      migrMask(e) = ( colorMask(e) && sizeMask(e) && targetMask(e) && edgeCutMask(e) );
    });
    return migrMask;
  }

  LIDs makePlan(std::string name, int size) {
    LIDs plan(name, size);
    Kokkos::parallel_for(size, KOKKOS_LAMBDA(const int v) {
      plan(v) = -1;
    });
    return plan;
  }

  /**
   * \brief the entry for a vertex is set to the destination process id
   */
  LIDs setVtxDestination(CSR cavs, LIDs migrMask, int numVerts, int tgtPeer) {
    assert(tgtPeer >= 0 && tgtPeer < PCU_Comm_Peers());
    const int self = PCU_Comm_Self();
    LIDs dest = makePlan("planNext", numVerts);
    auto cavs_off = cavs.off;
    auto cavs_items = cavs.items;
    Kokkos::parallel_for(cavs.n, KOKKOS_LAMBDA(const int e) {
      const int numAdjVerts = cavs_off(e+1)-cavs_off(e);
      const int firstVtx = cavs_off(e);
      for( int i=0; i<numAdjVerts; i++ ) {
        const int v = cavs_items(firstVtx+i);
        //using a conditional to keep the destination set
        //  to -1 if it is not being migrated
        if( migrMask(e) )
          dest(v) = tgtPeer;
      }
    });
    return dest;
  }

  /**
   * \brief return the weight of the vertices being migrated
   */
  wgt_t getVtxWeight(LIDs plan, WGTs vtxWeights,
      const wgt_t peerTgtW, const wgt_t tgtW) {
    const int nverts = plan.extent(0);
    wgt_t totWeight = 0;
    WGTs accW("accumulatedWeights", nverts);
    //inclusive scan over weights to find the saturation point for this peer or
    //the total
    Kokkos::parallel_scan(nverts,
      KOKKOS_LAMBDA (const int v, wgt_t& upd, const bool& final) {
      const wgt_t val = ( plan(v) >= 0 ) * vtxWeights(v);
      if (final) {
        accW(v) = upd;
      upd += val;
      }
    });
    //remove vertices from the plan if including them exceeds a weight limit
    //...this may cause poor cuts
    Kokkos::parallel_for(nverts, KOKKOS_LAMBDA(const int v) {
      const wgt_t w = accW(v);
      if( w > peerTgtW || w > tgtW )
        plan(v) = -1;
    });
    //sum the plan weight
    Kokkos::parallel_reduce(nverts,
      KOKKOS_LAMBDA(const int v, wgt_t& w) {
        w += ( plan(v) >= 0 ) * vtxWeights(v); }, totWeight);
    return totWeight;
  }

  /**
   * \brief return the weight of the edges that will be added to the part
   * boundary after migrating the cavity
   */
  wgt_t getEdgeWeight(CSR eoc, CSR pins,
      LIDs migrationMask, LIDs vtxOwner, 
      WGTs edgeWeights, int tgtPeer) {
    wgt_t totWeight = 0;
    auto eoc_off = eoc.off;
    auto eoc_items = eoc.items;
    auto pins_off = pins.off;
    auto pins_items = pins.items;
    Kokkos::parallel_reduce(eoc.n, KOKKOS_LAMBDA(const int e, wgt_t& w) {
      const int isMigrated = migrationMask(e);
      const int numCavEdges = eoc_off(e+1)-eoc_off(e);
      const int firstEdge = eoc_off(e);
      for(int edgeIdx = 0; edgeIdx < numCavEdges; edgeIdx++) {
        const int adjEdge = eoc_items(firstEdge+edgeIdx);
        const int pinDegree = pins_off(adjEdge+1)-pins_off(adjEdge);
        const int firstPin = pins_off(adjEdge);
        int residentOnPeer = 0;
        //The same vertex will be visited multiple times - consider forming list
        //of unique second adjacent vertices bounded by the cavity edge.
        //We have the following adjacencies,
        //edge1-verts1-edges2-verts2
        //where edge1 is the cavity edge that is cut by the part boundary, the set
        //of verts1 that will be migrated to peer, and the edges2 that will be on
        //the part boundary.  An edge in edges2 will already exist on peer if
        //one of its adjacent vertices is owned by peer.
        for(int pinIdx = 0; pinIdx < pinDegree; pinIdx++) {
          const int pin = pins_items(firstPin+pinIdx);
          const int owner = vtxOwner(pin); //TODO make vtxOwner random access
          residentOnPeer += (owner == tgtPeer);
        }
        //If the edge is not already resident on the peer (one of its adjacent
        //vertices is owned by the peer) then after migration the edge will
        //exist on the peer and its weight needs to be included.
        w += (!residentOnPeer && isMigrated) * edgeWeights(adjEdge);
      }
    }, totWeight);
    return totWeight;
  }

  LIDs buildVtxOwners(agi::PNgraph* pg) {
    const int M = pg->num_local_verts;
    const int G = pg->num_ghost_verts;
    LIDs vtxOwner("vtxOwner", M+G); //make this random access
    LIDs ghostOwners("ghostOwners", G);
    hostToDevice(ghostOwners,pg->owners);
    const int self = PCU_Comm_Self();
    Kokkos::parallel_for(M, KOKKOS_LAMBDA(const int v) {
      vtxOwner(v) = self;
    });
    Kokkos::parallel_for(G, KOKKOS_LAMBDA(const int v) {
      vtxOwner(v+M) = ghostOwners(v);
    });
    return vtxOwner;
  }

  void updatePlan(LIDs plan, LIDs planNext) {
    assert( plan.extent(0) == planNext.extent(0) );
    const int numPeers = PCU_Comm_Peers();
    const int n = plan.extent(0);
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int v) {
      if( planNext(v) != -1 )
        plan(v) = planNext(v);
    });
  }

  void buildMigrationPlan(agi::PNgraph* pg, LIDs plan, agi::Migration* migrPlan) {
    const int n = plan.extent(0);
    agi::lid_t* plan_h = new agi::lid_t[n];
    deviceToHost(plan,plan_h);
    for(int i=0; i<n; i++) {
      const int dest = plan_h[i];
      if( dest != -1 ) {
        if(dest <0 || dest >= PCU_Comm_Peers()) {
          printf("ERROR dest is not valid v dest %4d %4d\n", i, dest);
        }
        assert(dest >= 0 && dest <PCU_Comm_Peers());
        agi::GraphVertex* v = pg->getVertex(i);
        migrPlan->insert(std::make_pair(v,dest));
      }
    }
    delete [] plan_h;
  }

  wgt_t Selector::kkSelect(LIDs order, Targets* targets, agi::Migration* migrPlan,
                         wgt_t planW, unsigned int cavSize,int target_dimension) {
#ifdef KOKKOS_ENABLED
#if DEBUG_KK==1
    if( PCU_Comm_Self() == DEBUG_RANK ) {
      fprintf(stderr, "cavSize %3d\n", cavSize);
    }
#endif
    agi::etype t = target_dimension;
    agi::etype edgeType = t;
    if( edgeType == -1 )
      edgeType = 0;
    agi::PNgraph* pg = g->publicize();
    const int N = pg->num_local_edges[edgeType];
    const int M = pg->num_local_verts;
    WGTs vtxWeights("vtxWeights", M);
    hostToDevice(vtxWeights,pg->local_weights);
    WGTs edgeWeights("edgeWeights", N);
    hostToDevice(edgeWeights,pg->edge_weights[edgeType]);
    CSR pins("pins", N, pg->pin_degree_list[edgeType], pg->pin_list[edgeType]);
    LIDs vtxOwners = buildVtxOwners(pg);
    //create the owners mask - this can be moved outside the select function
    engpar::ColoringInput* inC = engpar::createBdryColoringInput(g, edgeType);
    agi::lid_t numColors;
    LIDs colors = engpar::EnGPar_KokkosColoring(inC, numColors);
    LIDs plan = makePlan("plan",M);

    const agi::lid_t maxPos = getMax(order);
    //printf("maxPos %d\n", maxPos);

    //for(agi::lid_t pos=0; pos
    //loop over colors
    for(agi::lid_t c=1; c<=numColors; c++) {
      if (planW > targets->total()) {
        //printf("sending enough total weight\n");
        break;
      }
      Targets::iterator tgt;
      for( tgt = targets->begin(); tgt != targets->end(); tgt++ ) {
        const int tgtPeer = tgt->first;
        if (planW > targets->total()) {
          //printf("sending enough total weight\n");
          break;
        }
#if DEBUG_KK==1
        if( PCU_Comm_Self() == DEBUG_RANK ) {
          fprintf(stderr, "color tgtPeer %3d %2d\n", c, tgtPeer);
        }
#endif
        const wgt_t tgtWeight = tgt->second;
        if( sending[tgtPeer] >= tgtWeight ) {
          //printf("sending enough weight to peer %d\n", tgtPeer);
          continue; //sent enough weight to this peer
        }
        //build all the cavities and peers
        CSR cavs("cavities",N);
        CSR peers("cavityPeers",N);
        CSR eoc("edgesOfCavity",N);
        getCavitiesAndPeers(edgeType,plan,vtxOwners,cavs,peers,eoc);
        //parallel for each cavity: compute cavity color mask
        LIDs colorMask = buildColorMask(colors,c);
        //parallel for each cavity: compute cavity size mask
        LIDs sizeMask = buildSizeMask(cavs,cavSize);
        //parallel for each cavity: compute target mask
        LIDs targetMask = buildTargetMask(peers,tgtPeer);
        //parallel for each cavity: compute edge cut growth mask
        LIDs edgeCutMask = buildEdgeCutMask(peers,tgtPeer,in->limitEdgeCutGrowth);
        //parallel for each cavity: logical and masks to determine if cavity is in plan
        LIDs migrationMask = buildMigrationMask(cavs,colorMask,sizeMask,targetMask,edgeCutMask);
        //parallel for each cavity: logical and masks to determine if cavity is in plan
        LIDs planNext = setVtxDestination(cavs,migrationMask,M,tgtPeer);
        //compute plan weight
        wgt_t w = 0;
        if( target_dimension == -1 ) { //vertices
          const wgt_t peerCap = tgtWeight-sending[tgtPeer];
          const wgt_t totCap = targets->total()-planW;
          w = getVtxWeight(planNext,vtxWeights,peerCap,totCap);
          if( w != 0 ) {
            printf("peer %d peerTgtW %f tgtW %f totWeight %f\n",
                tgtPeer, peerCap, totCap, w);
          }
        } else {
          w = getEdgeWeight(eoc,pins,migrationMask,vtxOwners,
              edgeWeights,tgtPeer);
          PCU_Debug_Print("peer %d totWeight %f\n", tgtPeer, w);
        }
        sending[tgtPeer] += w;
        planW += w;
        updatePlan(plan,planNext);
      }
    }
    //build Migration object from plan array
    buildMigrationPlan(pg,plan,migrPlan);
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
      if (planW > targets->total()) {
        //printf("sending enough total weight\n");
        break;
      }
      //Create Cavity and peers
      Cavity& cav = cavities[q->get(itr)];
      bool sent = false;
      if (cav.size() < cavSize) { //If the cavity is small enough
        Peers peers;
        getCavityPeers(g,plan,cav,peers);
        Peers::iterator pitr;
        for (pitr = peers.begin(); pitr != peers.end(); ++pitr) {
          part_t peer = *pitr;
          //if( targets->has(peer) && sending[peer] >= targets->get(peer) )
          //  printf("sending enough weight to peer %d\n", peer);
          if (targets->has(peer) && // Targeting this neighbor
              sending[peer]<targets->get(peer) && //Havent sent too much weight to this peer
              (in->limitEdgeCutGrowth <= 0 ||
               edgeCutGrowth(g, cav, peer) < in->limitEdgeCutGrowth)) {
                //add cavity to plan
                wgt_t w = addCavity(g,cav,peer,plan,target_dimension);
#if DEBUG_KK==1
                if( PCU_Comm_Self() == DEBUG_RANK ) {
                  printf("sending edge cav.size() peer %4d %3zu %3d\n", g->localID(q->get(itr)), cav.size(), peer);
                }
#endif
                planW+=w;
                sending[peer]+=w;
                sent=true;
                break;
          }
        }
      }
      if (!sent)
        q->addElement(q->get(itr));
      else
	cavities.erase(q->get(itr));
    }
    updateCavities(plan);
    return planW;
  }

  std::vector<torch::jit::IValue> getModelInput(agi::Ngraph* g, Cavity& cav, std::unordered_map<int,part_t>& idxToPeer) {
    // get cav edges and map to idx
    std::unordered_map<agi::GraphEdge*,int> cavEdgeToIdx(cav.size());
    for (const auto& graphVtx : cav) {
      agi::GraphEdge* e;
      agi::EdgeIterator* eitr = g->edges(graphVtx);
      while ((e=g->iterate(eitr))) {
        if (cavEdgeToIdx.find(e) == cavEdgeToIdx.end()) {
          cavEdgeToIdx[e] = cavEdgeToIdx.size();
        }
      }
      g->destroy(eitr);
    }
    if (cavEdgeToIdx.empty()) {
      return std::vector<torch::jit::IValue>();
    }

    // make adjacency matrix
    const int numNodes = cavEdgeToIdx.size();
    std::vector<int> adj(numNodes*numNodes, 0);
    for (const auto& graphVtx : cav) {
      agi::GraphEdge* e;
      agi::EdgeIterator* eitr = g->edges(graphVtx);
      while ((e=g->iterate(eitr))) {
        agi::GraphEdge* eB;
        agi::EdgeIterator* eitrB = g->edges(graphVtx);
        while ((eB=g->iterate(eitrB))) {
          if (eB != e) {
            adj[cavEdgeToIdx[e]*numNodes+cavEdgeToIdx[eB]]=1;
            adj[cavEdgeToIdx[eB]*numNodes+cavEdgeToIdx[e]]=1;
          }
        }
        g->destroy(eitrB);
      }
      g->destroy(eitr);
    }

    const int edgeCount = std::accumulate(adj.begin(), adj.end(), 0);

    // make edge list tensor from adjacency matrix
    at::Tensor edgeListTensor = at::empty({2, edgeCount}, at::kLong);
    int edgeListIdx = 0;
    for (int i = 0; i < numNodes; i++) {
      for (int j = 0; j < numNodes; j++) {
        const int idx = i*numNodes + j;
        if (adj[idx]) {
          edgeListTensor[0][edgeListIdx] = i;
          edgeListTensor[1][edgeListIdx] = j;
          edgeListIdx++;
        }
      }
    }

    // get processes of each vertex
    std::unordered_map<part_t, int> peerToIdx;
    std::vector<std::pair<int, int>> nodePeerPairs;
    for (const auto& kv : cavEdgeToIdx) {
      agi::Peers peers;
      g->getResidence(kv.first, peers);
      const int nodeIdx = kv.second;
      for (const part_t& peer : peers) {
        int peerIdx;
        const auto peerIdxItr = peerToIdx.find(peer);
        if (peerIdxItr == peerToIdx.end()) {
          peerIdx = peerToIdx.size();
          peerToIdx[peer] = peerIdx;
          idxToPeer[peerIdx] = peer;
        } else {
          peerIdx = peerIdxItr->second;
        }
        nodePeerPairs.push_back({nodeIdx, peerIdx});
      }
    }

    const int numProcs = peerToIdx.size();
    at::Tensor x = at::zeros({numNodes, numProcs}, at::kFloat);
    for (const auto& pair : nodePeerPairs) {
      x[pair.first][pair.second] = 1;
    }

    at::Tensor nodeSums = at::sum(x, 0, true).div(numNodes).expand({numNodes, -1});
    at::Tensor procSums = at::sum(x, 1, true).div(numProcs).expand({-1, numProcs});
    x = at::stack({x, nodeSums, procSums}, 2);

    std::vector<torch::jit::IValue> input = {x, edgeListTensor};
    return input;
  }

  wgt_t Selector::torchSelect(Targets* targets, agi::Migration* plan,
                              wgt_t planW, unsigned int cavSize,int target_dimension,
                              torch::jit::script::Module& model) {
    q->startIteration();
    std::unordered_map<int, part_t> idxToPeer;
    for (Queue::iterator itr = q->begin();itr!=q->end();itr++) {
      if (planW > targets->total()) {
        break;
      }
      Cavity& cav = cavities[q->get(itr)];
      bool sent = false;
      if (cav.size() < cavSize) {
        std::vector<torch::jit::IValue> input = getModelInput(g, cav, idxToPeer);
        if (!input.empty()) {
          at::Tensor output = model.forward(input).toTensor();
          int bestIdx = at::argmax(output).item<int>();
          part_t peer = idxToPeer[bestIdx];
          if (targets->has(peer) && // Targeting this neighbor
              sending[peer]<targets->get(peer) && //Havent sent too much weight to this peer
              (in->limitEdgeCutGrowth <= 0 ||
               edgeCutGrowth(g, cav, peer) < in->limitEdgeCutGrowth)) {
            //add cavity to plan
            wgt_t w = addCavity(g,cav,peer,plan,target_dimension);
            planW+=w;
            sending[peer]+=w;
            sent=true;
          }
        }
      }
      if (!sent)
        q->addElement(q->get(itr));
      else
	cavities.erase(q->get(itr));
    }
    updateCavities(plan);
    return planW;
  }


  void Selector::selectDisconnected(agi::Migration* plan, int target_dimension) {
    if (in->minConnectivity <= 1)
      return;
    q->startIteration();
    Queue::iterator itr;
    for (itr = q->begin();itr!=q->end();itr++) {
      //Create Cavity and peers
      Cavity& cav = cavities[q->get(itr)];
      Peers peers;
      getCavityPeers(g,plan,cav,peers);
      if (isPartiallyConnected(g,in->connectivityType,in->minConnectivity,plan,cav)) {
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

  void Selector::calculatePlanWeights(agi::Migration* plan,
				      std::unordered_map<int,double>& vtx_weight,
				      PeerEdgeSet* peerEdges,
				      std::unordered_set<part_t>& neighbors) {
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

  void Selector::sendPlanWeight(std::unordered_map<int,double>& vtx_weight,
				PeerEdgeSet* peerEdges, std::unordered_set<part_t>& neighbors) {
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
        for (unsigned int i = 0;i < completed_dimensions->size(); ++i)
          if ((*in).ws[i] > avail[i])
            hasSpace=false;
        accept[nbor] = new double[completed_dimensions->size()];
        if( hasSpace ) {
          for (unsigned int i = 0; i < completed_dimensions->size(); ++i)
            isAvail = (avail[i] -= accept[nbor][i] = (*in).ws[i])>0;
        } else {
          for (unsigned int i = 0; i < completed_dimensions->size(); ++i)
            isAvail = (avail[i] -= accept[nbor][i] = avail[i]) > 0;
        }
      } else {
        accept[nbor] = new double[completed_dimensions->size()];
        for (unsigned int i = 0; i < completed_dimensions->size(); ++i)
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
    for(pitr = plan->begin(); pitr != plan->end(); ++pitr) {
      agi::GraphVertex* v = *pitr;
      int dest = plan->get(v);
      bool isSpace=true;
      EdgeSet* tmpEdges = new EdgeSet[completed_dimensions->size()];
      
      for (unsigned int i = 0; i < completed_dimensions->size(); ++i) {
        if (completed_dimensions->at(i)!=-1) {
          tempInsertInteriorEdges(v, dest, tmpEdges[i],
                                  completed_dimensions->at(i),
                                  peerEdges[i][dest]);
          if(weight(peerEdges[i][dest]) + weight(tmpEdges[i])>(*capacity)[dest][i])
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

    constructCavities();
  }

  Selector* makeSelector(DiffusiveInput* in,Queue* q,
                         std::vector<int>* cd,
                         std::vector<double>* cw ) {
    return new Selector(in,q,cd,cw);
  }
}
