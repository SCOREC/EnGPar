#include "engpar_selector.h"
#include <engpar_metrics.h>
#include <set>
#include <PCU.h>
#include <agiMigration.h>

//Adapted from https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
namespace std {
  template <>
  struct hash<engpar::Connection> {
    std::size_t operator()(const engpar::Connection& k) const {
      return ((std::hash<agi::GraphVertex*>()(k.first)
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

  wgt_t Selector::select(Targets* targets, agi::Migration* plan,
                         wgt_t planW, unsigned int cavSize,int target_dimension) {
    q->startIteration();
    Queue::iterator itr;
    for (itr = q->begin();itr!=q->end();itr++) {
      if (planW > targets->total()) break;
      //Create Cavity and peers
      Cavity& cav = cavities[q->get(itr)];
      Peers peers;
      getCavityPeers(g,plan,cav,peers);
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
