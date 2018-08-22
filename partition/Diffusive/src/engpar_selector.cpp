#include "engpar_selector.h"
#include <engpar_metrics.h>
#include "engpar_version.h"
#include <set>
#include <PCU.h>
#include <agiMigration.h>

#include <iostream>
#include <fstream>
#include <sstream>

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

  void writeCavity(agi::Ngraph* g, Cavity& cav, part_t dest) {
    std::stringstream ss;
    ss << "cavities_" << PCU_Comm_Self() << ".txt";
    static std::ofstream outCav(ss.str());
    static int calls = 0;

    if (!calls) {
      outCav << "#engpar hash: " << engpar_version() << "\n";
      calls++;
    }

    outCav << "#num vertices\n";
    outCav << "1 " << cav.size() << "\n";
    // create adj matrix using the edges of the cavity
    // set<edges> edgesOfCavity
    // for vertex v in cavity 
    //   for each edge i in g->edges(v)
    //     edgesOfCavity.insert(i) 
    std::set<agi::GraphEdge*> edgesOfCavity;
    Cavity::iterator itr;
    for (itr = cav.begin();itr!=cav.end();itr++) {
      agi::GraphEdge* e;
      agi::EdgeIterator* eitr = g->edges(*itr);
      while ((e=g->iterate(eitr)))
        edgesOfCavity.insert(e);
    }
    // A = new int(size(edgesOfCavity)*size(edgesOfCavity))
    int numCavEdges = edgesOfCavity.size();
    int* A = new int[numCavEdges*numCavEdges];
    std::fill_n(A,numCavEdges*numCavEdges,0);
    // create a map from edge id to cavity idx
    std::map<agi::GraphEdge*,int> edgeToIdx;
    std::set<agi::GraphEdge*>::iterator sitr;
    int i = 0;
    for (sitr = edgesOfCavity.begin(); sitr != edgesOfCavity.end(); sitr++) {
      edgeToIdx[*sitr] = i++;
    }
    // for each cavity vertex v
    //   for each edge i in g->edges(v)
    //     for each edge j in g->edges(v)
    //       if i != j:
    //         A(i,j) = A(j,i) = 1
    for (itr = cav.begin();itr!=cav.end();itr++) {
      agi::GraphEdge* e;
      agi::EdgeIterator* eitr = g->edges(*itr);
      while ((e=g->iterate(eitr))) {
        agi::GraphEdge* eB;
        agi::EdgeIterator* eitrB = g->edges(*itr);
        while ((eB=g->iterate(eitrB))) {
          if( eB != e ) {
            A[edgeToIdx[e]*numCavEdges+edgeToIdx[eB]]=1;
            A[edgeToIdx[eB]*numCavEdges+edgeToIdx[e]]=1;
          }
        }
        g->destroy(eitrB);
      }
    }
    outCav << "#adj matrix\n";
    ss.str(""); //empty the stream
    for(int i = 0; i<numCavEdges; i++) {
      ss << numCavEdges;
      for(int j = 0; j<numCavEdges; j++) {
        ss << " " << A[i*numCavEdges+j];
      }
      ss << "\n";
    }
    outCav << ss.str();
    delete [] A;
    
    // for each edge in the cavity find the parts that 
    //  have a copy of it 
    // for edge i in edgesOfCavity
    //   Peers peers;
    //   getResidence(i,peers)
    //   write(i,peers)
    ss.str(""); //empty the stream
    for (sitr = edgesOfCavity.begin(); sitr != edgesOfCavity.end(); sitr++) {
      agi::Peers p;
      agi::GraphEdge* e = *sitr;
      g->getResidence(e,p);
      ss << edgeToIdx[*sitr];
      agi::Peers::iterator pitr;
      for (pitr = p.begin();pitr!=p.end();pitr++)
        ss << " " << *pitr;
      ss << "\n";
    }
    outCav << "#process ids\n";
    outCav << ss.str();

    outCav << "#destination process id\n";
    outCav << "1 " << dest << "\n";
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
      //For each peer of cavity
      Peers::iterator pitr;
      bool sent = false;
      for (pitr = peers.begin();pitr!=peers.end();pitr++) {
        part_t peer = *pitr;
        if (targets->has(peer) &&
            sending[*pitr]<targets->get(peer) &&
            cav.size()< cavSize) {
          writeCavity(g,cav,peer);
          //  addCavity to plan
          wgt_t w = addCavity(g,cav,peer,plan,target_dimension);
          planW+=w;
          sending[peer]+=w;
          sent=true;
          break;
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


  Selector* makeSelector(DiffusiveInput* in,Queue* q,
                         std::vector<int>* cd,
                         std::vector<double>* cw ) {
    return new Selector(in,q,cd,cw);
  }
}
