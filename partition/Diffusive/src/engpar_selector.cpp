#include "engpar_selector.h"
#include "../engpar.h"
#include <set>
#include <PCU.h>
namespace engpar {

  void getCavity(agi::Ngraph* g, agi::GraphEdge* edge, agi::Migration* plan,
                 Cavity& cav, Peers& peers) {
    agi::PinIterator* pitr = g->pins(edge);
    agi::lid_t deg = g->degree(edge);
    agi::GraphVertex* vtx;
    for (agi::lid_t i =0;i<deg;i++) {
      vtx = g->iterate(pitr);
      if (g->owner(vtx)==PCU_Comm_Self()) {
        if(!plan->has(vtx)) {
          cav.push_back(vtx);
        }
      }
      else
        peers.insert(g->owner(vtx));
    }
    g->destroy(pitr);
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

  wgt_t Selector::select(Targets* targets, agi::Migration* plan,
                         wgt_t planW, unsigned int cavSize,int target_dimension) {
    Queue::iterator itr;
    for (itr = q->begin();itr!=q->end();itr++) {
      //Create Cavity and peers
      Cavity cav;
      Peers peers;
      getCavity(g,*itr,plan,cav,peers);
      //For each peer of cavity
      Peers::iterator itr;
      for (itr = peers.begin();itr!=peers.end();itr++) {
        part_t peer = *itr;
        if (targets->has(peer) &&
            sending[*itr]<targets->get(peer) &&
            cav.size()< cavSize) {
          //  addCavity to plan
          wgt_t w = addCavity(g,cav,peer,plan,target_dimension);
          planW+=w;
          sending[peer]+=w;
        }
      }
    }
    return planW;
  
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
  typedef std::set<Migr,CompareMigr> MigrComm;

  //return  map<int neighbor, pair<double vtxW, double edgeW> > where
  //  neighbor is a neighbor's part id
  //  vtxW is the vtx weight capacity of neighbor
  //  edgeW is the edge weight capacity of neighbor
  Midd* Selector::trim(Targets*, agi::Migration* plan) {
    //compute the weight of the vertices and edges being sent to each peer
    PeerEdgeSet* peerEdges = new PeerEdgeSet[completed_dimensions->size()];
    agi::Migration::iterator itr;
    for(itr = plan->begin();itr!=plan->end();itr++) {
      agi::GraphVertex* vtx = *itr;
      const int dest = plan->get(vtx);
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        insertInteriorEdges(vtx, dest, peerEdges[i][dest],
                            completed_dimensions->at(i));
      }
    }
    //send vtx and edge weight
    PCU_Comm_Begin();
    PeerEdgeSet::iterator sitr;
    for (sitr = peerEdges[0].begin();sitr!=peerEdges[0].end();sitr++) {
      const int dest = sitr->first;
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        double w = weight(peerEdges[i][dest]);
        PCU_COMM_PACK(dest, w);
      }
    }
    delete [] peerEdges;
    PCU_Comm_Send();
    MigrComm incoming;
    double w;
    while (PCU_Comm_Receive()) {
      Migr migr(PCU_Comm_Sender());
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        PCU_COMM_UNPACK(w);
        migr.addWeight(w);
      }
      incoming.insert(migr);        
    }
    Midd accept;
    Ws avail = new double[completed_dimensions->size()];
    bool isAvail=true;
    for (unsigned int i=0;i<completed_dimensions->size();i++) {
      double totW = getWeight(g,completed_dimensions->at(i));
      avail[i] = completed_weights->at(i) - totW;
      if (avail[i]<0)
        isAvail=false;
    }

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
    PCU_Barrier();
    delete [] avail;
    PCU_Comm_Begin();
    Midd::iterator acc; 
    for (acc=accept.begin();acc!=accept.end();acc++) {
      for (unsigned int i=0;i<completed_dimensions->size();i++) 
        PCU_COMM_PACK(acc->first, acc->second[i]);
      delete [] acc->second;
    }
    PCU_Comm_Send();
    Midd* capacity = new Midd;
    while (PCU_Comm_Receive()) {
      int nbor = PCU_Comm_Sender();
      capacity->insert(std::make_pair(nbor,
                                      new double[completed_dimensions->size()]));
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        PCU_COMM_UNPACK(w);
        (*capacity)[nbor][i] = w;
      }
    }
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
    agi::Migration::iterator pitr;
    for(pitr = plan->begin();pitr!=plan->end();pitr++) {
      agi::GraphVertex* v = *pitr;
      int dest = plan->get(v);
      bool isSpace=true;
      EdgeSet* tmpEdges = new EdgeSet[completed_dimensions->size()];
      
      for (unsigned int i=0;i<completed_dimensions->size();i++) {
        tempInsertInteriorEdges(v, dest, tmpEdges[i],
                                completed_dimensions->at(i),
                                peerEdges[i][dest]);
        if(weight(peerEdges[i][dest])+weight(tmpEdges[i])>(*capacity)[dest][i])
          isSpace=false;
      }
      if (isSpace) {
        keep.push_back(PlanPair(v,dest));
        for (unsigned int i=0;i<completed_dimensions->size();i++) 
          combineSets(peerEdges[i][dest],tmpEdges[i]);
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


  Selector* makeSelector(DiffusiveInput* in,Queue* q,
                         std::vector<int>* cd,
                         std::vector<double>* cw ) {
    return new Selector(in,q,cd,cw);
  }
}
