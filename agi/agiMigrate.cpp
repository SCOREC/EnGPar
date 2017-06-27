#include "ngraph.h"
#include <PCU.h>
#include <unordered_set>
#include <vector>
#include <engpar_support.h>
namespace agi {
  typedef std::unordered_set<GraphVertex*> VertexVector;
  //TODO: Make a vector by using a "tag" on the edges to detect added or not
  typedef std::unordered_set<GraphEdge*> EdgeVector;
  void Ngraph::updateGhostOwners(Migration* plan) {
    PCU_Comm_Begin();
    Migration::iterator itr;
    for (itr = plan->begin();itr!=plan->end();itr++) {
      GraphVertex* v = itr->first;
      part_t toSend = itr->second;
      GraphIterator* gitr = adjacent(v);
      GraphVertex* other;
      while ((other=iterate(gitr))) {
        part_t own = owner(other);
        if (own!=PCU_Comm_Self()) {
          gid_t gid = globalID(v);
          PCU_COMM_PACK(own,gid);
          PCU_COMM_PACK(own,toSend);
        }
      }
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      gid_t gid;
      part_t own;
      PCU_COMM_UNPACK(gid);
      PCU_COMM_UNPACK(own);
      lid_t lid = vtx_mapping[gid];
      if (lid>=num_local_verts) {
        assert(lid>=num_local_verts);
        owners[lid-num_local_verts]=own;
      }
    }
  }
  /*
    Finds the vertices and edges that need to be communicated
    Also stores the local vertices and edges that will be owned
    on this part after migration
  */
  void getAffected(Ngraph* g, Migration* plan, VertexVector& verts,
                   EdgeVector* edges) {
    verts.reserve(plan->size());
    Migration::iterator itr;
    for (itr = plan->begin();itr!=plan->end();itr++) {
      GraphVertex* v = itr->first;
      part_t toSend = itr->second;
      if (toSend!=PCU_Comm_Self()) {
        verts.insert(v);
      }
    }

    //For each edge type
    for (etype t = 0;t<g->numEdgeTypes();t++) {
      //edges[t].reserve(verts.size());
      VertexVector::iterator itr;
      //loop through the vertices being sent
      for (itr = verts.begin();itr!=verts.end();itr++) {
        GraphIterator* gitr =g->adjacent(*itr);
        GraphEdge* e;
        GraphEdge* old=NULL;
        GraphVertex* v;
        while ((v = g->iterate(gitr))) {
          e=g->edge(gitr);
          if (old==NULL||e!=old) 
            edges[t].insert(e);
        }
      }      
    } 
  }
  
  void addVertices(Ngraph* g, std::vector<gid_t>& ownedVerts,
                   VertexVector& verts, std::vector<wgt_t>& weights,
                   std::vector<part_t>& originalOwners) {
    VertexIterator* vitr = g->begin();
    GraphVertex* v;
    while ((v = g->iterate(vitr))) {
      if (verts.find(v)==verts.end()) {
        ownedVerts.push_back(g->globalID(v));
        weights.push_back(g->weight(v));
        originalOwners.push_back(g->originalOwner(v));
      }
    }
  }
  void addEdges(Ngraph* g, Migration* plan, std::vector<gid_t>& ownedEdges,
                std::vector<wgt_t>& edgeWeights,
                std::vector<lid_t>& degrees,std::vector<gid_t>& pins,
                std::unordered_map<gid_t,part_t>& ghost_owners,
                std::unordered_set<gid_t>& addedEdges,
                etype t) {
    EdgeIterator* itr = g->begin(t);
    GraphEdge* e;
    while ((e = g->iterate(itr))) {
      agi::PinIterator* pitr = g->pins(e);
      agi::GraphVertex* end;
      lid_t i=0;
      while ((end=g->iterate(pitr))) {
        if (plan->find(end)!=plan->end())
          break;
        i++;
      }
      g->destroy(pitr);
      if (!end) {
        //Add the Edge
        gid_t gid = g->globalID(e);
        ownedEdges.push_back(gid);
        edgeWeights.push_back(g->weight(e));
        degrees.push_back(g->degree(e));
        addedEdges.insert(gid);
        pitr=g->pins(e);
        while ((end=g->iterate(pitr))) {
          gid_t other_gid = g->globalID(end);
          pins.push_back(other_gid);
          part_t owner = g->owner(end);
          if (plan->find(end)!=plan->end())
            owner = plan->find(end)->second;
          if (owner!=PCU_Comm_Self()) 
            ghost_owners[other_gid] = owner;
          
        }
      }
    }
    g->destroy(itr);
  }
  
  void Ngraph::sendVertex(GraphVertex* vtx, part_t toSend) {
    gid_t gid = globalID(vtx);
    wgt_t w = weight(vtx);
    part_t old_owner = originalOwner(vtx);
    PCU_COMM_PACK(toSend,gid);
    PCU_COMM_PACK(toSend,w);
    PCU_COMM_PACK(toSend,old_owner);
  }

  void Ngraph::recvVertex(std::vector<gid_t>& recv,
                          std::vector<wgt_t>& wgts,
                          std::vector<part_t>& old_owners) {
    gid_t gid;
    PCU_COMM_UNPACK(gid);
    wgt_t w;
    PCU_COMM_UNPACK(w);
    part_t old_owner;
    PCU_COMM_UNPACK(old_owner);
    recv.push_back(gid);
    wgts.push_back(w);
    old_owners.push_back(old_owner);
  }

  void Ngraph::sendEdges(Migration* plan, EdgeVector& affectedEdges) {
    EdgeVector::iterator eitr;
    for (eitr = affectedEdges.begin();eitr!=affectedEdges.end();
         eitr++) {
      GraphEdge* e = *eitr;
      gid_t* pin;
      part_t* pin_owners;
      lid_t deg=0;
      gid_t id;
      wgt_t eweight = weight(e);
      std::unordered_set<part_t> residence;
      if (isHyperGraph) {
        id = globalID(e);
        pin = new gid_t[degree(e)];
        pin_owners = new part_t[degree(e)];
        agi::PinIterator* pitr = this->pins(e);
        agi::GraphVertex* vtx;
        for (lid_t i=0;i<degree(e);i++) {
          vtx = iterate(pitr);
          part_t o = owner(vtx);

          if (plan->find(vtx)!=plan->end())
            o= plan->find(vtx)->second;
          pin_owners[deg]=o;
          pin[deg++] = globalID(vtx);
          residence.insert(o);
        }
        destroy(pitr);
      }
      else {
        id = localID(e);
        pin = new gid_t[2];
        pin_owners = new part_t[2];
        GraphVertex* source = u(e);
        pin_owners[deg] = owner(source);
        pin[deg++] = globalID(source);
        GraphVertex* dest = v(e);
        part_t o = owner(dest);
        if (plan->find(dest)!=plan->end())
          o= plan->find(dest)->second;
        pin_owners[deg] = o;
        pin[deg++] = globalID(dest);
        Migration::iterator mitr = plan->find(source);
        pin_owners[0] = mitr->second;
        assert(mitr!=plan->end());
        residence.insert(mitr->second);
      }
      std::unordered_set<part_t>::iterator sitr;
      for (sitr=residence.begin();sitr!=residence.end();sitr++) {
        PCU_COMM_PACK(*sitr,id);
        PCU_COMM_PACK(*sitr,eweight);
        PCU_COMM_PACK(*sitr,deg);
        PCU_Comm_Pack(*sitr,pin,deg*sizeof(gid_t));
        PCU_Comm_Pack(*sitr,pin_owners,deg*sizeof(part_t));

      }
      delete [] pin;
      delete [] pin_owners;
    }

  }

  void Ngraph::recvEdges(std::unordered_set<gid_t>& addedEdges,
                         std::vector<gid_t>& ownedEdges, std::vector<wgt_t>& edgeWeights,
                         std::vector<lid_t>& degrees, std::vector<gid_t>& pins,
                         std::unordered_map<gid_t,part_t>& ghost_owners,
                         etype t) {
    while (PCU_Comm_Receive()) {
      gid_t id;
      wgt_t eweight;
      lid_t deg;
      PCU_COMM_UNPACK(id);
      PCU_COMM_UNPACK(eweight);
      PCU_COMM_UNPACK(deg);
      gid_t* pin = new gid_t[deg];
      part_t* pin_owners = new part_t[deg];
      PCU_Comm_Unpack(pin,deg*sizeof(gid_t));
      PCU_Comm_Unpack(pin_owners,deg*sizeof(part_t));
      if (isHyperGraph) {
        if (addedEdges.find(id)!=addedEdges.end())
          continue;
      }
      addedEdges.insert(id);
      edge_mapping[t][id]=0;
      ownedEdges.push_back(id);
      edgeWeights.push_back(eweight);
      degrees.push_back(deg);
      for (lid_t i=0;i<deg;i++) {
        pins.push_back(pin[i]);
        if (pin_owners[i]!=PCU_Comm_Self())
          ghost_owners[pin[i]]=pin_owners[i];
      }
      delete [] pin;
      delete [] pin_owners;
    }

  }
  void Ngraph::migrate(Migration* plan) {
    double times[2];
    times[0]=PCU_Time();
    double mig[4];
    mig[0]=PCU_Time();
    int nt = num_types;
    updateGhostOwners(plan);
    VertexVector affectedVerts;
    EdgeVector* affectedEdges = new EdgeVector[nt];
    std::vector<gid_t> ownedVerts;
    std::vector<wgt_t> vertWeights;
    std::vector<part_t> old_owners;
    std::vector<gid_t>* ownedEdges = new std::vector<gid_t>[nt];
    std::vector<wgt_t>* edgeWeights = new std::vector<wgt_t>[nt];
    std::vector<lid_t>* degrees = new std::vector<lid_t>[nt];
    std::vector<gid_t>* pins = new std::vector<gid_t>[nt];
    std::unordered_map<gid_t,part_t> ghost_owners;
    ownedVerts.reserve(num_local_verts);
    vertWeights.reserve(num_local_verts);
    mig[0] = PCU_Time()-mig[0];
    mig[1]=PCU_Time();
    getAffected(this,plan,affectedVerts,affectedEdges);
    mig[1] = PCU_Time()-mig[1];
    mig[2] = PCU_Time();
    addVertices(this,ownedVerts,affectedVerts,vertWeights,old_owners);
    mig[2] = PCU_Time()-mig[2];
    mig[3] = PCU_Time();
    std::unordered_set<gid_t>* addedEdges = new std::unordered_set<gid_t>[nt];
    for (etype i=0;i<nt;i++)
      addEdges(this,plan,ownedEdges[i],edgeWeights[i],degrees[i],pins[i],
               ghost_owners,addedEdges[i],i);
    mig[3] = PCU_Time()-mig[3];
    //Send and recieve vertices
    Migration::iterator itr;
    PCU_Comm_Begin();
    for (itr = plan->begin();itr!=plan->end();itr++) {
      sendVertex(itr->first,itr->second);
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      recvVertex(ownedVerts,vertWeights,old_owners);
    }
    //Send and receive edges of each type
    for (etype t = 0;t < nt;t++) {
      PCU_Comm_Begin();
      sendEdges(plan,affectedEdges[t]);
      PCU_Comm_Send();
      recvEdges(addedEdges[t],ownedEdges[t],edgeWeights[t],degrees[t],
                pins[t],ghost_owners,t);
    }
    delete [] addedEdges;
    times[0]=PCU_Time()-times[0];
    times[1]=PCU_Time();
    //Construct the migrated graph
    constructVerts(isHyperGraph,ownedVerts,vertWeights);
    for (int i=0;i<nt;i++) {
      constructEdges(ownedEdges[i],degrees[i],pins[i]);
      setEdgeWeights(edgeWeights[i],i);
    }
    constructGhosts(ghost_owners);
    setOriginalOwners(old_owners);
    times[1] = PCU_Time()-times[1];

    PCU_Max_Doubles(times,2);
    PCU_Max_Doubles(mig,4);
    // if (!PCU_Comm_Self()) {
    //   printf("Create Vars %f, getAffected %f, addVerts %f, addEdges %f\n",
    //          mig[0],mig[1],mig[2],mig[3]);
    //   printf("Time to migrate %f, Time to construct %f\n",times[0],times[1]);
    // }
    
    delete [] affectedEdges;
    delete [] ownedEdges;
    delete [] edgeWeights;
    delete [] degrees;
    delete [] pins;
    delete plan;
  }

  PartitionMap* Ngraph::getPartition() {
    if (EnGPar_Is_Log_Open()) {
      char message[20];
      sprintf(message,"getPartition\n");
      EnGPar_Log_Function(message);
    }

    PCU_Comm_Begin();
    VertexIterator* itr = begin();
    GraphVertex* vtx;
    while ((vtx = iterate(itr))) {
      gid_t v = globalID(vtx);
      PCU_COMM_PACK(originalOwner(vtx),v);
    }
    PCU_Comm_Send();
    PartitionMap* map = new PartitionMap;
    while (PCU_Comm_Receive()) {
      gid_t v;
      PCU_COMM_UNPACK(v);
      map->insert(std::make_pair(v,PCU_Comm_Sender()));
    }
    if (EnGPar_Is_Log_Open()) {
      EnGPar_End_Function();
    }
    PCU_Barrier();
    return map;
  }


}
