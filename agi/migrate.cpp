#include "ngraph.h"
#include <PCU.h>
#include <set>
#include <vector>
namespace agi {

  typedef std::vector<GraphVertex*> VertexVector;
  //TODO: Make a vector by using a "tag" on the edges to detect added or not
  typedef std::set<GraphEdge*> EdgeVector;

  void getAffected(Ngraph* g, Migration* plan, VertexVector& verts,
		   EdgeVector* edges) {
    verts.reserve(plan->size());
    Migration::iterator itr;
    for (itr = plan->begin();itr!=plan->end();itr++) {
      GraphVertex* v = itr->first;
      part_t toSend = itr->second;
      if (toSend!=PCU_Comm_Self())
	verts.push_back(v);
    }
    for (etype t = 0;t<g->numEdgeTypes();t++) {
      //edges[t].reserve(verts.size());
      for (lid_t i = 0;i<verts.size();i++) {
	agi::EdgeIterator* itr = g->edges(verts[i]);
	agi::GraphEdge* e;
	while ((e = g->iterate(itr))) {
	  edges[t].insert(e);
	}
      }
    } 
  }

  void getSenders(Ngraph* g, VertexVector& verts, EdgeVector* edges,
		  VertexVector& vSenders,EdgeVector* eSenders) {
    /*    vSenders.reserve(verts.size());
    for (int i=0;i<verts.size();i++) {
      if (
    }
    */
  }
  void Ngraph::sendVertex(GraphVertex* vtx, part_t toSend) {
    gid_t gid = globalID(vtx);
    PCU_COMM_PACK(toSend,gid);
  }

  void Ngraph::recvVertex(std::vector<gid_t>& recv) {
    gid_t gid;
    PCU_COMM_UNPACK(gid);
    recv.push_back(gid);
  }
  
  void Ngraph::migrate(Migration* plan) {
    VertexVector affectedVerts;
    EdgeVector* affectedEdges = new EdgeVector[num_types];
    getAffected(this,plan,affectedVerts,affectedEdges);
    //Note: vertSenders = affectedVerts
    //VertexVector vertSenders;
    //Will work on if needed
    //EdgeVector* edgeSenders = new EdgeVector[num_types];
    //getSenders(this,affectedVerts,affectedEdges,vertSenders,edgeSenders);
    std::vector<gid_t> my_vertices;
    my_vertices.reserve(num_local_verts);
    
    Migration::iterator itr;
    PCU_Comm_Begin();
    //Send vertices
    for (itr = plan->begin();itr!=plan->end();itr++) {
      //printf("%d sending %p to %d\n",PCU_Comm_Self(),itr->first,itr->second);
      sendVertex(itr->first,itr->second);
      owners[localID(itr->first)] = itr->second;
      map_t::iterator mitr = vtx_mapping.find(globalID(itr->first));
      vtx_mapping.erase(mitr);
    }
    map_t::iterator mitr = vtx_mapping.begin();
    for (;mitr!=vtx_mapping.end();mitr++)
      my_vertices.push_back(mitr->first);
    PCU_Comm_Send();
    //Recieve vertices
    while (PCU_Comm_Receive()) {
      //printf("%d getting vertex\n",PCU_Comm_Self());
      recvVertex(my_vertices);
    }

    //Send Edges
    EdgeVector::iterator eitr = affectedEdges[0].begin();
    for (;eitr !=affectedEdges[0].end();eitr++) {
      //Get remotes of edge
      //send to each remote copy
      //sendEdge(*eitr,
    }
    //Recieve Edges
    

    //delete [] edgeSenders;
    delete [] affectedEdges;

    /*
    std::vector<gid_t> recv;
    std::vector<gid_t> recv_edges;

    //Reconstruct data
    lid_t num_added=recv.size();
    lid_t num_removed=sent.size();
    lid_t edge_added=recv_edges.size();
    lid_t edge_removed=0;
  
    lid_t new_verts = num_local_verts-num_removed+num_added;
    lid_t new_edges = num_local_edges-edge_removed+edge_added;
    lid_t* new_deg_list = new lid_t[new_local_verts];
    */
  }


}
