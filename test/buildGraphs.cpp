#include "buildGraphs.h"
#include <set>
#include <PCU.h>
#include <ngraph.h>

//Builds a ring of vertices
//  where each part gets 4 continuous vertices
agi::Ngraph* buildGraph() {
  agi::Ngraph* graph  = agi::createEmptyGraph();
  agi::lid_t local_verts = 4;
  agi::gid_t global_verts = 4*PCU_Comm_Peers();
  std::vector<agi::gid_t> verts;
  std::unordered_map<agi::gid_t,agi::part_t> owners;
  std::vector<agi::gid_t> edges;
  std::vector<agi::lid_t> degrees;
  std::vector<agi::gid_t> pins;
  for (agi::gid_t i=0;i<local_verts;i++)
    verts.push_back(local_verts*PCU_Comm_Self()+i);
  for (agi::gid_t i=0;i<local_verts;i++) {
    agi::gid_t e = local_verts*PCU_Comm_Self()+i;
    edges.push_back(e*2);
    degrees.push_back(2);
    pins.push_back(e);
    pins.push_back((e+1)%global_verts);
    if (i==local_verts-1&&PCU_Comm_Peers()>1) {
      owners[(e+1)%global_verts] = (PCU_Comm_Self()+1)%PCU_Comm_Peers();
      
    }
    else {
      edges.push_back(e*2+1);
      degrees.push_back(2);
      pins.push_back((e+1)%global_verts);
      owners[e] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();
      pins.push_back(e);
    }
  }
  if (PCU_Comm_Peers()>1) {
    agi::gid_t e = (local_verts*PCU_Comm_Self()+global_verts-1)%global_verts;
    edges.push_back(e*2);
    degrees.push_back(2);
    pins.push_back((e+1)%global_verts);
    pins.push_back(e);
  }
  //  owners[first_vert] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();
  //owners[last_vert] = (PCU_Comm_Self()+1)%PCU_Comm_Peers();
  std::vector<agi::wgt_t> weights;
  graph->constructGraph(false,verts,weights,edges,degrees,pins,owners);
  graph->setEdgeWeights(weights,0);
  return graph;
}

agi::Ngraph* buildHyperGraph() {
  agi::Ngraph* graph  = agi::createEmptyGraph();
  agi::lid_t local_verts = 4;
  agi::gid_t start_vert = local_verts*PCU_Comm_Self();
  agi::gid_t end_vert = local_verts*(PCU_Comm_Self()+1)-1;
  agi::lid_t local_edges = 3;
  agi::gid_t global_verts = local_verts*PCU_Comm_Peers();
  agi::gid_t global_edges = local_edges*PCU_Comm_Peers();
  std::vector<agi::gid_t> verts;
  std::unordered_map<agi::gid_t,agi::part_t> ghost_owners;
  std::vector<agi::gid_t> edges;
  std::vector<agi::lid_t> degrees;
  std::vector<agi::gid_t> pins;
  for (agi::gid_t i=0;i<local_verts;i++) 
    verts.push_back(local_verts*PCU_Comm_Self()+i);
   
  for (agi::gid_t i=0;i<local_edges;i++)
    edges.push_back(local_edges*PCU_Comm_Self()+i);
  degrees.push_back(4);
  degrees.push_back(2);
  degrees.push_back(4);
  for (agi::gid_t i=0;i<local_verts;i++)
    pins.push_back(verts[i]);
  pins.push_back(verts[1]);
  pins.push_back(verts[3]);
  for (agi::gid_t i=0;i<local_verts;i++) {
    gid_t v = (verts[i]+2)%global_verts;
    pins.push_back(v);
    
    if (v<start_vert||v>end_vert) {
      ghost_owners[v] = (PCU_Comm_Self()+1)%PCU_Comm_Peers();
    }
  }
  //add the edge backward from the first vertex of this part
  //  to the last of the previous
  if (PCU_Comm_Peers()>1) {
    edges.push_back((edges[0]+global_edges-1)%global_edges);
    degrees.push_back(4);
    for (agi::gid_t i=0;i<local_verts;i++) {
      gid_t v = (verts[i]+global_verts-2)%global_verts;
      pins.push_back(v);
      if (v<start_vert||v>end_vert)
        ghost_owners[v] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();

    }

  }
  std::vector<agi::wgt_t> weights;
  graph->constructGraph(true,verts,weights,edges,degrees,pins,ghost_owners);
  graph->setEdgeWeights(weights,0);
  return graph;
}

agi::Ngraph* buildGraphParts() {
  agi::Ngraph* graph  = agi::createEmptyGraph();
  agi::lid_t local_verts = 4;
  agi::gid_t global_verts = 4*PCU_Comm_Peers();
  std::vector<agi::gid_t> verts;
  std::unordered_map<agi::gid_t,agi::part_t> owners;
  std::vector<agi::gid_t> edges;
  std::vector<agi::lid_t> degrees;
  std::vector<agi::gid_t> pins;
  for (agi::gid_t i=0;i<local_verts;i++)
    verts.push_back(local_verts*PCU_Comm_Self()+i);
  std::vector<agi::wgt_t> weights;
  graph->constructVerts(false,verts,weights);

  for (agi::gid_t i=0;i<local_verts;i++) {
    agi::gid_t e = local_verts*PCU_Comm_Self()+i;
    edges.push_back(e*2);
    degrees.push_back(2);
    pins.push_back(e);
    pins.push_back((e+1)%global_verts);
    if (i==local_verts-1&&PCU_Comm_Peers()>1) {
      owners[(e+1)%global_verts] = (PCU_Comm_Self()+1)%PCU_Comm_Peers();
      
    }
    else {
      edges.push_back(e*2+1);
      degrees.push_back(2);
      pins.push_back((e+1)%global_verts);
      owners[e] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();
      pins.push_back(e);
    }
  }
  if (PCU_Comm_Peers()>1) {
    agi::gid_t e = (local_verts*PCU_Comm_Self()+global_verts-1)%global_verts;
    edges.push_back(e*2);
    degrees.push_back(2);
    pins.push_back((e+1)%global_verts);
    pins.push_back(e);
  }
  agi::etype t =graph->constructEdges(edges,degrees,pins);
  graph->setEdgeWeights(weights,t);

  std::vector<agi::gid_t> edges2;
  std::vector<agi::lid_t> degrees2;
  std::vector<agi::gid_t> pins2;
  agi::lid_t svert = PCU_Comm_Self()*local_verts;
  for (int i=0;i<4;i++) {
    degrees2.push_back(2);
    edges2.push_back(i);
    pins2.push_back(svert+i);
    pins2.push_back(svert+(i+2)%4);
  }
  if (PCU_Comm_Self()) {
    degrees2.push_back(2);
    edges2.push_back(4);
    pins2.push_back(svert+1);
    pins2.push_back(1);
    
  }
  agi::etype t2 = graph->constructEdges(edges2,degrees2,pins2);
  graph->setEdgeWeights(weights,t2);
  graph->constructGhosts(owners);
  return graph;
}

agi::Ngraph* buildHyperGraphParts() {
  agi::Ngraph* graph  = agi::createEmptyGraph();
  agi::lid_t local_verts = 4;
  agi::gid_t start_vert = local_verts*PCU_Comm_Self();
  agi::gid_t end_vert = local_verts*(PCU_Comm_Self()+1)-1;
  agi::lid_t local_edges = 3;
  agi::gid_t global_verts = local_verts*PCU_Comm_Peers();
  agi::gid_t global_edges = local_edges*PCU_Comm_Peers();
  std::vector<agi::gid_t> verts;
  std::unordered_map<agi::gid_t,agi::part_t> ghost_owners;
  std::vector<agi::gid_t> edges;
  std::vector<agi::lid_t> degrees;
  std::vector<agi::gid_t> pins;
  for (agi::gid_t i=0;i<local_verts;i++) 
    verts.push_back(local_verts*PCU_Comm_Self()+i);
   
  for (agi::gid_t i=0;i<local_edges;i++)
    edges.push_back(local_edges*PCU_Comm_Self()+i);
  degrees.push_back(4);
  degrees.push_back(2);
  degrees.push_back(4);
  for (agi::gid_t i=0;i<local_verts;i++)
    pins.push_back(verts[i]);
  pins.push_back(verts[1]);
  pins.push_back(verts[3]);
  for (agi::gid_t i=0;i<local_verts;i++) {
    gid_t v = (verts[i]+2)%global_verts;
    pins.push_back(v);
    
    if (v<start_vert||v>end_vert) {
      ghost_owners[v] = (PCU_Comm_Self()+1)%PCU_Comm_Peers();
    }
  }
  //add the edge backward from the first vertex of this part
  //  to the last of the previous
  if (PCU_Comm_Peers()>1) {
    edges.push_back((edges[0]+global_edges-1)%global_edges);
    degrees.push_back(4);
    for (agi::gid_t i=0;i<local_verts;i++) {
      gid_t v = (verts[i]+global_verts-2)%global_verts;
      pins.push_back(v);
      if (v<start_vert||v>end_vert)
        ghost_owners[v] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();

    }

  }
  std::vector<agi::wgt_t> weights;
  graph->constructVerts(true,verts,weights);
  agi::etype t = graph->constructEdges(edges,degrees,pins);
  graph->constructGhosts(ghost_owners);
  graph->setEdgeWeights(weights,t);
  return graph;
}


agi::Ngraph* buildHyperGraphLine() {
  agi::Ngraph* g  = agi::createEmptyGraph();
  agi::lid_t local_verts = 4;
  std::vector<agi::gid_t> verts;
  std::unordered_map<agi::gid_t,agi::part_t> owners;
  std::vector<agi::gid_t> edges;
  std::vector<agi::lid_t> degrees;
  std::vector<agi::gid_t> pins;
  if (PCU_Comm_Self()) {
    edges.push_back(local_verts*PCU_Comm_Self()-1);
    degrees.push_back(2);
    pins.push_back(local_verts*PCU_Comm_Self()-1);
    pins.push_back(local_verts*PCU_Comm_Self());
    owners[local_verts*PCU_Comm_Self()-1]=PCU_Comm_Self()-1;
    printf("Ghost vertex: %lu\n",local_verts*PCU_Comm_Self()-1);
  }
  for (agi::gid_t i=0;i<local_verts;i++) {
    int vert = local_verts*PCU_Comm_Self()+i;
    verts.push_back(vert);
    if (i==local_verts-1&&PCU_Comm_Self()==PCU_Comm_Peers()-1)
      continue;
    edges.push_back(vert);
    degrees.push_back(2);
    pins.push_back(vert);
    pins.push_back(vert+1);
    if (i==local_verts-1) {
      owners[vert+1]=PCU_Comm_Self()+1;
      printf("Ghost vertex: %d\n",vert+1);
    }
  }
  std::vector<agi::wgt_t> weights;
  g->constructGraph(true,verts,weights,edges,degrees,pins,owners);
  g->setEdgeWeights(weights,0);
  return g;
}
