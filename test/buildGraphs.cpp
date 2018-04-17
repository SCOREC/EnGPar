#include "buildGraphs.h"
#include <set>
#include <PCU.h>
#include <ngraph.h>

//Builds a traditional ring graph with an undirected edge between each
//sequentially numbered pair of vertices.
//Each part gets 4 sequentially numbered vertices.
//Traditional graphs only have directed edges.  Thus, to construct the ring we
//create one forward and one backward edge between each pair of sequentially
//numbered pair of vertices.  The source vertex of the pins forming an edge
//must be locally owned (i.e., exists in the vertex list on that part).  In
//addition, edges in a traditional graph have global ids, but, those ids are not
//required to be unique or consistent for a given vertex pair for correct EnGPar
//operation.
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
      pins.push_back(e);
    }
  }
  if (PCU_Comm_Peers()>1) {
    agi::gid_t e = (local_verts*PCU_Comm_Self()+global_verts-1)%global_verts;
    edges.push_back(e*2);
    degrees.push_back(2);
    pins.push_back((e+1)%global_verts);
    pins.push_back(e);
    owners[e] = (PCU_Comm_Self()+PCU_Comm_Peers()-1) % PCU_Comm_Peers();
  }
  std::vector<agi::wgt_t> weights;
  graph->constructGraph(false,verts,weights,edges,degrees,pins,owners);
  return graph;
}

agi::Ngraph* buildHyperGraph() {
  agi::Ngraph* graph = agi::createEmptyGraph();
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
  graph->constructVerts(false,verts.size(),&verts[0],&weights[0]);

  for (agi::gid_t i=0;i<local_verts;i++) {
    agi::gid_t e = local_verts*PCU_Comm_Self()+i;
    edges.push_back(e*2);
    degrees.push_back(2);
    pins.push_back(e);
    pins.push_back((e+1)%global_verts);
    if (i==local_verts-1&&PCU_Comm_Peers()>1) {
      if ((PCU_Comm_Self()+1)%PCU_Comm_Peers()!=PCU_Comm_Self()) {
        owners[(e+1)%global_verts] = (PCU_Comm_Self()+1)%PCU_Comm_Peers();
      }
    }
    else {
      edges.push_back(e*2+1);
      degrees.push_back(2);
      pins.push_back((e+1)%global_verts);
      if ((PCU_Comm_Self()+1)%PCU_Comm_Peers()!=PCU_Comm_Self()) {
        owners[e] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();
      }
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
  graph->constructEdges(edges.size(),&edges[0],
                                      &degrees[0],&pins[0]);

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
  graph->constructEdges(edges2.size(),&edges2[0],
                                        &degrees2[0],&pins2[0]);
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
  graph->constructEdges(edges,degrees,pins,weights);
  graph->constructGhosts(ghost_owners);

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
    printf("Ghost vertex: %d\n",local_verts*PCU_Comm_Self()-1);
  }
  for (agi::lid_t i=0;i<local_verts;i++) {
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
  return g;
}

agi::Ngraph* buildRequirementsGraph() {
  agi::Ngraph* g = agi::createEmptyGraph();
  std::unordered_map<agi::gid_t,agi::part_t> owners;
  std::vector<agi::gid_t> edges;
  std::vector<agi::lid_t> degrees;
  std::vector<agi::gid_t> pins;
  std::vector<agi::gid_t> verts;

  verts.push_back(1);
  verts.push_back(3);
  verts.push_back(5);
  verts.push_back(9);

  edges.push_back(3);
  edges.push_back(9);
  edges.push_back(27);

  degrees.push_back(0);
  degrees.push_back(4);
  degrees.push_back(2);

  pins.push_back(1);
  pins.push_back(3);
  pins.push_back(5);
  pins.push_back(9);
  pins.push_back(1);
  pins.push_back(1);
  std::vector<agi::wgt_t> weights;
  g->constructGraph(true,verts,weights,edges,degrees,pins,owners);
  return g;

}

agi::Ngraph* buildDisconnected2Graph() {
  //The graph is 2 parts where one has two vertices that have edges to 4 different components in the other part
  //The core component is a line of 5 vertices, then there are two components with depth 2. One is a tree and the other is a line. The last component is a line of 1
  //The first part also has an island of vertices completely disconnected from the rest of the graph
  assert(PCU_Comm_Peers()==2);

  if (!PCU_Comm_Self()) {
    printf("Distance Queue of part 0 for this graph should be:\n");
    printf("400\n");
    printf("200\n");
    printf("300\n");
    printf("100\n");
  }
  agi::Ngraph* g = agi::createEmptyGraph();
  std::unordered_map<agi::gid_t,agi::part_t> owners;
  std::vector<agi::gid_t> edges;
  std::vector<agi::lid_t> degrees;
  std::vector<agi::gid_t> pins;
  std::vector<agi::gid_t> verts;
  agi::coord_t* cs;
  if (PCU_Comm_Self()) {
    cs = new agi::coord_t[5];
    verts.push_back(-1);
    for (int i=0;i<3;i++)
      cs[0][i] = 0;
    cs[0][1]=1;
    verts.push_back(0);
    for (int i=0;i<3;i++)
      cs[1][i] = 0;
    edges.push_back(0);
    degrees.push_back(2);
    pins.push_back(-1);
    pins.push_back(0);
    edges.push_back(100);
    edges.push_back(200);
    edges.push_back(300);
    edges.push_back(400);
    degrees.push_back(2);
    degrees.push_back(2);
    degrees.push_back(2);
    degrees.push_back(2);
    pins.push_back(0);
    pins.push_back(100);
    pins.push_back(0);
    pins.push_back(200);
    pins.push_back(0);
    pins.push_back(300);
    pins.push_back(0);
    pins.push_back(400);
    owners[100]=0;
    owners[200]=0;
    owners[300]=0;
    owners[400]=0;

    //Island of vertices
    verts.push_back(1);
    cs[2][0]=1;cs[2][1]=2;cs[2][2]=0;
    verts.push_back(2);
    cs[3][0]=1.5;cs[3][1]=2.5;cs[3][2]=0;
    verts.push_back(3);
    cs[4][0]=1.5;cs[4][1]=1.5;cs[4][2]=0;
    edges.push_back(2);
    degrees.push_back(3);
    pins.push_back(1);
    pins.push_back(2);
    pins.push_back(3);
  }
  else {
    cs = new agi::coord_t[17];
    //Core Component, 5 vertices in a straight line
    for (int i=100;i<105;i++) {
      verts.push_back(i);
      cs[i-100][0]=-.5*(i-99);cs[i-100][1]=0; cs[i-100][2]=0;
      edges.push_back(i);
      degrees.push_back(2);
      pins.push_back(i);
      if (i==100)
        pins.push_back(0);
      else
        pins.push_back(i-1);
    }
    owners[0]=1;
    
    //Second Line component, 3 vertices in a straight line
    for (int i=200;i<203;i++) {
      verts.push_back(i);
      int index = i-200+5;
      cs[index][0]=.5*(i-199);cs[index][1]=0; cs[index][2]=0;
      edges.push_back(i);
      degrees.push_back(2);
      pins.push_back(i);
      if (i==200)
        pins.push_back(0);
      else
        pins.push_back(i-1);
    }
    
    /*Third component, tree
     *  root is connected to part 1
     *  depth=3
     *  each non leaf vertex has 2 children
     */
    
    verts.push_back(300);
    int index=8;
    cs[index][0] = 0; cs[index][1]=-.5; cs[index++][2]=0;
    edges.push_back(300);
    degrees.push_back(2);
    pins.push_back(300);
    pins.push_back(0);
    for (int i=0;i<3;i++) {
      int base = 300+i*2;
      verts.push_back(base+1);
      verts.push_back(base+2);
      edges.push_back(300+i+1);
      degrees.push_back(3);
      pins.push_back(300+i);
      pins.push_back(base+1);
      pins.push_back(base+2);
    }
    cs[index][0] = -.5; cs[index][1]=-1; cs[index++][2]=0;
    cs[index][0] = .5; cs[index][1]=-1; cs[index++][2]=0;
    cs[index][0] = -.75; cs[index][1]=-2; cs[index++][2]=0;
    cs[index][0] = -.25; cs[index][1]=-2; cs[index++][2]=0;
    cs[index][0] = .25; cs[index][1]=-2; cs[index++][2]=0;
    cs[index][0] = .75; cs[index][1]=-2; cs[index++][2]=0;

    
    //Last component, line of 1 vertex
    for (int i=400;i<401;i++) {
      verts.push_back(i);
      int ind = i-400+15;
      cs[ind][1]=-.5*(i-399);cs[ind][1]=.5*(i-399); cs[ind][2]=0;
      edges.push_back(i);
      degrees.push_back(2);
      pins.push_back(i);
      if (i==400)
        pins.push_back(0);
      else
        pins.push_back(i-1);
    }
    
  }
  std::vector<agi::wgt_t> weights;
  g->constructGraph(true,verts,weights,edges,degrees,pins,owners);
  g->setCoords(cs);
  return g;
 
}

agi::Ngraph* buildEmptyGraph() {
  agi::Ngraph* g = agi::createEmptyGraph();
  std::unordered_map<agi::gid_t,agi::part_t> owners;
  std::vector<agi::gid_t> edges;
  std::vector<agi::lid_t> degrees;
  std::vector<agi::gid_t> pins;
  std::vector<agi::gid_t> verts;

  if (PCU_Comm_Self()==1) {
    verts.push_back(0);
    if (PCU_Comm_Peers()>2) {
      edges.push_back(1);
      pins.push_back(0);
      pins.push_back(1);
      degrees.push_back(2);
      owners[1]=2;
    }
  }
  if (PCU_Comm_Self()==2) {
    for (int i=0;i<5;i++)
      verts.push_back(i);
    edges.push_back(1);
    pins.push_back(0);
    pins.push_back(1);
    degrees.push_back(2);
    owners[0]=1;
  }
  
  std::vector<agi::wgt_t> weights;
  g->constructGraph(false,verts,weights,edges,degrees,pins,owners);
  return g;  
}
