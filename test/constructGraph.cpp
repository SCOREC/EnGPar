#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
typedef agi::lid_t lid_t;

void buildGraph();
void buildHyperGraph();
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  buildGraph();

  PCU_Barrier();
  
  buildHyperGraph();


  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n");
  

  EnGPar_Finalize();
  MPI_Finalize();
}

//Builds a ring of vertices
//  where each part gets 4 continuous vertices
void buildGraph() {
  if (!PCU_Comm_Self())
    printf("Building Regular Graph\n");
  agi::Ngraph* graph  = new agi::Ngraph;
  agi::lid_t local_verts = 4;
  agi::gid_t global_verts = 4*PCU_Comm_Peers();
  std::vector<agi::gid_t> verts;
  std::vector<agi::gid_t> edges;
  std::vector<lid_t> degrees;
  std::vector<agi::gid_t> pins;
  for (agi::gid_t i=0;i<local_verts;i++)
    verts.push_back(local_verts*PCU_Comm_Self()+i);
  for (agi::gid_t i=0;i<local_verts;i++) {
    agi::gid_t e = local_verts*PCU_Comm_Self()+i;
    edges.push_back(e);
    degrees.push_back(2);
    pins.push_back(e);
    pins.push_back((e+1)%global_verts);
  }
  //add the edge backward from th first vertex of this part
  //  to the last of the previous
  if (PCU_Comm_Peers()>1) {
    agi::gid_t e = (local_verts*PCU_Comm_Self()-1+global_verts)%global_verts;
    edges.push_back(e);
    degrees.push_back(2);
    pins.push_back(e);
    pins.push_back((e+1)%global_verts);
  }
  graph->constructGraph(verts,edges,degrees,pins);

  assert(graph->numLocalVtxs()==local_verts);
  if (PCU_Comm_Peers()>1) {
    assert(graph->numGhostVtxs()==2);
  }
  assert(graph->numLocalEdges()==local_verts*2);
  assert(graph->numEdgeTypes()==1);
  assert(!graph->isHyper());
  agi::destroyGraph(graph);

}

void buildHyperGraph() {
  if (!PCU_Comm_Self())
    printf("Building HyperGraph\n");
  agi::Ngraph* graph  = new agi::Ngraph;
  agi::lid_t local_verts = 4;
  agi::lid_t local_edges = 3;
  agi::gid_t global_verts = local_verts*PCU_Comm_Peers();
  agi::gid_t global_edges = local_edges*PCU_Comm_Peers();
  std::vector<agi::gid_t> verts;
  std::vector<agi::gid_t> edges;
  std::vector<lid_t> degrees;
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
  for (agi::gid_t i=0;i<local_verts;i++)
    pins.push_back((verts[i]+2)%global_verts);
  //add the edge backward from the first vertex of this part
  //  to the last of the previous
  if (PCU_Comm_Peers()>1) {
    edges.push_back((edges[0]-1+global_edges)%global_edges);
    degrees.push_back(4);
    for (agi::gid_t i=0;i<local_verts;i++)
      pins.push_back((verts[i]-2+global_verts)%global_verts);

  }
  graph->constructGraph(verts,edges,degrees,pins);

  assert(graph->numLocalVtxs()==local_verts);
  if (PCU_Comm_Peers()>1) {
    assert(graph->numGhostVtxs()==4);
  }
  assert(graph->numLocalEdges()==local_edges+(PCU_Comm_Peers()>1));
  assert(graph->numEdgeTypes()==1);
  assert(graph->isHyper());
  assert(graph->numLocalPins()==pins.size()-4*(PCU_Comm_Peers()>1));
  
  agi::destroyGraph(graph);
}
