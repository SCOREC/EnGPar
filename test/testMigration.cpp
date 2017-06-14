#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include <set>
agi::Ngraph* buildGraph();
agi::Ngraph* buildHyperGraph();
int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  agi::Ngraph* g = buildGraph();
  PCU_Barrier();
  agi::Migration* plan = new agi::Migration;
  agi::GraphVertex* v;
  agi::VertexIterator* itr=g->begin();
  while ((v = g->iterate(itr))) {
    agi::GraphIterator* gitr = g->adjacent(v);
    agi::GraphVertex* other;
    while ((other=g->iterate(gitr))) {
      int owner = g->owner(other);
      if (owner>PCU_Comm_Self()) {
	plan->insert(std::make_pair(v,owner));
	break;
      }
    }
  }
  g->setOriginalOwners();
  g->migrate(plan);

  agi::destroyGraph(g);
  
  PCU_Barrier();

  g = buildHyperGraph();
  PCU_Barrier();
  plan = new agi::Migration;

  itr=g->begin();
  v = g->iterate(itr);
  (*plan)[v] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();
  g->setOriginalOwners();
  g->migrate(plan);
  
  agi::destroyGraph(g);
  PCU_Barrier();
  
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n");
  

  EnGPar_Finalize();
  MPI_Finalize();

}

//Builds a ring of vertices
//  where each part gets 4 continuous vertices
agi::Ngraph* buildGraph() {
    if (!PCU_Comm_Self())
    printf("Building Regular Graph\n");
  agi::Ngraph* graph  = new agi::Ngraph;
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
    if (PCU_Comm_Peers()>1&&i==local_verts-1) {
      owners.insert(std::make_pair((e+1)%global_verts,
				   (PCU_Comm_Self()+1)%PCU_Comm_Peers()));
    }
    else {
      edges.push_back(e*2+1);
      degrees.push_back(2);
      pins.push_back((e+1)%global_verts);
      pins.push_back(e);
    }
  }
  //add the edge backward from the first vertex of this part
  //  to the last of the previous
  if (PCU_Comm_Peers()>1) {
    agi::gid_t e = (local_verts*PCU_Comm_Self()+global_verts-1)%global_verts;
    edges.push_back(e);
    degrees.push_back(2);
    pins.push_back((e+1)%global_verts);
    pins.push_back(e);
    owners[e] = (PCU_Comm_Self()+PCU_Comm_Peers()-1)%PCU_Comm_Peers();
  }
  std::vector<agi::wgt_t> weights;
  graph->constructGraph(false,verts,weights,edges,degrees,pins,owners);
  graph->setEdgeWeights(weights,0);
  return graph;
}

agi::Ngraph* buildHyperGraph() {
  if (!PCU_Comm_Self())
    printf("Building HyperGraph\n");
  agi::Ngraph* graph  = new agi::Ngraph;
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
