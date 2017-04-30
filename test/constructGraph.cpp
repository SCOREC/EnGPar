#include <ngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include <set>
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
  graph->constructGraph(false,verts,edges,degrees,pins,owners);

  assert(graph->numLocalVtxs()==local_verts);
  if (PCU_Comm_Peers()>1) {
    assert(graph->numGhostVtxs()==2);
  }
  else
    assert(graph->numGhostVtxs()==0);
  assert(graph->numLocalEdges()==local_verts*2);
  assert(graph->numEdgeTypes()==1);
  assert(!graph->isHyper());

  std::set<agi::gid_t> vs;
  for (int i = -1 ; i<(int)local_verts+1;i++)  {
    agi::gid_t vert = (PCU_Comm_Self()*local_verts+global_verts+i)%global_verts;
    vs.insert(vert);
  }
  agi::GraphVertex* v;
  agi::VertexIterator* vitr = graph->begin();
  while ((v = graph->iterate(vitr))) {
    assert(vs.find(graph->globalID(v))!=vs.end());
    assert(graph->localID(v)<graph->numTotalVtxs());
    
    if (graph->localID(v)>graph->numLocalVtxs())
      assert(graph->owner(v)!=PCU_Comm_Self());
    agi::GraphVertex* other;
    agi::GraphIterator* gitr = graph->adjacent(v);
    while ((other = graph->iterate(gitr))) {
      assert(vs.find(graph->globalID(other))!=vs.end());
    }
  }
  agi::destroyGraph(graph);

}

void buildHyperGraph() {
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
  graph->constructGraph(true,verts,edges,degrees,pins,ghost_owners);

  assert(graph->numLocalVtxs()==local_verts);
  if (PCU_Comm_Peers()>1) {
    assert(graph->numGhostVtxs()==4);
  }
  assert(graph->numLocalEdges()==local_edges+(PCU_Comm_Peers()>1));
  assert(graph->numEdgeTypes()==1);
  assert(graph->isHyper());
  assert(graph->numLocalPins()==pins.size());

  agi::lid_t vert_degs[4] = {2,3,2,3};
  agi::lid_t edge_degs[4] = {4,2,4,4};
  agi::lid_t ghost_degs[4] = {0,0,2,2};
  agi::gid_t pin[14] = {0,1,2,3,1,3,2,3,4,5,6,7,0,1};
  agi::VertexIterator* itr = graph->begin();
  agi::GraphVertex* vtx;
  lid_t i=0;
  while ((vtx = graph->iterate(itr))) {
    assert(graph->localID(vtx)<graph->numLocalVtxs());
    assert(graph->globalID(vtx)<graph->numGlobalVtxs());
    assert(graph->degree(vtx)==vert_degs[i]);
    lid_t count = 0;
    agi::EdgeIterator* eitr = graph->edges(vtx);
    agi::GraphEdge* e;
    while ((e = graph->iterate(eitr))) {
      count++;
    }
    assert(count==vert_degs[i]);
    graph->destroy(eitr);
    i++;
  }
  assert(i==graph->numLocalVtxs());
  i=0;

  agi::GraphEdge* e;
  agi::EdgeIterator* eitr = graph->begin(0);
  int k=0;
  while ((e = graph->iterate(eitr))) {
    agi::lid_t deg = graph->degree(e);
    assert(deg==edge_degs[i]);
    lid_t count =0;
    lid_t ghost =0;
    agi::PinIterator* pitr = graph->pins(e);
    for (lid_t j=0;j<deg;j++) {
      vtx = graph->iterate(pitr);
      if (PCU_Comm_Peers()>1||pin[k]<graph->numLocalVtxs())
      assert(graph->localID(vtx)==pin[k++]);
      assert(graph->localID(vtx)<graph->numTotalVtxs());
      assert(graph->globalID(vtx)<graph->numGlobalVtxs());
      count++;
      if (graph->localID(vtx)>=graph->numLocalVtxs()) {
	assert(PCU_Comm_Peers()>1);
	assert(graph->owner(vtx)!=PCU_Comm_Self());
	ghost++;
      }
    }
    assert(ghost==ghost_degs[i]*(PCU_Comm_Peers()>1));
    assert(count==edge_degs[i]);
    i++;
  }
  graph->destroy(eitr);
  agi::destroyGraph(graph);
}
