#include <apfGraph.h>
#include <apfMesh2.h>
#include <cassert>
#include <PCU.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apf.h>
#include <cstdlib>
#include <iostream>
#include <stdint.h>

void testSizes(apf::Mesh* m,agi::apfGraph& g,int primary,int* seconds,int n);
void testVertices(apf::Mesh* m,agi::apfGraph& g);
void testEdges(apf::Mesh* m,agi::apfGraph& g,int primary, int* seconds,int n);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <primary_dimension> <secondary_dimension>\n", argv[0]);
    MPI_Finalize();
    assert(false);
  }

  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  //Check dimension too high
  //This doesn't get caught since primary is assumed to be 3 for now
  /*
  try {
    agi::apfGraph g(m,5,0);
    throw "FAIL\n";
  }
  catch(int e) {
    assert(e==1);
    printf("Error caught successfully\n");
  } 
  */ 
  //check dimension too low
  try {
    agi::apfGraph g(m,3,-1);
    throw "FAIL\n";
  }
  catch(int e) {
    assert(e==1);
    if (!PCU_Comm_Self())
      printf("Error caught successfully\n");
  }
  //check primary==secondary
  try {
    agi::apfGraph g(m,3,3);
    throw "FAIL!\n";
  }
  catch(int e) {
    assert(e==2);
    if (!PCU_Comm_Self())
      printf("Error caught successfully\n");
  }

  
  int primary = atoi(argv[3]);
  int second = atoi(argv[4]);
  if (!PCU_Comm_Self())
    printf("\nConstructing Graph with Vertex Dim: %d, and Edge Dim: %d\n",
           primary,second);
  agi::apfGraph g(m,primary,second);
  int secondaries1[1] = {second};
  testSizes(m,g,primary,secondaries1,1);
  testVertices(m,g);
  testEdges(m,g,primary,secondaries1,1);
  
  PCU_Barrier();
  
  //Test of multiple edge types, vertices=0,faces=1
  int secondaries[2] = {0,2};
  if (!PCU_Comm_Self())
    printf("\nConstructing Graph with Vertex Dim: %d, and Edge Dims: %d,%d\n",
           primary,secondaries[0],secondaries[1]);

  agi::apfGraph g2(m,primary,secondaries,2);
  MPI_Barrier(MPI_COMM_WORLD);
  testSizes(m,g2,primary,secondaries,2);
  testVertices(m,g2);
  testEdges(m,g2,primary,secondaries,2);

  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  
  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}

size_t getNumPins(apf::Mesh* m,int primary,int second) {
  apf::MeshEntity* ent;
  apf::MeshIterator* itr = m->begin(second);
  size_t num_pins=0;
  while ((ent = m->iterate(itr))) {
    apf::Adjacent adj;
    m->getAdjacent(ent,primary,adj);
    num_pins+=adj.size();
  }
  return num_pins;
}

size_t getNumNaiveEdges(apf::Mesh* m,int primary,int second) {
  apf::MeshEntity* ent;
  apf::MeshIterator* itr = m->begin(primary);
  size_t num_edges=0;
  while ((ent = m->iterate(itr))) {
    apf::Downward down;
    int nd = m->getDownward(ent,second,down);
    for (int i=0;i<nd;i++) {
      apf::Adjacent adj;
      m->getAdjacent(down[i],primary,adj);
      num_edges+=adj.size()-1;
    }
  }
  return num_edges;
}

void testSizes(apf::Mesh* m,agi::apfGraph& g,int primary, int* seconds,int n) {
  if (!PCU_Comm_Self())
    printf("Checking Sizes\n");
  size_t num_verts = countOwned(m,primary);
  size_t global_verts = countOwned(m,primary);
  PCU_Add_Long(global_verts);
  assert(g.numLocalVtxs()==num_verts);
  assert(g.numGlobalVtxs()==global_verts);

  for (int i=0;i<n;i++) {
    size_t num_edges = m->count(seconds[i]);
    size_t global_edges = countOwned(m,seconds[i]);
    assert(g.numLocalEdges(i)==num_edges);
    assert(g.numGlobalEdges(i)==global_edges);
    if (PCU_Comm_Peers()==1)
      assert(getNumPins(m,primary,seconds[i])==g.numLocalPins(i));
    else
      assert(getNumPins(m,primary,seconds[i])<=g.numLocalPins(i));
  }
}

void testVertices(apf::Mesh* m,agi::apfGraph& g) {
  if (!PCU_Comm_Self())
    printf("Iterating over vertices\n");
  //Test iterating through vertices
  agi::VertexIterator* gitr = g.begin();
  agi::GraphVertex* vtx=NULL;
  size_t i=0;
  while (vtx = g.iterate(gitr)) {
    i++;
    assert(g.weight(vtx)==1.0);
    assert(g.degree(vtx,0)>0);
    assert(i<=g.numLocalVtxs());
  }
  assert(i==g.numLocalVtxs());

}
void testEdges(apf::Mesh* m,agi::apfGraph& g,int primary,int* seconds,int n) {
  if (!PCU_Comm_Self())
    printf("Iterating over edges & pins\n");
  //Test iterating through edges & pins on vertices
  agi::VertexIterator* gitr = g.begin();
  agi::GraphVertex* vtx=NULL;
  for (int i=0;i<n;i++) {
    int num_pins = getNumNaiveEdges(m,primary,seconds[i]);
    int tot_pins=0;
    int ghost_pins=0;
    int other_count=0;
    gitr=g.begin();
    while (vtx = g.iterate(gitr)) {
      int deg = g.degree(vtx,i);
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g.edges(vtx,i);
      for (int j=0;j<deg;j++) {
        edge = g.iterate(eitr);
        assert(g.weight(edge)>0);
        int np = g.degree(edge);
        agi::GraphVertex* out;
        agi::PinIterator* pitr = g.pins(edge);
        other_count+=np-1;
        for (int k=0;k<np;k++) {
          out = g.iterate(pitr);
          if (g.owner(out)!=PCU_Comm_Self()) {
            ghost_pins++;
	  }
          if (g.isEqual(out,vtx))
              continue;
          tot_pins++;
        }
      }
      g.destroy(eitr);
    }
    assert(tot_pins-ghost_pins==num_pins);
  }

}
