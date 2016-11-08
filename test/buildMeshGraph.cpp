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
    printf("Error caught successfully\n");
  }
  //check primary==secondary
  try {
    agi::apfGraph g(m,3,3);
    throw "FAIL!\n";
  }
  catch(int e) {
    assert(e==2);
    printf("Error caught successfully\n");
  }

  
  int primary = atoi(argv[3]);
  int second = atoi(argv[4]);
  printf("Constructing Graph with Vertices: %d, and edges: %d\n",primary,second);
  agi::apfGraph g(m,primary,second);
  int secondaries1[1] = {second};
  testSizes(m,g,primary,secondaries1,1);
  testVertices(m,g);
  testEdges(m,g,primary,secondaries1,1);
  

  

  //Test of multiple edge types, vertices=0,faces=1
  
  int secondaries[2] = {0,2};
  printf("Constructing Graph with Vertices: %d, and edges: %d,%d\n",primary,secondaries[0],secondaries[1]);
  
  agi::apfGraph g2(m,primary,secondaries,2);
  
  testSizes(m,g2,primary,secondaries,2);
  testVertices(m,g2);
  testEdges(m,g2,primary,secondaries,2);
  
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();

  printf("All tests passed\n");
  return 0;
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
  printf("Checking Sizes\n");
  size_t num_verts = countOwned(m,primary);
  assert(g.numVtxs()==num_verts);
  for (int i=0;i<n;i++)
    assert(getNumNaiveEdges(m,primary,seconds[i])>=g.numEdges(i));
}

void testVertices(apf::Mesh* m,agi::apfGraph& g) {
  printf("Iterating over vertices\n");
  //Test iterating through vertices
  agi::VertexIterator* gitr = g.begin();
  agi::GraphVertex* vtx=NULL;
  size_t i=0;
  while (vtx = g.iterate(gitr)) {
    i++;
    assert(g.weight(vtx)==1.0);
    assert(g.degree(vtx,0)>0);
    assert(i<=g.numVtxs());
  }
  assert(i==g.numVtxs());

}

void testEdges(apf::Mesh* m,agi::apfGraph& g,int primary,int* seconds,int n) {
  printf("Iterating over edges & pins\n");
  //Test iterating through edges & pins on vertices
  agi::VertexIterator* gitr = g.begin();
  agi::GraphVertex* vtx=NULL;
  for (int i=0;i<n;i++) {
    int num_edges = getNumNaiveEdges(m,primary,seconds[i]);
    double tot_edges=0;
    int other_count=0;
    gitr=g.begin();
    while (vtx = g.iterate(gitr)) {
      int deg = g.degree(vtx,i);
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g.edges(vtx,i);
      for (int j=0;j<deg;j++) {
        edge = g.iterate(eitr);
        int np = g.degree(edge);
        agi::GraphVertex* out;
        agi::PinIterator* pitr = g.pins(edge);
        other_count+=np-1;
        for (int k=0;k<np;k++) {
          out = g.iterate(pitr);
          if (g.isEqual(out,vtx))
              continue;
          tot_edges++;
        }
      }
      g.destroy(eitr);
    }
    //NOTE: with current edge weight implementation the sum of naive edges = sum of edge weights
    assert(tot_edges==num_edges);
  }

}
