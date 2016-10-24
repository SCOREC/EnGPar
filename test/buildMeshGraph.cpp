#include <apfGraph.h>
#include <apfMesh2.h>
#include <cassert>
#include <PCU.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apf.h>
#include <cstdlib>
#include <iostream>
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
  try {
    agi::apfGraph g(m,5,0);
  }
  catch(int e) {
    assert(e==1);
    printf("Error caught successfully\n");
  }  
  //check dimension too low
  try {
    agi::apfGraph g(m,3,-1);
  }
  catch(int e) {
    assert(e==1);
    printf("Error caught successfully\n");
  }
  //check primary==secondary
  try {
    agi::apfGraph g(m,3,3);
  }
  catch(int e) {
    assert(e==2);
    printf("Error caught successfully\n");
  }
  
  int primary = atoi(argv[3]);
  int second = atoi(argv[4]);
  agi::apfGraph g(m,primary,second);

  int num_verts = countOwned(m,primary);
  assert(g.numVtxs()==num_verts);

  apf::MeshEntity* ent;
  apf::MeshIterator* itr = m->begin(primary);

  int num_edges=0;
  while ((ent = m->iterate(itr))) {
    apf::Downward down;
    int nd = m->getDownward(ent,second,down);
    for (int i=0;i<nd;i++) {
      apf::Adjacent adj;
      m->getAdjacent(down[i],primary,adj);
      num_edges+=adj.size()-1;
    }
  }
  assert(num_edges==g.numEdges());
  
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();

  printf("All tests passed\n");
  return 0;
}
