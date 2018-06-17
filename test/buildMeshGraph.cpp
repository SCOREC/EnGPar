#include <apfGraph.h>
#include <apfMesh2.h>
#include <cassert>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apf.h>
#include <cstdlib>
#include <stdint.h>
#include <apfNumbering.h>
#include <cstring>
#include <engpar_support.h>

void testSizes(apf::Mesh* m,agi::Ngraph* g,int primary,int* seconds,int n);
void testIds(apf::Mesh* m,agi::Ngraph* g,int primary,int* seconds,int n);
void testVertices(apf::Mesh* m,agi::Ngraph* g);
void testEdges(apf::Mesh* m,agi::Ngraph* g,int primary, int* seconds,int n);
void testGhosts(apf::Mesh* m,agi::Ngraph* g,int primary, int* seconds,int n);

void testGraph(apf::Mesh* m,agi::Ngraph* g,int primary, int* seconds,int n) {
  testSizes(m,g,primary,seconds,n);
  testIds(m,g,primary,seconds,n);
  testVertices(m,g);
  testEdges(m,g,primary,seconds,n);
  testGhosts(m,g,primary,seconds,n);
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  EnGPar_Initialize();
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <primary_dimension> <secondary_dimension>\n", argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  //Check dimension too high
  //This doesn't get caught since primary is assumed to be 3 for now
  try {
    agi::Ngraph* g = agi::createAPFGraph(m,5,0);
    agi::destroyGraph(g);
    throw "FAIL\n";
  }
  catch(int e) {
    assert(e==1);
    printf("Error caught successfully\n");
  } 
   
  //check dimension too low
  try {
    agi::Ngraph* g = agi::createAPFGraph(m,3,-1);
    agi::destroyGraph(g);
    throw "FAIL\n";
  }
  catch(int e) {
    assert(e==1);
    if (!PCU_Comm_Self())
      printf("Error caught successfully\n");
  }
  //check primary==secondary
  try {
    agi::Ngraph* g = agi::createAPFGraph(m,3,3);
    agi::destroyGraph(g);
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
  agi::Ngraph* g = agi::createAPFGraph(m,primary,second);
  int secondaries1[1] = {second};
  testGraph(m,g,primary,secondaries1,1);

  agi::checkValidity(g);
  agi::destroyGraph(g);

  PCU_Barrier();
  //Test of multiple edge types, vertices=0,faces=1
  try {
    int secondaries[2] = {0,2};
    if (!PCU_Comm_Self())
      printf("\nConstructing Graph with Vertex Dim: %d, and Edge Dims: %d,%d\n",
             primary,secondaries[0],secondaries[1]);

    agi::Ngraph* g2 = agi::createAPFGraph(m,primary,secondaries,2);
    MPI_Barrier(MPI_COMM_WORLD);
    testGraph(m,g2,primary,secondaries,2);

    agi::checkValidity(g2);
    agi::destroyGraph(g2);
  }
  catch(int e) {
    if (primary==0||primary==2) {
      assert(e==2);
      if (!PCU_Comm_Self())
        printf("Error caught successfully\n");
    }
  }
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  EnGPar_Finalize();
  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}

agi::lid_t getNumPins(apf::Mesh* m,int primary,int second) {
  apf::MeshEntity* ent;
  apf::MeshIterator* itr = m->begin(second);
  size_t num_pins=0;
  while ((ent = m->iterate(itr))) {
    apf::Adjacent adj;
    m->getAdjacent(ent,primary,adj);
    num_pins+=adj.size();
  }
  m->end(itr);
  return num_pins;
}

agi::lid_t getNumNaiveEdges(apf::Mesh* m,int primary,int second) {
  apf::MeshEntity* ent;
  apf::MeshIterator* itr = m->begin(primary);
  agi::lid_t num_edges=0;
  while ((ent = m->iterate(itr))) {
    apf::Adjacent adj1;
    m->getAdjacent(ent,second,adj1);
    for (size_t i=0;i<adj1.size();i++) {
      apf::Adjacent adj;
      m->getAdjacent(adj1[i],primary,adj);
      num_edges+=adj.size()-1;
    }
  }
  m->end(itr);

  return num_edges;
}

void testSizes(apf::Mesh* m,agi::Ngraph* g,int primary, int* seconds,int n) {
  if (!PCU_Comm_Self())
    printf("Checking Sizes\n");
  agi::lid_t num_verts = countOwned(m,primary);
  agi::lid_t global_verts = countOwned(m,primary);
  global_verts = PCU_Add_Long(global_verts);
  assert(g->numLocalVtxs()==num_verts);
  assert(g->numGlobalVtxs()==global_verts);

  for (int i=0;i<n;i++) {
    agi::lid_t num_edges = m->count(seconds[i]);
    agi::lid_t global_edges = countOwned(m,seconds[i]);
    global_edges = PCU_Add_Long(global_edges);
    assert(g->numLocalEdges(i)==num_edges);
    assert(g->numGlobalEdges(i)==global_edges);
    if (PCU_Comm_Peers()==1)
      assert(getNumPins(m,primary,seconds[i])==g->numLocalPins(i));
    else
      assert(getNumPins(m,primary,seconds[i])<=g->numLocalPins(i));
  }
}

void testIds(apf::Mesh* m,agi::Ngraph* g,int primary,int* seconds,int n) {
  //Test the ids of the graph vertices
  if (!PCU_Comm_Self())
    printf("Testing Ids\n");
  agi::VertexIterator* gitr = g->begin();
  agi::GraphVertex* vtx = NULL;
  apf::MeshIterator* mitr = m->begin(primary);
  apf::MeshEntity* ent = NULL;
  apf::GlobalNumbering* id_nums=NULL;
  id_nums = m->findGlobalNumbering("primary_ids_global");
  assert(id_nums);
  while ((vtx = g->iterate(gitr)) && (ent = m->iterate(mitr))) {
    assert(g->globalID(vtx)==(gid_t)apf::getNumber(id_nums,ent,0,0));
  }
  m->end(mitr);
  for (int i=0;i<n;i++) {
    agi::EdgeIterator* eitr = g->begin(i);
    char name[100];
    sprintf(name,"secondary_ids%d_global",i);
    id_nums = m->findGlobalNumbering(name);
    assert(id_nums);
    mitr = m->begin(seconds[i]);
    agi::GraphEdge* edge;
    while ((edge = g->iterate(eitr)) && (ent = m->iterate(mitr))) {
      assert(g->globalID(edge)==(gid_t)apf::getNumber(id_nums,ent,0,0));
    }
    m->end(mitr);
    g->destroy(eitr);
  }
}

void testVertices(apf::Mesh*,agi::Ngraph* g) {
  if (!PCU_Comm_Self())
    printf("Iterating over vertices\n");
  //Test iterating through vertices
  agi::VertexIterator* gitr = g->begin();
  agi::GraphVertex* vtx=NULL;
  agi::lid_t i=0;
  while ((vtx = g->iterate(gitr))) {
    i++;
    assert(g->weight(vtx)==1.0);
    g->coord(vtx);
    assert(g->degree(vtx,0)>0);
    assert(i<=g->numLocalVtxs());
  }
  assert(i==g->numLocalVtxs());

}

void testEdges(apf::Mesh* m,agi::Ngraph* g,int primary,int* seconds,int n) {
  if (!PCU_Comm_Self())
    printf("Iterating over edges & pins\n");
  //Test iterating through edges & pins on vertices
  agi::VertexIterator* gitr = g->begin();
  agi::GraphVertex* vtx=NULL;
  for (int i=0;i<n;i++) {
    agi::lid_t num_pins = getNumNaiveEdges(m,primary,seconds[i]);
    int tot_pins=0;
    int ghost_pins=0;
    int other_count=0;
    gitr=g->begin();
    while ((vtx = g->iterate(gitr))) {
      int deg = g->degree(vtx,i);
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g->edges(vtx,i);
      for (int j=0;j<deg;j++) {
        edge = g->iterate(eitr);
        assert(g->weight(edge)>0);
        int np = g->degree(edge);
        agi::GraphVertex* out;
        agi::PinIterator* pitr = g->pins(edge);
        other_count+=np-1;
        for (int k=0;k<np;k++) {
          out = g->iterate(pitr);
          if (g->owner(out)!=PCU_Comm_Self()) {
            ghost_pins++;
          }
          if (g->isEqual(out,vtx))
              continue;
          tot_pins++;
        }
        g->destroy(pitr);
      }
      g->destroy(eitr);
    }
    assert(tot_pins-ghost_pins==num_pins);
  }

}

void testGhosts(apf::Mesh* m,agi::Ngraph* g, int primary, int* seconds, int n) {
  if (!PCU_Comm_Self())
    printf("Iterating to check ghost counts\n");
  //For each edge type
  for (int i=0;i<n;i++) {
    agi::GraphEdge* edge;
    agi::EdgeIterator* eitr = g->begin(i);

    apf::MeshIterator* mitr = m->begin(seconds[i]);
    apf::MeshEntity* ent;
    //Iterate over graph hyperedges
    while ((edge = g->iterate(eitr))) {
      //Iterate the mesh entity
      ent = m->iterate(mitr);
      //Iterate over the graph pins
      agi::GraphVertex* pin;
      agi::PinIterator* pitr = g->pins(edge);
      agi::lid_t deg = g->degree(edge);
      agi::lid_t owned_pins=0,ghost_pins=0;
      for (agi::lid_t i=0;i<deg;i++) {
        pin = g->iterate(pitr);
        //Get the owner of the graph vertex
        agi::part_t owner = g->owner(pin);
        if (PCU_Comm_Self()==owner)
          owned_pins++;
        else
          ghost_pins++;
      }
      g->destroy(pitr);
      //Get the number of adjacent mesh primaries
      apf::Adjacent adj;
      m->getAdjacent(ent,primary,adj);
      assert(owned_pins==(agi::lid_t)adj.getSize());
      assert(owned_pins+ghost_pins == deg);
    }
    m->end(mitr);
    g->destroy(eitr);
  }
}
