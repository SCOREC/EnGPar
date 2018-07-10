#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosGraph_graph_color.hpp>
#include <KokkosKernels_Handle.hpp>
#include <binGraph.h>
#include <mpi.h>
#include <stdio.h>
#include <cstring>
#include "engpar_support.h"
#include "engpar_input.h"
#include "engpar_kokkosColoring.h"


bool Check_Directed(agi::Ngraph* g, agi::etype t=0) {

  agi::PNgraph* pg = g->publicize();

  const agi::lid_t numverts = pg->num_local_verts; 

  const agi::lid_t* degree_list = pg->degree_list[t];

  const agi::lid_t* edge_list = pg->edge_list[t];

  for (int i=0; i<numverts; ++i) { 
    for (int j=degree_list[i]; j<degree_list[i+1]; ++j) {
      bool back_edge = false; 
      for (int k=degree_list[edge_list[j]]; k<degree_list[1+edge_list[j]]; ++k) {
        if (i == edge_list[k])
          back_edge = true;
      }
      if (!back_edge) {
        printf ("Directed graphs are not supported\n");
        return false;
      }
    }
  }
  return true;
}


int main(int argc, char* argv[]) {

  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  if ( argc != 2 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file>",argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  Kokkos::initialize(argc,argv);
  
  agi::Ngraph* g = agi::createBinGraph(argv[1]);
  //assert(Check_Directed(g));

  // Call EnGPar graph coloring
  double t0 = PCU_Time(); 
  engpar::ColoringInput* in = engpar::createColoringInput(g, VTX_TYPE);
  agi::lid_t** colors = new agi::lid_t*[1];
  engpar::EnGPar_KokkosColoring(in, colors); 
  printf ("Coloring time: %f\n", PCU_Time()-t0);

  // Assign colors to graph vertices
  agi::GraphTag* tag = g->createIntTag(-1);
  agi::GraphVertex* v;
  agi::VertexIterator* vitr = g->begin();
  int i=0;
  while ((v=g->iterate(vitr))) {
    g->setIntTag(tag, v, (*colors)[i++]);
  }

  // Check that the coloring is valid
  agi::GraphEdge* e;
  agi::EdgeIterator* eitr = g->begin(0);
  int conflicts = 0;
  while ((e=g->iterate(eitr))) {
    int u = g->getIntTag(tag, g->u(e));
    int v = g->getIntTag(tag, g->v(e));
    if (u==v) ++conflicts;
  }
  g->destroy(eitr);
  assert(conflicts == 0);

  destroyGraph(g);
  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n"); 
  Kokkos::finalize();
  EnGPar_Finalize();
  MPI_Finalize();

  return 0;

}
