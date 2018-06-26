#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <cstring>
#include <engpar_support.h>
#include <fstream>
#include <binGraph.h>

void attachCoords(agi::Ngraph* g, char* filename);

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  EnGPar_Initialize();
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <ebin graph> <coords_file> <save prefix>\n", argv[0]);
    EnGPar_Finalize();
    MPI_Finalize();
    assert(false);
  }

  if (PCU_Comm_Peers()>1) {
    EnGPar_Error_Message("This tool is only built for serial graphs\n");
    assert(false);
  }
    
  agi::Ngraph* g = agi::createBinGraph(argv[1]);

  attachCoords(g,argv[2]);
  
  g->saveToFile(argv[3]);

  agi::destroyGraph(g);
  
  EnGPar_Finalize();
  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}


void attachCoords(agi::Ngraph* g, char* filename) {
  agi::coord_t* cs = new agi::coord_t[g->numLocalVtxs()];

  std::ifstream in_str(filename);
  if (!in_str) {
    EnGPar_Error_Message("%s could not be opened",filename);
  }

  double x,y,z;
  int index = 0;
  while (in_str>>x>>y>>z) {
    cs[index][0] = x;
    cs[index][1] = y;
    cs[index][2] = z;
    index++;
  }

  g->setCoords(cs);

  agi::VertexIterator* vitr = g->begin();
  agi::GraphVertex* v;

  while ((v = g->iterate(vitr))) {
    const agi::coord_t& c = g->coord(v);
    printf(" %f %f %f\n",c[0],c[1],c[2]);
  }
  delete [] cs;
}
