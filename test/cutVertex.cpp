#include <binGraph.h>
#include <mpi.h>
#include <PCU.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <ZoltanCutVertex.h>

typedef agi::lid_t lid_t;
void partitionInfo(agi::binGraph* g) {
  lid_t total_verts = g->numLocalVtxs();
  lid_t global_verts = g->numGlobalVtxs();
  lid_t edges = g->numLocalEdges()+g->numLocalEdges(SPLIT_TYPE);
  lid_t min = PCU_Min_Int(total_verts);
  lid_t max = PCU_Max_Int(total_verts);
  long tot = PCU_Add_Long(total_verts);
  double avg = ((double)tot)/PCU_Comm_Peers();
  double imb = max/avg;
  double inc = ((double)(tot-global_verts))/global_verts*100;
  if (!PCU_Comm_Self()) 
    printf("Vertices: Min %d Max %d Tot %lu Inc %1.4f Avg %1.4f Imb %1.3f\n",min,max,tot,inc,avg,imb);
  min = PCU_Min_Int(edges);
  max = PCU_Max_Int(edges);
  tot = PCU_Add_Long(edges);
  avg = ((double)tot)/PCU_Comm_Peers();
  imb = max/avg;
  if (!PCU_Comm_Self()) 
    printf("Edges: Min %d Max %d Tot %lu Avg %1.4f Imb %1.3f\n",min,max,tot,avg,imb);
  
  lid_t edge_cut =0;
  agi::EdgeIterator* eitr = g->begin(0);
  for (int i=0;i<g->numLocalEdges();i++) {
    agi::GraphEdge* e = g->iterate(eitr);
    if (g->owner(g->v(e))!=PCU_Comm_Self())
      edge_cut++;
  }
  edge_cut+=g->numLocalEdges(SPLIT_TYPE);
  edge_cut = PCU_Add_Long(edge_cut);
  if (!PCU_Comm_Self())
    printf("Edge Cut: %d\n",edge_cut);
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Debug_Open();
  if ( argc != 2&&argc!=3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <binary_graph_file> [vertex_partition_file]",argv[0]);
    PCU_Comm_Free();
    MPI_Finalize();
    assert(false);
  }
  agi::binGraph* g;
  zagi::ZoltanCutVertex* ptn;
  int self = PCU_Comm_Self();
  int num_parts = PCU_Comm_Peers();
  PCU_Switch_Comm(MPI_COMM_SELF);
  if (!self) {
    g = new agi::binGraph(argv[1]);
    ptn = new zagi::ZoltanCutVertex(g,num_parts);
    ptn->run();
  }
  else
    g = new agi::binGraph();
  agi::EdgePartitionMap map;
  if (!self) {
    ptn->createPtn(map);
    delete ptn;
  }

  PCU_Switch_Comm(MPI_COMM_WORLD);
  g->migrate(map);

  partitionInfo(g);

  delete g;

  if (!PCU_Comm_Self())
    printf("Block Partitioning\n");

  g = new agi::binGraph(argv[1]);

  partitionInfo(g);

  delete g;

  if (argc==3) {
  if (!PCU_Comm_Self())
    printf("Partitioned: %s\n",argv[2]);
    g = new agi::binGraph(argv[1],argv[2]);
    partitionInfo(g);
    delete g;
  }

  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");

  PCU_Comm_Free();
  MPI_Finalize();

}
