#include <ngraph.h>
#include <pngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include "buildGraphs.h"
#include "sellCSigma.h"

void writeGraphArrays(agi::Ngraph* g, const char* name) {
  agi::PNgraph* pg = g->publicize();
  printf("\n---- start %s ----\n", name);
  printf("isSellCSigma  %d\n", pg->isSellCSigma);
  printf("isHyperGraph  %d\n", pg->isHyperGraph);
  printf("num verts %d\n", pg->num_local_verts);
  printf("num types %d\n", pg->num_types);
  if( pg->isSellCSigma ) {
    printf("num vtx-to-edge chunks  %d\n", pg->num_vtx_chunks);
    printf("chunk size  %d\n", pg->chunk_size);
  }
  for (agi::etype t=0;t<pg->num_types;t++) {
    printf("type %d\n", t);
    printf("  num edges %d\n", pg->num_local_edges[t]);
    printf("  num pins  %d\n", pg->num_local_pins[t]);
    if( pg->isSellCSigma ) {
      // SellCSigma
      printf("    degree_list[t][0] %d\n", pg->degree_list[t][0]);
      for(agi::lid_t chunk = 0; chunk < pg->num_vtx_chunks; chunk++) {
        printf("    degree_list[t][%d] %d\n", chunk+1, pg->degree_list[t][chunk+1]);
      }
      agi::lid_t chunkStart = 0;
      // loop over chunks
      for(agi::lid_t chunk = 0; chunk < pg->num_vtx_chunks; chunk++) {
        agi::lid_t maxChunkDeg = pg->degree_list[t][chunk+1] - pg->degree_list[t][chunk];
        printf("    chunkStart %d maxChunkDegree %d\n", chunkStart, maxChunkDeg);
        agi::lid_t chunkSize = pg->chunk_size*maxChunkDeg;
        // write the chunk
        printf("    chunk values\n");
        for (agi::lid_t i = chunkStart; i < chunkStart+chunkSize; i++)
          printf("      %d\n", pg->edge_list[t][i]);
        // loop over chunk 'rows' 
        for (agi::lid_t i = 0; i < pg->chunk_size; i++) {
          agi::lid_t vtx = chunkStart+i;
          printf("    vtx %d\n", i);
          for (agi::lid_t j = vtx;
                          j < vtx+(maxChunkDeg*pg->chunk_size);
                          j += pg->chunk_size) {
            printf("      edge_list[t][%d] %d\n", j, pg->edge_list[t][j]);
          }
        }
        chunkStart+=chunkSize; // the number of entries in the chunk
      }
    } else {
      // CSR
      printf("    degree_list[t][0] %d\n", pg->degree_list[t][0]);
      for(agi::lid_t vtx = 0; vtx < pg->num_local_verts; vtx++) {
        printf("    degree_list[t][%d] %d\n", vtx+1, pg->degree_list[t][vtx+1]);
      }
      for(agi::lid_t vtx = 0; vtx < pg->num_local_verts; vtx++) {
        printf("    vtx %d\n", vtx);
        for(agi::lid_t j = pg->degree_list[t][vtx]; j < pg->degree_list[t][vtx+1]; j++) {
          printf("      edge_list[t][%d] %d\n", j, pg->edge_list[t][j]);
        }
      }
    }
  }

  printf("---- end %s ----\n", name);
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();

  agi::Ngraph* g;
  if (argc==1)
    if (PCU_Comm_Peers()==2)
      g = buildDisconnected2Graph();
    else
      g=buildHyperGraphLine();
  else {
    g = agi::createEmptyGraph();
    g->loadFromFile(argv[1]);
  }
  PCU_Barrier();

  agi::lid_t C = 4;
  agi::lid_t sigma = 1;
  agi::Ngraph* scg = ssg::convertFromAGI(g,C,sigma);

  writeGraphArrays(g,"csr");
  writeGraphArrays(scg,"scg");

  agi::destroyGraph(g);

  agi::destroyGraph(scg);

  PCU_Barrier();
  if (!PCU_Comm_Self()) 
    printf("All Tests Passed\n"); 

  EnGPar_Finalize();
  MPI_Finalize();

  return 0;
}
