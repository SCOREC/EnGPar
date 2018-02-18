#include <ngraph.h>
#include <pngraph.h>
#include <PCU.h>
#include <engpar_support.h>
#include "buildGraphs.h"
#include "sellCSigma.h"

/* @brief write the vertex to hyperedge graph
 * @param t (in) the hyperedge type
 * @param pg (in) the publicized SCG
 **/
void writeVtxToEdgeScg(agi::lid_t t, agi::PNgraph* pg) {
  assert(pg->isSellCSigma);
  if(!pg->isSellCSigma) 
    return;

  agi::lid_t chunkStart = 0;
  // loop over chunks
  for(agi::lid_t chunk = 0; chunk < pg->num_vtx_chunks; chunk++) {
    agi::lid_t maxChunkDeg = pg->degree_list[t][chunk+1] - pg->degree_list[t][chunk];
    printf("    vtx->edge chunkStart %ld maxChunkDegree %ld\n", chunkStart, maxChunkDeg);
    agi::lid_t chunkSize = pg->chunk_size*maxChunkDeg;
    // write the chunk
    printf("    vtx->edge chunk values\n");
    for (agi::lid_t i = chunkStart; i < chunkStart+chunkSize; i++)
      printf("      %ld\n", pg->edge_list[t][i]);
    // loop over chunk 'rows' 
    for (agi::lid_t i = 0; i < pg->chunk_size; i++) {
      agi::lid_t vtx = chunkStart+i;
      if (chunk*pg->chunk_size+i >= pg->num_local_verts)
        printf("    vtx %ld gid:-1\n", i);
      else
        printf("    vtx %ld gid:%ld\n", i,pg->local_unmap[chunk*pg->chunk_size+i]);
      for (agi::lid_t j = vtx;
          j < vtx+(maxChunkDeg*pg->chunk_size);
          j += pg->chunk_size) {
        if (pg->edge_list[t][j]==-1)
          printf("      edge_list[t][%ld] gid:-1\n", j);
        else
          printf("      edge_list[t][%ld] gid:%ld\n", j,
              pg->edge_unmap[t][pg->edge_list[t][j]]);
      }
    }
    chunkStart+=chunkSize; // the number of entries in the chunk
  }
}

/* @brief write the hyperedge to vertex graph
 * @param t (in) the hyperedge type
 * @param pg (in) the publicized SCG
 **/
void writeEdgeToVtxScg(agi::lid_t t, agi::PNgraph* pg) {
  assert(pg->isSellCSigma && pg->isHyperGraph);
  if(!pg->isSellCSigma || !pg->isHyperGraph) {
    printf("skipping edge->vtx; not an scg hypergraph\n");
    return;
  }
  const agi::lid_t C = pg->chunk_size;
  // loop over chunks
  agi::lid_t chunkStart = 0; // the first index in the current chunk
  for(agi::lid_t chunk = 0; chunk < pg->num_edge_chunks[t]; chunk++) {
    agi::lid_t maxChunkDeg = pg->pin_degree_list[t][chunk+1] - pg->pin_degree_list[t][chunk];
    printf("    edge->vtx chunkStart %ld maxChunkDegree %ld\n", chunkStart, maxChunkDeg);
    agi::lid_t chunkSize = C*maxChunkDeg; // the number of entries in the chunk
    // write the chunk
    printf("    edge->vtx chunk values\n");
    for (agi::lid_t i = chunkStart; i < chunkStart+chunkSize; i++)
      printf("      %ld\n", pg->pin_list[t][i]);
    // loop over chunk 'rows'
    for (agi::lid_t i = 0; i < C; i++) {
      // the starting index of vertices adjacent to the ith edge in the chunk
      agi::lid_t edge = chunkStart+i;
      if (chunk*C+i >= pg->num_local_edges[t])
        printf("    edge %ld gid:-1\n", i);
      else
        printf("    edge %ld gid:%ld\n", i,pg->edge_unmap[t][chunk*C+i]);
      for (agi::lid_t j = edge;
          j < edge+(maxChunkDeg*C);
          j += C) {
        if (pg->edge_list[t][j]==-1)
          printf("      pin_list[t][%ld] gid:-1\n", j);
        else
          printf("      pin_list[t][%ld] gid:%ld\n", j,
              pg->local_unmap[pg->pin_list[t][j]]);
      }
    }
    chunkStart+=chunkSize; // set the starting chunk index for the next iteration
  }
}

void writeGraphArrays(agi::Ngraph* g, const char* name) {
  agi::PNgraph* pg = g->publicize();
  printf("\n---- start %s ----\n", name);
  printf("isSellCSigma  %d\n", pg->isSellCSigma);
  printf("isHyperGraph  %d\n", pg->isHyperGraph);
  printf("num verts %ld\n", pg->num_local_verts);
  printf("num types %d\n", pg->num_types);
  if( pg->isSellCSigma ) {
    printf("num vtx-to-edge chunks  %d\n", pg->num_vtx_chunks);
    printf("chunk size  %d\n", pg->chunk_size);
  }
  for (agi::etype t=0;t<pg->num_types;t++) {
    printf("type %d\n", t);
    printf("  num edges %ld\n", pg->num_local_edges[t]);
    printf("  num pins  %ld\n", pg->num_local_pins[t]);
    if( pg->isSellCSigma ) {
      // SellCSigma
      printf("    degree_list[t][0] %ld\n", pg->degree_list[t][0]);
      for(agi::lid_t chunk = 0; chunk < pg->num_vtx_chunks; chunk++) {
        printf("    degree_list[t][%ld] %ld\n", chunk+1, pg->degree_list[t][chunk+1]);
      }
      writeVtxToEdgeScg(t,pg);
      writeEdgeToVtxScg(t,pg);
    } else {
      // CSR
      printf("    degree_list[t][0] %ld\n", pg->degree_list[t][0]);
      for(agi::lid_t vtx = 0; vtx < pg->num_local_verts; vtx++) {
        printf("    degree_list[t][%ld] %ld\n", vtx+1, pg->degree_list[t][vtx+1]);
      }
      for(agi::lid_t vtx = 0; vtx < pg->num_local_verts; vtx++) {
        printf("    vtx lid:%ld gid:%ld\n", vtx,pg->local_unmap[vtx]);
        for(agi::lid_t j = pg->degree_list[t][vtx]; j < pg->degree_list[t][vtx+1]; j++) {
          printf("      edge_list[t][%ld] gid:%ld\n", j, pg->edge_unmap[t][pg->edge_list[t][j]]);
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
  agi::lid_t sigma = g->numLocalVtxs();
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
