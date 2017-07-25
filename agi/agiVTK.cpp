#include "ngraph.h"
#include <PCU.h>
#include "agi.h"
namespace agi {

  void writePoints(Ngraph* g,FILE* f) {
    fprintf(f,"<Points>\n");
    fprintf(f,"<DataArray NumberOfComponents=\"3\"");
    fprintf(f," type=\"Float64\" format=\"ascii\" Name=\"Points\">\n");
    //Write the position of each vertex
    agi::GraphVertex* v;
    agi::VertexIterator* vitr = g->begin();
    int i=0;
    while ((v=g->iterate(vitr))) {
      const coord_t& c = g->coord(v);
      fprintf(f," %f %f %f\n",c[0],c[1],c[2]);
      i++;
    }
    printf("%d\n",i);
    //Write the position of each edge
    agi::GraphEdge* e;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((e=g->iterate(eitr))) {
      coord_t c = {0,0,0};
      agi::PinIterator* pitr = g->pins(e);
      while ((v = g->iterate(pitr))) {
        const coord_t& cv = g->coord(v);
        c[0]+=cv[0];c[1]+=cv[1];c[2]+=cv[2];
      }
      g->destroy(pitr);
      c[0]/=g->degree(e);c[1]/=g->degree(e);c[2]/=g->degree(e);
      fprintf(f," %f %f %f\n",c[0],c[1],c[2]);
      i++;
    }
    g->destroy(eitr);
    fprintf(f,"\n</DataArray>\n</Points>\n");
  }

  void writeCells(Ngraph* g, FILE* f) {
    fprintf(f,"<Cells>\n");
    PNgraph* pg = g->publicize();
    fprintf(f,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (lid_t i=0;i<pg->num_local_edges[0];i++) {
      for (lid_t j=pg->pin_degree_list[0][i];j<pg->pin_degree_list[0][i+1];j++)
        fprintf(f,"%lu %lu\n",i+pg->num_local_verts,pg->pin_list[0][j]);
    }
    int off=0;
    fprintf(f,"\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for (lid_t i=0;i<pg->num_local_edges[0];i++) {
      for (lid_t j=pg->pin_degree_list[0][i];j<pg->pin_degree_list[0][i+1];j++)
        fprintf(f," %d\n",off+=2);
    }
    fprintf(f,"\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (lid_t i=0;i<pg->num_local_edges[0];i++) {
      for (lid_t j=pg->pin_degree_list[0][i];j<pg->pin_degree_list[0][i+1];j++)
        fprintf(f," %d\n",3);
    }
    
    fprintf(f,"\n</DataArray>\n</Cells>\n");
  }
  void writePointData(Ngraph* g, FILE* f) {
    fprintf(f,"<PointData>\n");
    fprintf(f,"<DataArray type=\"Int32\" Name=\"VorE\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    agi::GraphVertex* v;
    agi::VertexIterator* vitr = g->begin();
    while ((v=g->iterate(vitr))) {
      fprintf(f,"0\n");
    }
    //Write the position of each edge
    agi::GraphEdge* e;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((e=g->iterate(eitr))) {
      fprintf(f,"1\n");
    }
    g->destroy(eitr);
    fprintf(f,"\n</DataArray>\n");
    fprintf(f,"\n</PointData>\n");
  }
  void writeVTK(Ngraph* g, const char* prefix) {
    char filename[256];
    sprintf(filename,"%s.vtu",prefix);
    FILE* f = fopen(filename,"w");
    fprintf(f,"<VTKFile type=\"UnstructuredGrid\">\n");
    fprintf(f,"<UnstructuredGrid>\n");
    fprintf(f,"<Piece NumberOfPoints=\"%lu\" NumberOfCells=\"%lu\">\n",
            g->numLocalVtxs()+g->numLocalEdges(),g->numLocalPins());
    //writeCellData
    writePoints(g,f);
    writeCells(g,f);
    writePointData(g,f);
    fprintf(f,"</Piece>\n");
    fprintf(f,"</UnstructuredGrid>\n");
    fprintf(f,"</VTKFile>\n");
  }
}
