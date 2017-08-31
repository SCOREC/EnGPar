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
    while ((v=g->iterate(vitr))) {
      const coord_t& c = g->coord(v);
      fprintf(f," %f %f %f\n",c[0],c[1],c[2]);
    }
    //Write the position of each edge
    agi::GraphEdge* e;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((e=g->iterate(eitr))) {
      coord_t c = {0,0,0};
      agi::PinIterator* pitr = g->pins(e);
      int ghosts=0;
      while ((v = g->iterate(pitr))) {
        if (g->owner(v)!=PCU_Comm_Self()){
          ghosts++;
          continue;
        }
        const coord_t& cv = g->coord(v);
        c[0]+=cv[0];c[1]+=cv[1];c[2]+=cv[2];
      }
      g->destroy(pitr);
      lid_t deg = g->degree(e)-ghosts;
      for (int i=0;i<3;i++) {
        if (deg==1)
          c[i]+=.1;
        c[i]/=deg;
      }
      fprintf(f," %f %f %f\n",c[0],c[1],c[2]);
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
        fprintf(f,"%ld %ld\n",i+pg->num_local_verts,pg->pin_list[0][j]);
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
  void writePointData(Ngraph* g, FILE* f, GraphTag* tag,etype t) {
    fprintf(f,"<PointData>\n");
    fprintf(f,"<DataArray type=\"Int32\" Name=\"VorE\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    agi::VertexIterator* vitr = g->begin();
    while (g->iterate(vitr)) {
      fprintf(f,"0\n");
    }
    agi::EdgeIterator* eitr = g->begin(0);
    while (g->iterate(eitr)) {
      fprintf(f,"1\n");
    }
    g->destroy(eitr);
    fprintf(f,"</DataArray>\n");
    if (tag!=NULL&&t!=NO_TYPE) {
      fprintf(f,"<DataArray type=\"Int32\" Name=\"Tag\" NumberOfComponents=\"1\" format=\"ascii\">\n");
      agi::GraphVertex* v;
      agi::VertexIterator* vitr = g->begin();
      while ((v=g->iterate(vitr))) {
        if (t==VTX_TYPE)
          fprintf(f,"%d\n",g->getIntTag(tag,v));
        else
          fprintf(f,"-1\n");
      }
      agi::GraphEdge* e;
      agi::EdgeIterator* eitr = g->begin(0);
      while ((e=g->iterate(eitr))) {
        if (t==0)
          fprintf(f,"%d\n",g->getIntTag(tag,e));
        else
          fprintf(f,"-1\n");
      }
      g->destroy(eitr);
      fprintf(f,"</DataArray>\n");
    }
    fprintf(f,"</PointData>\n");
  }
  void writeVTK(Ngraph* g, const char* prefix,GraphTag* tag,etype t) {
    char filename[256];
    sprintf(filename,"%s%d.vtu",prefix,PCU_Comm_Self());
    FILE* f = fopen(filename,"w");
    fprintf(f,"<VTKFile type=\"UnstructuredGrid\">\n");
    fprintf(f,"<UnstructuredGrid>\n");
    fprintf(f,"<Piece NumberOfPoints=\"%ld\" NumberOfCells=\"%ld\">\n",
            g->numLocalVtxs()+g->numLocalEdges(),g->numLocalPins());
    writePoints(g,f);
    writeCells(g,f);
    writePointData(g,f,tag,t);
    fprintf(f,"</Piece>\n");
    fprintf(f,"</UnstructuredGrid>\n");
    fprintf(f,"</VTKFile>\n");
    fclose(f);
  }
}
