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
    //Gather the coordinates of ghost vertices
    coord_t* ghost_cs = new coord_t[g->numGhostVtxs()];
    GhostIterator* g_itr = g->beginGhosts();
    PCU_Comm_Begin();
    while ((v = g->iterate(g_itr))) {
      part_t owner = g->owner(v);
      gid_t gid = g->globalID(v);
      PCU_COMM_PACK(owner,gid);
    }
    PCU_Comm_Send();
    std::vector<std::pair<gid_t,part_t> > requests;
    while (PCU_Comm_Receive()) {
      gid_t gid;
      PCU_COMM_UNPACK(gid);
      requests.push_back(std::make_pair(gid,PCU_Comm_Sender()));
    }

    PCU_Comm_Begin();
    for (size_t i=0;i<requests.size();i++) {
      agi::GraphVertex* my_v = g->findGID(requests[i].first);
      double vals[3];
      const coord_t& c = g->coord(my_v);
      for (int i=0;i<3;i++) {
        vals[i] = c[i];
      }
      PCU_COMM_PACK(requests[i].second,requests[i].first);
      PCU_Comm_Pack(requests[i].second,vals,sizeof(double)*3);
    }
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      gid_t gid;
      PCU_COMM_UNPACK(gid);
      double vals[3];
      PCU_Comm_Unpack(vals,sizeof(double)*3);
      agi::GraphVertex* my_v = g->findGID(gid);
      assert(g->owner(my_v)!=PCU_Comm_Self());
      lid_t lid = g->localID(my_v)-g->numLocalVtxs();
      for (int i=0;i<3;i++)
        ghost_cs[lid][i] = vals[i];
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
          lid_t lid = g->localID(v)-g->numLocalVtxs();
          for (int i=0;i<3;i++)
            c[i]+=ghost_cs[lid][i];
        }
        else {
          const coord_t& cv = g->coord(v);
          c[0]+=cv[0];c[1]+=cv[1];c[2]+=cv[2];
        }
      }
      g->destroy(pitr);
      lid_t deg = g->degree(e)-ghosts;
      for (int i=0;i<3;i++) {
        c[i]/=deg;
      }
      fprintf(f," %f %f %f\n",c[0],c[1],c[2]);
    }
    g->destroy(eitr);
    fprintf(f,"\n</DataArray>\n</Points>\n");
    delete [] ghost_cs;
  }

  void writeCells(Ngraph* g, FILE* f) {
    fprintf(f,"<Cells>\n");
    fprintf(f,"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    GraphEdge* e;
    EdgeIterator* eitr = g->begin(0);
    while ((e = g->iterate(eitr))) {
      GraphVertex* v;
      PinIterator* pitr = g->pins(e);
      while ((v = g->iterate(pitr))) {
        if (g->owner(v)==PCU_Comm_Self())
          fprintf(f,"%d %d\n",g->localID(e)+g->numLocalVtxs(),g->localID(v));
      }
      g->destroy(pitr);
    }
    g->destroy(eitr);
    int off=0;
    fprintf(f,"\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    eitr = g->begin(0);
    while ((e = g->iterate(eitr))) {
      GraphVertex* v;
      PinIterator* pitr = g->pins(e);
      while ((v = g->iterate(pitr))) {
        if (g->owner(v)==PCU_Comm_Self())
          fprintf(f," %d\n",off+=2);
      }
      g->destroy(pitr);
    }
    g->destroy(eitr);
    fprintf(f,"\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    eitr = g->begin(0);
    while ((e = g->iterate(eitr))) {
      GraphVertex* v;
      PinIterator* pitr = g->pins(e);
      while ((v = g->iterate(pitr))) {
        if (g->owner(v)==PCU_Comm_Self())
          fprintf(f," %d\n",3);
      }
      g->destroy(pitr);
    }
    g->destroy(eitr);
    
    fprintf(f,"\n</DataArray>\n</Cells>\n");
  }
  void writePointData(Ngraph* g, FILE* f, GraphTag* tag,etype t) {
    fprintf(f,"<PointData>\n");
    fprintf(f,"<DataArray type=\"Int32\" Name=\"Part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    agi::GraphVertex* v;
    agi::VertexIterator* vitr = g->begin();
    while ((v=g->iterate(vitr))) {
      fprintf(f,"%d\n",PCU_Comm_Self());
    }
    agi::GraphEdge* e;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((e=g->iterate(eitr))) {
      fprintf(f,"%d\n",PCU_Comm_Self());
    }
    g->destroy(eitr);
    fprintf(f,"</DataArray>\n");

    fprintf(f,"<DataArray type=\"Int32\" Name=\"VorE\" NumberOfComponents=\"1\" format=\"ascii\">\n");

    vitr = g->begin();
    while ((v=g->iterate(vitr))) {
      fprintf(f,"0\n");
    }
    eitr = g->begin(0);
    while ((e=g->iterate(eitr))) {
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
  
  gid_t getNumCells(Ngraph* g) {
    gid_t total = 0;
    EdgeIterator* eitr = g->begin(0);
    GraphEdge* edge;
    while ((edge = g->iterate(eitr))) {
      PinIterator* pitr = g->pins(edge);
      GraphVertex* v;
      while ((v = g->iterate(pitr))) {
        total+=(g->owner(v)==PCU_Comm_Self());
      }
      g->destroy(pitr);
    }
    g->destroy(eitr);
    return total;
  }

  void writePPoints(FILE* pf) {
    fprintf(pf,"<PPoints>\n");
    fprintf(pf,"<PDataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\" Name=\"Points\"/>\n");
    fprintf(pf,"</PPoints>\n");

  }
  void writePPointData(FILE* pf, GraphTag* tag, etype t) {
    fprintf(pf,"<PPointData>\n");
    fprintf(pf,"<PDataArray type=\"Int32\" Name=\"Part\" NumberOfComponents=\"1\" format=\"ascii\"/>\n");
    fprintf(pf,"<PDataArray type=\"Int32\" Name=\"VorE\" NumberOfComponents=\"1\" format=\"ascii\"/>\n");
    if (tag!=NULL&&t!=NO_TYPE) {
      fprintf(pf,"<PDataArray type=\"Int32\" Name=\"Tag\" NumberOfComponents=\"1\" format=\"ascii\"/>\n");
    }
    fprintf(pf,"</PPointData>\n");
  }

  void writePSources(FILE* pf, const char* prefix) {
    char filename[256];
    for (int i=0;i<PCU_Comm_Peers();i++) {
      sprintf(filename,"%s%d.vtu",prefix,i);
      fprintf(pf,"<Piece Source=\"%s\"/>\n",filename);
    }
  }
  
  void writePVTU(FILE* pf,const char* prefix, GraphTag* tag, etype t) {
      fprintf(pf,"<VTKFile type=\"PUnstructuredGrid\">\n");
      fprintf(pf,"<PUnstructuredGrid GhostLevel=\"0\">\n");
      writePPoints(pf);
      writePPointData(pf,tag,t);
      writePSources(pf,prefix);
      fprintf(pf,"</PUnstructuredGrid>\n");
      fprintf(pf,"</VTKFile>\n");

  }
  void writeVTK(Ngraph* g, const char* prefix,GraphTag* tag,etype t) {
    char filename[256];
    sprintf(filename,"%s%d.vtu",prefix,PCU_Comm_Self());
    FILE* f = fopen(filename,"w");
    fprintf(f,"<VTKFile type=\"UnstructuredGrid\">\n");
    fprintf(f,"<UnstructuredGrid>\n");
    gid_t numCells = getNumCells(g);
    fprintf(f,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%ld\">\n",
            g->numLocalVtxs()+g->numLocalEdges(),numCells);
    writePoints(g,f);
    writeCells(g,f);
    writePointData(g,f,tag,t);
    fprintf(f,"</Piece>\n");
    fprintf(f,"</UnstructuredGrid>\n");
    fprintf(f,"</VTKFile>\n");
    fclose(f);

    if (!PCU_Comm_Self()) {
      char pfilename[256];
      sprintf(pfilename,"%s.pvtu",prefix);
      FILE* pf = fopen(pfilename,"w");
      writePVTU(pf,prefix,tag,t);
      fclose(pf);
    }
  }
}
