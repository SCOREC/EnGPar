#include <engpar_support.h>
#include "ngraph.h"
#include <stdio.h>
#include <PCU.h>
#include <cstring>
#include <sys/stat.h>
#include <pcu_io.h>
namespace agi {

  /* bgd file format, field = < description - [ field {field...} | type=[valid values] ] - count >
   * <isHyperGraph - unsigned=[0(false)|1(true)] - 1 >
   * <number of vertices - unsigned=[>0] - 1 >
   * <vertices -
   *   <global id - unsigned=[>0] - 1 >
   *   <weight - double=[>0] - 1 >
   *   -
   *   |V|
   * >
   * <hasCoordinates - unsigned=[0(false)|1(true)] - 1 >
   * { <vertex coordinates -
   *     <xyz - double=[real] - 3 > -
   *     |V|
   *   >
   * }
   * <number of edge types - unsigned=[>0] - 1 >
   * <edge types -
   *   <type - unsigned=[>0] - 1 >
   *   <count - unsigned=[>0] - 1 >
   *   <edges
   *     <global id - unsigned=[>0] - 1 >
   *     <weight - double=[>0] - 1 >
   *     <degree - unsinged=[>0] - 1 >
   *     <pins-to-vertices -
   *       <vertex global id - unsigned=[>0] - 1 > -
   *       |degree(edge[i])|
   *     > -
   *     |E(type)|
   *   > -
   *   number of edge types
   * >
   * <number of ghost vertices - unsigned=[>=0] - 1 >
   * <ghost vertices -
   *   <vertex global id - unsigned=[>0] - 1 >
   *   <owning process - unsigned=[>=0] - 1 >
   *   -
   *   |ghosts|
   * >
   */

  void writeHeader(struct pcu_file* f, Ngraph* g) {
    unsigned int isH = g->isHyper();
    PCU_WRITE_UNSIGNED(f,isH);
  }
  void writeVertices(struct pcu_file* f, Ngraph* g) {
    unsigned int numVtxs = g->numLocalVtxs();
    PCU_WRITE_UNSIGNED(f,numVtxs);
    agi::VertexIterator* itr = g->begin();
    agi::GraphVertex* vtx;
    while ((vtx = g->iterate(itr))) {
      unsigned int gid =  g->globalID(vtx);
      wgt_t w = g->weight(vtx);
      PCU_WRITE_UNSIGNED(f,gid);
      pcu_write_doubles(f,&w,1);          
    }
    //Write the coordinates if they exist
    if (g->hasCoords()) {
      unsigned int one = 1;
      PCU_WRITE_UNSIGNED(f,one);
      agi::VertexIterator* itr = g->begin();
      agi::GraphVertex* vtx;
      while ((vtx = g->iterate(itr))) {
        const coord_t& c = g->coord(vtx);
        coord_t copy;
        for (int i = 0; i < 3; i++)
          copy[i] = c[i];
        pcu_write_doubles(f,copy,3);
      }
    }
    else {
      unsigned int zero = 0;
      PCU_WRITE_UNSIGNED(f,zero);
    }
  }
  void writeEdges(struct pcu_file* f,Ngraph* g,unsigned int t,
                  std::unordered_map<gid_t,part_t>& owns) {
    PCU_WRITE_UNSIGNED(f,t);
    unsigned int numEdges = g->numLocalEdges(t);
    PCU_WRITE_UNSIGNED(f,numEdges);
    agi::EdgeIterator* eitr = g->begin(t);
    agi::GraphEdge* edge;
    while ((edge = g->iterate(eitr))) {
      unsigned int gid = g->globalID(edge);
      wgt_t w = g->weight(edge);
      unsigned int deg = g->degree(edge);
      unsigned int* ps = new unsigned int[deg];
      agi::PinIterator* pitr = g->pins(edge);
      for (unsigned int i = 0;i<deg;i++) {
        agi::GraphVertex* v = g->iterate(pitr);
        ps[i] = g->globalID(v);
        if (g->owner(v)!=PCU_Comm_Self())
          owns[ps[i]] = g->owner(v);
      }
      g->destroy(pitr);
      PCU_WRITE_UNSIGNED(f,gid);
      pcu_write_doubles(f,&w,1);
      PCU_WRITE_UNSIGNED(f,deg);
      pcu_write_unsigneds(f,ps,deg);
      delete [] ps;
    }
    g->destroy(eitr);
  }
  void writeGhosts(struct pcu_file* f, const std::unordered_map<gid_t,part_t>& owns) {
    unsigned int numOwners = owns.size();
    PCU_WRITE_UNSIGNED(f,numOwners);
    std::unordered_map<gid_t,part_t>::const_iterator itr;
    for (itr=owns.begin();itr!=owns.end();itr++) {
      unsigned int first = itr->first;
      PCU_WRITE_UNSIGNED(f,first);
      unsigned int s = itr->second;
      PCU_WRITE_UNSIGNED(f,s);
    }
  }

  //Recursive mkdir function modified from here:
  //http://nion.modprobe.de/blog/archives/357-Recursive-directory-creation.html
  void mkdir_r(const char *dir) {
    char tmp[256];
    char *p = NULL;
    size_t len;
    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);
    if(tmp[len - 1] == '/')
      tmp[len - 1] = 0;
    for(p = tmp + 1; *p; p++)
      if(*p == '/') {
        *p = 0;
        mkdir(tmp, S_IRWXU);
        *p = '/';
      }
  }
  void Ngraph::saveToFile(const char* prefix) {
    char filename[256];
    sprintf(filename,"%s_%d.bgd",prefix,PCU_Comm_Self());
    mkdir_r(filename);
    struct pcu_file* file = pcu_fopen(filename,true,false);
    if (!file) {
      EnGPar_Error_Message("Could not open file for saving: %s#.bgd\n",prefix);
      throw 1;
    }
    writeHeader(file,this); 
    writeVertices(file,this);
    unsigned int nt = numEdgeTypes();
    PCU_WRITE_UNSIGNED(file,nt);
    std::unordered_map<gid_t,part_t> owns;
    for (unsigned int t=0;t<nt;t++)
      writeEdges(file,this,t,owns);
    writeGhosts(file,owns);
    pcu_fclose(file);
  }

  bool readHeader(struct pcu_file* f) {
    unsigned int isHG;
    PCU_READ_UNSIGNED(f,isHG);
    return isHG;
  }
  void readVertices(struct pcu_file* f, std::vector<gid_t>& verts,
                    std::vector<wgt_t>& weights,coord_t*& cs) {
    unsigned int nv;
    PCU_READ_UNSIGNED(f,nv);
    for (unsigned int i=0;i<nv;i++) {
      unsigned int gid;
      wgt_t w;
      PCU_READ_UNSIGNED(f,gid);
      pcu_read_doubles(f,&w,1);
      verts.push_back(gid);
      weights.push_back(w);
    }
    unsigned int hasC;
    PCU_READ_UNSIGNED(f,hasC);
    if (hasC==1) {
      cs = new coord_t[nv];
      for (unsigned int i=0;i<nv;i++) {
        pcu_read_doubles(f,cs[i],3);
      }
    }
    
  }
  void readEdges(struct pcu_file* f,std::vector<gid_t>& edges,
                 std::vector<wgt_t>& eWeights, std::vector<lid_t>& degs,
                 std::vector<gid_t>& pins2v, unsigned int& t){
    PCU_READ_UNSIGNED(f,t);
    unsigned int ne;
    PCU_READ_UNSIGNED(f,ne);
    for (unsigned int i=0;i<ne;i++) {
      unsigned int gid;
      wgt_t w;
      unsigned int deg=0;
      
      PCU_READ_UNSIGNED(f,gid);
      pcu_read_doubles(f,&w,1);
      edges.push_back(gid);
      eWeights.push_back(w);
      PCU_READ_UNSIGNED(f,deg);
      degs.push_back(deg);
      unsigned int pin;
      for (unsigned int j=0;j<deg;j++) {
        PCU_READ_UNSIGNED(f,pin);
        pins2v.push_back(pin);
      }
    }
    
  }
  void readGhosts(struct pcu_file* f, std::unordered_map<gid_t,part_t>& owns) {
    unsigned int num_gs;
    PCU_READ_UNSIGNED(f,num_gs);
    unsigned int v;
    unsigned int o;
    for (unsigned int i =0;i<num_gs;i++) {
      PCU_READ_UNSIGNED(f,v);
      PCU_READ_UNSIGNED(f,o);
      owns[v] = o;
    }
  }

  void Ngraph::loadFromFile(const char* prefix) {
    char filename[256];
    sprintf(filename,"%s_%d.bgd",prefix,PCU_Comm_Self());
    struct pcu_file* file = pcu_fopen(filename,false,false);
    if (!file) {
      EnGPar_Error_Message("Could not open file for loading: %s#.bgd",prefix);
      throw 1;
    }
    bool isHG = readHeader(file);
    if (isHG) {}
    std::vector<gid_t> verts;
    std::vector<wgt_t> weights;
    coord_t* cs = NULL;
    readVertices(file,verts,weights,cs);
    constructVerts(isHG,verts,weights);
    if (cs!=NULL)
      setCoords(cs);
    delete [] cs;
    unsigned int nt;
    PCU_READ_UNSIGNED(file,nt);
    std::unordered_map<gid_t,part_t> owns;
    for (unsigned int t=0;t<nt;t++) {
      std::vector<gid_t> es;
      std::vector<wgt_t> eWeights;
      std::vector<lid_t> degs;
      std::vector<gid_t> pins2v;
      readEdges(file,es,eWeights,degs,pins2v,t);
      constructEdges(es,degs,pins2v,eWeights);
    }
    readGhosts(file,owns);
    constructGhosts(owns);
    pcu_fclose(file);

    if (EnGPar_Check_Verbosity(0));
      //TODO: Add printStats(); to print the loaded graph
  }
}
