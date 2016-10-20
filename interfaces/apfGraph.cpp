#include "apfGraph.h"
#include <stdio.h>
//TODO: may want to replace with something more light weight in future
#include <vector>
#include <apfNumbering.h>

namespace agi {

APFGraph::APFGraph(apf::Mesh* mesh,int primary_dimension,
                   int secondary_dimension) : Ngraph() {
  //TODO: Error checking on primary and secondary dims
  //TODO: create vertex and edge weights
  //TODO: copy over coordinates
  m = mesh;
  num_verts = m->count(primary_dimension);
  if (num_verts==0) {
    fprintf(stderr,"[ERROR] Mesh is empty, exiting\n");
    throw 1;
  }

  //Create a global numbering on the mesh over primary_dimension
  apf::Numbering* numbers = apf::numberOwnedDimension(m,"primary_ids",
                                                      primary_dimension);
  apf::GlobalNumbering* global_nums = apf::makeGlobal(numbers);
  apf::synchronize(global_nums);
  
  num_edges = 0;
  out_degree_list = new int[num_verts];
  std::vector<int> out_verts;
  apf::MeshEntity* ent;
  apf::MeshIterator* itr= m->begin(primary_dimension);
  int vert=0;
  while ((ent = m->iterate(itr))) {
    if (!m->isOwned(ent))
      continue;
    apf::Adjacent adj;
    m->getAdjacent(ent,secondary_dimension,adj);
    int degree = adj.getSize();
    num_edges+=degree;
    out_degree_list[vert++] =degree;
    for (int i=0;i<degree;i++) {
      out_verts.push_back(getNumber(global_nums,adj[i],0));
    }
  }
  out_vertices = new int[num_edges];
  std::copy(out_verts.begin(),out_verts.end(),out_vertices);

}
  

}
