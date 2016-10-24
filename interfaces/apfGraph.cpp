#include "apfGraph.h"
#include <stdio.h>
//TODO: may want to replace with something more light weight in future
#include <vector>
#include <apfNumbering.h>

namespace agi {

apfGraph::apfGraph(apf::Mesh* mesh,int primary_dimension,
                   int secondary_dimension) : Ngraph() {
  //Error checking on primary and secondary dims
  int dim = mesh->getDimension();
  if (primary_dimension>dim ||
      primary_dimension<0 ||
      secondary_dimension>dim ||
      secondary_dimension<0) {
    fprintf(stderr,"[ERROR] primary or secondary dimensions are invalid\n");
    throw 1;
  }
  
  if (primary_dimension==secondary_dimension) {
    fprintf(stderr,"[ERROR] primary and secondary dimensions cannot be equal\n");
    throw 2;
  }
  
  //TODO: create vertex and edge weights
  //TODO: copy over coordinates
  m = mesh;
  num_verts = m->count(primary_dimension);
  if (num_verts==0) {
    fprintf(stderr,"[ERROR] Mesh is empty, exiting\n");
    throw 3;
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
    int nents = adj.getSize();
    for (int i=0;i<nents;i++) {
      apf::MeshEntity* bridge = adj[i];
      apf::Adjacent neighbors;
      m->getAdjacent(bridge,primary_dimension,neighbors);
      int npri = neighbors.size();
      for (int j=0;j<npri;j++) {
        apf::MeshEntity* neighbor = neighbors[j];
        if (neighbor!=ent) {
          out_verts.push_back(getNumber(global_nums,neighbor,0));
        }
      }

    }
    out_degree_list[vert]=out_verts.size()-out_degree_list[vert-1];
    vert++;
  }
  num_edges = out_verts.size();
  out_vertices = new int[num_edges];
  std::copy(out_verts.begin(),out_verts.end(),out_vertices);

}
  

}
