#include "agi.h"
#include <cstdlib>
namespace agi {

Ngraph::Ngraph(bool isUndirected) {
  isUnd = isUndirected;
  num_verts=num_edges=0;
  weights=out_weights=NULL;
  out_degree_list=out_vertices=NULL;
}

Ngraph::~Ngraph() {
  if (weights)
    delete [] weights;
  if (out_degree_list)
    delete [] out_degree_list;
  if (out_vertices)
    delete [] out_vertices;
  if (out_weights)
    delete [] out_weights;
}

  
//Protected functions
void Ngraph::create_csr(int nv, int ne, int* srcs,
                          int* dsts, int* wgts) {
  num_verts = nv;
  num_edges = ne;
  weights = new int[num_verts];
  out_vertices = new int[num_edges];
  out_weights = new int[num_edges];
  out_degree_list = new int[num_verts+1];

  for (int i = 0; i < num_edges; ++i)
    out_vertices[i] = 0;
  for (int i = 0; i < num_edges; ++i)
    out_weights[i] = 0;
  for (int i = 0; i < num_verts+1; ++i)
    out_degree_list[i] = 0;

  int* temp_counts = new int[num_verts];
  for (int i = 0; i < num_verts; ++i)
    temp_counts[i] = 0;
  for (int i = 0; i < num_edges; ++i)
    ++temp_counts[srcs[i]];
  for (int i = 0; i < num_verts; ++i)
    out_degree_list[i+1] = out_degree_list[i] + temp_counts[i];
  std::copy(out_degree_list, out_degree_list + num_verts, temp_counts);
  for (int i = 0; i < num_edges; ++i) {
    out_vertices[temp_counts[srcs[i]]] = dsts[i];
    out_weights[temp_counts[srcs[i]]++] = wgts[i];
  }
  
  delete [] temp_counts;

}

}
