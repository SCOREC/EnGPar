#include "sellCSigma.h"
#include <cstring>
namespace ssg {

  agi::Ngraph* convertFromAGI(agi::Ngraph* old, lid_t C, lid_t sigma) {
    agi::Ngraph* g = new SellCSigma(C,sigma);
    g->isHyperGraph = old->isHyperGraph;
    g->num_types = old->num_types;
    g->num_local_verts = old->num_local_verts;
    g->num_global_verts = old->num_global_verts;
    g->num_ghost_verts = old->num_ghost_verts;

    //Get Edge Counts
    for (etype t=0;t<num_types;t++) {
      g->num_local_edges[t] = old->num_local_edges[t];
      g->num_global_edges[t] = old->num_global_edges[t];
      g->num_local_pins[t] = old->num_local_pins[t];
      g->num_global_pins[t] = old->num_global_pins[t];
    }

    //Apply sigma sorting
    etype sort_type = 0;
    //Sort each block of sigma vertices
    //temporary container for degree_list (pair of <deg,local id>
    std::pair<lid_t,lid_t>* degree_list = new std::pair<lid_t,lid_t>[g->num_local_verts];
    agi::GraphVertex* v;
    agi::VertexIterator* vitr = old->begin();
    int count=0;
    while ((v = old->iterate(vitr))) {
      degree_list[count].first = old->degree(v,sort_type);
      degree_list[count++].second = old->localID(v);
    }

    if (sigma>1) {
      lid_t i;
      for (i =0;i<g->num_local_verts-sigma,i+=sigma)
        std::sort(degree_list+i,degree_list+i+sigma,operator>);
      std::sort(degree_list+i,degree_list+num_local_verts,operator>);
    }

    //Set gid->lid and lid->gid containers
    g->local_unmap = new gid_t[num_local_verts];
    for (lid_t i=0;i<g->num_local_verts;i++) {
      g->local_unmap[i] = old->local_unmap[degree_list[i].second];
      g->vtx_mapping[g->local_unmap[i]] = i;
    }
    
    //TODO: Reorder theses
    /*
    if (old->original_owners)
      memcpy(g->original_owners,old->original_owners,g->num_local_verts);
    if (old->local_weights)
      memcpy(g->local_weights,old->local_weights,g->num_local_verts);
    if (old->local_coords)
      memcpy(g->local_coords,old->local_coords,g->num_local_verts);
    */
    
    for (etype t=0;t<num_types;t++) {
      //TODO: Sigma sort edges

      //TODO: Reorder these
      if (old->edge_weights[t])
        memcpy(g->edge_weights[t],old->edge_weights[t],g->num_local_edges[t]);

      //Get the degrees of each block of C vertices
      lid_t numBlocks = g->num_local_verts / C + (g->num_local_verts % C != 0);
      g->degree_list[t] = new lid_t[numBlocks+1];
      lid_t total = 0;
      g->degree_list[t]=0;
      for (lid_t i =0;i<numBlocks,i++) {
        g->degree_list[t][i+1] = 0;
        for (lid_t j = i * C; j < (i + 1) * C && j < g->num_local_verts; j++)
          g->degree_list[t][i+1] = std::max(g->degree_list[t][i+1],degree_list[j].first);
        total=(g->degree_list[t][i+1]+=total);
      }
      
      //build the padded edge_list
      g->edge_list[t] = new lid_t[total*C];
      int index =0;
      for (lid_t i =0;i<numBlocks,i++) {
        for (lid_t deg = g->degree_list[t][i];deg < g->degree_list[t][i+1];deg++) {
          for (lid_t j = i * C; j < (i + 1) * C; j++) {
            //if there is a edge at this location
            if (j<g->num_local_verts && deg<degree_list[j].first) {
              lid_t old_lid = degree_list[j].second;
              lid_t ind = old->degree_list[old_lid]+deg;
              //NOTE: this isnt quite right yet, but easier to test this way
              g->edge_list[t][index++] = old->edge_list[ind];
            }
            //else add padding
            else
              g->edge_list[t][index++] = -1;
          }
        }
      }
      delete [] degree_list;
    }
    return g;
  }
}
