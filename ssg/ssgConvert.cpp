#include "sellCSigma.h"
#include <cstring>
#include <algorithm>

namespace ssg {
  SellCSigma::SellCSigma(agi::Ngraph* g, lid_t c, lid_t sig) : agi::Ngraph() {
    isSellCSigma = true;
    chunk_size=C=c;
    sigma = sig;

    PNgraph* old = g->publicize();
    
    isHyperGraph = old->isHyperGraph;
    num_types = old->num_types;
    num_local_verts = old->num_local_verts;
    num_global_verts = old->num_global_verts;
    num_ghost_verts = old->num_ghost_verts;

    //Get Edge Counts
    for (etype t=0;t<num_types;t++) {
      num_local_edges[t] = old->num_local_edges[t];
      num_global_edges[t] = old->num_global_edges[t];
      num_local_pins[t] = old->num_local_pins[t];
      num_global_pins[t] = old->num_global_pins[t];
    }

    //Apply sigma sorting
    etype sort_type = 0;
    //Sort each block of sigma vertices
    typedef std::pair<lid_t,lid_t> pair_t;
    //temporary container for degree_list (pair of <deg,local id>)
    pair_t* temp_degree_list = new pair_t[num_local_verts];
    agi::GraphVertex* v;
    agi::VertexIterator* vitr = g->begin();
    int count=0;
    while ((v = g->iterate(vitr))) {
      temp_degree_list[count].first = g->degree(v,sort_type);
      temp_degree_list[count++].second = g->localID(v);
    }

    if (sigma>1) {
      lid_t i;
      for (i =0;i<num_local_verts-sigma;i+=sigma)
        std::sort(temp_degree_list+i,temp_degree_list+i+sigma,std::greater<pair_t>());
      std::sort(temp_degree_list+i,temp_degree_list+num_local_verts,std::greater<pair_t>());
    }

    //Set gid->lid and lid->gid containers
    local_unmap = new gid_t[num_local_verts];
    for (lid_t i=0;i<num_local_verts;i++) {
      local_unmap[i] = old->local_unmap[temp_degree_list[i].second];
      vtx_mapping[local_unmap[i]] = i;
    }
    
    //TODO: Reorder these and determine if necessary
    /*
    if (old->original_owners)
      memcpy(original_owners,old->original_owners,num_local_verts);
    if (old->local_weights)
      memcpy(local_weights,old->local_weights,num_local_verts);
    if (old->local_coords)
      memcpy(local_coords,old->local_coords,num_local_verts);
    */
    num_vtx_chunks = num_local_verts / C + (num_local_verts % C != 0);
    
    for (etype t=0;t<num_types;t++) {
      //TODO: Sigma sort edges

      //TODO: Reorder these
      if (old->edge_weights[t]) {
        edge_weights[t] = new agi::wgt_t[num_local_edges[t]];
        memcpy(edge_weights[t],old->edge_weights[t],num_local_edges[t]);
      }

      //Get the degrees of each block of C vertices
      degree_list[t] = new lid_t[num_vtx_chunks+1];
      lid_t total = 0;
      degree_list[t][0] = 0;
      for (lid_t i =0;i<num_vtx_chunks;i++) {
        degree_list[t][i+1] = 0;
        for (lid_t j = i * C; j < (i + 1) * C && j < num_local_verts; j++)
          degree_list[t][i+1] = std::max(degree_list[t][i+1],temp_degree_list[j].first);
        total=(degree_list[t][i+1]+=total);
      }
      
      //build the padded edge_list
      edge_list[t] = new lid_t[total*C];
      int index =0;
      for (lid_t i =0;i<num_vtx_chunks;i++) {
        for (lid_t deg = 0;deg < degree_list[t][i+1]-degree_list[t][i];deg++) {
          for (lid_t j = i * C; j < (i + 1) * C; j++) {
            //if there is a edge at this location
            if (j<num_local_verts && deg<temp_degree_list[j].first) {
              lid_t old_lid = temp_degree_list[j].second;
              lid_t ind = old->degree_list[t][old_lid]+deg;
              //NOTE: this isnt quite right yet, but easier to test this way
              edge_list[t][index++] = old->edge_list[t][ind];
            }
            //else add padding
            else
              edge_list[t][index++] = -1;
          }
        }
      }
      delete [] temp_degree_list;
    }
  }
  agi::Ngraph* convertFromAGI(agi::Ngraph* old, lid_t C, lid_t sigma) {
    return new SellCSigma(old,C,sigma);
  }

}
