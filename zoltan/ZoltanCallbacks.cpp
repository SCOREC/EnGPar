#include "ZoltanCallbacks.h"
#include <iostream>
#include <ngraph.h>
namespace zagi {

  int num_obj(void* data,int* ierr) {
    agi::Ngraph*  g = static_cast<agi::Ngraph*>(data);
    *ierr = 0;
    return g->numLocalEdges();
  }

  void obj_list(void* data,int num_gid_entries, int num_lid_entries,ZOLTAN_ID_PTR global_ids,
                ZOLTAN_ID_PTR local_ids, int wgt_dim,float* wgts,int* ierr) {
    agi::Ngraph*  g = static_cast<agi::Ngraph*>(data);

    wgt_dim=0;
    wgts=NULL;
    for (agi::gid_t i=0;i<g->numLocalEdges();i++) {
      local_ids[i] = i;
      global_ids[i] = i;
    }
    *ierr = 0;
  }

  void hg_size(void* data,int* num_lists, int* num_pins,int* format,int* ierr) {
    agi::Ngraph*  g = static_cast<agi::Ngraph*>(data);
    *num_lists = g->numLocalVtxs();
    *num_pins = g->numLocalEdges();
    *format = ZOLTAN_COMPRESSED_EDGE;
    *ierr = 0;
  }

  void hg_data(void* data,int num_gid_entries, int num_vtx_edge,int num_pins,int format, 
               ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr,ZOLTAN_ID_PTR pin_GID,int* ierr) {
    agi::Ngraph*  g = static_cast<agi::Ngraph*>(data);
    agi::VertexIterator* itr = g->begin();
    agi::GraphVertex* vtx;
    int i=0;
    vtxedge_ptr[0] =0;
    int k=0;
    while ((vtx = g->iterate(itr))) {
      vtxedge_GID[i] = g->globalID(vtx);
      agi::lid_t deg = g->degree(vtx);
      if (i<g->numLocalVtxs())
        vtxedge_ptr[i+1] = vtxedge_ptr[i]+deg;
      for (agi::lid_t j=0;j<deg;j++) {
        pin_GID[vtxedge_ptr[i]+j]=vtxedge_ptr[i]+j;
        assert(vtxedge_ptr[i]+j<num_pins);
        assert(vtxedge_ptr[i]+j<vtxedge_ptr[i+1]);
        ++k;
      }
      ++i;
    }
    *ierr = 0;
  }

  void hg_ew_size(void* data, int* num_edges, int* ierr) {
    agi::Ngraph*  g = static_cast<agi::Ngraph*>(data);
    *num_edges = g->numLocalVtxs();
    *ierr = 0;
  }

  void hg_ew(void* data,int num_gid_entries,int num_lid_entries,int num_edges, int edge_weight_dim,
             ZOLTAN_ID_PTR edge_GID,ZOLTAN_ID_PTR edge_LID,float* edge_weights,int* ierr) {
    agi::Ngraph*  g = static_cast<agi::Ngraph*>(data);
    agi::VertexIterator* itr = g->begin();
    agi::GraphVertex* vtx;
    int i=0;
    while ((vtx = g->iterate(itr))) {
      edge_GID[i] = g->globalID(vtx);
      edge_LID[i] = i;
      edge_weights[i] = 1.0/g->degree(vtx);
    }
    *ierr=0;
  }

}
