#ifndef __ZOLTAN_CALLBACKS_AGI__
#define __ZOLTAN_CALLBACKS_AGI__

#include <zoltan.h>

namespace zagi {

  int num_obj(void* data,int* ierr);

  void obj_list(void* data,int num_gid_entries, int num_lid_entries,ZOLTAN_ID_PTR global_ids,
		ZOLTAN_ID_PTR local_lids, int wgt_dim,float* wgts,int* ierr);

  void hg_size(void* data,int* num_lists, int* num_pins,int* format,int* ierr);

  void hg_data(void* data,int num_gid_entries, int num_vtx_edge,int num_pins,int format, 
	       ZOLTAN_ID_PTR vtxedge_GID, int* vtxedge_ptr,ZOLTAN_ID_PTR pin_GID,int* ierr);

  void hg_ew_size(void* data, int* num_edges, int* ierr);

  void hg_ew(void* data,int num_gid_entries,int num_lid_entries,int num_edges, int edge_weight_dim,
	     ZOLTAN_ID_PTR edge_GID,ZOLTAN_ID_PTR edge_LID,float* edge_weights,int* ierr);
}

#endif
