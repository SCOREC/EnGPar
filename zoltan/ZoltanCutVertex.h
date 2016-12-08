#ifndef __ZOLTAN_CUT_VERTEX__
#define __ZOLTAN_CUT_VERTEX__

#include "ZoltanCallbacks.h"
#include <ngraph.h>
#include <zoltan.h>
namespace zagi {

class ZoltanCutVertex {
 public:
  ZoltanCutVertex(agi::Ngraph* g,int num_parts);
  ~ZoltanCutVertex();

  void run();
 private:
  agi::Ngraph* g;
  Zoltan_Struct* ztn;
  ZOLTAN_ID_PTR import_gids;
  ZOLTAN_ID_PTR import_lids; /* Pointers to nodes imported */
  int *import_procs; /* Proc IDs of procs owning nodes to be imported.*/
  ZOLTAN_ID_PTR export_gids; /* Global ids of nodes exported */
  ZOLTAN_ID_PTR export_lids; /* Pointers to nodes exported */
  int *export_procs;
  int num_imported; /* Number of nodes to be imported. */
  int num_exported; /* Number of nodes to be exported. */
  int *import_to_part;
  int *export_to_part;
  int changes; 
  int lidSz;
  int gidSz;

};

}
#endif
