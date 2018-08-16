#include "engpar_zoltan.h"
#include <PCU.h>
#include <agiMigration.h>
#include <ngraph.h>
#include <stdexcept>
#include <engpar_support.h>
#include "engpar_split_input.h"
#ifdef HAS_ZOLTAN
#include <zoltan.h>
#endif

namespace engpar {
#ifdef HAS_ZOLTAN
  void setParameters(SplitInput* input, struct Zoltan_Struct* zz, int target_parts);
  void setCallbacks(SplitInput* input, struct Zoltan_Struct* zz);
  
  agi::Migration* EnGPar_Zoltan(SplitInput* input, int target_parts, bool isLocal) {
    float version;
    int ret = Zoltan_Initialize(0, NULL, &version);
    if (ret != ZOLTAN_OK) {
      if (!PCU_Comm_Self())
        EnGPar_Error_Message("Zoltan failed to initialize\n");
      throw std::runtime_error("Zoltan failed to initialize");
    }

    MPI_Comm comm;
    if (isLocal)
      comm = MPI_COMM_SELF;
    else
      comm = PCU_Get_Comm();

    struct Zoltan_Struct* zz = Zoltan_Create(comm);
    //TODO: make some of these parameters changable in the input
    Zoltan_Set_Param(zz,"LB_METHOD", "HYPERGRAPH");
    Zoltan_Set_Param(zz,"HYPERGRAPH_PACKAGE", "PHG");
    Zoltan_Set_Param(zz,"LB_APPROACH", "PARTITION");

    setParameters(input, zz, target_parts);
    setCallbacks(input, zz);

    int changes = 0;
    int ngids = 1, nlids = 1;
    int num_import = 0;
    ZOLTAN_ID_PTR import_gids = NULL, import_lids = NULL;
    int* import_procs = NULL, *import_parts = NULL;
    int num_export = 0;
    ZOLTAN_ID_PTR export_gids = NULL, export_lids = NULL;
    int* export_procs = NULL, *export_parts = NULL;

    ret = Zoltan_LB_Partition(zz, &changes, &ngids, &nlids,
                              &num_import, &import_gids, &import_lids,
                              &import_procs, &import_parts,
                              &num_export, &export_gids, &export_lids,
                              &export_procs, &export_parts);
    if (ret != ZOLTAN_OK) {
      if (!PCU_Comm_Self())
        EnGPar_Error_Message("Zoltan_LB_Partition failed\n");
      throw std::runtime_error("Zoltan_LB_Partition failed");
    }
 
    agi::Migration* plan = new agi::Migration(input->g);

    agi::PNgraph* pg = input->g->publicize();
    for (int i = 0; i < num_export; ++i) {
      agi::GraphVertex* v = pg->getVertex(export_lids[i]);
      plan->insert(std::make_pair(v,export_procs[i]));
    }

    Zoltan_LB_Free_Part(&import_gids, &import_lids, &import_procs, &import_parts);
    Zoltan_LB_Free_Part(&export_gids, &export_lids, &export_procs, &export_parts);
    Zoltan_Destroy(&zz);
    
    return plan;
  }

  void setParameters(SplitInput* input, struct Zoltan_Struct* zz, int target_parts) {
    Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","1");
    Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
    char param[100];
    sprintf(param,"%d",target_parts);
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS",param);
    Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS","1");
    Zoltan_Set_Param(zz,"RETURN_LISTS","EXPORT");
    sprintf(param,"%f",input->tolerance);
    Zoltan_Set_Param(zz,"IMBALANCE_TOL",param);
    if (!EnGPar_Check_Verbosity(2))
      Zoltan_Set_Param(zz,"DEBUG_LEVEL","0");
  }

  /******  Callback functions  ******/
  //ZOLTAN_NUM_OBJ_FN
  int numVertices(void* data, int *ierr) {
    agi::Ngraph* g = static_cast<agi::Ngraph*>(data);
    *ierr = ZOLTAN_OK;
    return g->numLocalVtxs();
  }

  //ZOLTAN_OBJ_LIST_FN
  void getVertices(void* data, int ngid, int nlid,
                   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                   int nwgts, float* wgts, int* ierr) {
    agi::Ngraph* g = static_cast<agi::Ngraph*>(data);

    agi::GraphVertex* v;
    agi::VertexIterator* vitr = g->begin();
    while ((v = g->iterate(vitr))) {
      agi::lid_t l = g->localID(v);
      lids[l*nlid] = l;
      gids[l*ngid] = g->globalID(v);
      wgts[l*nwgts] = g->weight(v);
    }
    *ierr = ZOLTAN_OK;
  }

  //ZOLTAN_HG_SIZE_CS_FN
  void numHyperedges(void* data, int* numVtxs, int* numPins, int* format, int* ierr) {
    agi::Ngraph* g = static_cast<agi::Ngraph*>(data);
    *numVtxs = g->numLocalVtxs();
    //TODO: this number is inaccurate in parallel
    *numPins = g->numLocalPins();
    *format = ZOLTAN_COMPRESSED_VERTEX;
    *ierr = ZOLTAN_OK;
  }

  //ZOLTAN_HG_CS_FN
  void getHyperedges(void* data, int ngid, int, int, int,
                     ZOLTAN_ID_PTR vtx_gids, int* offsets, ZOLTAN_ID_PTR edge_gids, int* ierr) {
    agi::Ngraph* g = static_cast<agi::Ngraph*>(data);
    agi::GraphVertex* v;
    agi::VertexIterator* vitr = g->begin();
    int index = 0;
    offsets[0] = 0;
    while ((v = g->iterate(vitr))) {
      agi::lid_t lid = g->localID(v);
      vtx_gids[lid*ngid] = g->globalID(v);
      if (lid < g->numLocalVtxs()-1)
        offsets[lid+1] = offsets[lid] + g->degree(v); 
      agi::GraphEdge* e;
      agi::EdgeIterator* eitr = g->edges(v);
      while ((e = g->iterate(eitr))) {
        edge_gids[index++] = g->globalID(e);
      }
      g->destroy(eitr);
    }
    *ierr = ZOLTAN_OK;
  }

  //ZOLTAN_HG_SIZE_EDGE_WTS_FN
  void numHyperedgeWeights(void *data, int* num_edges, int* ierr) {
    agi::Ngraph* g = static_cast<agi::Ngraph*>(data);
    *num_edges = g->numLocalEdges();
    *ierr = ZOLTAN_OK;
  }
  
  //ZOLTAN_HG_EDGE_ETS_FN
  void getHyperedgeWeights(void* data, int ngids, int nlids, int, int nweights,
                           ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, float* wgts, int* ierr) {
    agi::Ngraph* g = static_cast<agi::Ngraph*>(data);
    agi::GraphEdge* e;
    agi::EdgeIterator* eitr = g->begin(0);
    while ((e = g->iterate(eitr))) {
      agi::lid_t lid = g->localID(e);
      lids[nlids*lid] = lid;
      gids[ngids*lid] = g->globalID(e);
      wgts[nweights*lid] = g->weight(e);
    }
    g->destroy(eitr);
    *ierr = ZOLTAN_OK;
  }
  
  //TODO: add graph callbacks for when using a traditional graph

  void setCallbacks(SplitInput* input, struct Zoltan_Struct* zz) {
    Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE,
                  (void (*)())numVertices,(void*) input->g);
    Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE,
                  (void (*)())getVertices,(void*) input->g);
    Zoltan_Set_Fn(zz, ZOLTAN_HG_SIZE_CS_FN_TYPE,
                  (void (*)())numHyperedges,(void*) input->g);
    Zoltan_Set_Fn(zz, ZOLTAN_HG_CS_FN_TYPE,
                  (void (*)())getHyperedges,(void*) input->g);
    Zoltan_Set_Fn(zz, ZOLTAN_HG_SIZE_EDGE_WTS_FN_TYPE,
                  (void (*)())numHyperedgeWeights,(void*) input->g);
    Zoltan_Set_Fn(zz, ZOLTAN_HG_EDGE_WTS_FN_TYPE,
                  (void (*)())getHyperedgeWeights,(void*) input->g);
  }
#else
  agi::Migration* EnGPar_Zoltan(SplitInput*, int, bool) {
    throw std::runtime_error("ZOLTAN not found\n");
  }
#endif
}
