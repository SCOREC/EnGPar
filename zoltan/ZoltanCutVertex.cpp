#include "ZoltanCutVertex.h"
#include <iostream>
namespace zagi {

  ZoltanCutVertex::ZoltanCutVertex(agi::Ngraph* graph,int num_parts) {
    g = graph;
    float ver;
    int ret = Zoltan_Initialize(0,0,&ver);
    if (ZOLTAN_OK != ret) {
      fprintf(stderr, "ERROR: Zoltan initialization failed\n");
      exit(1);
    }
    ztn =Zoltan_Create(MPI_COMM_SELF);
    import_gids = NULL;
    import_lids = NULL;
    import_procs = NULL;
    export_gids = NULL;
    export_lids = NULL;
    export_procs = NULL;
    num_imported = 0;
    num_exported = 0;
    import_to_part = NULL;
    export_to_part = NULL;
    changes=0;
    lidSz=1;
    gidSz=1;

    char paramStr[128];

  

    std::string lbMethod = "HYPERGRAPH";
    Zoltan_Set_Param(ztn,"LB_METHOD",lbMethod.c_str());
    lbMethod = "PARTITION";
    Zoltan_Set_Param(ztn,"LB_APPROACH",lbMethod.c_str());

    sprintf(paramStr,"%d",num_parts);
    Zoltan_Set_Param(ztn, "NUM_GLOBAL_PARTS", paramStr);
    Zoltan_Set_Param(ztn, "NUM_LOCAL_PARTS", paramStr);
    sprintf(paramStr,"PARTS");
    Zoltan_Set_Param(ztn, "RETURN_LISTS", paramStr);
    Zoltan_Set_Fn(ztn,ZOLTAN_NUM_OBJ_FN_TYPE,(void (*)())num_obj,(void*)(g));
    Zoltan_Set_Fn(ztn,ZOLTAN_OBJ_LIST_FN_TYPE,(void (*)())obj_list,(void*)(g));
    Zoltan_Set_Fn(ztn,ZOLTAN_HG_SIZE_CS_FN_TYPE,(void (*)())hg_size,(void*)(g));
    Zoltan_Set_Fn(ztn,ZOLTAN_HG_CS_FN_TYPE,(void (*)())hg_data,(void*)(g));
  }

  ZoltanCutVertex::~ZoltanCutVertex() {
    Zoltan_LB_Free_Part(&import_gids, &import_lids, &import_procs, &import_to_part);
    Zoltan_LB_Free_Part(&export_gids, &export_lids, &export_procs, &export_to_part);
    Zoltan_Destroy(&ztn);
    
  }
  void ZoltanCutVertex::run() {
    int ret = Zoltan_LB_Partition(ztn,&changes,&gidSz,&lidSz,&num_imported, 
				  &import_gids,&import_lids, &import_procs,
				  &import_to_part,&num_exported,&export_gids,
				  &export_lids,&export_procs,&export_to_part);
    
    if( ZOLTAN_OK != ret ) {
      fprintf(stderr, "ERROR Zoltan partitioning failed\n");
      exit(EXIT_FAILURE);
    }

  }
  void ZoltanCutVertex::createPtn(agi::EdgePartitionMap& ptn) {

    for (int i=0;i<num_exported;i++) {
      ptn[export_lids[i]] =  export_to_part[i];
    }
  }
}
