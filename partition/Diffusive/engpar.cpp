#include "engpar.h"
#include <PCU.h>
#include "src/engpar_sides.h"

namespace engpar {
  wgt_t getWeight(agi::Ngraph* g, int dimension, bool countGhosts) {
    wgt_t w =0;
    if (dimension==-1) {
      
      //calculate the total weight of the vertices
      agi::GraphVertex* vtx;
      agi::VertexIterator* vitr = g->begin();
      while ((vtx = g->iterate(vitr))) 
        w+=g->weight(vtx);
      if (countGhosts)
	w+=g->numGhostVtxs();
    }
    else {
      agi::GraphEdge* edge;
      agi::EdgeIterator* eitr = g->begin(dimension);
      while ((edge = g->iterate(eitr))) 
        w+=g->weight(edge);
      g->destroy(eitr);
    }
    return w;
  }

  double EnGPar_Get_Imbalance(wgt_t my_val) {
    wgt_t max,total;
    double avg;
    MPI_Datatype type = MPI_DOUBLE;
    MPI_Allreduce(&my_val,&max,1,type,MPI_MAX,PCU_Get_Comm());
    MPI_Allreduce(&my_val,&total,1,type,MPI_SUM,PCU_Get_Comm());
    avg = total*1.0/PCU_Comm_Peers();
    return max*1.0/avg;
  }

  void printImbalances(agi::Ngraph* g) {
    for (agi::etype type = -1;type<g->numEdgeTypes();type++) {
      agi::wgt_t w = getWeight(g,type);
      double imb = EnGPar_Get_Imbalance(w);
      if (!PCU_Comm_Self())
        printf("%2.4f ",imb);
    }
    if (!PCU_Comm_Self())
      printf("\n");
  }
  
  void printMaxMinAvgImb(agi::lid_t my_val,std::string prefix) {
    agi::gid_t in = my_val;
    agi::gid_t max,min,total;
    double avg,imb;
    MPI_Datatype type = MPI_UNSIGNED_LONG;
    MPI_Allreduce(&in,&max,1,type,MPI_MAX,PCU_Get_Comm());
    MPI_Allreduce(&in,&min,1,type,MPI_MIN,PCU_Get_Comm());
    MPI_Allreduce(&in,&total,1,type,MPI_SUM,PCU_Get_Comm());
    avg = total*1.0/PCU_Comm_Peers();
    imb = max*1.0/avg;
    if (!PCU_Comm_Self())
      printf("%s <max, min, avg, imb> = <%lu %lu %f %f>\n",
             prefix.c_str(),max,min,avg,imb);
  }
  void printMaxMinAvgImb(agi::wgt_t my_val,std::string prefix) {
    agi::wgt_t max,min,total;
    double avg,imb;
    MPI_Datatype type = MPI_DOUBLE;
    MPI_Allreduce(&my_val,&max,1,type,MPI_MAX,PCU_Get_Comm());
    MPI_Allreduce(&my_val,&min,1,type,MPI_MIN,PCU_Get_Comm());
    MPI_Allreduce(&my_val,&total,1,type,MPI_SUM,PCU_Get_Comm());
    avg = total*1.0/PCU_Comm_Peers();
    imb = max*1.0/avg;
    if (!PCU_Comm_Self())
      printf("%s <max, min, avg, imb> = <%f %f %f %f>\n",
             prefix.c_str(),max,min,avg,imb);
  }

  
  void evaluatePartition(agi::Ngraph* g,std::string) {
    agi::wgt_t* my_vals = new agi::wgt_t[2+g->numEdgeTypes()+3];
    agi::wgt_t* min = new agi::wgt_t[2+g->numEdgeTypes()+3];
    agi::wgt_t* max = new agi::wgt_t[2+g->numEdgeTypes()+3];
    agi::wgt_t* total = new agi::wgt_t[2+g->numEdgeTypes()+3];
    agi::wgt_t* avg = new agi::wgt_t[2+g->numEdgeTypes()+3];

    //Disconnected comp
    my_vals[0] = 0;

    //Neighbors & Part Boundary
    DiffusiveInput* in = static_cast<DiffusiveInput*>(createDiffusiveInput(g,0));
    Sides* sides = makeSides(in);
    my_vals[1] = sides->size();
    my_vals[2] = sides->total();
    delete in;

    //Vertex Imbalance
    my_vals[3] = getWeight(g,-1);
    my_vals[4] = getWeight(g,-1,true);

    //Edge type imbalance
    for (agi::etype t = 0;t<g->numEdgeTypes();t++) {
      my_vals[t+5] = getWeight(g,t);
    }


    //Empty parts
    int empty = my_vals[3]==0;

    MPI_Reduce(my_vals,max,5+g->numEdgeTypes(),MPI_DOUBLE,MPI_MAX,0,PCU_Get_Comm());
    MPI_Reduce(my_vals,min,5+g->numEdgeTypes(),MPI_DOUBLE,MPI_MIN,0,PCU_Get_Comm());
    MPI_Reduce(my_vals,total,5+g->numEdgeTypes(),MPI_DOUBLE,MPI_SUM,0,PCU_Get_Comm());
    empty = PCU_Add_Int(empty);

    if (!PCU_Comm_Self()) {
      for (int i=0;i<5+g->numEdgeTypes();i++) {
        avg[i] = total[i]/PCU_Comm_Peers();
      }
      printf("ENGPAR STATUS: Empty Parts: %d\n"
             "ENGPAR STATUS: Disconnected Components: <min,max,avg,imb> %.3f %.3f %.3f %.3f\n"
             "ENGPAR STATUS: Neighbors: <min,max,avg,imb> %.3f %.3f %.3f %.3f\n"
             "ENGPAR STATUS: Edge Cut: <min,max,avg,imb> %.3f %.3f %.3f %.3f\n"
             "ENGPAR STATUS: Local Vertex Imbalance: <min,max,avg,imb> %.3f %.3f %.3f %.3f\n"
             "ENGPAR STATUS: Total Vertex Imbalance: <min,max,avg,imb> %.3f %.3f %.3f %.3f\n",
             empty,max[0],min[0],avg[0],total[0]/avg[0],max[1],min[1],avg[1],total[1]/avg[1],
             max[2],min[2],avg[2],total[2]/avg[2],max[3],min[3],avg[3],total[3]/avg[3],
             max[4],min[4],avg[4],total[4]/avg[4]);
      for (int i=0;i<g->numEdgeTypes();i++) {
        char edge_name[15];
        sprintf(edge_name,"Edges type %d",i);
        printf("ENGPAR STATUS: %s: <min,max,avg,imb> %.3f %.3f %.3f %.3f\n",edge_name,
               max[5+i],min[5+i],avg[5+i],total[5+i]/avg[5+i]);
      }
    }
  }
}
