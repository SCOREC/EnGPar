#include "engpar.h"
#include <PCU.h>
#include "Diffusive/src/engpar_sides.h"
#include <engpar_support.h>

namespace engpar {
  wgt_t getWeight(agi::Ngraph* g, int dimension, bool countGhosts) {
    wgt_t w =0;
    if (dimension==-1) {
      
      //calculate the total weight of the vertices
      agi::GraphVertex* vtx;
      agi::VertexIterator* vitr = g->begin();
      if (!countGhosts)
        while ((vtx = g->iterate(vitr))) 
          w+=g->weight(vtx);
      else
	w=g->numTotalVtxs();
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
    char imbalances[100];
    int n=0;
    for (agi::etype type = -1;type<g->numEdgeTypes();type++) {
      agi::wgt_t w = getWeight(g,type);
      double imb = EnGPar_Get_Imbalance(w);
      if (!PCU_Comm_Self())
        n+=sprintf(imbalances+n, "%2.4f ", imb);
    }
    if (!PCU_Comm_Self())
      printf("%s\n",imbalances);
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
      EnGPar_Status_Message("%s <max, min, avg, imb> = <%lu %lu %f %f>\n",
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
      EnGPar_Status_Message("%s <max, min, avg, imb> = <%f %f %f %f>\n",
             prefix.c_str(),max,min,avg,imb);
  }

  void bfs(agi::Ngraph* g, agi::GraphVertex* v, agi::GraphTag* visited) {
    agi::GraphVertex** queue = new agi::GraphVertex*[g->numLocalVtxs()];
    queue[0] = v;
    g->setIntTag(visited,v,1);
    agi::lid_t size=1;
    agi::lid_t i=0;
    while (i<size) {
      agi::GraphVertex* u = queue[i++];
      agi::GraphVertex* other;
      agi::GraphIterator* gitr = g->adjacent(u);
      while ((other = g->iterate(gitr))) {
        if (g->owner(other)!=PCU_Comm_Self())
          continue;
        if (g->getIntTag(visited,other)==0) {
          queue[size++] = other;
          g->setIntTag(visited,other,1);
        }
      }
      g->destroy(gitr);
    }
    delete [] queue;
  }

  int countDisconnected(agi::Ngraph* g) {
    int dis = 0;

    agi::GraphTag* visited = g->createIntTag();
    agi::VertexIterator* vitr = g->begin();
    agi::GraphVertex* vtx;
    while ((vtx = g->iterate(vitr))) {
      g->setIntTag(visited,vtx,0);
    }

    vitr = g->begin();
    while ((vtx = g->iterate(vitr))) {
      if (g->getIntTag(visited,vtx)==0) {
        bfs(g,vtx,visited);
        dis++;
      }
    }

    g->destroyTag(visited);
    return dis-1;
  }

  
  void evaluatePartition(agi::Ngraph* g,const char* name) {
    agi::wgt_t* my_vals = new agi::wgt_t[2+g->numEdgeTypes()+3];
    agi::wgt_t* min = new agi::wgt_t[2+g->numEdgeTypes()+3];
    agi::wgt_t* max = new agi::wgt_t[2+g->numEdgeTypes()+3];
    agi::wgt_t* total = new agi::wgt_t[2+g->numEdgeTypes()+3];
    agi::wgt_t* avg = new agi::wgt_t[2+g->numEdgeTypes()+3];

    //Disconnected comp
    my_vals[0] = countDisconnected(g);

    //Neighbors & Part Boundary
    DiffusiveInput* in = static_cast<DiffusiveInput*>(createDiffusiveInput(g,0));
    Sides* sides = makeSides(in);
    my_vals[1] = sides->size();
    my_vals[2] = sides->total();
    delete in;
    delete sides;
    
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
    delete [] my_vals;
    if (!PCU_Comm_Self()) {
      for (int i=0;i<5+g->numEdgeTypes();i++) {
        avg[i] = total[i]/PCU_Comm_Peers();
      }
      char prefix[100];
      sprintf(prefix,"PARTITION %s",name);
      EnGPar_Status_Message("%s: Empty Parts: %d\n",prefix,empty);
      EnGPar_Status_Message("%s: Disconnected Components: <max,tot> %.3f %.3f\n",
                            prefix,max[0],total[0]);
      EnGPar_Status_Message("%s: Neighbors: <max,min,avg,imb> %.3f %.3f %.3f %.3f\n",
                            prefix,max[1],min[1],avg[1],max[1]/avg[1]);
      EnGPar_Status_Message("%s: Edge Cut: <max,min,avg,imb> %.3f %.3f %.3f %.3f\n",
                            prefix,max[2],min[2],avg[2],max[2]/avg[2]);
      EnGPar_Status_Message("%s: Local Vertex Imbalance: <max,min,avg,imb> %.3f %.3f %.3f %.3f\n",
                            prefix,max[3],min[3],avg[3],max[3]/avg[3]);
      EnGPar_Status_Message("%s: Total Vertex Imbalance: <max,min,avg,imb> %.3f %.3f %.3f %.3f\n",
                            prefix,max[4],min[4],avg[4],max[4]/avg[4]);
      for (int i=0;i<g->numEdgeTypes();i++) {
        char edge_name[30];
        sprintf(edge_name,"%s: Edges type %d",prefix,i);
        EnGPar_Status_Message("%s: <max,min,avg,imb> %.3f %.3f %.3f %.3f\n",edge_name,
               max[5+i],min[5+i],avg[5+i],max[5+i]/avg[5+i]);
      }
    }
    
  delete [] max;
  delete [] min;
  delete [] total;
  delete [] avg;
  }
}
