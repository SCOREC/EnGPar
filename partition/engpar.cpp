#include "engpar.h"
#include <PCU.h>
#include "Diffusive/src/engpar_sides.h"
#include <engpar_support.h>
#include <engpar_metrics.h>

namespace engpar {

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
    int num_vals = 4 + g->numEdgeTypes() * 2;
    agi::wgt_t* my_vals = new agi::wgt_t[num_vals];
    agi::wgt_t* min = new agi::wgt_t[num_vals];
    agi::wgt_t* max = new agi::wgt_t[num_vals];
    agi::wgt_t* total = new agi::wgt_t[num_vals];
    agi::wgt_t* avg = new agi::wgt_t[num_vals];

    //Disconnected comp
    my_vals[0] = countDisconnected(g);

    //Neighbors
    DiffusiveInput* in = static_cast<DiffusiveInput*>(createDiffusiveInput(g,0));
    Sides* sides = makeSides(in);
    my_vals[1] = sides->size();
    delete in;
    delete sides;
    
    //Vertex Imbalance
    my_vals[2] = getWeight(g,-1);
    my_vals[3] = getWeight(g,-1,true);

    for (agi::etype t = 0;t<g->numEdgeTypes();t++) {
      //Edge Imbalance
      my_vals[2*t+4] = getWeight(g,t);

      //Edge Cut
      DiffusiveInput* in = static_cast<DiffusiveInput*>(createDiffusiveInput(g,0));
      in->sides_edge_type = t;
      Sides* sides = makeSides(in);
      my_vals[2*t+5] = sides->total();
      delete in;
      delete sides;

    }


    //Empty parts
    int empty = my_vals[3]==0;

    MPI_Reduce(my_vals,max,num_vals,MPI_DOUBLE,MPI_MAX,0,PCU_Get_Comm());
    MPI_Reduce(my_vals,min,num_vals,MPI_DOUBLE,MPI_MIN,0,PCU_Get_Comm());
    MPI_Reduce(my_vals,total,num_vals,MPI_DOUBLE,MPI_SUM,0,PCU_Get_Comm());
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
      EnGPar_Status_Message("%s: Local Vertex: <max,min,avg,imb> %.3f %.3f %.3f %.3f\n",
                            prefix,max[2],min[2],avg[2],max[2]/avg[2]);
      EnGPar_Status_Message("%s: Total Vertex: <max,min,avg,imb> %.3f %.3f %.3f %.3f\n",
                            prefix,max[3],min[3],avg[3],max[3]/avg[3]);
      for (int i=0;i<g->numEdgeTypes();i++) {
        char edge_name[30];
        sprintf(edge_name,"%s: Edges type %d",prefix,i);
        EnGPar_Status_Message("%s: <max,min,avg,imb> %.3f %.3f %.3f %.3f\n",edge_name,
                              max[4+i*2],min[4+i*2],avg[4+i*2],max[4+i*2]/avg[4+i*2]);
        EnGPar_Status_Message("%s Cut: <max,tot> %.3f %.3f\n", edge_name,
                              max[5+i*2], total[5+i*2]);

      }
    }
    
    delete [] max;
    delete [] min;
    delete [] total;
    delete [] avg;
  }
}
