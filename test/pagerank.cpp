#include <binGraph.h>
#include <mpi.h>
#include <PCU.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <ZoltanCutVertex.h>

typedef agi::lid_t lid_t;
void partitionInfo(agi::binGraph* g) {
  lid_t total_verts = g->numTotalVtxs();
  lid_t global_verts = g->numGlobalVtxs();
  lid_t edges = g->numLocalEdges()+g->numLocalEdges(SPLIT_TYPE);
  lid_t min = PCU_Min_Int(total_verts);
  lid_t max = PCU_Max_Int(total_verts);
  long tot = PCU_Add_Long(total_verts);
  double avg = ((double)tot)/PCU_Comm_Peers();
  double imb = max/avg;
  double inc = ((double)(tot-global_verts))/global_verts*100;
  if (!PCU_Comm_Self()) 
    printf("Vertices: Min %d Max %d Tot %ld Inc %1.4f Avg %1.4f Imb %1.3f\n",min,max,tot,inc,avg,imb);
  min = PCU_Min_Int(edges);
  max = PCU_Max_Int(edges);
  tot = PCU_Add_Long(edges);
  avg = ((double)tot)/PCU_Comm_Peers();
  imb = max/avg;
  if (!PCU_Comm_Self()) 
    printf("Edges: Min %d Max %d Tot %ld Avg %1.4f Imb %1.3f\n",min,max,tot,avg,imb);
  
  lid_t edge_cut =0;
  agi::EdgeIterator* eitr = g->begin(0);
  for (int i=0;i<g->numLocalEdges();i++) {
    agi::GraphEdge* e = g->iterate(eitr);
    if (g->owner(g->v(e))!=PCU_Comm_Self())
      edge_cut++;
  }
  edge_cut+=g->numLocalEdges(SPLIT_TYPE);
  edge_cut = PCU_Add_Long(edge_cut);
  if (!PCU_Comm_Self())
    printf("Edge Cut: %d\n",edge_cut);
}

#define NUM_ITERS 5
#define DAMPING_FACTOR 0.85

int pagerank_verify(agi::binGraph* g,double* pageranks,std::unordered_map<gid_t,int>& extra_edges) {
  double* global_pr = new double[g->numGlobalVtxs()];
  for (gid_t i =0;i<g->numGlobalVtxs();++i)
    global_pr[i]=-1.0;
  agi::GraphVertex* vtx;
  agi::VertexIterator* itr = g->begin();
  while ((vtx = g->iterate(itr))) {
    lid_t deg = g->degree(vtx)+extra_edges[g->globalID(vtx)];
    if (deg!=0)
      global_pr[g->globalID(vtx)] = pageranks[g->localID(vtx)]*deg;
    else
      global_pr[g->globalID(vtx)] = pageranks[g->localID(vtx)]*g->numGlobalVtxs();
  }
  PCU_Max_Doubles(global_pr,g->numGlobalVtxs());

  if (!PCU_Comm_Self()) {
    double pr_sum=0.0;
    for (gid_t i=0;i<g->numGlobalVtxs();i++)
      pr_sum+=global_pr[i];
    printf("Pageranks sum %9.6lf\n",pr_sum);
  }
  return 0;
}
void do_pagerank(agi::binGraph* g) {
  double* pageranks = new double[g->numTotalVtxs()];
  double* pageranks_next = new double[g->numTotalVtxs()];
  double sum_noouts =0.0;
  double sum_noouts_next =0.0;
  agi::EdgeIterator* eitr = g->begin(SPLIT_TYPE);
  agi::GraphEdge* edge;
  std::unordered_map<gid_t,int> extra_edges;
  PCU_Comm_Begin();
  for (lid_t i =0;i<g->numLocalEdges(SPLIT_TYPE);++i) {
    edge = g->iterate(eitr);
    agi::GraphVertex* vtx = g->v(edge);
    agi::GraphVertex* mine = g->find(vtx);
    gid_t gid = g->globalID(mine);
    lid_t deg = g->degree(mine);
    assert(g->globalID(mine)==g->globalID(vtx));
    PCU_COMM_PACK(g->owner(vtx),gid);
    PCU_COMM_PACK(g->owner(vtx),deg);
  }
  g->destroy(eitr);
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    gid_t gid;
    PCU_COMM_UNPACK(gid);
    lid_t deg;
    PCU_COMM_UNPACK(deg);
    extra_edges[gid]=deg;
  }
  agi::GraphVertex* vtx;
  agi::VertexIterator* itr = g->begin();
  while ((vtx = g->iterate(itr))) {
    lid_t lid = g->localID(vtx);
    gid_t gid = g->globalID(vtx);
    lid_t deg = g->degree(vtx)+extra_edges[gid];
    pageranks[lid] = 1.0/g->numGlobalVtxs();
    if (deg>0)
      pageranks[lid] /= deg;
    else {
      pageranks[lid] /= g->numGlobalVtxs();
      sum_noouts_next +=pageranks[lid];
    }
  }
  for (lid_t i = g->numLocalVtxs();i<g->numTotalVtxs();i++) 
    pageranks[i] = 1.0/g->numGlobalVtxs()/g->numGlobalVtxs();
  for (lid_t i = 0;i<g->numTotalVtxs();i++)
    pageranks_next[i] = pageranks[i];

  sum_noouts = PCU_Add_Double(sum_noouts_next);
  sum_noouts_next=0.0;

  for (int iter =0;iter<NUM_ITERS;++iter) {
    for (lid_t i =0;i<g->numLocalVtxs();i++)
      pageranks_next[i]=0.0;
    itr = g->begin();
    while ((vtx = g->iterate(itr))) {
      eitr = g->edges(vtx);
      lid_t lid = g->localID(vtx);
      lid_t deg = g->degree(vtx);
      for (lid_t i =0;i<deg;i++) {
        edge = g->iterate(eitr);
        pageranks_next[g->localID(g->v(edge))]+=pageranks[lid];
      }
      g->destroy(eitr);
    }
    itr = g->begin();
    while ((vtx = g->iterate(itr))) {
      lid_t lid = g->localID(vtx);
      pageranks_next[lid] += sum_noouts/(double)g->numGlobalVtxs();
      pageranks_next[lid]*=DAMPING_FACTOR;
      pageranks_next[lid]+=(1.0-DAMPING_FACTOR)/g->numGlobalVtxs();
      agi::gid_t gid = g->globalID(vtx);
      lid_t deg = g->degree(vtx) + extra_edges[gid];
      if (deg>0)
        pageranks_next[lid] /= deg;
      else {
        pageranks_next[lid] /= g->numGlobalVtxs();
        sum_noouts_next +=pageranks[lid];
      }
    }
    eitr=g->begin(SPLIT_TYPE);
    PCU_Comm_Begin();
    for (lid_t i =0;i<g->numLocalEdges(SPLIT_TYPE);++i) {
      edge = g->iterate(eitr);
      agi::GraphVertex* vtx = g->v(edge);
      agi::GraphVertex* mine = g->find(vtx);
      assert(g->globalID(mine)==g->globalID(vtx));
      gid_t gid = g->globalID(mine);
      
      PCU_COMM_PACK(g->owner(vtx),gid);
      PCU_COMM_PACK(g->owner(vtx),pageranks[g->localID(mine)]);
    }
    g->destroy(eitr);
    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      gid_t gid;
      PCU_COMM_UNPACK(gid);
      double pr;
      PCU_COMM_UNPACK(pr);
      pageranks_next[g->localID(g->findGID(gid))]+=pr;
    }
    sum_noouts = PCU_Add_Double(sum_noouts_next);
    sum_noouts_next=0.0;
    double* temp = pageranks;
    pageranks=pageranks_next;
    pageranks_next = temp;
  }
  pagerank_verify(g,pageranks,extra_edges);
}

int main(int argc,char* argv[]) {

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if (argc!=2 &&argc!=3) {
    if (!PCU_Comm_Self())
      printf("Usage: %s <binary_graph_file> [vertex_partition_file]/cut",argv[0]);
    PCU_Comm_Free();
    MPI_Finalize();
    assert(false);
  }
  agi::binGraph* g;
  if (argc==2)
    g = new agi::binGraph(argv[1]);
  else if (std::string(argv[2])=="cut") {
    int self = PCU_Comm_Self();
    int num_parts = PCU_Comm_Peers();
    PCU_Switch_Comm(MPI_COMM_SELF);
    zagi::ZoltanCutVertex* ptn;
    if (!self) {
      g = new agi::binGraph(argv[1]);
      ptn = new zagi::ZoltanCutVertex(g,num_parts);
      ptn->run();
    }
    else
      g = new agi::binGraph();
    agi::EdgePartitionMap map;
    if (!self) {
      ptn->createPtn(map);
      delete ptn;
    }
    PCU_Switch_Comm(MPI_COMM_WORLD);
    g->migrate(map);

  }
  else 
    g=new agi::binGraph(argv[1],argv[2]);

  partitionInfo(g);


  //Run Pagerank
  do_pagerank(g);

  delete g;


  PCU_Barrier();
  if (!PCU_Comm_Self())
    printf("\nAll tests passed\n");
  
  PCU_Comm_Free();
  MPI_Finalize();
}
