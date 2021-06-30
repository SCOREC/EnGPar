#include <apfGraph.h>
#include <apfMesh2.h>
#include <cassert>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apf.h>
#include <cstdlib>
#include <stdint.h>
#include <apfNumbering.h>
#include <cstring>
#include <engpar_support.h>
#include <engpar.h>
#include <engpar_input.h>
#include <binGraph.h>
#include <parma.h>
#include <engpar_diffusive_input.h>

void applyDOFWeights(apf::Mesh2* mesh, agi::Ngraph* graph, int order);
void printImbalance(agi::Ngraph* graph);

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  EnGPar_Initialize();
  if (argc != 4 && argc != 5) {
    if (!PCU_Comm_Self()) {
      fprintf(stderr, "Usage: %s <model> <mesh> <0=separate, 1=single edge type> [dof order]\n",
              argv[0]);
    }
    EnGPar_Finalize();
    PCU_Comm_Free();
    MPI_Finalize();
    assert(false);
  }

  //Load mesh
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);


  int method=atoi(argv[3]);
  int order = 1;
  if(argc >= 5)
    order = atoi(argv[4]);

  agi::Ngraph* g = NULL;
  if (method == 0) {
    if (!PCU_Comm_Self())
      fprintf(stderr, "Creating N-graph with 2 hyperedge types, one for mesh vertex, one for mesh edges\n");
    //Create graph via built in apf->ngraph function
    int types[2] = {0,1};
    g = agi::createAPFGraph(m, "twoedgestypes", 3, types, 2);
    //Apply weights to each type
    applyDOFWeights(m, g, order);
  }
  else if (method == 1) {
    if (!PCU_Comm_Self())
      fprintf(stderr, "Creating N-graph with 1 hyperedge type, one for dof holders\n");
    //Create ngraph with 1 hyperedge type for dof holders and apply weights directly
    g = agi::createSerendipityGraph(m, "dofgraph", order);
  }

  //Check to make sure graph is constructed
  if (g == NULL) {
    if (!PCU_Comm_Self())
      fprintf(stderr, "[ERROR] No graph was constructed\n");
    return 1;
  }

  //Print Imbalances
  printImbalance(g);

  //Run EnGPar
  double step_factor=0.2;
  double tol = 1.05;
  engpar::DiffusiveInput* input = engpar::createDiffusiveInput(g, step_factor);
  //Hyperedge types have first priority
  for (agi::etype t = 0; t< g->numEdgeTypes(); ++t) {
    input->addPriority(t, tol);
  }
  //Graph vertices have last priority
  input->addPriority(-1, tol);

  int verbosity = 1;
  engpar::balance(input, verbosity);

  printImbalance(g);

  agi::destroyGraph(g);

  m->destroyNative();
  apf::destroyMesh(m);

  EnGPar_Finalize();
  PCU_Comm_Free();
  MPI_Finalize();

  return 0;

}

void applyDOFWeights(apf::Mesh2* mesh, agi::Ngraph* graph, int order) {
  graph->setEdgeWeights(1, 0);
  graph->setEdgeWeights(order - 1, 1);
}

void printImbalance(agi::Ngraph* graph) {
  engpar::evaluatePartition(graph);
  agi::wgt_t dofs = 0;
  for (int t = 0; t < graph->numEdgeTypes(); ++t) {
    agi::EdgeIterator* eitr = graph->begin(t);
    agi::GraphEdge* edge;
    while ((edge = graph->iterate(eitr)) != NULL) {
      dofs += graph->weight(edge);
    }
    graph->destroy(eitr);
  }
  agi::wgt_t avg, max, min, imb;
  max = PCU_Max_Double(dofs);
  min = PCU_Min_Double(dofs);
  avg = PCU_Add_Double(dofs) / PCU_Comm_Peers();
  imb = max / avg;
  if (!PCU_Comm_Self())
    fprintf(stderr, "DOF <max, min, avg, imb>: %f %f %f %f\n", max, min, avg, imb);
}
