#include "cengpar.h"
#include "engpar_support.h"
#include "ngraph.h"
#include "engpar.h"

void cengpar_initialize() {
  EnGPar_Initialize();
}

void cengpar_finalize() {
  EnGPar_Finalize();
}

void cengpar_setftncommunicator(MPI_Fint fcomm) {
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  PCU_Switch_Comm(comm);
}

ngraph cengpar_createEmptyGraph() {
  agi::Ngraph* ng = agi::createEmptyGraph();
  return (ngraph)ng;
}

void cengpar_constructVerts(ngraph g, bool isHg,
    agi::gid_t* verts, agi::wgt_t* weights, int nverts) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  std::vector<agi::gid_t> v(verts, verts + nverts);
  std::vector<agi::wgt_t> w(weights, weights + nverts);
  ng->constructVerts(isHg,v,w);
}

void cengpar_constructEdges(ngraph g, agi::gid_t* edges,
    agi::lid_t* degs, agi::gid_t* pins, int nedges, int npins) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  std::vector<agi::gid_t> e(edges, edges + nedges);
  std::vector<agi::lid_t> d(degs, degs + nedges);
  std::vector<agi::gid_t> p(pins, pins + npins);
  ng->constructEdges(e,d,p);
}

void cengpar_constructGhosts(ngraph g, agi::gid_t* verts, agi::part_t* owners,
    int nghosts) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  std::unordered_map<agi::gid_t,agi::part_t> ghosts;
  for(int i=0; i<nghosts; i++)
    ghosts[verts[i]] = owners[i];
  ng->constructGhosts(ghosts);
}

void cengpar_checkValidity(ngraph g) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  agi::checkValidity(ng);
}

void cengpar_destroyGraph(ngraph g) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  agi::destroyGraph(ng);
}

void cengpar_balanceVertices(ngraph g, double tol, double stepfactor, int verbosity) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  agi::Balancer* b = engpar::makeVtxBalancer(ng, stepfactor, verbosity);
  b->balance(tol);
}
