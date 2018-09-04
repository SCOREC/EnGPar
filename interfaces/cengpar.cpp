#include "cengpar.h"
#include "engpar_support.h"
#include "ngraph.h"
#include "engpar.h"
#include <engpar_split_input.h>
#include <engpar_diffusive_input.h>

void cengpar_initialize() {
  EnGPar_Initialize();
}

void cengpar_finalize() {
  EnGPar_Finalize();
}

void cengpar_setftncommunicator(MPI_Fint fcomm) {
  MPI_Comm comm = MPI_Comm_f2c(fcomm);
  EnGPar_Switch_Comm(comm);
}

ngraph cengpar_createEmptyGraph() {
  agi::Ngraph* ng = agi::createEmptyGraph();
  return (ngraph)ng;
}


engparInput cengpar_createLocalSplitInput(ngraph g, MPI_Fint smallComm, MPI_Fint largeComm,
    bool isOrig, int splitFactor, double tol, agi::etype t, agi::part_t* ranks) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  engpar::Input* input = engpar::createLocalSplitInput(ng,smallComm,largeComm,
                                                       isOrig,splitFactor,tol,ranks,t);
  return (engparInput)input;
}

engparInput cengpar_createGlobalSplitInput(ngraph g, MPI_Fint smallComm, MPI_Fint largeComm,
                                      bool isOrig, double tol, agi::etype t) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  engpar::Input* input = engpar::createGlobalSplitInput(
      ng,smallComm,largeComm,isOrig,tol,t);
  return (engparInput)input;
}

void cengpar_loadFromFile(ngraph g, const char fileName[]) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  ng->loadFromFile(fileName);
}

void cengpar_saveToFile(ngraph g, const char fileName[]) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  ng->saveToFile(fileName);
}

void cengpar_evaluatePartition(ngraph g) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  engpar::evaluatePartition(ng);
}

void cengpar_split(engparInput in, const char method[]) {
  engpar::Input* input = (engpar::Input*)in;
  std::map<std::string,int> methods;
  methods["GLOBAL_PARMETIS"] = engpar::GLOBAL_PARMETIS;
  methods["LOCAL_PARMETIS"] = engpar::LOCAL_PARMETIS;
  assert( methods.count(method) );
  engpar::split(input, (engpar::SPLIT_METHOD)methods[method]);
}

void cengpar_constructVerts(ngraph g, bool isHg,
    agi::gid_t* verts, agi::wgt_t* weights, int nverts) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  std::vector<agi::gid_t> v(verts, verts + nverts);
  std::vector<agi::wgt_t> w(weights, weights + nverts);
  ng->constructVerts(isHg,v,w);
}

agi::etype cengpar_constructEdges(ngraph g, agi::gid_t* edges,
    agi::lid_t* degs, agi::wgt_t* weights,
    agi::gid_t* pins, int nedges, int npins) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  std::vector<agi::gid_t> e(edges, edges + nedges);
  std::vector<agi::lid_t> d(degs, degs + nedges);
  std::vector<agi::wgt_t> w(weights, weights + nedges);
  std::vector<agi::gid_t> p(pins, pins + npins);
  return ng->constructEdges(e,d,p,w);
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

void cengpar_getPartition(ngraph g, agi::gid_t* verts,
    agi::part_t* parts, int nverts) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  assert(ng);
  agi::PartitionMap* ptn = ng->getPartition();
  assert( (int)ptn->size() == nverts );
  agi::PartitionMap::iterator itr;
  int i=0;
  for (itr=ptn->begin();itr!=ptn->end();itr++) {
    verts[i] = itr->first;
    parts[i] = itr->second;
    i++;
  }
  delete ptn;
}

engparInput cengpar_createDiffusiveInput(ngraph g, double stepfactor) {
  agi::Ngraph* ng = (agi::Ngraph*)g;
  engpar::Input* input = engpar::createDiffusiveInput(ng,stepfactor);
  return (engparInput)input;
}

void cengpar_addPriority(engparInput in, agi::etype t, double tol) {
  engpar::Input* input = (engpar::Input*)in;
  input->addPriority(t,tol);
}

void cengpar_balanceGhosts(engparInput in, bool enable) {
  engpar::Input* input = (engpar::Input*)in;
  input->balanceGhosts(enable);
}

void cengpar_balance(engparInput in, int verbosity) {
  engpar::Input* input = (engpar::Input*)in;
  engpar::balance(input,verbosity);
}
