#include <ngraph.h>
#include <engpar_support.h>
#include <iostream>
#include <PCU.h>
#include <vector>
#include "buildGraphs.h"
#include "gatherGraphs.h"
#include "TestingSuite.h"
#include "../partition/Diffusive/engpar_diffusive_input.h"
#include "../partition/Diffusive/src/engpar_queue.h"
#include <dirent.h>
#include <agiMigration.h>
int switchComm();
int tagGraph();
int tagHyperGraph();
int testAeroDQs();

int traverseEdges(agi::Ngraph*);
int testAdjacent(agi::Ngraph* g);
int compareTraversal(agi::Ngraph* g);
int testDistanceQueue(agi::Ngraph* g);
int testMigration(agi::Ngraph* g);
int testRepartition(agi::Ngraph* g);
int testVEVAdjacency(agi::Ngraph* g);
int testEVEAdjacency(agi::Ngraph* g);

int main(int argc, char* argv[]) {
  int trial = -1;
  if (argc>1)
    trial = atoi(argv[1]);
  MPI_Init(&argc,&argv);
  EnGPar_Initialize();
  EnGPar_Set_Verbosity(-1);
  //Create the testing suite
  TestingSuite suite(argv[0]);
    
  //Gather specific tests that have more fine grain checks
  if (PCU_Comm_Peers()>=2)
    suite.addFineTest("Switch Comm", switchComm);
  suite.addFineTest("Tag Graph", tagGraph);
  suite.addFineTest("Tag HyperGraph", tagHyperGraph);
  if (PCU_Comm_Peers()==1)
    suite.addFineTest("Build Aero Distance Queues",testAeroDQs);
  
  //Gather the graphs for the general tests
  gatherBuildGraphs(suite);
  gatherEBINGraphs(suite);
  gatherBGDGraphs(suite);
  
  //Gather general tests that run on the graphs collected prior to this
  suite.addGeneralTest("Traverse Edges",traverseEdges);
  suite.addGeneralTest("Traverse Adjacency",testAdjacent);
  suite.addGeneralTest("Compare Traversals",compareTraversal);
  suite.addGeneralTest("Build Distance Queue", testDistanceQueue);
  if (PCU_Comm_Peers()>1) {
    suite.addGeneralTest("Migration",testMigration);
    suite.addGeneralTest("Repartition", testRepartition);
  }
  suite.addGeneralTest("VEV adjacency", testVEVAdjacency);
  suite.addGeneralTest("EVE adjacency", testEVEAdjacency);
  
  //Run the tests and get the number of failures
  int ierr = suite.runTests(trial);

  suite.deleteTestGraphs();
  EnGPar_Finalize();
  MPI_Finalize();
  return ierr;
}

int switchComm() {
  MPI_Comm new_comm;
  int new_size = PCU_Comm_Peers()/2;
  int old_self = PCU_Comm_Self();
  MPI_Comm_split(MPI_COMM_WORLD,old_self/new_size,
                 old_self%new_size,&new_comm);

  PCU_Switch_Comm(new_comm);
  
  if (old_self<new_size) {
    agi::Ngraph* g = buildGraph();
    agi::VertexIterator* itr = g->begin();
    agi::GraphVertex* vert;
    while ((vert = g->iterate(itr)));
    agi::destroyGraph(g);
  }
  else {
    agi::Ngraph* g = buildHyperGraph();
    agi::EdgeIterator* itr = g->begin(0);
    agi::GraphEdge* edge;
    while ((edge = g->iterate(itr)));
    g->destroy(itr);  
    agi::destroyGraph(g);
  }
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&new_comm);
  return 0;
}

int tagGraph() {
  agi::Ngraph* g = buildGraphParts();

  agi::GraphTag* tag = g->createIntTag();
  agi::VertexIterator* itr= g->begin();
  agi::GraphVertex* vtx;

  int i=0;
  while ((vtx = g->iterate(itr))) {
    g->setIntTag(tag,vtx,i);
    i++;
  }

  itr = g->begin();
  i=0;
  while ((vtx = g->iterate(itr))) {
    if (g->getIntTag(tag,vtx)!=i)
      return 1;
    i++;
  }

  g->destroyTag(tag);

  tag = g->createLongTag(0);
  agi::GraphEdge* edge;
  agi::EdgeIterator* eitr = g->begin(0);
  while ((edge=g->iterate(eitr))) {
    g->setDoubleTag(tag,edge,g->degree(edge));
  }
  g->destroy(eitr);
  eitr = g->begin(0);
  while ((edge=g->iterate(eitr))) {
    if (g->getDoubleTag(tag,edge)!=2)
      return 2;
  }
  g->destroy(eitr);
  g->destroyTag(tag);

  tag = g->createLongTag(1);

  eitr = g->begin(1);
  while ((edge=g->iterate(eitr))) {
    g->setDoubleTag(tag,edge,g->degree(edge));
  }
  g->destroy(eitr);
  eitr = g->begin(1);
  while ((edge=g->iterate(eitr))) {
    if (g->getDoubleTag(tag,edge)!=2)
      return 3;
  }
  g->destroy(eitr);
  g->destroyTag(tag);
  
  agi::destroyGraph(g);
  return 0;
}

int tagHyperGraph() {
  agi::Ngraph* g = buildHyperGraphParts();

  agi::GraphTag* tag = g->createIntTag();
  agi::VertexIterator* itr= g->begin();
  agi::GraphVertex* vtx;

  int i=0;
  while ((vtx = g->iterate(itr))) {
    g->setIntTag(tag,vtx,i);
    i++;
  }

  itr = g->begin();
  i=0;
  while ((vtx = g->iterate(itr))) {
    if (g->getIntTag(tag,vtx)!=i)
      return 1;
    i++;
  }

  g->destroyTag(tag);

  tag = g->createLongTag(0);
  agi::GraphEdge* edge;
  agi::EdgeIterator* eitr = g->begin(0);
  while ((edge=g->iterate(eitr))) {
    g->setDoubleTag(tag,edge,g->degree(edge));
  }
  g->destroy(eitr);
  eitr = g->begin(0);
  while ((edge=g->iterate(eitr))) {
    if (g->getDoubleTag(tag,edge)!=g->degree(edge))
      return 2;
  }
  g->destroy(eitr);
  g->destroyTag(tag);

  agi::destroyGraph(g);
  return 0;
}


int traverseEdges(agi::Ngraph* g) {
  if (g->isHyper())
    return 0;
  agi::EdgeIterator* eitr = g->begin(0);
  agi::GraphEdge* edge =NULL;
  agi::lid_t numEdges=0;
  std::map<agi::GraphVertex*,agi::lid_t> edgesPerVertex;
  while ((edge = g->iterate(eitr))) {
    numEdges++;
    agi::GraphVertex* u = g->u(edge);
    edgesPerVertex[u]++;
    agi::GraphVertex* v = g->v(edge);
    agi::EdgeIterator* eitr2 = g->edges(u);
    agi::GraphEdge* e2 = NULL;
    while ((e2 = g->iterate(eitr2))) {
      if (e2==edge) {
        if (!g->isEqual(v,g->v(e2)))
          return 1;
        if (!g->isEqual(u,g->u(e2)))
          return 2;
      }
    }
    g->destroy(eitr2);
  }
  if (numEdges!=g->numLocalEdges())
    return 3;
  std::map<agi::GraphVertex*,agi::lid_t>::iterator itr;
  for (itr=edgesPerVertex.begin();itr!=edgesPerVertex.end();itr++) {
    if (itr->second!=g->degree(itr->first))
      return 4;
  }
  g->destroy(eitr);
  return 0;
}

int testAdjacent(agi::Ngraph* g) {
  for (agi::etype t = 0;t<g->numEdgeTypes();t++) {
    agi::VertexIterator* vitr = g->begin();
    agi::GraphVertex* vtx = NULL;
    agi::lid_t num_pins =0;
    agi::lid_t num_edges=0;
    while ((vtx = g->iterate(vitr))) {
      agi::GraphIterator* gitr = g->adjacent(vtx,t);
      agi::GraphVertex* other = NULL;
      agi::GraphEdge* edge = NULL;
      while ((other = g->iterate(gitr))) {
        if (!other)
          return 1;
        edge = g->edge(gitr);
        if (!edge)
          return 2;
        if (g->isEqual(vtx,other))
          num_pins++;
        num_edges++;      
      }
      g->destroy(gitr);
    }
    if (g->isHyper()) {
      if (PCU_Comm_Peers()==1)
        if (num_pins!=g->numLocalPins(t))
          return 3;
    }
    else {
      if (num_edges!=g->numLocalEdges(t))
        return 4;
    }
  }
  return 0;
}

int compareTraversal(agi::Ngraph* g) {
  for (agi::etype t = 0;t<g->numEdgeTypes();t++) {
    agi::VertexIterator* vitr = g->begin();
    agi::GraphVertex* vtx = NULL;
    std::vector<agi::GraphEdge*> edges;
    std::vector<agi::GraphVertex*> vtxs;
    while ((vtx = g->iterate(vitr))) {
      agi::GraphIterator* gitr = g->adjacent(vtx,t);
      agi::GraphVertex* other = NULL;
      while ((other = g->iterate(gitr))) {
        edges.push_back(g->edge(gitr));
        vtxs.push_back(other);
      }
      g->destroy(gitr);
    }
    vitr = g->begin();
    agi::lid_t i=0;
    while ((vtx = g->iterate(vitr))) {
      agi::EdgeIterator* eitr = g->edges(vtx,t);
      agi::GraphVertex* other;
      agi::GraphEdge* edge;
      while ((edge = g->iterate(eitr))) {
        if (edges[i]!=edge)
          return 1;
        if (g->isHyper()) {
          agi::PinIterator* pitr = g->pins(edge);
          for (agi::lid_t j=0;j<g->degree(edge);j++) {
            other = g->iterate(pitr);
            if (!g->isEqual(other,vtxs[i]))
              return 5;
            if (edges[i]!=edge)
              return 6;
            i++;
          }
          g->destroy(pitr);
        }
        else {
          if (g->v(edge)!=vtxs[i])
            return 2;
          i++;
        }
      }
      g->destroy(eitr);
    }
    if (i!=(int)vtxs.size())
      return 3;
    if (i!=(int)edges.size())
      return 4;
  }
  return 0;
}

int testAeroDQs() {

  agi::Ngraph* g = agi::createEmptyGraph();
  char aeroDir[1024];
  char* filePlace = aeroDir +sprintf(aeroDir,"%s/aero1Belm/",ENGPAR_GRAPHS);
  
  struct dirent *pDirent;
  DIR *pDir;
  pDir = opendir (aeroDir);
  
  if (pDir == NULL) {
    EnGPar_Error_Message("Cannot open directory '%s'\n", aeroDir);
    return 1;
  }

  while ((pDirent = readdir(pDir)) != NULL) {
    if (pDirent->d_name[0]=='.')
      continue;
    int n = sprintf(filePlace,"%s",pDirent->d_name);
    filePlace[n-6] = '\0';
    g->loadFromFile(aeroDir);
    testDistanceQueue(g);
    agi::destroyGraph(g);
    break;
  }
  closedir (pDir);
  return 0;
}

int testDistanceQueue(agi::Ngraph* g) {
  engpar::Input* input = engpar::createDiffusiveInput(g,0);
  engpar::DiffusiveInput* inp = static_cast<engpar::DiffusiveInput*>(input);
  engpar::Queue* q = engpar::createDistanceQueue(inp);
  delete q;
  delete input;
  return 0;
}

int testMigration(agi::Ngraph* g) {
  agi::Migration* plan = new agi::Migration(g);
  agi::GraphVertex* v;
  agi::VertexIterator* itr=g->begin();
  while ((v = g->iterate(itr))) {
    agi::GraphIterator* gitr = g->adjacent(v);
    agi::GraphVertex* other;
    while ((other=g->iterate(gitr))) {
      int owner = g->owner(other);
      if (owner>PCU_Comm_Self()) {
        plan->insert(std::make_pair(v,owner));
        break;
      }
    }
    g->destroy(gitr);
  }
  g->setOriginalOwners();
  g->migrate(plan);

  agi::checkValidity(g);
  return 0;
}

int testRepartition(agi::Ngraph* g) {
  agi::part_t* partition = new agi::part_t[g->numLocalVtxs()];

  for (int i=0;i<g->numLocalVtxs()/2;i++)
    partition[i] = PCU_Comm_Self();

  for (int i=g->numLocalVtxs()/2;i<g->numLocalVtxs();i++)
    partition[i] = PCU_Comm_Self();
  g->setOriginalOwners();
  g->repartition(partition);
  delete [] partition;

  checkValidity(g);

  return 0;
}

int testVEVAdjacency(agi::Ngraph* g) {
  g->create_vev_adjacency(0,false);
  
  g->create_vev_adjacency(0,true);

}

int testEVEAdjacency(agi::Ngraph* g) {
  g->create_eve_adjacency(0,false);

  g->create_eve_adjacency(0,true);
}



