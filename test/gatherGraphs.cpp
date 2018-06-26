#include "gatherGraphs.h"
#include "buildGraphs.h"
#include <PCU.h>
#include <binGraph.h>


void gatherBuildGraphs(TestingSuite& suite) {
  suite.addTestGraph("Built Graph" , buildGraph());
  suite.addTestGraph("Built HyperGraph", buildHyperGraph());
  suite.addTestGraph("Built Graph Parts", buildGraphParts());
  suite.addTestGraph("Built HyperGraph Parts", buildHyperGraphParts());
  if (PCU_Comm_Peers()==1)
    suite.addTestGraph("Built Edge Cases", buildRequirementsGraph());
  suite.addTestGraph("Built Empty Graph", buildEmptyGraph());
  if (PCU_Comm_Peers()==2) {
    suite.addTestGraph("Built Disconnected Graph" , buildDisconnected2Graph());
    suite.addTestGraph("Built Unbalanced Line HyperGraph", buildUnbalancedLineHG());
  }
}

void gatherEBINGraphs(TestingSuite& suite) {
  char file[1024];
  #ifndef ENGPAR_BIG_ENDIAN
  if (PCU_Comm_Peers() <= 2) {
    //Ring graph
    sprintf(file,"%s/ring.ebin",ENGPAR_GRAPHS);
    suite.addTestGraph("Ring Binary Graph", agi::createBinGraph(file));

    //Tree graph
    sprintf(file,"%s/tree.ebin",ENGPAR_GRAPHS);
    suite.addTestGraph("Tree Binary Graph", agi::createBinGraph(file));
  }
  //Gnutella graph
  sprintf(file,"%s/gnutella.ebin",ENGPAR_GRAPHS);
  suite.addTestGraph("Gnutella Binary Graph", agi::createBinGraph(file));
  #endif
}

void gatherBGDGraphs(TestingSuite& suite) {
  //Cube tests (comm size = 1 or 4)
  char file[1024];
  file[0] = '\0';
  if (PCU_Comm_Peers() == 1) {
    sprintf(file,"%s/cube/cube",ENGPAR_GRAPHS);
  }
  else if (PCU_Comm_Peers()==4) {
    sprintf(file,"%s/cube/4/",ENGPAR_GRAPHS);
  }
  if (file[0]!='\0') {
    agi::Ngraph* g = agi::createEmptyGraph();
    g->loadFromFile(file);
    suite.addTestGraph("Cube elm-vtx Mesh Graph", g);
  }

  if (PCU_Comm_Peers() == 4) {
    sprintf(file,"%s/torus/4/",ENGPAR_GRAPHS);
    agi::Ngraph* g = agi::createEmptyGraph();
    g->loadFromFile(file);
    suite.addTestGraph("Torus elm-vtx Mesh Graph", g);
    sprintf(file,"%s/torus/4_01/",ENGPAR_GRAPHS);
    g = agi::createEmptyGraph();
    g->loadFromFile(file);
    suite.addTestGraph("Torus elm-vtx/edge Mesh Graph", g);
  }
}

