#include "gatherGraphs.h"
#include "buildGraphs.h"
#include <PCU.h>
#include <binGraph.h>


void gatherBuildGraphs(std::vector<agi::Ngraph*>& graphs) {
  graphs.push_back(buildGraph());
  graphs.push_back(buildHyperGraph());
  graphs.push_back(buildGraphParts());
  graphs.push_back(buildHyperGraphParts());
  graphs.push_back(buildRequirementsGraph());
  graphs.push_back(buildEmptyGraph());
  if (PCU_Comm_Peers()==2) {
    graphs.push_back(buildDisconnected2Graph());
    graphs.push_back(buildUnbalancedLineHG());
  }
}

void gatherEBINGraphs(std::vector<agi::Ngraph*>& graphs) {
  char file[1024];
  //Ring graph
  sprintf(file,"%s/ring.ebin",ENGPAR_GRAPHS);
  graphs.push_back(agi::createBinGraph(file));

  //Tree graph
  sprintf(file,"%s/tree.ebin",ENGPAR_GRAPHS);
  graphs.push_back(agi::createBinGraph(file));

  //Gnutella graph
  sprintf(file,"%s/gnutella.ebin",ENGPAR_GRAPHS);
  graphs.push_back(agi::createBinGraph(file));
}

void gatherBGDGraphs(std::vector<agi::Ngraph*>& graphs) {
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
    graphs.push_back(g);
  }

  if (PCU_Comm_Peers() == 4) {
    sprintf(file,"%s/torus/4/",ENGPAR_GRAPHS);
    agi::Ngraph* g = agi::createEmptyGraph();
    g->loadFromFile(file);
    graphs.push_back(g);
    sprintf(file,"%s/torus/4_01/",ENGPAR_GRAPHS);
    g = agi::createEmptyGraph();
    g->loadFromFile(file);
    graphs.push_back(g);
  }
}

