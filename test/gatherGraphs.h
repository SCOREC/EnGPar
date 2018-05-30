#ifndef GATHER_GRAPHS_H__
#define GATHER_GRAPHS_H__

#include <vector>
#include <ngraph.h>
void gatherBuildGraphs(std::vector<agi::Ngraph*>&);

void gatherEBINGraphs(std::vector<agi::Ngraph*>&);

void gatherBGDGraphs(std::vector<agi::Ngraph*>&);

#endif
