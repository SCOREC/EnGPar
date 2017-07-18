//This file is for type definitions used throughout the graph interface

#ifndef _AGI_H__
#define _AGI_H__

#include <stdint.h>
#include <apfVector.h>
#include <set>
#include "agi_types.h"
//ID types
//TODO: discuss what these should be?
namespace agi {
  
  typedef AGI_GID_T gid_t;
  typedef AGI_LID_T lid_t;
  typedef AGI_WGT_T wgt_t;
  typedef AGI_PART_T part_t;
  typedef double coord_t[3];

  typedef std::set<part_t> Peers;
  
  class GraphVertex;
  typedef std::unordered_map<gid_t,part_t> PartitionMap;
  typedef std::unordered_map<GraphVertex*,part_t> Migration;
  typedef std::unordered_map<lid_t,part_t> VertexPartitionMap;
  typedef std::unordered_map<lid_t,part_t> EdgePartitionMap;

//Definitions for edge types
//static size of edge types
#define MAX_TYPES 7
typedef int etype;
#define VTX_TYPE -1
//predefined edge type for split vtx
#define SPLIT_TYPE MAX_TYPES-1 
}
#endif
