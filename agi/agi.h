//This file is for type definitions used throughout the graph iterface

#ifndef _AGI_H__
#define _AGI_H__

#include <stdint.h>
#include <apfVector.h>
//ID types
//TODO: discuss what these should be?
namespace agi {
  
  typedef uint64_t gid_t;
  typedef uint64_t lid_t;
  typedef double wgt_t;
  typedef int32_t part_t;
  typedef apf::Vector3 coord_t;

  class GraphVertex;
  typedef std::unordered_map<GraphVertex*,part_t> Migration;
  typedef std::unordered_map<lid_t,part_t> VertexPartitionMap;
  typedef std::unordered_map<lid_t,part_t> EdgePartitionMap;

//Definitions for edge types
//static size of edge types
#define MAX_TYPES 10 
typedef int etype; 
//predefined edge type for split vtx
#define SPLIT_TYPE 9 

}
#endif
