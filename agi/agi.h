//This file is for type definitions used throughout the graph iterface

#ifndef _AGI_H__
#define _AGI_H__

#include <stdint.h>
//ID types
//TODO: discuss what these should be?
namespace agi {
  
typedef uint64_t gid_t;
typedef uint64_t lid_t;
typedef double wgt_t;
typedef int32_t part_t;
 
//Definitions for edge types
#define MAX_TYPES 10 //static size of edge types
typedef int etype; 
#define SPLIT_TYPE 9 //predefined edge type for split vtx

}
#endif
