# SCOREC EnGPar #

EnGPar is a set of C/C++ libraries for partitioning relational data structures using a parallel graph structure to represent the data.

### Repository Layout ###

* AGI: APIs for the graph structure, Ngraph.
* interfaces: Extensions of the Ngraph for certain data structures.
* zoltan: Partitioning methods utilizing zoltan
* partition: Diffusive load balancing routines
* test: A set of tests to utilize the features of EnGPar

### Getting started ###

* Dependencies: CMake, MPI, PCU ([can be grabbed from SCOREC Core](https://github.com/SCOREC/core))
* Configuration: two example configuration files are provided in the EnGPar directory minimal_config.sh and config.sh. The first has the minimum requirements to setup a directory while the second lists the specific flags to turn on and off the various portions of EnGPar.
* Testing: After building the code, the test directory will have several tests and standalone executables that use the different functionalities of EnGPar.

### Contact Developers ###
* Gerrett Diamond <diamog@rpi.edu>