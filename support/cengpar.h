#ifndef CENGPAR
#define CENGPAR

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

void cengpar_initialize();
void cengpar_finalize();
void cengpar_setftncommunicator(MPI_Fint fcomm);

#ifdef __cplusplus
}
#endif

#endif
