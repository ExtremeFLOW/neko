#ifndef __COMM_H
#define __COMM_H

/**
 * Wrapper for Neko's MPI related data
 */
#include <mpi.h>

/**
 * MPI communicator
 * @note Same name as for the Fortran equivalent
 */
extern MPI_Comm NEKO_COMM;


#endif
