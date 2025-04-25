/**
 * Wrapper for Neko's MPI related data
 */
#include <comm/comm.h>

/**
 * MPI communicator
 * @note Same name as for the Fortran equivalent
 */
MPI_Comm NEKO_COMM;

/**
 * Setup Neko's communicator for C/C++
 */
void neko_comm_wrapper_init(MPI_Fint fcomm) {
  NEKO_COMM = MPI_Comm_f2c(fcomm);
}
