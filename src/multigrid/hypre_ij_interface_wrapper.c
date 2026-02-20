/*
 *  C wrappers for IJ interface
 *  Removin the issue of figureing out the data type of the MPI
 *  communicator from within Fortran
 *
 */

#include <mpi.h>
#include <HYPRE.h>
#include "HYPRE_parcsr_ls.h"

#include <comm/comm.h>

int HYPRE_IJMatrixCreate_wrapper(int ilower, int iupper, int jlower, int jupper, HYPRE_IJMatrix *matrix)
{
  int rcode;

  rcode = HYPRE_IJMatrixCreate(NEKO_COMM, ilower, iupper, jlower, jupper, matrix);

  return rcode;
}

int HYPRE_IJVectorCreate_wrapper(int jlower, int jupper, HYPRE_IJVector *vector)
{
  int rcode;

  rcode = HYPRE_IJVectorCreate(NEKO_COMM, jlower, jupper, vector);

  return rcode;
}

// TODO: Sticking this here for now. Move to an appropriate place.
int HYPRE_init_wrapper()
{
  int rcode;
  rcode = HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
  rcode = HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
  HYPRE_SetSpGemmUseVendor(0);
  return rcode;
}
