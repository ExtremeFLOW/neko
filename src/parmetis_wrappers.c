/*

  Fortran wrappers for ParMETIS (from CUBE)

  Removing the issue of figuring out the data type of the MPI
  communicator from within Fortran

*/

#include <mpi.h>
#include <parmetis.h>

#ifdef HAVE_PARMETIS
#if REALTYPEWIDTH == 64
#error "ParMETIS built with 64bit support"
#endif
#endif

/*!
  Fortran wrapper for ParMETIS PartGeom
*/
int ParMETIS_V3_PartGeom_wrapper(int *vtxdist, int *ndims,
				 float *xyz, int *part)
{
  int rcode;
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  rcode = METIS_ERROR;

#ifdef HAVE_PARMETIS

  rcode = ParMETIS_V3_PartGeom(vtxdist, ndims, xyz, part, &comm);

#endif

  return rcode;
  
}
