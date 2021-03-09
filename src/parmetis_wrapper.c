/*

  Fortran wrappers for ParMETIS (from CUBE)

  Removing the issue of figuring out the data type of the MPI
  communicator from within Fortran

*/

#include <mpi.h>
#include <parmetis.h>


/*!
  Fortran wrapper for ParMETIS PartGeom
*/
int ParMETIS_V3_PartGeom_wrapper(idx_t *vtxdist, idx_t *ndims,
				 real_t *xyz, idx_t *part)
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

/*!
  Fortran wrapper for ParMETIS PartMeshKway
*/
int ParMETIS_V3_PartMeshKway_wrapper(idx_t *elmdist, idx_t *eptr, idx_t *eind,
				     idx_t *elmwgt, idx_t *wgtflag,
				     idx_t *numflag, idx_t *ncon,
				     idx_t *ncommonnodes, idx_t *nparts,
				     real_t *tpwgts, real_t *ubvec,
				     idx_t *options, idx_t *edgecut, idx_t *part)
{

  int rcode;
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  rcode = METIS_ERROR;

#ifdef HAVE_PARMETIS
  
  rcode = ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt, wgtflag,
				   numflag, ncon, ncommonnodes, nparts,tpwgts,
				   ubvec, options, edgecut, part, &comm);
#endif

  return rcode;

}
