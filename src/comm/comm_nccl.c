
#include <comm/comm_nccl.h>
#include <comm/comm.h>

#if defined(HAVE_NCCL) || defined(HAVE_RCCL)

/**
 * NCCL communicator
 */
ncclComm_t NEKO_COMM_NCCL;

/** 
 * Setup Neko's NCCL communicator
 */
void neko_comm_nccl_init(){
  int pe_rank, pe_size;
  MPI_Comm_rank(NEKO_COMM, &pe_rank);
  MPI_Comm_size(NEKO_COMM, &pe_size);

  ncclUniqueId id;
  
  if (pe_rank == 0) {
    ncclGetUniqueId(&id);
  }
  MPI_Bcast((void *) &id, sizeof(id), MPI_BYTE, 0, MPI_COMM_WORLD);
  ncclCommInitRank(&NEKO_COMM_NCCL, pe_size, id, pe_rank);
    
}

void neko_comm_nccl_finalize(){

  ncclCommDestroy(NEKO_COMM_NCCL);
  
}

#endif


