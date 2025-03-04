
#ifdef HAVE_NVSHMEM

#include <nvshmemx.h>
#include <comm/comm.h>

extern "C" {

  /**
   * Setup Neko's nvshmem communicator
   */  
  void neko_comm_nvshmem_init(){

    nvshmemx_init_attr_t attr;
    attr.mpi_comm = &NEKO_COMM;      
    if (nvshmemx_init_status() == NVSHMEM_STATUS_NOT_INITIALIZED)
    {
      nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attr);
    }
  }

  void neko_comm_nvshmem_finalize() {
    nvshmem_finalize();
  }
    
  
}

#endif
