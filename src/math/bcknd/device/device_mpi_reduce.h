#ifndef __MATH_DEVICE_MPI_REDUCE_H__
#define __MATH_DEVICE_MPI_REDUCE_H__

/**
 * C wrapper for MPI calls, until we integrate with NCCL/RCCL
 */

extern "C" {
  void device_mpi_allreduce(void *buf_d, void *buf, int count,
                            size_t nbytes, int op);

  void device_mpi_allreduce_inplace(void *buf_d, int count,
                                    size_t nbytes, int op);
}

#endif // __MATH_DEVICE_MPI_REDUCE_H__
