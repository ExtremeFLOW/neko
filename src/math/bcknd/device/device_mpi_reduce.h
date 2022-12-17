
/**
 * C wrapper for MPI calls, until we integrate with NCCL/RCCL
 * @note We use @c MPI_COMM_WORLD which @e should be equal to @c NEKO_COMM.
 */

extern "C" {
  void device_mpi_allreduce(void *buf_d, void *buf, int count, size_t nbytes);
}
