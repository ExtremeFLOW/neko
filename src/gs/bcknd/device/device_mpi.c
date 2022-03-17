#include <stdlib.h>
#include <mpi.h>

void device_mpi_init_request(void **req_out) {
  MPI_Request *req = malloc(sizeof(MPI_Request));
  *req_out = req;
}

void device_mpi_free_request(void *req) {
  free(req);
}

void device_mpi_isend(void *buf_d, int *nbytes, int *rank, void *req) {
  MPI_Isend(buf_d, *nbytes, MPI_BYTE, *rank, 0, MPI_COMM_WORLD, req);
}

void device_mpi_irecv(void *buf_d, int *nbytes, int *rank, void *req) {
  MPI_Irecv(buf_d, *nbytes, MPI_BYTE, *rank, 0, MPI_COMM_WORLD, req);
}

int device_mpi_test(void *req) {
  int flag = 0;
  //MPI_Wait(req, MPI_STATUS_IGNORE);
  //flag = 1;
  MPI_Test(req, &flag, MPI_STATUS_IGNORE);
  return flag;
}
