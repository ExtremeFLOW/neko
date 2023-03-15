/*
 Copyright (c) 2022, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * C wrapper for MPI calls, since passing device pointers does not work in the
 * Fortran MPI interface.
 * @note We use @c MPI_COMM_WORLD which @e should be equal to @c NEKO_COMM.
 */

#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

void device_mpi_init_reqs(int n, void **reqs_out) {
  MPI_Request *reqs = malloc(n * sizeof(MPI_Request));
  *reqs_out = reqs;
}

void device_mpi_free_reqs(void **reqs) {
  free(*reqs);
  *reqs = NULL;
}

void device_mpi_isend(void *buf_d, int offset, int nbytes, int rank,
		      void *vreqs, int i) {
  MPI_Request *reqs = vreqs;
  int tid = omp_get_thread_num();
  MPI_Isend(buf_d+offset, nbytes, MPI_BYTE, rank, tid, MPI_COMM_WORLD, &reqs[i-1]);
}

void device_mpi_irecv(void *buf_d, int offset, int nbytes, int rank,
		      void *vreqs, int i) {
  MPI_Request *reqs = vreqs;
  int tid = omp_get_thread_num();
  MPI_Irecv(buf_d+offset, nbytes, MPI_BYTE, rank, tid, MPI_COMM_WORLD, &reqs[i-1]);
}

int device_mpi_test(void *vreqs, int i) {
  MPI_Request *reqs = vreqs;
  int flag = 0;
  MPI_Test(&reqs[i-1], &flag, MPI_STATUS_IGNORE);
  return flag;
}

void device_mpi_waitall(int n, void *vreqs) {
  MPI_Request *reqs = vreqs;
  MPI_Waitall(n, reqs, MPI_STATUSES_IGNORE);
}

int device_mpi_waitany(int n, void *vreqs, int *i) {
  MPI_Request *reqs = vreqs;
  int j;
  MPI_Waitany(n, reqs, &j, MPI_STATUSES_IGNORE);
  if (j == MPI_UNDEFINED) {
    return 0;
  } else {
    *i = j + 1;
    return 1;
  }
}
