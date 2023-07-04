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
 * C wrapper for MPI calls, until we integrate with NCCL/RCCL
 * @note We use @c MPI_COMM_WORLD which @e should be equal to @c NEKO_COMM.
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "device_mpi_op.h"

void device_mpi_allreduce(void *buf_d, void *buf, int count, int nbytes, int op) {

  if (nbytes == sizeof(float)) {
    if (op == DEVICE_MPI_SUM)
      MPI_Allreduce(buf_d, buf, count, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    else if (op == DEVICE_MPI_MAX)
      MPI_Allreduce(buf_d, buf, count, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    else {
      fprintf(stderr, __FILE__ ": Invalid reduction op)\n");
      exit(1);
    }
  }
  else if (nbytes == sizeof(double)) {
    if (op == DEVICE_MPI_SUM)
      MPI_Allreduce(buf_d, buf, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    else if (op == DEVICE_MPI_MAX)
      MPI_Allreduce(buf_d, buf, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    else {
      fprintf(stderr, __FILE__ ": Invalid reduction op)\n");
      exit(1);
    }
  }
  else {
    fprintf(stderr, __FILE__ ": Invalid data type)\n");
    exit(1);
  }
}

void device_mpi_allreduce_inplace(void *buf_d, int count, int nbytes, int op) {

  if (nbytes == sizeof(float)) {
    if (op == DEVICE_MPI_SUM)
      MPI_Allreduce(MPI_IN_PLACE, buf_d, count,
                    MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    else if (op == DEVICE_MPI_MAX)
      MPI_Allreduce(MPI_IN_PLACE, buf_d, count,
                    MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    else {
      fprintf(stderr, __FILE__ ": Invalid reduction op)\n");
      exit(1);
    }
  }
  else if (nbytes == sizeof(double)) {
    if (op == DEVICE_MPI_SUM)
      MPI_Allreduce(MPI_IN_PLACE, buf_d, count,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    else if (op == DEVICE_MPI_MAX)
      MPI_Allreduce(MPI_IN_PLACE, buf_d, count,
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  else {
    fprintf(stderr, __FILE__ ": Invalid data type)\n");
    exit(1);
  }
}

