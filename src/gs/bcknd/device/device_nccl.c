/*
 Copyright (c) 2025, The Neko Authors
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
 * C wrapper for NCCL point-to-point calls,
 */

#include <stdlib.h>
#include <stdio.h>
#include <comm/comm_nccl.h>

void device_nccl_sendrecv(void *sbuf_d, int soffset, int scount, int srank,
                          void *rbuf_d, int roffset, int rcount, int rrank,
                          int nbytes, void *stream) {

#if defined(HAVE_NCCL) || defined(HAVE_RCCL)
  if (nbytes == sizeof(float)) {
    ncclGroupStart();
    ncclSend(sbuf_d+soffset, scount, ncclFloat, srank, NEKO_COMM_NCCL, stream);
    ncclRecv(rbuf_d+roffset, rcount, ncclFloat, rrank, NEKO_COMM_NCCL, stream);
    ncclGroupEnd();
  }
  else if (nbytes == sizeof(double)) {
    ncclGroupStart();
    ncclSend(sbuf_d+soffset, scount, ncclFloat64, srank, NEKO_COMM_NCCL, stream);
    ncclRecv(rbuf_d+roffset, rcount, ncclFloat64, rrank, NEKO_COMM_NCCL, stream);
    ncclGroupEnd();
  }
  else {
    fprintf(stderr, __FILE__ ": Invalid data type)\n");
    exit(1);
  }
#endif
}
