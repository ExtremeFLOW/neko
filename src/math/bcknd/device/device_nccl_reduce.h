#ifndef __MATH_DEVICE_NCCL_REDUCE_H__
#define __MATH_DEVICE_NCCL_REDUCE_H__

/**
 * C wrapper for NCCL reduction calls
 */

extern "C" {
void device_nccl_allreduce(void *sbuf_d, void *rbuf_d, int count,
                        int nbytes, int op, void *stream);  
}

#endif
