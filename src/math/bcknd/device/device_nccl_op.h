#ifndef __MATH_DEVICE_NCCL_OP_H__
#define __MATH_DEVICE_NCCL_OP_H__
/** 
 * Dummy defintion of NCCL reduce operations
 * @note we can't put this in device_nccl_reduce.h and include it in device_nccl_reduce.c 
 * due to the C++ interface
 */

#define DEVICE_NCCL_SUM 0
#define DEVICE_NCCL_MAX 1

#endif 
