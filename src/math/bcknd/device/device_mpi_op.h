#ifndef __MATH_DEVICE_MPI_OP_H__
#define __MATH_DEVICE_MPI_OP_H__
/** 
 * Dummy defintion of device MPI reduce operations
 * @note we can't put this in device_mpi_reduce.h and include it in device_mpi_reduce.c 
 * due to the C++ interface
 */

#define DEVICE_MPI_SUM 0
#define DEVICE_MPI_MAX 1

#endif // __MATH_DEVICE_MPI_OP_H__
