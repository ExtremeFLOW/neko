#ifndef __COMM_NCCL_H
#define __COMM_NCCL_H

#ifdef HAVE_NCCL
#include <nccl.h>

/** 
 * NCCL communicator
 */
extern ncclComm_t NEKO_COMM_NCCL;


#endif

#endif
