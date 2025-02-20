#ifndef __COMM_NCCL_H
#define __COMM_NCCL_H

#if defined(HAVE_NCCL) || defined(HAVE_RCCL)

#ifdef HAVE_NCCL
#include <nccl.h>
#endif

#ifdef HAVE_RCCL
#include <rccl.h>
#endif

/** 
 * NCCL communicator
 */
extern ncclComm_t NEKO_COMM_NCCL;


#endif

#endif
