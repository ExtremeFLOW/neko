/** Wrapper interface of p4est and sc libraries
 *
 */
#ifndef NEKO_AMR_H_
#define NEKO_AMR_H_

/* Numbers designating the type and action for refinement  */
#define AMR_RM_NONE    0  /* No refinement action */
#define AMR_RM_H_REF   1  /* Refine by splitting element */
#define AMR_RM_H_CRS  -1  /* Coarsen by merging element */
#define AMR_RM_P_REF   2  /* Refine by rising polynomial order in element */
#define AMR_RM_P_CRS  -2  /* Coarsen by lowering polynomial order in element */

#endif /* NEKO_AMR_H_ */
