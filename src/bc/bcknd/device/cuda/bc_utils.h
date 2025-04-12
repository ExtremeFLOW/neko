
#ifndef __BC_UTILS_KERNEL__
#define __BC_UTILS_KERNEL__

/**
 * Computes the linear index for area and normal arrays
 * @note Fortran indexing input, C indexing output
 */

 #define coef_normal_area_idx(i, j, k, l, lx, nf) \
 (((i) + (lx) * (((j) - 1) + (lx) * (((k) - 1) + (nf) * (((l) - 1))))) - 1)

/**
* Device function to compute i,j,k,e indices from a linear index
* @note Assumes idx is a Fortran index
*/
__device__
void nonlinear_index(const int idx, const int lx, int *index) {
 const int idx2 = idx -1;
 index[3] = idx2/(lx * lx * lx) ;
 index[2] = (idx2 - (lx*lx*lx)*index[3])/(lx * lx);
 index[1] = (idx2 - (lx*lx*lx)*index[3] - (lx*lx) * index[2]) / lx;
 index[0] = (idx2 - (lx*lx*lx)*index[3] - (lx*lx) * index[2]) - lx*index[1];
 index[0]++;
 index[1]++;
 index[2]++;
 index[3]++;
}

#endif