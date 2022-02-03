
/**
 * Device kernel for schwarz extrude
 * This can probably be done i a better way...
 * We sum the "shell" of a1 that is l1 steps in scaled with f1
 * with the shell of a2 that is l2 steps in and scale with f3.
 * Right now we just thorw away all arrays that are not
 * on the first face of dimension (nx-2)(nx-2)
 * It should be noted that a1, a2 are often the same array.
 * If l1,l2 are not the same or if one is not 0 this might lead to a race.
 */
template< typename T>
__global__ void schwarz_extrude_kernel(T * a1,
                                       const int l1,
                                       const T f1,
                                       T * a2,
                                       const int l2,
                                       const T f2,
                                       const int nx) {

  const int idx = threadIdx.x;
  const int el = blockIdx.x*nx*nx*nx;

  for(int ijk = idx; ijk<nx*nx*nx; ijk+=blockDim.x){
     int jk = ijk/nx;
     int i = ijk - nx*jk;
     int k = jk/nx;
     int j = jk -k*nx;
     if(j>0 && j< nx-1 && k > 0 && k < nx -1){
       int idx1 = i + j*nx + k*nx*nx + el;
       if(i == l1){
         int idx2 = l2 + j*nx + k*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
       if(i == nx-1-l1){
         int idx2 = nx-1-l2 + j*nx + k*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
     }
     if( i > 0 && i < nx-1 && k > 0 && k < nx -1){
       int idx1 = i + j*nx + k*nx*nx + el;
       if(j == l1){
         int idx2 = i + l2*nx + k*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
       if(j == nx-1-l1){
         int idx2 = i + (nx-1-l2)*nx + k*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
     }
     if( i > 0 && i < nx-1 && j>0 && j< nx-1 ){
       int idx1 = i + j*nx + k*nx*nx + el;
       if(k == l1){
         int idx2 = i + j*nx + l2*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
       if(k == nx-1-l1){
         int idx2 = i + j*nx + (nx-l2-1)*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
     }
  }    
}

/**
 * Device kernel for schwarz extrude
 */
template< typename T>
__global__ void schwarz_toext3d_kernel(T * __restrict__ a,
                                       T *__restrict__ b,
                                       const int nx) {

  const int idx = threadIdx.x;
  const int nx2 = nx+2;
  const int el2 = blockIdx.x*nx2*nx2*nx2;
  const int el = blockIdx.x*nx*nx*nx;
  for(int i = idx; i<nx2*nx2*nx2; i+=blockDim.x){
    a[i+el2] = 0.0;
  }
  for(int ijk = idx; ijk<nx*nx*nx; ijk+=blockDim.x){
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    a[(i+1)+(j+1)*nx2+(k+1)*nx2*nx2+el2] = b[ijk+el];
  }
}

template< typename T>
__global__ void schwarz_toreg3d_kernel(T * __restrict__ b,
                                       T *__restrict__ a,
                                       const int nx) {

  const int idx = threadIdx.x;
  const int nx2 = nx+2;
  const int el2 = blockIdx.x*nx2*nx2*nx2;
  const int el = blockIdx.x*nx*nx*nx;
  for(int ijk = idx; ijk<nx*nx*nx; ijk+=blockDim.x){
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    b[ijk+el] = a[(i+1)+(j+1)*nx2+(k+1)*nx2*nx2+el2];
  }
}
