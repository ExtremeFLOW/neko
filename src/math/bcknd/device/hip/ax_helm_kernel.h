/*
 Copyright (c) 2021-2022, The Neko Authors
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
 * Device kernel for Ax helm
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void ax_helm_kernel(T * __restrict__ w,            
                               const T * __restrict__ u,      
                               const T * __restrict__ dx,     
                               const T * __restrict__ dy,     
                               const T * __restrict__ dz,     
                               const T * __restrict__ dxt,    
                               const T * __restrict__ dyt,    
                               const T * __restrict__ dzt,    
                               const T * __restrict__ h1,     
                               const T * __restrict__ g11,    
                               const T * __restrict__ g22,    
                               const T * __restrict__ g33,    
                               const T * __restrict__ g12,    
                               const T * __restrict__ g13,    
                               const T * __restrict__ g23) {  
                                                                               
  __shared__ T shdx[LX*LX];                                                    
  __shared__ T shdy[LX*LX];                                                    
  __shared__ T shdzt[LX*LX];                                                   
                                                                               
  __shared__ T shdxt[LX*LX];                                                   
  __shared__ T shdyt[LX*LX];                                                   
  __shared__ T shdz[LX*LX];                                                    
                                                                               
  __shared__ T shu[LX*LX*LX];                                                  
  __shared__ T shur[LX*LX*LX];                                                 
  __shared__ T shus[LX*LX*LX];                                                 
  __shared__ T shut[LX*LX*LX];                                                 
                                                                               
  int l,i,j,k,n;                                                               
                                                                               
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1)/CHUNKS + 1;                           
                                                                               
  if (iii<LX*LX) {                                                             
    shdx[iii] = dx[iii];                                                       
    shdy[iii] = dy[iii];                                                       
    shdz[iii] = dz[iii];                                                       
  }                                                                            
  i = iii;                                                                     
  while (i < LX * LX * LX){                                                    
    shu[i] = u[i+e*LX*LX*LX];                                                  
    i = i + CHUNKS;                                                            
  }                                                                            
                                                                               
  __syncthreads();
                                                                               
  if (iii<LX*LX){                                                              
    shdxt[iii] = dxt[iii];                                                     
    shdyt[iii] = dyt[iii];                                                     
    shdzt[iii] = dzt[iii];                                                     
  }                                                                            
                                                                               
  for (n=0; n<nchunks; n++){                                                   
    const int ijk = iii+n*CHUNKS;                                              
    const int jk = ijk/LX;                                                     
    i = ijk-jk*LX;                                                             
    k = jk/LX;                                                                 
    j = jk-k*LX;                                                               
    if (i<LX && j<LX && k<LX && ijk < LX*LX*LX){        
      T rtmp = 0.0;                                                         
      T stmp = 0.0;                                                            
      T ttmp = 0.0;                                                            
      for (l = 0; l<LX; l++){                                                  
        rtmp = rtmp + shdx[i+l*LX] * shu[l+j*LX+k*LX*LX];                      
        stmp = stmp + shdy[j+l*LX] * shu[i+l*LX+k*LX*LX];                      
        ttmp = ttmp + shdz[k+l*LX] * shu[i+j*LX+l*LX*LX];                      
      }                                                                        
      shur[ijk] = h1[ijk+e*LX*LX*LX]                                           
                * (g11[ijk+e*LX*LX*LX] * rtmp                                  
                   + g12[ijk+e*LX*LX*LX] * stmp                                
                   + g13[ijk+e*LX*LX*LX] * ttmp);                              
      shus[ijk] = h1[ijk+e*LX*LX*LX]                                           
                * (g12[ijk+e*LX*LX*LX] * rtmp                                  
                   + g22[ijk+e*LX*LX*LX] * stmp                                
                   + g23[ijk+e*LX*LX*LX] * ttmp);                              
      shut[ijk] = h1[ijk+e*LX*LX*LX]                                           
                * (g13[ijk+e*LX*LX*LX] * rtmp                                  
                   + g23[ijk+e*LX*LX*LX] * stmp                                
                   + g33[ijk+e*LX*LX*LX] * ttmp);                              
    }                                                                          
  }                                                                            
                                                                               
  __syncthreads();
                                                                               
  for (n=0; n<nchunks; n++){                                                   
    const int ijk = iii+n*CHUNKS;                                              
    const int jk = ijk/LX;                                                     
    i = ijk-jk*LX;                                                             
    k = jk/LX;                                                                 
    j = jk-k*LX;                                                               
    if (i<LX && j<LX && k<LX && ijk <LX*LX*LX){                 
      T wijke = 0.0;                                                        
      for (l = 0; l<LX; l++){                                                  
        wijke = wijke                                                          
              + shdxt[i+l*LX] * shur[l+j*LX+k*LX*LX]                           
              + shdyt[j+l*LX] * shus[i+l*LX+k*LX*LX]                           
              + shdzt[k+l*LX] * shut[i+j*LX+l*LX*LX];                          
      }                                                                        
      w[ijk+e*LX*LX*LX] = wijke;                                               
    }                                                                          
  }                                                                            
}

