! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC. 
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF 
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract 
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE, 
! LLC nor any of their employees, makes any warranty, 
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process, 
! or services by trade name, trademark, manufacturer or otherwise does 
! not necessarily constitute or imply its endorsement, recommendation, 
! or favoring by the United States Government or UCHICAGO ARGONNE LLC. 
! The views and opinions of authors expressed 
! herein do not necessarily state or reflect those of the United States 
! Government or UCHICAGO ARGONNE, LLC, and shall 
! not be used for advertising or product endorsement purposes.
!
!> Fast diagonalization methods from NEKTON
module fast3d
  use num_types
  use speclib
  use math
  implicit none

  public :: fd_weights_full, semhat, setup_intp

contains

  !> Compute finite-difference stencil weights for evaluating derivatives up to
  !! order \f$m\f$ at a point.
  !! @details
  !! This routine comes from the Appendix C of
  !! "A Practical Guide to Pseudospectral Methods" by B. Fornberg,
  !! Cambridge University Press, 1996.
  !!
  !! Given gridpoints \f$ x_0, x_1, \dots x_n \f$ and some point \f$\xi\f$
  !! (not necessarily a grid point!) find weights \f$ c_{j, k} \f$, such that
  !! the expansions 
  !! \f$ \frac{d^k f}{d x^k}|_{x=\xi} \approx \sum_{j=0}^n c_{j,k} f(x_j)\f$,
  !! \f$k=0, \dots m\f$ are optimal.
  !! Note that finite-difference stencils are exactly such type of expansions.
  !! For the derivation of the algorithm, refer to 3.1 in the reference above.
  !!
  !! @note - Setting \f$m=0\f$ makes is this a polynomial interpolation routine.
  !! It is the fastest such routine possible for a single interpolation point,
  !! according to the above reference.
  !! @note - The name `_full` refers to the fact that we use the values \f$f(x_j)\f$
  !! at all available nodes \f$x\f$ to construct the expansion. So we always
  !! get the finite difference stencil of maximum order possible.
  !!
  !! @warning The calculation of the wieghts is numerically stable.
  !! But applying the weights to a function can be ill-conditioned in the case
  !! of high-order derivatives. 
  !!
  !! @param xi Point at which the approximations are to be accurate
  !! @param x The coordinates for the grid points
  !! @param[in] n The size of `x` is `n + 1`
  !! @param[in] m Highest order of derivative to be approximated
  !! @param c The stencil weights. Row j corresponds to weight of \f$f(x_j)\f$
  !! and column k to the kth derivative.
  subroutine fd_weights_full(xi, x, n, m, c)
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(kind=rp), intent(in) :: x(0:n)
    real(kind=rp), intent(out) :: c(0:n,0:m)
    real(kind=rp), intent(in) :: xi
    real(kind=rp) :: c1, c2, c3, c4, c5
    integer :: i, j, k, mn

    c1 = 1d0
    c4 = x(0) - xi

    do k = 0, m
       do j = 0, n
          c(j,k) = 0d0
       end do
    end do

    c(0,0) = 1d0

    do i = 1, n                                                              
       mn = min(i,m)
       c2 = 1d0  
       c5 = c4                                                       
       c4 = x(i) - xi
       do j = 0, i - 1                                                  
          c3 = x(i) - x(j)
          c2 = c2 * c3                                                    
          do k = mn, 1, -1                    
             c(i,k) = c1 * (k * c(i-1,k-1) - c5 * c(i-1,k)) / c2
          end do
          c(i,0) = -c1 * c5 * c(i-1,0) / c2
          do k = mn, 1, -1                           
             c(j,k) = (c4 * c(j,k) - k * c(j,k-1)) / c3
          end do
          c(j,0) = c4 * c(j,0) / c3
       end do
       c1 = c2
    end do
    
  end subroutine fd_weights_full
     
  
!>  Generate matrices for single element, 1D operators:
!!        a    = Laplacian
!!        b    = diagonal mass matrix
!!        c    = convection operator b*d
!!        d    = derivative matrix
!!        dgll = derivative matrix,    mapping from pressure nodes to velocity
!!        jgll = interpolation matrix, mapping from pressure nodes to velocity
!!        z    = GLL points
!!
!!        zgl  = GL points
!!        bgl  = diagonal mass matrix on GL
!!        dgl  = derivative matrix,    mapping from velocity nodes to pressure
!!        jgl  = interpolation matrix, mapping from velocity nodes to pressure
!!
!!        n    = polynomial degree (velocity space)
!!        w    = work array of size 2*n+2
!!
!!     Currently, this is set up for pressure nodes on the interior GLL pts.
!!
!!
  subroutine semhat(a, b, c, d, z, dgll, jgll, bgl, zgl, dgl, jgl, n, w)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: a(0:n,0:n)
    real(kind=rp), intent(inout) :: b(0:n)
    real(kind=rp), intent(inout) :: c(0:n,0:n)
    real(kind=rp), intent(inout) :: d(0:n,0:n)
    real(kind=rp), intent(inout) :: z(0:n)
    real(kind=rp), intent(inout) :: dgll(0:n,1:n-1),jgll(0:n,1:n-1)
    real(kind=rp), intent(inout) :: bgl(1:n-1)
    real(kind=rp), intent(inout) :: zgl(1:n-1)
    real(kind=rp), intent(inout) :: dgl(1:n-1,0:n)
    real(kind=rp), intent(inout) :: jgl(1:n-1,0:n)
    real(kind=rp), intent(inout) :: w(0:2*n+1)
    integer :: np, nm, n2, i, j, k
    np = n+1
    nm = n-1
    n2 = n-2
    call zwgll(z, b, np)
    do i = 0,n
       call fd_weights_full(z(i), z, n, 1, w)
       do j = 0,n
          d(i,j) = w(j+np)                   !  Derivative matrix
       end do
    end do

    if (n.eq.1) return                       !  No interpolation for n=1

    do i = 0,n
       call fd_weights_full(z(i), z(1), n2, 1, w(1))
       do j = 1,nm
          jgll(i,j) = w(j   )                  !  Interpolation matrix
          dgll(i,j) = w(j+nm)                  !  Derivative    matrix
       end do
    end do
    call rzero(a, np*np)
    do j = 0,n
    do i = 0,n
       do k = 0,n
          a(i,j) = a(i,j) + d(k,i)*b(k)*d(k,j)
       end do
       c(i,j) = b(i)*d(i,j)
    end do
    end do
    call zwgl(zgl, bgl, nm)
    do i = 1,n-1
       call fd_weights_full(zgl(i), z, n, 1, w)
       do j = 0,n
          jgl(i,j) = w(j   )                  !  Interpolation matrix
          dgl(i,j) = w(j+np)                  !  Derivative    matrix
       end do
    end do
  end subroutine semhat

  !> Computes interpolation between points zf, zc
  !! The interpolation vectors are stored in jh, jht
  !! nf, nc is the number of points in the spaces
  !! derivate specifies if we want the derivative interpolation instead
  !! derivate = 1 gives the first derivative etc.
  subroutine setup_intp(jh, jht, zf, zc, nf, nc, derivate)
    integer, intent(in) :: nf, nc, derivate
    real(kind=rp), intent(inout) :: jh(nf,nc), zf(nf), zc(nc), jht(nc,nf)
    real(kind=rp) ::  w(nc,0:derivate)
    integer :: i, j
    do i = 1, nf
       call fd_weights_full(zf(i), zc, nc-1, derivate, w)
       do j = 1, nc
          jh(i,j) = w(j, derivate)
          jht(j,i) = w(j, derivate)
       end do
    end do
  end subroutine setup_intp
end module fast3d
