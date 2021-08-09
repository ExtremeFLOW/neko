!> Fast diagonalization methods from NEKTON
module fast3d
  use num_types
  use speclib
  use math
  implicit none  

contains

  !> Evaluates the derivative based on all points in the stencils  
  !! @details
  !! This set of routines comes from the appendix of                       
  !! A Practical Guide to Pseudospectral Methods, B. Fornberg              
  !! Cambridge Univ. Press, 1996.
  subroutine fd_weights_full(xx, x, n, m, c)
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(kind=rp), intent(in) :: x(0:n)
    real(kind=rp), intent(out) :: c(0:n,0:m)
    real(kind=rp), intent(in) :: xx
    real(kind=rp) :: c1, c2, c3, c4, c5
    integer :: i, j, k, mn

    c1 = 1d0
    c4 = x(0) - xx

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
       c4 = x(i) - xx
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

end module fast3d
