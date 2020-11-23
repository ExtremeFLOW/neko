!> Fast diagonalization methods from NEKTON
module fast3d
  use num_types
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
    real(kind=dp), intent(in) :: x(0:n)
    real(kind=dp), intent(out) :: c(0:n,0:m)
    real(kind=dp), intent(in) :: xx
    real(kind=dp) :: c1, c2, c3, c4, c5
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
 
end module fast3d
