! Set mask for Dirichlet conditions
subroutine set_mask(msk)
  use field
  use mesh
  implicit none
  
  type(field_t), target, intent(inout) :: msk
  type(mesh_t), pointer :: msh
  integer :: i, j, k, l
  integer :: lx, ly, lz
  msh => msk%msh
  lx = msk%Xh%lx
  ly = msk%Xh%ly
  lz = msk%Xh%lz

  do i = 1, msh%nelv
     if (msh%facet_neigh(1, i) .eq. 0) then
        do l = 1, lz
           do k = 1, ly
              msk%x(1, k, l, i) = 0d0
           end do
        end do
     end if

     if (msh%facet_neigh(2, i) .eq. 0) then
        do l = 1, lz
           do k = 1, ly
              msk%x(lx, k, l, i) = 0d0
           end do
        end do
     end if

     if (msh%facet_neigh(3, i) .eq. 0) then
        do l = 1, lz
           do j = 1, lx
              msk%x(j, 1, l, i) = 0d0
           end do
        end do
     end if

     if (msh%facet_neigh(4, i) .eq. 0) then
        do l = 1, lz
           do j = 1, lx
              msk%x(j, ly, l, i) = 0d0
           end do
        end do
     end if

     if (msh%facet_neigh(5, i) .eq. 0) then
        do k = 1, ly
           do j = 1, lx
              msk%x(j, k, 1, i) = 0d0
           end do
        end do
     end if

     if (msh%facet_neigh(6, i) .eq. 0) then
        do k = 1, ly
           do j = 1, lx
              msk%x(j, k, lz, i) = 0d0
           end do
        end do
     end if
  end do

end subroutine set_mask

! Inverse of counting matrix
subroutine set_multiplicity(c, n, gs_h)
  use gather_scatter
  use num_types
  implicit none
  
  real(kind=dp), intent(inout), dimension(n) :: c
  type(gs_t), intent(inout) :: gs_h
  integer, intent(inout) :: n
  integer :: i

  call rone(c, n)
  call gs_op(gs_h, c, n, GS_OP_ADD)

  do i = 1, n
     c(i) = 1d0 / c(i)
  end do

end subroutine set_multiplicity

! Setup rhs
subroutine set_f(f, c, n, gs_h)
  use gather_scatter
  use num_types
  implicit none
  
  real(kind=dp), intent(inout), dimension(n) :: f
  real(kind=dp), intent(inout), dimension(n) :: c
  integer,  intent(inout) :: n
  type(gs_t), intent(inout) :: gs_h
  real(kind=dp) :: arg
  integer :: i

  do i = 1, n
     arg = 1d9 * (i * i)
     arg = 1d9 * cos(arg)
     f(i) = sin(arg)
  end do

  call gs_op(gs_h, f, n, GS_OP_ADD)
  call col2(f,c,n)

end subroutine set_f

subroutine setup_g(g, w, lx, ly, lz, n)
  use num_types
  implicit none
  
  real(kind=dp), intent(inout), dimension(6, lx, ly, lz, n) :: g
  real(kind=sp), intent(inout), dimension(lx) :: w
  integer, intent(in) :: lx, ly, lz, n
  integer :: i, j, k, l

  g = 0d0
  
  do i = 1, n
     do l = 1, lz
        do k = 1, ly
           do j = 1, lx
              g(1, j, k, l, i) = dble(w(j) * w(k) * w(l))
              g(4, j, k, l, i) = dble(w(j) * w(k) * w(l))
              g(6, j, k, l, i) = dble(w(j) * w(k) * w(l))
           end do
        end do
     end do
  end do

end subroutine setup_g
  
subroutine semhat(w, d, z, n)
  use num_types
  implicit none
  
  real(kind=sp), intent(inout), dimension(0:2*n) :: w
  real(kind=sp), intent(inout), dimension(0:n,0:n) :: d
  real(kind=sp), intent(inout), dimension(0:n) :: z
  integer, intent(inout) :: n
  integer :: i, j, np

  w = 0.
  d = 0.

  np = n + 1
  do i = 0, n
     call fd_weights_full(z(i), z, n, 1, w)
     do j = 0, n
        d(i, j) = w(j+np)  
     end do
  end do

  
  
end subroutine semhat

subroutine fd_weights_full(xx,x,n,m,c)
  use num_types
  implicit none
!
!     This routine evaluates the derivative based on all points
!     in the stencils.  It is more memory efficient than "fd_weights"
!
!     This set of routines comes from the appendix of 
!     A Practical Guide to Pseudospectral Methods, B. Fornberg
!     Cambridge Univ. Press, 1996.   (pff)
!
!    Input parameters:
!       xx -- point at wich the approximations are to be accurate
!       x  -- array of x-ordinates:   x(0:n)
!       n  -- polynomial degree of interpolant (# of points := n+1)
!       m  -- highest order of derivative to be approxxmated at xi
!
!     Output:
!       c  -- set of coefficients c(0:n,0:m).
!             c(j,k) is to be applied at x(j) when
!             the kth derivative is approxxmated by a 
!             stencil extending over x(0),x(1),...x(n).
!
!
  real(kind=sp), intent(inout) :: xx
  real(kind=sp), intent(inout) :: x(0:n)
  real(kind=sp), intent(inout) :: c(0:n, 0:m)
  integer, intent(inout) :: n
  integer :: m
  real(kind=sp) :: c1, c2, c3, c4, c5
  integer :: mn, i, j, k

  c1       = 1.
  c4       = x(0) - xx
  do k=0,m
     do j=0,n
        c(j,k) = 0.
     enddo
  enddo
  c(0,0) = 1.
 
  do i=1,n
     mn = min(i,m)
     c2 = 1.
     c5 = c4
     c4 = x(i)-xx
     do j=0,i-1
        c3 = x(i)-x(j)
        c2 = c2*c3
        do k=mn,1,-1
           c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
        enddo
        c(i,0) = -c1*c5*c(i-1,0)/c2
        do k=mn,1,-1
           c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            enddo
            c(j,0) = c4*c(j,0)/c3
         enddo
         c1 = c2
      enddo
      return
    end subroutine fd_weights_full
    
