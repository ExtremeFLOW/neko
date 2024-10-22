! Set mask for Dirichlet conditions
subroutine set_mask(msk, msh, lx, ly, lz, n)
  use utils
  use mesh
  implicit none
  
  integer, intent(inout) :: msk(0:n)
  type(mesh_t), intent(in) :: msh
  integer, intent(in) :: lx
  integer, intent(in) :: ly
  integer, intent(in) :: lz
  integer, intent(in) :: n
  integer :: i, j, k, l, msk_c

  msk_c = 0
  do i = 1, msh%nelv
     if (msh%facet_neigh(1, i) .eq. 0) then
        do l = 1, lz
           do k = 1, ly
              msk_c = msk_c + 1
              msk(msk_c) = linear_index(1,k,l,i,lx,ly,lz)
           end do
        end do
     end if

     if (msh%facet_neigh(2, i) .eq. 0) then
        do l = 1, lz
           do k = 1, ly
              msk_c = msk_c + 1
              msk(msk_c) = linear_index(lx,k,l,i,lx,ly,lz)
           end do
        end do
     end if

     if (msh%facet_neigh(3, i) .eq. 0) then
        do l = 1, lz
           do j = 1, lx
              msk_c = msk_c + 1
              msk(msk_c) = linear_index(j,1,l,i,lx,ly,lz)
           end do
        end do
     end if

     if (msh%facet_neigh(4, i) .eq. 0) then
        do l = 1, lz
           do j = 1, lx
              msk_c = msk_c + 1
              msk(msk_c) = linear_index(j,ly,l,i,lx,ly,lz)
           end do
        end do
     end if

     if (msh%facet_neigh(5, i) .eq. 0) then
        do k = 1, ly
           do j = 1, lx
              msk_c = msk_c + 1
              msk(msk_c) = linear_index(j,k,1,i,lx,ly,lz)
           end do
        end do
     end if

     if (msh%facet_neigh(6, i) .eq. 0) then
        do k = 1, ly
           do j = 1, lx
              msk_c = msk_c + 1
              msk(msk_c) = linear_index(j,k,lz,i,lx,ly,lz)
           end do
        end do
     end if
  end do
  msk(0) = msk_c
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
  real(kind=dp), intent(inout), dimension(lx) :: w
  integer, intent(in) :: lx, ly, lz, n
  integer :: i, j, k, l

  g = 0d0
  
  do i = 1, n
     do l = 1, lz
        do k = 1, ly
           do j = 1, lx
              g(1, j, k, l, i) = w(j) * w(k) * w(l)
              g(4, j, k, l, i) = w(j) * w(k) * w(l)
              g(6, j, k, l, i) = w(j) * w(k) * w(l)
           end do
        end do
     end do
  end do

end subroutine setup_g
  
