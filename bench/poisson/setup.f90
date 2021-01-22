! Set Dirichlet conditions
subroutine set_bc(bc_, msh)
  use neko
  implicit none
  
  type(mesh_t), intent(in) :: msh
  type(dirichlet_t), intent(inout) :: bc_
  integer :: i

  do i = 1, msh%nelv
     if (msh%facet_neigh(1, i) .eq. 0) then
       call bc_%mark_facet(1, i) 
     end if
     if (msh%facet_neigh(2, i) .eq. 0) then
       call bc_%mark_facet(2, i) 
     end if
     if (msh%facet_neigh(3, i) .eq. 0) then
       call bc_%mark_facet(3, i) 
     end if
     if (msh%facet_neigh(4, i) .eq. 0) then
       call bc_%mark_facet(4, i) 
     end if
     if (msh%facet_neigh(5, i) .eq. 0) then
       call bc_%mark_facet(5, i) 
     end if
     if (msh%facet_neigh(6, i) .eq. 0) then
       call bc_%mark_facet(6, i) 
     end if
  enddo
end subroutine set_bc

! Setup rhs
subroutine set_f(f, c, dm, n, gs_h)
  use gather_scatter
  use num_types
  use utils
  implicit none
  
  real(kind=dp), intent(inout), dimension(n) :: f
  real(kind=dp), intent(inout), dimension(n) :: c
  type(dofmap_t), intent(in) :: dm
  integer,  intent(inout) :: n
  type(gs_t), intent(inout) :: gs_h
  real(kind=dp) :: arg, dx, dy, dz
  integer :: i, idx(4)

  do i = 1, n
     idx = nonlinear_index(i, dm%Xh%lx, dm%Xh%ly, dm%Xh%lz)
     dx = dm%x(idx(1), idx(2), idx(3), idx(4)) - 4;
     dy = dm%y(idx(1), idx(2), idx(3), idx(4)) - 4;
     dz = dm%z(idx(1), idx(2), idx(3), idx(4)) - 4;
     f(i) = 500.0*exp(-(dx**2 + dy**2 + dz**2)/2);
  end do

  call gs_op(gs_h, f, n, GS_OP_ADD)
  call col2(f,c,n)
end subroutine set_f

