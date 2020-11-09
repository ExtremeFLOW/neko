!> Defines a Matrix-vector product
module ax_product
  use gather_scatter
  use space
  use field
  use mesh
  use num_types
  implicit none

  !> Base type for a matrix-vector product providing \f$ Ax \f$
  type, abstract :: ax_t
   contains
     procedure(ax_compute), pass(this), deferred :: compute
  end type ax_t

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w vector of length @a n
  !! @param z vector of length @a n
  !! @param g geometric factors
  !! @param gs_h gather-scatter handle
  !! @param msh mesh
  !! @param Xh function space \f$ X_h \f$
  !! @param n integer, size of vectors
  abstract interface
     subroutine ax_compute(this, w, z, g, gs_h, msh, Xh, n)
       import space_t
       import mesh_t
       import gs_t
       import ax_t
       import dp
       implicit none
       class(ax_t), intent(inout) :: this
       type(gs_t), intent(inout) :: gs_h
       type(mesh_t), intent(inout) :: msh
       type(space_t), intent(inout) :: Xh
       integer, intent(inout) :: n
       real(kind=dp), dimension(n), intent(inout) :: w
       real(kind=dp), dimension(n), intent(inout) :: z
       real(kind=dp), intent(inout) :: g(6, Xh%lx, Xh%ly, Xh%lz)
     end subroutine ax_compute
  end interface
  
end module ax_product
