!> Defines a Matrix-vector product
module ax_product
  use gather_scatter
  use num_types
  use space
  use field
  use mesh
  use bc
  implicit none

  !> Base type for a matrix-vector product providing \f$ Ax \f$
  type, abstract :: ax_t
   contains
     procedure(ax_compute), nopass, deferred :: compute
  end type ax_t

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w vector of length @a n
  !! @param z vector of length @a n
  !! @param g geometric factors
  !! @param msh mesh
  !! @param Xh function space \f$ X_h \f$
  !! @param n integer, size of vectors
  abstract interface
  subroutine ax_compute(w, u, dof, Xh, n)
       import bc_list_t
       import space_t
       import dofmap_t
       import gs_t
       import ax_t
       import dp
       implicit none
       type(dofmap_t), intent(inout) :: dof
       type(space_t), intent(inout) :: Xh
       integer, intent(inout) :: n
       real(kind=dp), dimension(n), intent(inout) :: w
       real(kind=dp), dimension(n), intent(inout) :: u
     end subroutine ax_compute
  end interface
  
end module ax_product
