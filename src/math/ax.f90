!> Defines a Matrix-vector product
module ax_product
  use mxm_wrapper
  use num_types
  use coefs
  use space
  use mmb_field
  use field
  use mesh
  implicit none

  !> Base type for a matrix-vector product providing \f$ Ax \f$
  type, abstract :: ax_t
   contains
     procedure(ax_compute), nopass, deferred :: compute
  end type ax_t

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w vector of size @a (lx,ly,lz,nelv)
  !! @param z vector of size @a (lx,ly,lz,nelv)
  !! @param coef Coefficients
  !! @param msh mesh
  !! @param Xh function space \f$ X_h \f$
  abstract interface
  subroutine ax_compute(w, u, coef, msh, Xh)
       import space_t
       import mesh_t
       import coef_t
       import ax_t
       import rp
       implicit none
       type(space_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh       
       type(coef_t), intent(inout) :: coef
       real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
       real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
     end subroutine ax_compute
  end interface


  !> Base type for a matrix-vector product providing \f$ Ax \f$
  type, abstract :: ax_mmb_t
   contains
     procedure(ax_mmb_compute), nopass, deferred :: compute
  end type ax_mmb_t

  !> Abstract interface for computing\f$ Ax \f$ inside a Krylov method
  !!
  !! @param w vector of size @a (lx,ly,lz,nelv)
  !! @param z vector of size @a (lx,ly,lz,nelv)
  !! @param coef Coefficients
  !! @param msh mesh
  !! @param Xh function space \f$ X_h \f$
  abstract interface
  subroutine ax_mmb_compute(w, u, coef, msh, Xh)
       import space_t
       import mesh_t
       import coef_t
       import ax_t
       import mmb_field_t
       import rp
       implicit none
       type(space_t), intent(inout) :: Xh
       type(mesh_t), intent(inout) :: msh       
       type(coef_t), intent(inout) :: coef
       type(mmb_field_t), intent(inout) :: w, u
     end subroutine ax_mmb_compute
  end interface
  
end module ax_product
