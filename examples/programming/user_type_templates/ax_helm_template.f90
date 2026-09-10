!> Template for a user-defined Helmholtz matrix-vector product.
module ax_helm_template
  use ax_product, only : ax_t, ax_helm_allocate, register_ax_helm
  use ax_helm, only : ax_helm_t
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use num_types, only : rp
  use space, only : space_t
  implicit none
  private

  !> A scalar Helmholtz-like operator.
  !! Extending ax_helm_t supplies compute_vector by applying compute to each
  !! component independently.
  type, extends(ax_helm_t) :: ax_helm_template_t
   contains
     procedure, pass(this) :: compute => ax_helm_template_compute
  end type ax_helm_template_t

  public :: ax_helm_template_register_types

contains

  !> Apply the user-defined Helmholtz matrix-vector product.
  !! @param w Result vector of size @a (lx,ly,lz,nelv).
  !! @param u Input vector of size @a (lx,ly,lz,nelv).
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  subroutine ax_helm_template_compute(this, w, u, coef, msh, Xh)
    class(ax_helm_template_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    ! TODO: Replace the identity operation with the desired Ax product.
    w = u
  end subroutine ax_helm_template_compute

  !> Register the user-defined Helmholtz operator.
  subroutine ax_helm_template_register_types()
    procedure(ax_helm_allocate), pointer :: allocator
    allocator => ax_helm_template_allocate
    call register_ax_helm('ax_helm_template', allocator)
  end subroutine ax_helm_template_register_types

  !> Allocate the user-defined Helmholtz operator.
  !! @param obj Operator to allocate.
  subroutine ax_helm_template_allocate(obj)
    class(ax_t), allocatable, intent(inout) :: obj
    allocate(ax_helm_template_t :: obj)
  end subroutine ax_helm_template_allocate
end module ax_helm_template
