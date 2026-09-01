!> Template for a user-defined Krylov solver.
module krylov_solver_template
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use field, only : field_t
  use gather_scatter, only : gs_t
  use krylov, only : ksp_t, ksp_monitor_t, krylov_allocate, register_krylov
  use num_types, only : rp
  use precon, only : pc_t
  use scalar_bc_projector, only : scalar_bc_projector_t
  use utils, only : neko_error
  use vector_bc_projector, only : vector_bc_projector_t, &
       vector_bc_projector_components
  implicit none
  private

  type, extends(ksp_t) :: krylov_solver_template_t
   contains
     procedure :: init => krylov_solver_template_init
     procedure :: free => krylov_solver_template_free
     procedure :: solve => krylov_solver_template_solve
     procedure :: solve_coupled => krylov_solver_template_solve_coupled
  end type krylov_solver_template_t

  public :: krylov_solver_template_register_types

contains

  subroutine krylov_solver_template_init(this, n, max_iter, M, rel_tol, &
       abs_tol, monitor)
    class(krylov_solver_template_t), target, intent(inout) :: this
    integer, intent(in) :: n, max_iter
    class(pc_t), target, optional, intent(in) :: M
    real(kind=rp), optional, intent(in) :: rel_tol, abs_tol
    logical, optional, intent(in) :: monitor

    call this%free()
    if (present(M)) call this%set_pc(M)
    if (present(rel_tol) .and. present(abs_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol)
    else if (present(rel_tol)) then
       call this%ksp_init(max_iter, rel_tol = rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol)
    else if (present(monitor)) then
       call this%ksp_init(max_iter, monitor = monitor)
    else
       call this%ksp_init(max_iter)
    end if
    ! TODO: Allocate work arrays of length n.
  end subroutine krylov_solver_template_init

  subroutine krylov_solver_template_free(this)
    class(krylov_solver_template_t), intent(inout) :: this
    ! TODO: Deallocate work arrays and other derived-type state.
    call this%ksp_free()
  end subroutine krylov_solver_template_free

  function krylov_solver_template_solve(this, Ax, x, f, n, coef, bc_projector, &
       gs_h, niter) result(ksp_results)
    class(krylov_solver_template_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: f(n)
    type(coef_t), intent(inout) :: coef
    class(scalar_bc_projector_t), intent(inout) :: bc_projector
    type(gs_t), intent(inout) :: gs_h
    integer, optional, intent(in) :: niter
    type(ksp_monitor_t) :: ksp_results

    ! TODO: Implement the Krylov iteration and return its monitor data.
    call neko_error('krylov_solver_template: solve is not implemented')
  end function krylov_solver_template_solve

  function krylov_solver_template_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, bc_projector, gs_h, niter) result(ksp_results)
    class(krylov_solver_template_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x, y, z
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: fx(n), fy(n), fz(n)
    type(coef_t), intent(inout) :: coef
    class(vector_bc_projector_t), intent(inout) :: bc_projector
    type(gs_t), intent(inout) :: gs_h
    integer, optional, intent(in) :: niter
    type(ksp_monitor_t) :: ksp_results(3)
    type(scalar_bc_projector_t), pointer :: bc_x, bc_y, bc_z

    ! A valid default when a solver applies independently to each component.
    call vector_bc_projector_components(bc_projector, bc_x, bc_y, bc_z)
    ksp_results(1) = this%solve(Ax, x, fx, n, coef, bc_x, gs_h, niter)
    ksp_results(2) = this%solve(Ax, y, fy, n, coef, bc_y, gs_h, niter)
    ksp_results(3) = this%solve(Ax, z, fz, n, coef, bc_z, gs_h, niter)
  end function krylov_solver_template_solve_coupled

  subroutine krylov_solver_template_register_types()
    procedure(krylov_allocate), pointer :: allocator
    allocator => krylov_solver_template_allocate
    call register_krylov('krylov_solver_template', allocator)
  end subroutine krylov_solver_template_register_types

  subroutine krylov_solver_template_allocate(obj)
    class(ksp_t), allocatable, intent(inout) :: obj
    allocate(krylov_solver_template_t :: obj)
  end subroutine krylov_solver_template_allocate
end module krylov_solver_template
