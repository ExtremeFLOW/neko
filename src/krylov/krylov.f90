! Copyright (c) 2020-2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the base abstract type for Krylov solvers plus helper types.
module krylov
  use gather_scatter, only : gs_t, GS_OP_ADD
  use ax_product, only : ax_t
  use num_types, only: rp, c_rp
  use precon, only : pc_t
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use field, only : field_t
  use utils, only : neko_error, neko_warning
  use bc, only : bc_list_t
  use identity, only : ident_t
  use device_identity, only : device_ident_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use logger, only : neko_log, LOG_SIZE
  implicit none
  private

  integer, public, parameter :: KSP_MAX_ITER = 1e3       !< Maximum number of iters.
  real(kind=rp), public, parameter :: KSP_ABS_TOL = 1d-9 !< Absolut tolerance
  real(kind=rp), public, parameter :: KSP_REL_TOL = 1d-9 !< Relative tolerance

  !> Type for storing initial and final residuals in a Krylov solver.
  type, public :: ksp_monitor_t
     !> Iteration number.
     integer :: iter
     !> Initial residual.
     real(kind=rp) :: res_start
     !> FInal residual
     real(kind=rp) :: res_final
  end type ksp_monitor_t

  !> Base abstract type for a canonical Krylov method, solving \f$ Ax = f \f$.
  type, public, abstract :: ksp_t
     class(pc_t), pointer :: M => null() !< Preconditioner
     real(kind=rp) :: rel_tol            !< Relative tolerance
     real(kind=rp) :: abs_tol            !< Absolute tolerance
     integer :: max_iter                 !< Maximum number of iterations
     class(pc_t), allocatable :: M_ident !< Internal preconditioner (Identity)
     logical :: monitor                  !< Turn on/off monitoring
   contains
     !> Base type constructor.
     procedure, pass(this) :: ksp_init => krylov_init
     !> Base type destructor.
     procedure, pass(this) :: ksp_free => krylov_free
     !> Set preconditioner.
     procedure, pass(this) :: set_pc => krylov_set_pc
     !> Solve the system.
     procedure(ksp_method), pass(this), deferred :: solve
     !> Solve the system (coupled version).
     procedure(ksp_method_coupled), pass(this), deferred :: solve_coupled
     !> Monitor start
     procedure, pass(this) :: monitor_start => krylov_monitor_start
     !> Monitor stop
     procedure, pass(this) :: monitor_stop => krylov_monitor_stop
     !> Monitor iteration
     procedure, pass(this) :: monitor_iter => krylov_monitor_iter
     !> Destructor.
     procedure(ksp_t_free), pass(this), deferred :: free
  end type ksp_t


  !> Abstract interface for a Krylov method's solve routine
  !!
  !! @param x field to solve for
  !! @param f right hand side
  !! @param n integer, size of vectors
  !! @param coef Coefficients
  !! @param blst list of  boundary conditions
  !! @param gs_h Gather-scatter handle
  !! @param niter iteration trip count
  abstract interface
     function ksp_method(this, Ax, x, f, n, coef, blst, gs_h, niter) &
          result(ksp_results)
       import :: bc_list_t
       import :: field_t
       import :: ksp_t
       import :: coef_t
       import :: gs_t
       import :: ax_t
       import :: ksp_monitor_t
       import rp
       implicit none
       class(ksp_t), intent(inout) :: this
       class(ax_t), intent(inout) :: Ax
       type(field_t), intent(inout) :: x
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: f
       type(coef_t), intent(inout) :: coef
       type(bc_list_t), intent(inout) :: blst
       type(gs_t), intent(inout) :: gs_h
       integer, optional, intent(in) :: niter
       type(ksp_monitor_t) :: ksp_results
     end function ksp_method
  end interface

  !> Abstract interface for a Krylov method's coupled solve routine
  !!
  !! @param x field to solve for
  !! @param y field to solve for
  !! @param z field to solve for
  !! @param fx right hand side
  !! @param fy right hand side
  !! @param fz right hand side
  !! @param n integer, size of vectors
  !! @param coef Coefficients
  !! @param blst list of  boundary conditions
  !! @param gs_h Gather-scatter handle
  !! @param niter iteration trip count
  abstract interface
     function ksp_method_coupled(this, Ax, x, y, z, fx, fy, fz, &
          n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
       import :: bc_list_t
       import :: field_t
       import :: ksp_t
       import :: coef_t
       import :: gs_t
       import :: ax_t
       import :: ksp_monitor_t
       import rp
       implicit none
       class(ksp_t), intent(inout) :: this
       class(ax_t), intent(inout) :: Ax
       type(field_t), intent(inout) :: x
       type(field_t), intent(inout) :: y
       type(field_t), intent(inout) :: z
       integer, intent(in) :: n
       real(kind=rp), dimension(n), intent(inout) :: fx
       real(kind=rp), dimension(n), intent(inout) :: fy
       real(kind=rp), dimension(n), intent(inout) :: fz
       type(coef_t), intent(inout) :: coef
       type(bc_list_t), intent(inout) :: blstx
       type(bc_list_t), intent(inout) :: blsty
       type(bc_list_t), intent(inout) :: blstz
       type(gs_t), intent(inout) :: gs_h
       integer, optional, intent(in) :: niter
       type(ksp_monitor_t), dimension(3) :: ksp_results
     end function ksp_method_coupled
  end interface

  !> Abstract interface for deallocating a Krylov method
  abstract interface
     subroutine ksp_t_free(this)
       import :: ksp_t
       class(ksp_t), intent(inout) :: this
     end subroutine ksp_t_free
  end interface

  interface
     !> Factory for Krylov solvers. Both creates and initializes the object.
     !! @param object The object to be allocated.
     !! @param n Size of the vectors the solver operates on.
     !! @param type_name The name of the solver type.
     !! @param max_iter The maximum number of iterations
     !! @param abstol The absolute tolerance, optional.
     !! @param M The preconditioner, optional.
     !! @param monitor Enable/disable monitoring, optional.
     module subroutine krylov_solver_factory(object, n, type_name, &
          max_iter, abstol, M, monitor)
       class(ksp_t), allocatable, target, intent(inout) :: object
       integer, intent(in), value :: n
       character(len=*), intent(in) :: type_name
       integer, intent(in) :: max_iter
       real(kind=rp), optional :: abstol
       class(pc_t), optional, intent(inout), target :: M
       logical, optional, intent(in) :: monitor
     end subroutine krylov_solver_factory

     !> Destroy an iterative Krylov type_name
     module subroutine krylov_solver_destroy(object)
       class(ksp_t), allocatable, intent(inout) :: object
     end subroutine krylov_solver_destroy
  end interface

  public :: krylov_solver_factory, krylov_solver_destroy
contains

  !> Constructor for the base type.
  !! @param max_iter Maximum number of iterations.
  !! @param rel_tol Relative tolarance for converence.
  !! @param rel_tol Absolute tolarance for converence.
  !! @param M The preconditioner.
  subroutine krylov_init(this, max_iter, rel_tol, abs_tol, M, monitor)
    class(ksp_t), target, intent(inout) :: this
    integer, intent(in) :: max_iter
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    class(pc_t), optional, target, intent(in) :: M
    logical, optional, intent(in) :: monitor

    call krylov_free(this)

    if (present(rel_tol)) then
       this%rel_tol = rel_tol
    else
       this%rel_tol = KSP_REL_TOL
    end if

    if (present(abs_tol)) then
       this%abs_tol = abs_tol
    else
       this%abs_tol = KSP_ABS_TOL
    end if

    this%max_iter = max_iter

    if (present(M)) then
       this%M => M
    else
       if (.not. associated(this%M)) then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             allocate(device_ident_t::this%M_ident)
          else
             allocate(ident_t::this%M_ident)
          end if
          this%M => this%M_ident
       end if
    end if

    if (present(monitor)) then
       this%monitor = monitor
    else
       this%monitor = .false.
    end if

  end subroutine krylov_init

  !> Deallocate a Krylov solver
  subroutine krylov_free(this)
    class(ksp_t), intent(inout) :: this

    !> @todo add calls to destroy precon. if necessary

  end subroutine krylov_free

  !> Setup a Krylov solver's preconditioner.
  !! @param M The preconditioner.
  subroutine krylov_set_pc(this, M)
    class(ksp_t), intent(inout) :: this
    class(pc_t), target, intent(in) :: M

    if (associated(this%M)) then
       select type (pc => this%M)
       type is (ident_t)
       type is (device_ident_t)
       class default
          call neko_error('Preconditioner already defined')
       end select
    end if

    this%M => M

  end subroutine krylov_set_pc

  !> Monitor start
  subroutine krylov_monitor_start(this, name)
    class(ksp_t), intent(in) :: this
    character(len=*) :: name
    character(len=LOG_SIZE) :: log_buf
    
    if (this%monitor) then
       write(log_buf, '(A)') 'Krylov monitor (' // trim(name) // ')'
       call neko_log%section(trim(log_buf))
       call neko_log%newline()
       call neko_log%begin()
       write(log_buf, '(A)') ' Iter.       Residual'
       call neko_log%message(log_buf)
       write(log_buf, '(A)') '---------------------'
       call neko_log%message(log_buf)
    end if
  end subroutine krylov_monitor_start

  !> Monitor stop
  subroutine krylov_monitor_stop(this)
    class(ksp_t), intent(in) :: this

    if (this%monitor) then
       call neko_log%end()
       call neko_log%end_section()
       call neko_log%newline()
    end if
  end subroutine krylov_monitor_stop

  
  !> Monitor iteration
  subroutine krylov_monitor_iter(this, iter, rnorm)
    class(ksp_t), intent(in) :: this
    integer, intent(in) :: iter
    real(kind=rp), intent(in) :: rnorm
    character(len=LOG_SIZE) :: log_buf

    if (this%monitor) then
       write(log_buf, '(I6,E15.7)') iter, rnorm
       call neko_log%message(log_buf)
    end if
    
  end subroutine krylov_monitor_iter

end module krylov
