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
  use precon,  only : pc_t
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use field, only : field_t
  use utils, only : neko_error, neko_warning
  use bc, only : bc_list_t, bc_list_apply_vector, bc_list_apply_scalar, &
                 bc_list_apply
  use identity, only : ident_t
  use device_identity, only : device_ident_t
  use neko_config
  implicit none
  private

  integer, public, parameter :: KSP_MAX_ITER = 1e4 !< Maximum number of iters.
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
     class(pc_t), allocatable :: M_ident !< Internal preconditioner (Identity)
   contains
     procedure, pass(this) :: ksp_init => krylov_init
     procedure, pass(this) :: ksp_free => krylov_free
     procedure, pass(this) :: set_pc => krylov_set_pc
     procedure(ksp_method), pass(this), deferred :: solve
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
     function ksp_method(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
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

  !> Abstract interface for deallocating a Krylov method
  abstract interface
     subroutine ksp_t_free(this)
       import :: ksp_t
       class(ksp_t), intent(inout) :: this
     end subroutine ksp_t_free
  end interface
  
contains

  !> Create a krylov solver
  !! @param rel_tol Relative tolarance for converence.
  !! @param rel_tol Absolute tolarance for converence.
  !! @param M The preconditioner.
  subroutine krylov_init(this, rel_tol, abs_tol, M)    
    class(ksp_t), target, intent(inout) :: this
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    class(pc_t), optional, target, intent(in) :: M
    
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
       select type(pc => this%M)
       type is (ident_t)
       type is (device_ident_t)
       class default
          call neko_error('Preconditioner already defined')
       end select
    end if
    
    this%M => M
    
  end subroutine krylov_set_pc
  
end module krylov
