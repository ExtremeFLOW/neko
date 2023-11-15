! Copyright (c) 2021-2022, The Neko Authors
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
module krylov_fctry
  use cg, only : cg_t
  use cg_sx, only : sx_cg_t
  use cg_device, only : cg_device_t
  use cacg, only : cacg_t
  use pipecg, only : pipecg_t
  use pipecg_sx, only : sx_pipecg_t
  use pipecg_device, only : pipecg_device_t
  use bicgstab, only : bicgstab_t
  use gmres, only : gmres_t
  use gmres_sx, only : sx_gmres_t
  use gmres_device, only : gmres_device_t
  use num_Types, only : rp
  use krylov, only : ksp_t, ksp_monitor_t
  use precon, only : pc_t
  use utils, only : neko_error
  use neko_config
  implicit none

  public :: krylov_solver_factory, krylov_solver_destroy

contains

  !> Initialize an interative Krylov solver
  subroutine krylov_solver_factory(ksp, n, solver, abstol, M)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=*), intent(in) :: solver
    real(kind=rp), optional :: abstol
    class(pc_t), optional, intent(inout), target :: M
 
    if (allocated(ksp)) then
       call krylov_solver_destroy(ksp)
       deallocate(ksp)
    end if
    if (trim(solver) .eq. 'cg') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_cg_t::ksp)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(cg_device_t::ksp)
       else
          allocate(cg_t::ksp)
       end if
    else if (trim(solver) .eq. 'pipecg') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_pipecg_t::ksp)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          if (NEKO_BCKND_OPENCL .eq. 1) then
             call neko_error('PipeCG not supported for OpenCL')
          end if
          allocate(pipecg_device_t::ksp)
       else
          allocate(pipecg_t::ksp)
       end if
    else if (trim(solver) .eq. 'cacg') then
       allocate(cacg_t::ksp)
    else if (trim(solver) .eq. 'gmres') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_gmres_t::ksp)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(gmres_device_t::ksp)
       else
          allocate(gmres_t::ksp)
       end if
    else if (trim(solver) .eq. 'bicgstab') then
       allocate(bicgstab_t::ksp)
    else
       call neko_error('Unknown Krylov solver '//trim(solver))
    end if

    if (present(abstol) .and. present(M)) then
       select type(kp => ksp)
       type is(cg_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(sx_cg_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(cg_device_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(pipecg_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(sx_pipecg_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(pipecg_device_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(cacg_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(gmres_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(sx_gmres_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(gmres_device_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(bicgstab_t)
          call kp%init(n, M = M, abs_tol = abstol)
       end select
    else if (present(abstol)) then
       select type(kp => ksp)
       type is(cg_t)
          call kp%init(n, abs_tol = abstol)
       type is(sx_cg_t)
          call kp%init(n, abs_tol = abstol)       
       type is(cg_device_t)
          call kp%init(n, abs_tol = abstol)       
       type is(pipecg_t)
          call kp%init(n, abs_tol = abstol)
       type is(sx_pipecg_t)
          call kp%init(n, abs_tol = abstol)
       type is (pipecg_device_t)
          call kp%init(n, abs_tol = abstol)
       type is(cacg_t)
          call kp%init(n, abs_tol = abstol)
       type is(gmres_t)
          call kp%init(n, abs_tol = abstol)
       type is(sx_gmres_t)
          call kp%init(n, abs_tol = abstol)
       type is(gmres_device_t)
          call kp%init(n, abs_tol = abstol)
       type is(bicgstab_t)
          call kp%init(n, abs_tol = abstol)
       end select
    else if (present(M)) then
       select type(kp => ksp)
       type is(cg_t)
          call kp%init(n, M = M)
       type is(sx_cg_t)
          call kp%init(n, M = M)       
       type is(cg_device_t)
          call kp%init(n, M = M)
       type is(pipecg_t)
          call kp%init(n, M = M)
       type is(sx_pipecg_t)
          call kp%init(n, M = M)
       type is (pipecg_device_t)
          call kp%init(n, M = M)
       type is(cacg_t)
          call kp%init(n, M = M)
       type is(gmres_t)
          call kp%init(n, M = M)
       type is(sx_gmres_t)
          call kp%init(n, M = M)
       type is(gmres_device_t)
          call kp%init(n, M = M)
       type is(bicgstab_t)
          call kp%init(n, M = M)
       end select
    else
       select type(kp => ksp)
       type is(cg_t)
          call kp%init(n)
       type is(sx_cg_t)
          call kp%init(n)       
       type is(cg_device_t)
          call kp%init(n)       
       type is(pipecg_t)
          call kp%init(n)
       type is(sx_pipecg_t)
          call kp%init(n)
       type is (pipecg_device_t)
          call kp%init(n)
       type is(cacg_t)
          call kp%init(n)
       type is(gmres_t)
          call kp%init(n)
       type is(sx_gmres_t)
          call kp%init(n)
       type is(gmres_device_t)
          call kp%init(n)
       type is(bicgstab_t)
          call kp%init(n)
       end select
    end if

  end subroutine krylov_solver_factory

  !> Destroy an interative Krylov solver
  subroutine krylov_solver_destroy(ksp)
    class(ksp_t), allocatable, intent(inout) :: ksp

    if (allocated(ksp)) then
       select type(kp => ksp)
       type is(cg_t)
          call kp%free()
       type is(sx_cg_t)
          call kp%free()
       type is(cg_device_t)
          call kp%free()       
       type is(pipecg_t)
          call kp%free()
       type is(sx_pipecg_t)
          call kp%free()
       type is (pipecg_device_t)
          call kp%free()
       type is(cacg_t)
          call kp%free()
       type is(gmres_t)
          call kp%free()
       type is(sx_gmres_t)
          call kp%free()
       type is(gmres_device_t)
          call kp%free()
       type is(bicgstab_t)
          call kp%free()
       end select
    end if
 
  end subroutine krylov_solver_destroy
    
end module krylov_fctry
  
