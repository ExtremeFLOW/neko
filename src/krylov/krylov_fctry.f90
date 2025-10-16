! Copyright (c) 2021-2025, The Neko Authors
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
submodule (krylov) krylov_fctry
  use cg, only : cg_t
  use cg_sx, only : sx_cg_t
  use cg_cpld, only : cg_cpld_t
  use cg_device, only : cg_device_t
  use cg_cpld_device, only : cg_cpld_device_t
  use cacg, only : cacg_t
  use pipecg, only : pipecg_t
  use pipecg_sx, only : sx_pipecg_t
  use pipecg_device, only : pipecg_device_t
  use fusedcg_device, only : fusedcg_device_t
  use fusedcg_cpld_device, only : fusedcg_cpld_device_t
  use bicgstab, only : bicgstab_t
  use gmres, only : gmres_t
  use cheby, only : cheby_t
  use cheby_device, only : cheby_device_t
  use gmres_sx, only : sx_gmres_t
  use gmres_device, only : gmres_device_t
  use num_Types, only : rp
  use precon, only : pc_t
  use utils, only : neko_type_error
  use neko_config, only : NEKO_BCKND_SX, NEKO_BCKND_OPENCL
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: KSP_KNOWN_TYPES(9) = [character(len=20) :: &
       "cg", &
       "pipecg", &
       "fused_cg", &
       "cacg", &
       "gmres", &
       "cheby", &
       "bicgstab", &
       "fused_coupled_cg", &
       "coupled_cg"]

contains

  !> Factory for Krylov solvers. Both creates and initializes the object.
  !! @param object The object to be allocated.
  !! @param n Size of the vectors the solver operates on.
  !! @param type_name The name of the solver type.
  !! @param max_iter The maximum number of iterations
  !! @param abstol The absolute tolerance, optional.
  !! @param M The preconditioner, optional.
  !! @param monitor Enable/disable residual history, optional.
  module subroutine krylov_solver_factory(object, n, type_name, &
       max_iter, abstol, M, monitor)
    class(ksp_t), allocatable, intent(inout) :: object
    integer, intent(in), value :: n
    character(len=*), intent(in) :: type_name
    integer, intent(in) :: max_iter
    real(kind=rp), optional :: abstol
    class(pc_t), optional, intent(in), target :: M
    logical, optional, intent(in) :: monitor

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ('cg')
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_cg_t::object)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(cg_device_t::object)
       else
          allocate(cg_t::object)
       end if

    case ('coupled_cg')
       if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(cg_cpld_device_t::object)
       else
          allocate(cg_cpld_t::object)
       end if

    case ('pipecg')
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_pipecg_t::object)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          if (NEKO_BCKND_OPENCL .eq. 1) then
             call neko_error('PipeCG not supported for OpenCL')
          end if
          allocate(pipecg_device_t::object)
       else
          allocate(pipecg_t::object)
       end if

    case ('fused_cg')
       if (NEKO_BCKND_DEVICE .eq. 1) then
          if (NEKO_BCKND_OPENCL .eq. 1) then
             call neko_error('FusedCG not supported for OpenCL')
          end if
          allocate(fusedcg_device_t::object)
       else
          call neko_error('FusedCG only supported for CUDA/HIP')
       end if

    case ('fused_coupled_cg')
       if (NEKO_BCKND_DEVICE .eq. 1) then
          if (NEKO_BCKND_OPENCL .eq. 1) then
             call neko_error('Coupled FusedCG not supported for OpenCL')
          end if
          allocate(fusedcg_cpld_device_t::object)
       else
          call neko_error('Coupled FusedCG only supported for CUDA/HIP')
       end if

    case ('cacg')
       allocate(cacg_t::object)

    case ('gmres')
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_gmres_t::object)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(gmres_device_t::object)
       else
          allocate(gmres_t::object)
       end if

    case ('cheby')
       if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(cheby_device_t::object)
       else
          allocate(cheby_t::object)
       end if

    case ('bicgstab')
       allocate(bicgstab_t::object)

    case default
       call neko_type_error('Krylov solver', type_name, KSP_KNOWN_TYPES)
    end select

    call object%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)

  end subroutine krylov_solver_factory

end submodule krylov_fctry
