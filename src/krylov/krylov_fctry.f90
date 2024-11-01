! Copyright (c) 2021-2024, The Neko Authors
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
  use utils, only : concat_string_array
  use neko_config, only : NEKO_BCKND_SX, NEKO_BCKND_OPENCL
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: KSP_KNOWN_TYPES(8) = [character(len=20) :: &
     "cg", &
     "pipecg", &
     "fusedcg", &
     "cacg", &
     "gmres", &
     "cheby", &
     "bicgstab", &
     "cpldcg"]

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
    class(ksp_t), allocatable, target, intent(inout) :: object
    integer, intent(in), value :: n
    character(len=*), intent(in) :: type_name
    integer, intent(in) :: max_iter
    real(kind=rp), optional :: abstol
    class(pc_t), optional, intent(inout), target :: M
    logical, optional, intent(in) :: monitor
    character(len=:), allocatable :: type_string

    if (allocated(object)) then
       call krylov_solver_destroy(object)
       deallocate(object)
    end if

    if (trim(type_name) .eq. 'cg') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_cg_t::object)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(cg_device_t::object)
       else
          allocate(cg_t::object)
       end if
    else if (trim(type_name) .eq. 'cpldcg') then
       allocate(cg_cpld_t::object)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call neko_error('Coupled CG only supported for CPU')
       end if
    else if (trim(type_name) .eq. 'pipecg') then
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
    else if (trim(type_name) .eq. 'fusedcg') then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          if (NEKO_BCKND_OPENCL .eq. 1) then
             call neko_error('FusedCG not supported for OpenCL')
          end if
          allocate(fusedcg_device_t::object)
       else
          call neko_error('FusedCG only supported for CUDA/HIP')
       end if
    else if (trim(type_name) .eq. 'fcpldcg') then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          if (NEKO_BCKND_OPENCL .eq. 1) then
             call neko_error('Coupled FusedCG not supported for OpenCL')
          end if
          allocate(fusedcg_cpld_device_t::object)
       else
          call neko_error('Coupled FusedCG only supported for CUDA/HIP')
       end if
    else if (trim(type_name) .eq. 'cacg') then
       allocate(cacg_t::object)
    else if (trim(type_name) .eq. 'gmres') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_gmres_t::object)
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(gmres_device_t::object)
       else
          allocate(gmres_t::object)
       end if
    else if (trim(type_name) .eq. 'cheby') then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          allocate(cheby_device_t::object)
       else
          allocate(cheby_t::object)
       end if
    else if (trim(type_name) .eq. 'bicgstab') then
       allocate(bicgstab_t::object)
    else
       type_string =  concat_string_array(KSP_KNOWN_TYPES,&
            NEW_LINE('A') // "-  ",  .true.)
       call neko_error("Unknown Krylov solver type: " &
                       // trim(type_name) // ". Known types are: " &
                       // type_string)
    end if

    ! This select type is in principle not necessary,but we have it due to
    ! issues with compilers, when it was not there. However, at some point we
    ! should check if we can get away with just having one obj%init statement.
    ! Same applies to the code in the "destroy" routine below.
    if (present(abstol) .and. present(M) .and. present(monitor)) then
       select type (obj => object)
       type is (cg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (sx_cg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (cg_cpld_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (cg_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (pipecg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (sx_pipecg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (pipecg_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (fusedcg_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (fusedcg_cpld_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (cacg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (gmres_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (sx_gmres_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (gmres_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (bicgstab_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (cheby_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       type is (cheby_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol, monitor = monitor)
       end select
    else if (present(abstol) .and. present(M)) then
       select type (obj => object)
       type is (cg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (sx_cg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (cg_cpld_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (cg_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (pipecg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (sx_pipecg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (pipecg_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (fusedcg_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (fusedcg_cpld_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (cacg_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (gmres_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (sx_gmres_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (gmres_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (bicgstab_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (cheby_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       type is (cheby_device_t)
          call obj%init(n, max_iter, M = M, abs_tol = abstol)
       end select
    else if (present(monitor) .and. present(M)) then
       select type (obj => object)
       type is (cg_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (sx_cg_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (cg_cpld_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (cg_device_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (pipecg_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (sx_pipecg_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (pipecg_device_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (fusedcg_device_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (fusedcg_cpld_device_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (cacg_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (gmres_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (sx_gmres_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (gmres_device_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (bicgstab_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (cheby_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       type is (cheby_device_t)
          call obj%init(n, max_iter, M = M, monitor = monitor)
       end select
    else if (present(abstol) .and. present(monitor)) then
          select type (obj => object)
       type is (cg_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (sx_cg_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (cg_cpld_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (cg_device_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (pipecg_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (sx_pipecg_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (pipecg_device_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (fusedcg_device_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (fusedcg_cpld_device_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (cacg_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (gmres_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (sx_gmres_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (gmres_device_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (bicgstab_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (cheby_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       type is (cheby_device_t)
          call obj%init(n, max_iter, monitor = monitor, abs_tol = abstol)
       end select
    else if (present(abstol)) then
       select type (obj => object)
       type is (cg_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (sx_cg_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (cg_cpld_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (cg_device_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (pipecg_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (sx_pipecg_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (pipecg_device_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (fusedcg_device_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (fusedcg_cpld_device_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (cacg_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (gmres_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (sx_gmres_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (gmres_device_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (bicgstab_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (cheby_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       type is (cheby_device_t)
          call obj%init(n, max_iter, abs_tol = abstol)
       end select
    else if (present(monitor)) then
       select type (obj => object)
       type is (cg_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (sx_cg_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (cg_cpld_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (cg_device_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (pipecg_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (sx_pipecg_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (pipecg_device_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (fusedcg_device_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (fusedcg_cpld_device_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (cacg_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (gmres_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (sx_gmres_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (gmres_device_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (bicgstab_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (cheby_t)
          call obj%init(n, max_iter, monitor = monitor)
       type is (cheby_device_t)
          call obj%init(n, max_iter, monitor = monitor)
       end select
    else if (present(M)) then
       select type (obj => object)
       type is (cg_t)
          call obj%init(n, max_iter, M = M)
       type is (sx_cg_t)
          call obj%init(n, max_iter, M = M)
       type is (cg_cpld_t)
          call obj%init(n, max_iter, M = M)
       type is (cg_device_t)
          call obj%init(n, max_iter, M = M)
       type is (pipecg_t)
          call obj%init(n, max_iter, M = M)
       type is (sx_pipecg_t)
          call obj%init(n, max_iter, M = M)
       type is (pipecg_device_t)
          call obj%init(n, max_iter, M = M)
       type is (fusedcg_device_t)
          call obj%init(n, max_iter, M = M)
       type is (fusedcg_cpld_device_t)
          call obj%init(n, max_iter, M = M)
       type is (cacg_t)
          call obj%init(n, max_iter, M = M)
       type is (gmres_t)
          call obj%init(n, max_iter, M = M)
       type is (sx_gmres_t)
          call obj%init(n, max_iter, M = M)
       type is (gmres_device_t)
          call obj%init(n, max_iter, M = M)
       type is (bicgstab_t)
          call obj%init(n, max_iter, M = M)
       type is (cheby_t)
          call obj%init(n, max_iter, M = M)
       type is (cheby_device_t)
          call obj%init(n, max_iter, M = M)
       end select
    else
       select type (obj => object)
       type is (cg_t)
          call obj%init(n, max_iter)
       type is (sx_cg_t)
          call obj%init(n, max_iter)
       type is (cg_cpld_t)
          call obj%init(n, max_iter)
       type is (cg_device_t)
          call obj%init(n, max_iter)
       type is (pipecg_t)
          call obj%init(n, max_iter)
       type is (sx_pipecg_t)
          call obj%init(n, max_iter)
       type is (pipecg_device_t)
          call obj%init(n, max_iter)
       type is (fusedcg_device_t)
          call obj%init(n, max_iter)
       type is (fusedcg_cpld_device_t)
          call obj%init(n, max_iter)
       type is (cacg_t)
          call obj%init(n, max_iter)
       type is (gmres_t)
          call obj%init(n, max_iter)
       type is (sx_gmres_t)
          call obj%init(n, max_iter)
       type is (gmres_device_t)
          call obj%init(n, max_iter)
       type is (bicgstab_t)
          call obj%init(n, max_iter)
       type is (cheby_t)
          call obj%init(n, max_iter)
       type is (cheby_device_t)
          call obj%init(n, max_iter)
       end select
    end if

  end subroutine krylov_solver_factory

  !> Destroy an iterative Krylov type_name
  module subroutine krylov_solver_destroy(object)
    class(ksp_t), allocatable, intent(inout) :: object

    if (allocated(object)) then
       select type (obj => object)
       type is (cg_t)
          call obj%free()
       type is (sx_cg_t)
          call obj%free()
       type is (cg_cpld_t)
          call obj%free()
       type is (cg_device_t)
          call obj%free()
       type is (pipecg_t)
          call obj%free()
       type is (sx_pipecg_t)
          call obj%free()
       type is (pipecg_device_t)
          call obj%free()
       type is (fusedcg_device_t)
          call obj%free()
       type is (fusedcg_cpld_device_t)
          call obj%free()
       type is (cacg_t)
          call obj%free()
       type is (gmres_t)
          call obj%free()
       type is (sx_gmres_t)
          call obj%free()
       type is (gmres_device_t)
          call obj%free()
       type is (bicgstab_t)
          call obj%free()
       type is (cheby_t)
          call obj%free()
       end select
    end if

  end subroutine krylov_solver_destroy

end submodule krylov_fctry

