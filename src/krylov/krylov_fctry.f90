module krylov_fctry
  use cg
  use cg_sx
  use cg_device
  use cacg
  use pipecg
  use pipecg_sx
  use bicgstab
  use gmres
  use gmres_sx
  use krylov
  use neko_config
  implicit none
  
contains

  !> Initialize an interative Krylov solver
  subroutine krylov_solver_factory(ksp, n, solver, abstol, M)
    class(ksp_t), allocatable, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=*) :: solver
    real(kind=rp), optional :: abstol
    class(pc_t), optional, intent(inout), target :: M
 
    if (allocated(ksp)) then
       call krylov_solver_destroy(ksp)
       deallocate(ksp)
    end if

    if (trim(solver) .eq. 'cg') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_cg_t::ksp)
       else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1)) then
          allocate(cg_device_t::ksp)
       else
          allocate(cg_t::ksp)
       end if
    else if (trim(solver) .eq. 'pipecg') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_pipecg_t::ksp)
       else
          allocate(pipecg_t::ksp)
       end if
    else if (trim(solver) .eq. 'cacg') then
       allocate(cacg_t::ksp)
    else if (trim(solver) .eq. 'gmres') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_gmres_t::ksp)
       else
          allocate(gmres_t::ksp)
       end if
    else if (trim(solver) .eq. 'bicgstab') then
       allocate(bicgstab_t::ksp)
    else
       call neko_error('Unknown Krylov solver')
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
       type is(cacg_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(gmres_t)
          call kp%init(n, M = M, abs_tol = abstol)
       type is(sx_gmres_t)
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
       type is(cacg_t)
          call kp%init(n, abs_tol = abstol)
       type is(gmres_t)
          call kp%init(n, abs_tol = abstol)
       type is(sx_gmres_t)
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
       type is(cacg_t)
          call kp%init(n, M = M)
       type is(gmres_t)
          call kp%init(n, M = M)
       type is(sx_gmres_t)
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
       type is(cacg_t)
          call kp%init(n)
       type is(gmres_t)
          call kp%init(n)
       type is(sx_gmres_t)
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
       type is(pipecg_t)
          call kp%free()
       type is(sx_pipecg_t)
          call kp%free()
       type is(cacg_t)
          call kp%free()
       type is(gmres_t)
          call kp%free()
       type is(sx_gmres_t)
          call kp%free()
       type is(bicgstab_t)
          call kp%free()
       end select

       call ksp%free()

    end if
 
  end subroutine krylov_solver_destroy
    
end module krylov_fctry
  
