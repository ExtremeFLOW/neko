module precon_fctry
  use precon
  use jacobi
  use sx_jacobi
  use device_jacobi
  use hsmg
  use utils
  use neko_config
  implicit none

contains
  
  !> Create a preconditioner
  subroutine precon_factory(pc, pctype)
    class(pc_t), allocatable, intent(inout) :: pc
    character(len=*) :: pctype

    if (allocated(pc)) then
       call precon_destroy(pc)
       deallocate(pc)
    end if

    if (trim(pctype) .eq. 'jacobi') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_jacobi_t::pc)
       else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
            (NEKO_BCKND_OPENCL .eq. 1)) then
          allocate(device_jacobi_t::pc)
       else
          allocate(jacobi_t::pc)
       end if
    else if (pctype(1:4) .eq. 'hsmg') then
       allocate(hsmg_t::pc)
    else
       call neko_error('Unknown preconditioner')
    end if
    
  end subroutine precon_factory

  !> Destroy a preconditioner
  subroutine precon_destroy(pc)
    class(pc_t), allocatable, intent(inout) :: pc

    if (allocated(pc)) then
       select type(pcp => pc)
       type is(jacobi_t)
          call pcp%free()
       type is(sx_jacobi_t)
          call pcp%free()
       type is(device_jacobi_t)
          call pcp%free()
       type is (hsmg_t)
          call pcp%free()
       end select                 
    end if
    
  end subroutine precon_destroy
  
end module precon_fctry
