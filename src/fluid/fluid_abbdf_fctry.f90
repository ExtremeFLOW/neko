module fluid_abbdf_fctry
  use fluid_abbdf
  use fluid_abbdf_cpu
  use fluid_abbdf_sx
  use fluid_abbdf_device
  use neko_config
  implicit none

contains

  subroutine fluid_sumab_fctry(sumab)
    class(fluid_sumab_t), allocatable, intent(inout) :: sumab

    if (allocated(sumab)) then
       deallocate(sumab)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(fluid_sumab_sx_t::sumab)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &         
         (NEKO_BCKND_OPENCL .eq. 1)) then
       allocate(fluid_sumab_device_t::sumab)
    else
       allocate(fluid_sumab_cpu_t::sumab)
    end if
    
  end subroutine fluid_sumab_fctry

  subroutine fluid_makeabf_fctry(makeabf)
    class(fluid_makeabf_t), allocatable, intent(inout) :: makeabf

    if (allocated(makeabf)) then
       deallocate(makeabf)
    end if
    
    if (NEKO_BCKND_SX .eq. 1) then
       allocate(fluid_makeabf_sx_t::makeabf)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       allocate(fluid_makeabf_device_t::makeabf)
    else
       allocate(fluid_makeabf_cpu_t::makeabf)
    end if
    
  end subroutine fluid_makeabf_fctry

  subroutine fluid_makebdf_fctry(makebdf)
    class(fluid_makebdf_t), allocatable, intent(inout) :: makebdf

    if (allocated(makebdf)) then
       deallocate(makebdf)
    end if

    if (NEKO_BCKND_SX .eq. 1) then
       allocate(fluid_makebdf_sx_t::makebdf)
    else if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
         (NEKO_BCKND_OPENCL .eq. 1)) then
       allocate(fluid_makebdf_device_t::makebdf)
    else       
       allocate(fluid_makebdf_cpu_t::makebdf)
    end if
    
  end subroutine fluid_makebdf_fctry
  
end module fluid_abbdf_fctry
