module fluid_abbdf_fctry
  use fluid_abbdf
  use fluid_abbdf_cpu
  implicit none

contains

  subroutine fluid_sumab_fctry(sumab)
    class(fluid_sumab_t), allocatable, intent(inout) :: sumab

    if (allocated(sumab)) then
       deallocate(sumab)
    end if

    allocate(fluid_sumab_cpu_t::sumab)
    
  end subroutine fluid_sumab_fctry

  subroutine fluid_makeabf_fctry(makeabf)
    class(fluid_makeabf_t), allocatable, intent(inout) :: makeabf

    if (allocated(makeabf)) then
       deallocate(makeabf)
    end if

    allocate(fluid_makeabf_cpu_t::makeabf)
    
  end subroutine fluid_makeabf_fctry

  subroutine fluid_makebdf_fctry(makebdf)
    class(fluid_makebdf_t), allocatable, intent(inout) :: makebdf

    if (allocated(makebdf)) then
       deallocate(makebdf)
    end if

    allocate(fluid_makebdf_cpu_t::makebdf)
    
  end subroutine fluid_makebdf_fctry
  
end module fluid_abbdf_fctry
