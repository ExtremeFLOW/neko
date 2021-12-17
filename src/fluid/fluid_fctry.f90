!> Factory for all fluid schemes
module fluid_fctry
  use device_fluid_plan4
  use fluid_plan1
  use fluid_plan4
  use fluid_method
  use neko_config
  use utils
  implicit none

contains

  !> Initialise a fluid scheme
  subroutine fluid_scheme_factory(fluid, fluid_scheme)
    class(fluid_scheme_t), intent(inout), allocatable :: fluid
    character(len=*) :: fluid_scheme

    if (trim(fluid_scheme) .eq. 'plan1') then
       allocate(fluid_plan1_t::fluid)
    else if (trim(fluid_scheme) .eq. 'plan4') then
       if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
            (NEKO_BCKND_OPENCL .eq. 1)) then
          allocate(device_fluid_plan4_t::fluid)
       else
          allocate(fluid_plan4_t::fluid)
       end if
    else
       call neko_error('Invalid fluid scheme')
    end if
    
  end subroutine fluid_scheme_factory

end module fluid_fctry
