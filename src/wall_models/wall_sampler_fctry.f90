! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
!> Factory for wall-model samplers.
module wall_sampler_fctry
  use json_module, only : json_file
  use json_utils, only : json_get
  use wall_sampler, only : wall_sampler_t
  use wall_gll_sampler, only : wall_gll_sampler_t
  use wall_distance_sampler, only : wall_distance_sampler_t
  use utils, only : neko_error
  implicit none
  private
  public :: wall_sampler_factory, wall_sampler_allocator

contains

  subroutine wall_sampler_factory(object, json)
    class(wall_sampler_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(json_file) :: sampling
    character(len=:), allocatable :: type_name
    integer :: h_index

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    if (.not. json%valid_path('sampling')) then
       call json_get(json, 'h_index', h_index)
       call wall_sampler_allocator(object, 'gll')
       select type (object)
       type is (wall_gll_sampler_t)
          call object%init_from_indices([h_index])
       end select
       return
    end if

    call json_get(json, 'sampling', sampling)
    call json_get(sampling, 'type', type_name)
    call wall_sampler_allocator(object, type_name)
    call object%init(sampling)
  end subroutine wall_sampler_factory

  subroutine wall_sampler_allocator(object, type_name)
    class(wall_sampler_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if
    select case (trim(type_name))
    case ('gll')
       allocate(wall_gll_sampler_t :: object)
    case ('distance')
       allocate(wall_distance_sampler_t :: object)
    case default
       call neko_error('Unknown wall sampler type: ' // trim(type_name))
    end select
  end subroutine wall_sampler_allocator

end module wall_sampler_fctry
