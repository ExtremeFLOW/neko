! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright
!     notice, this list of conditions and the following disclaimer in
!     the documentation and/or other materials provided with the
!     distribution.
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

  !> Wall sampler factory.
  !! @details If no `sampling` block is present, this routine falls back to
  !! legacy behaviour using `h_index` and the GLL sampler.
  !! @param object Allocatable wall sampler to create.
  !! @param json Case configuration data.
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
          call object%init_from_indices([h_index], .true.)
       end select
       return
    end if

    call json_get(json, 'sampling', sampling)
    call json_get(sampling, 'type', type_name)
    call wall_sampler_allocator(object, type_name)
    call object%init(sampling)
  end subroutine wall_sampler_factory

  !> Allocate a concrete wall sampler implementation by type name.
  !! @param object Allocatable wall sampler to allocate.
  !! @param type_name Sampler type identifier read from configuration.
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
