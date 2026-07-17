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
!> Implements `wall_distance_sampler_t`.
module wall_distance_sampler
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use vector, only : vector_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use wall_sampler, only : wall_sampler_t
  use user_intf, only : user_t
  use global_interpolation, only : global_interpolation_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use device, only : HOST_TO_DEVICE
  implicit none
  private

  !> Implements wall-normal sampling at prescribed physical distances.
  type, public, extends(wall_sampler_t) :: wall_distance_sampler_t
     !> Physical distances from the wall for each sampling point.
     real(kind=rp), allocatable :: distances(:)
     !> Sampling point locations with shape `(3, n_nodes * n_samples)`.
     real(kind=rp), allocatable :: xyz(:,:)
     !> Interpolator for sampling the solution field at the sampling points.
     type(global_interpolation_t) :: interpolator
   contains
     !> Partial constructor from JSON.
     procedure, pass(this) :: init => wall_distance_sampler_init
     !> Partial constructor from physical wall-normal distances.
     procedure, pass(this) :: init_from_distances => &
          wall_distance_sampler_init_from_distances
     !> Partial constructor for user-provided sampling distances.
     procedure, private, pass(this) :: init_from_user => &
          wall_distance_sampler_init_from_user
     !> Fully construct a sampler from its resolved components.
     procedure, pass(this) :: init_from_components => &
          wall_distance_sampler_init_from_components
     !> Finalize by computing sampling point locations and building the
     !! interpolator.
     procedure, pass(this) :: finalize => wall_distance_sampler_finalize
     !> Sample the solution field at the sampling points.
     procedure, pass(this) :: sample => wall_distance_sampler_sample
     !> Destructor.
     procedure, pass(this) :: free => wall_distance_sampler_free
  end type wall_distance_sampler_t

contains

  !> Partial constructor from JSON.
  !! @param json Sampler configuration data.
  subroutine wall_distance_sampler_init(this, json)
    class(wall_distance_sampler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: distances(:)
    real(kind=rp) :: distance
    integer :: n_samples, var_type
    logical :: found, output_h
    character(len=:), allocatable :: value

    call json_get_or_default(json, 'output_h', output_h, .true.)

    call json%info('value', found = found, var_type = var_type)
    if (.not. found) then
       call neko_error('Wall distance sampler requires a value')
    end if

    select case (var_type)
    case (3) ! json_array
       call json_get(json, 'value', distances)
       call this%init_from_distances(distances, output_h)
    case (5, 6) ! json_integer, json_real
       call json_get(json, 'value', distance)
       call this%init_from_distances([distance], output_h)
    case (7) ! json_string
       call json_get(json, 'value', value)
       if (trim(value) /= 'user') then
          call neko_error('Wall distance sampler value must be a real,' // &
               ' array, or user')
       end if
       call json_get(json, 'n_samples', n_samples)
       call this%init_from_user(n_samples, output_h)
    case default
       call neko_error('Wall distance sampler value must be a real, array,' // &
            ' or user')
    end select
  end subroutine wall_distance_sampler_init

  !> Partial constructor from physical wall-normal distances.
  !! @param distances Positive distances, one for each sampling point per
  !! wall node.
  !! @param output_h Whether to write the sampling-distance diagnostic.
  subroutine wall_distance_sampler_init_from_distances(this, distances, &
       output_h)
    class(wall_distance_sampler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: distances(:)
    logical, intent(in) :: output_h

    if (size(distances) < 1 .or. any(distances <= 0.0_rp)) then
       call neko_error('Wall sampling distances must be positive')
    end if
    this%distances = distances
    this%n_samples = size(distances)
    this%user_values = .false.
    this%output_h_enabled = output_h
  end subroutine wall_distance_sampler_init_from_distances

  !> Partial constructor for user-provided sampling distances.
  !! @param n_samples Number of samples per wall node.
  !! @param output_h Whether to write the sampling-distance diagnostic.
  subroutine wall_distance_sampler_init_from_user(this, n_samples, output_h)
    class(wall_distance_sampler_t), intent(inout) :: this
    integer, intent(in) :: n_samples
    logical, intent(in) :: output_h

    if (n_samples < 1) then
       call neko_error('Wall distance sampler n_samples must be positive')
    end if

    this%n_samples = n_samples
    this%user_values = .true.
    this%output_h_enabled = output_h
  end subroutine wall_distance_sampler_init_from_user

  !> Fully construct a distance sampler from its resolved components.
  !! @note No compatibility checks are performed, input data is assumed to be
  !! valid. Checking would require geometric information.
  !! @param coef SEM coefficients used to construct global interpolation.
  !! @param distances Physical distances, one for each sample per wall node.
  !! @param xyz Sampling locations with shape `(3, n_nodes * n_samples)`.
  !! @param h Wall-normal distances in sampler ordering.
  !! @param output_h Sampling-distance diagnostic output flag.
  subroutine wall_distance_sampler_init_from_components(this, coef, distances, &
       xyz, h, output_h)
    class(wall_distance_sampler_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: distances(:)
    real(kind=rp), intent(in) :: xyz(:,:)
    type(vector_t), intent(in) :: h
    logical, intent(in) :: output_h
    integer :: n_nodes, n_points

    if (size(distances) < 1 .or. any(distances <= 0.0_rp)) then
       call neko_error('Wall sampling distances must be positive')
    end if
    if (size(xyz, 1) /= 3) then
       call neko_error('Wall sampler coordinates must have three components')
    end if

    n_points = size(xyz, 2)
    if (n_points < 1 .or. mod(n_points, size(distances)) /= 0) then
       call neko_error('Wall sampler coordinates have an invalid size')
    end if

    n_nodes = n_points / size(distances)
    call this%free()
    call this%init_base(n_nodes, size(distances), h, output_h)
    this%distances = distances
    this%xyz = xyz
    call this%interpolator%init(coef%dof)
    call this%interpolator%find_points(this%xyz, n_points)
  end subroutine wall_distance_sampler_init_from_components

  !> Construct sampling locations and initialise global interpolation.
  !! @param coef SEM coefficients.
  !! @param msk Mask selecting local wall nodes.
  !! @param facet Facet index for each selected wall node.
  !! @param n_x X-component of wall-normal vectors at wall nodes.
  !! @param n_y Y-component of wall-normal vectors at wall nodes.
  !! @param n_z Z-component of wall-normal vectors at wall nodes.
  !! @param bc_name Name of the owning boundary condition for user values.
  !! @param user User interface providing the distance-sampling callback.
  subroutine wall_distance_sampler_finalize(this, coef, msk, facet, &
       n_x, n_y, n_z, bc_name, user)
    class(wall_distance_sampler_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(0:)
    integer, intent(in) :: facet(0:)
    type(vector_t), intent(in) :: n_x, n_y, n_z
    character(len=*), optional, intent(in) :: bc_name
    type(user_t), target, optional, intent(in) :: user
    integer :: i, j, p, linear
    real(kind=rp) :: distance
    real(kind=rp), allocatable :: user_distances(:,:)

    if (.not. this%user_values .and. .not. allocated(this%distances)) then
       call neko_error('Wall distance sampler has not been initialized')
    end if

    call this%h%free()
    call this%interpolator%free()

    if (allocated(this%xyz)) deallocate(this%xyz)

    this%n_nodes = msk(0)
    if (.not. this%user_values) this%n_samples = size(this%distances)

    if (this%user_values) then
       if (.not. present(bc_name) .or. .not. present(user)) then
          call neko_error('Wall distance user sampler has no ' // &
               'configured callback')
       end if

       allocate(user_distances(this%n_samples, this%n_nodes))
       call user%wall_sampling_distance(bc_name, msk(1:this%n_nodes), &
            user_distances)
       if (any(user_distances <= 0.0_rp)) then
          call neko_error('Wall sampling distances must be positive')
       end if
    end if
    allocate(this%xyz(3, this%n_nodes * this%n_samples))
    call this%h%init(this%n_nodes * this%n_samples)

    do i = 1, this%n_nodes
       linear = msk(i)
       do j = 1, this%n_samples
          p = (i - 1) * this%n_samples + j
          if (this%user_values) then
             distance = user_distances(j,i)
          else
             distance = this%distances(j)
          end if
          this%xyz(1,p) = coef%dof%x(linear,1,1,1) - &
               distance * n_x%x(i)
          this%xyz(2,p) = coef%dof%y(linear,1,1,1) - &
               distance * n_y%x(i)
          this%xyz(3,p) = coef%dof%z(linear,1,1,1) - &
               distance * n_z%x(i)
          this%h%x(p) = distance
       end do
    end do

    call this%interpolator%init(coef%dof)
    call this%interpolator%find_points(this%xyz, &
         this%n_nodes * this%n_samples)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%h%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if

    if (this%output_h_enabled .and. present(bc_name)) then
       call this%output_h(coef, msk, bc_name)
    end if
  end subroutine wall_distance_sampler_finalize

  !> Sample a field at the configured physical locations.
  !! @param field Field to be sampled.
  !! @param values Output sampled values in sampler ordering.
  subroutine wall_distance_sampler_sample(this, field, values)
    class(wall_distance_sampler_t), intent(inout) :: this
    type(field_t), intent(inout) :: field
    type(vector_t), intent(inout) :: values
    integer :: n

    n = this%n_nodes * this%n_samples
    if (values%size() /= n) then
       call neko_error('Wall sampler destination has an invalid size')
    end if

    call this%interpolator%evaluate(values%x, field%x, &
         NEKO_BCKND_DEVICE .eq. 0)
  end subroutine wall_distance_sampler_sample

  !> Release sampler resources.
  subroutine wall_distance_sampler_free(this)
    class(wall_distance_sampler_t), intent(inout) :: this

    call this%interpolator%free()
    call this%h%free()
    if (allocated(this%distances)) deallocate(this%distances)
    if (allocated(this%xyz)) deallocate(this%xyz)
    this%n_nodes = 0
    this%n_samples = 0
    this%user_values = .false.
    this%output_h_enabled = .true.
  end subroutine wall_distance_sampler_free

end module wall_distance_sampler
