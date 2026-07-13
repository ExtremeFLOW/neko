! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
!> Implements `wall_distance_sampler_t`.
module wall_distance_sampler
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use vector, only : vector_t
  use json_module, only : json_file
  use json_utils, only : json_get
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
     !> Sampling point locations.
     real(kind=rp), allocatable :: xyz(:,:)
     !> Interpolator for sampling the solution field at the sampling points.
     type(global_interpolation_t) :: interpolator
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => wall_distance_sampler_init
     !> Constructor from physical distances.
     procedure, pass(this) :: init_from_distances => &
          wall_distance_sampler_init_from_distances
     !> Finalize by computing sampling point locations and building the
     !! interpolator.
     procedure, pass(this) :: finalize => wall_distance_sampler_finalize
     !> Sample the solution field at the sampling points.
     procedure, pass(this) :: sample => wall_distance_sampler_sample
     !> Destructor.
     procedure, pass(this) :: free => wall_distance_sampler_free
  end type wall_distance_sampler_t

contains

  subroutine wall_distance_sampler_init(this, json)
    class(wall_distance_sampler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: distances(:)
    real(kind=rp) :: distance
    integer :: var_type
    logical :: found
    character(len=:), allocatable :: value

    call json%info('value', found = found, var_type = var_type)
    if (.not. found) then
       call neko_error('Wall distance sampler requires a value')
    end if

    select case (var_type)
    case (3) ! json_array
       call json_get(json, 'value', distances)
       call this%init_from_distances(distances)
    case (5, 6) ! json_integer, json_real
       call json_get(json, 'value', distance)
       call this%init_from_distances([distance])
    case (7) ! json_string
       call json_get(json, 'value', value)
       if (trim(value) /= 'user') then
          call neko_error('Wall distance sampler value must be a real,' // &
               ' array, or user')
       end if
       call json_get(json, 'n_samples', this%n_samples)
       if (this%n_samples < 1) then
          call neko_error('Wall distance sampler n_samples must be positive')
       end if
       this%user_values = .true.
    case default
       call neko_error('Wall distance sampler value must be a real, array,' // &
            ' or user')
    end select
  end subroutine wall_distance_sampler_init

  subroutine wall_distance_sampler_init_from_distances(this, distances)
    class(wall_distance_sampler_t), intent(inout) :: this
    real(kind=rp), intent(in) :: distances(:)

    if (size(distances) < 1 .or. any(distances <= 0.0_rp)) then
       call neko_error('Wall sampling distances must be positive')
    end if
    this%distances = distances
    this%n_samples = size(distances)
  end subroutine wall_distance_sampler_init_from_distances

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
  end subroutine wall_distance_sampler_finalize

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

  subroutine wall_distance_sampler_free(this)
    class(wall_distance_sampler_t), intent(inout) :: this

    call this%interpolator%free()
    call this%h%free()
    if (allocated(this%distances)) deallocate(this%distances)
    if (allocated(this%xyz)) deallocate(this%xyz)
    this%n_nodes = 0
    this%n_samples = 0
    this%user_values = .false.
  end subroutine wall_distance_sampler_free

end module wall_distance_sampler
