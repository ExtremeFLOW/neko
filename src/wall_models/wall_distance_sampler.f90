! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
!> Implements wall-normal sampling at prescribed physical distances.
module wall_distance_sampler
  use num_types, only : rp
  use field, only : field_t
  use coefs, only : coef_t
  use vector, only : vector_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use wall_sampler, only : wall_sampler_t
  use global_interpolation, only : global_interpolation_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use device, only : HOST_TO_DEVICE
  implicit none
  private

  type, public, extends(wall_sampler_t) :: wall_distance_sampler_t
     real(kind=rp), allocatable :: distances(:)
     real(kind=rp), allocatable :: xyz(:,:)
     type(global_interpolation_t) :: interpolator
   contains
     procedure, pass(this) :: init => wall_distance_sampler_init
     procedure, pass(this) :: init_from_distances => &
          wall_distance_sampler_init_from_distances
     procedure, pass(this) :: finalize => wall_distance_sampler_finalize
     procedure, pass(this) :: sample => wall_distance_sampler_sample
     procedure, pass(this) :: free => wall_distance_sampler_free
  end type wall_distance_sampler_t

contains

  subroutine wall_distance_sampler_init(this, json)
    class(wall_distance_sampler_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: distances(:)

    call json_get(json, 'distances', distances)
    call this%init_from_distances(distances)
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
       n_x, n_y, n_z)
    class(wall_distance_sampler_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(0:)
    integer, intent(in) :: facet(0:)
    type(vector_t), intent(in) :: n_x, n_y, n_z
    integer :: i, j, p, linear

    if (.not. allocated(this%distances)) then
       call neko_error('Wall distance sampler has not been initialized')
    end if

    call this%h%free()
    call this%interpolator%free()
    if (allocated(this%xyz)) deallocate(this%xyz)
    this%n_nodes = msk(0)
    this%n_samples = size(this%distances)
    allocate(this%xyz(3, this%n_nodes * this%n_samples))
    call this%h%init(this%n_nodes * this%n_samples)

    do i = 1, this%n_nodes
       linear = msk(i)
       do j = 1, this%n_samples
          p = (i - 1) * this%n_samples + j
          this%xyz(1,p) = coef%dof%x(linear,1,1,1) - &
               this%distances(j) * n_x%x(i)
          this%xyz(2,p) = coef%dof%y(linear,1,1,1) - &
               this%distances(j) * n_y%x(i)
          this%xyz(3,p) = coef%dof%z(linear,1,1,1) - &
               this%distances(j) * n_z%x(i)
          this%h%x(p) = this%distances(j)
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
  end subroutine wall_distance_sampler_free

end module wall_distance_sampler
