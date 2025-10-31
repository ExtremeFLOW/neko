! Copyright (c) 2022-2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
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
!> Defines a dong outflow condition
module dong_outflow
  use neko_config
  use dirichlet, only : dirichlet_t
  use device, only : device_memcpy, device_alloc, HOST_TO_DEVICE, device_free
  use num_types, only : rp, c_rp
  use bc, only : bc_t
  use field, only : field_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use utils, only : nonlinear_index
  use device_dong_outflow, only : device_dong_outflow_apply_scalar
  use field_registry, only : neko_field_registry
  use, intrinsic :: iso_c_binding, only : c_ptr, c_sizeof, c_null_ptr, &
       c_associated
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use utils, only : neko_error
  use time_state, only : time_state_t
  implicit none
  private

  !> Dong outflow condition
  !! Follows
  !! "A Convective-like Energy-Stable Open Boundary Condition for
  !! Simulations of Incompressible Flows"
  !! by S. Dong
  type, public, extends(bc_t) :: dong_outflow_t
     type(field_t), pointer :: u
     type(field_t), pointer :: v
     type(field_t), pointer :: w
     real(kind=rp) :: delta
     real(kind=rp) :: uinf
     type(c_ptr) :: normal_x_d = c_null_ptr
     type(c_ptr) :: normal_y_d = c_null_ptr
     type(c_ptr) :: normal_z_d = c_null_ptr
   contains
     procedure, pass(this) :: apply_scalar => dong_outflow_apply_scalar
     procedure, pass(this) :: apply_vector => dong_outflow_apply_vector
     procedure, pass(this) :: apply_scalar_dev => dong_outflow_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => dong_outflow_apply_vector_dev
     !> Constructor
     procedure, pass(this) :: init => dong_outflow_init
     !> Destructor.
     procedure, pass(this) :: free => dong_outflow_free
     !> Finalize.
     procedure, pass(this) :: finalize => dong_outflow_finalize
  end type dong_outflow_t

contains
  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine dong_outflow_init(this, coef, json)
    class(dong_outflow_t), target, intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    call this%free()
    call this%init_base(coef)

    call json_get_or_default(json, 'delta', this%delta, 0.01_rp)
    call json_get_or_default(json, 'velocity_scale', this%uinf, 1.0_rp)

  end subroutine dong_outflow_init

  !> Boundary condition apply for a generic Dirichlet condition
  !! to a vector @a x
  subroutine dong_outflow_apply_scalar(this, x, n, time, strong)
    class(dong_outflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    integer :: i, m, k, facet, idx(4)
    real(kind=rp) :: vn, S0, ux, uy, uz, normal_xyz(3)
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    !Im actually not sure what to do if one has two dong that share a corner.
    if (strong_) then
       m = this%msk(0)
       do i = 1, m
          k = this%msk(i)
          facet = this%facet(i)
          ux = this%u%x(k,1,1,1)
          uy = this%v%x(k,1,1,1)
          uz = this%w%x(k,1,1,1)
          idx = nonlinear_index(k, this%Xh%lx, this%Xh%lx, this%Xh%lx)
          normal_xyz = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), &
               facet)
          vn = ux*normal_xyz(1) + uy*normal_xyz(2) + uz*normal_xyz(3)
          S0 = 0.5_rp*(1.0_rp - tanh(vn / (this%uinf * this%delta)))

          x(k) = -0.5*(ux*ux+uy*uy+uz*uz)*S0
       end do
    end if
  end subroutine dong_outflow_apply_scalar

  !> Boundary condition apply for a generic Dirichlet condition
  !! to vectors @a x, @a y and @a z
  subroutine dong_outflow_apply_vector(this, x, y, z, n, time, strong)
    class(dong_outflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

  end subroutine dong_outflow_apply_vector

  !> Boundary condition apply for a generic Dirichlet condition
  !! to a vector @a x (device version)
  subroutine dong_outflow_apply_scalar_dev(this, x_d, time, strong, strm)
    class(dong_outflow_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_ .and. this%msk(0) .gt. 0) then
       call device_dong_outflow_apply_scalar(this%msk_d, x_d, &
            this%normal_x_d, this%normal_y_d, this%normal_z_d, &
            this%u%x_d, this%v%x_d, this%w%x_d, &
            this%uinf, this%delta, &
            this%msk(0), strm)
    end if

  end subroutine dong_outflow_apply_scalar_dev

  !> Boundary condition apply for a generic Dirichlet condition
  !! to vectors @a x, @a y and @a z (device version)
  subroutine dong_outflow_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(dong_outflow_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    !call device_dong_outflow_apply_vector(this%msk_d, x_d, y_d, z_d, &
    !                                   this%g, size(this%msk))

  end subroutine dong_outflow_apply_vector_dev

  !> Destructor
  subroutine dong_outflow_free(this)
    class(dong_outflow_t), target, intent(inout) :: this

    call this%free_base
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)

    if (c_associated(this%normal_x_d)) then
       call device_free(this%normal_x_d)
       this%normal_x_d = c_null_ptr
    end if

    if (c_associated(this%normal_y_d)) then
       call device_free(this%normal_y_d)
       this%normal_y_d = c_null_ptr
    end if

    if (c_associated(this%normal_z_d)) then
       call device_free(this%normal_z_d)
       this%normal_z_d = c_null_ptr
    end if

  end subroutine dong_outflow_free

  !> Finalize
  subroutine dong_outflow_finalize(this, only_facets)
    class(dong_outflow_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    real(kind=rp), allocatable :: temp_x(:)
    real(kind=rp), allocatable :: temp_y(:)
    real(kind=rp), allocatable :: temp_z(:)
    real(c_rp) :: dummy
    integer :: i, m, k, facet, idx(4)
    real(kind=rp) :: normal_xyz(3)

    if (present(only_facets)) then
       if (only_facets .eqv. .false.) then
          call neko_error("For dong_outflow_t, only_facets has to be true.")
       end if
    end if

    call this%finalize_base(.true.)

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")
    if ((NEKO_BCKND_DEVICE .eq. 1) .and. (this%msk(0) .gt. 0)) then
       call device_alloc(this%normal_x_d, c_sizeof(dummy)*this%msk(0))
       call device_alloc(this%normal_y_d, c_sizeof(dummy)*this%msk(0))
       call device_alloc(this%normal_z_d, c_sizeof(dummy)*this%msk(0))
       m = this%msk(0)
       allocate(temp_x(m))
       allocate(temp_y(m))
       allocate(temp_z(m))
       do i = 1, m
          k = this%msk(i)
          facet = this%facet(i)
          idx = nonlinear_index(k, this%Xh%lx, this%Xh%lx, this%Xh%lx)
          normal_xyz = &
               this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), facet)
          temp_x(i) = normal_xyz(1)
          temp_y(i) = normal_xyz(2)
          temp_z(i) = normal_xyz(3)
       end do
       call device_memcpy(temp_x, this%normal_x_d, m, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(temp_y, this%normal_y_d, m, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(temp_z, this%normal_z_d, m, HOST_TO_DEVICE, &
            sync = .true.)
       deallocate( temp_x, temp_y, temp_z)
    end if
  end subroutine dong_outflow_finalize

end module dong_outflow
