! Copyright (c) 2020-2025, The Neko Authors
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
!> Axis-aligned mixed Dirichlet-Neumann symmetry plane.
module symmetry_aligned
  use device_symmetry_aligned, only : device_symmetry_aligned_apply_vector
  use num_types, only : rp
  use bc, only : bc_t, BC_TYPES
  use tuple, only : tuple_i4_t
  use coefs, only : coef_t
  use json_module, only : json_file
  use zero_dirichlet, only : zero_dirichlet_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  use time_state, only : time_state_t
  implicit none
  private

  type, public, extends(bc_t) :: symmetry_aligned_t
     type(zero_dirichlet_t) :: bc_x
     type(zero_dirichlet_t) :: bc_y
     type(zero_dirichlet_t) :: bc_z
   contains
     procedure, pass(this) :: apply_scalar => symmetry_aligned_apply_scalar
     procedure, pass(this) :: apply_vector => symmetry_aligned_apply_vector
     procedure, pass(this) :: apply_scalar_dev => symmetry_aligned_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => symmetry_aligned_apply_vector_dev
     procedure, pass(this) :: init => symmetry_aligned_init
     procedure, pass(this) :: init_from_components => &
          symmetry_aligned_init_from_components
     procedure, pass(this) :: free => symmetry_aligned_free
     procedure, pass(this) :: get_normal_axis => symmetry_aligned_get_normal_axis
     procedure, pass(this) :: finalize => symmetry_aligned_finalize
  end type symmetry_aligned_t

contains

  subroutine symmetry_aligned_init(this, coef, json)
    class(symmetry_aligned_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%init_from_components(coef)
  end subroutine symmetry_aligned_init

  subroutine symmetry_aligned_init_from_components(this, coef)
    class(symmetry_aligned_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef

    call this%free()

    call this%init_base(coef)
    this%constraints = (/ .true., .false., .false. /)
    this%bc_type = BC_TYPES%MIXED_CONSTRAINS_NORMAL
    call this%bc_x%init_from_components(this%coef)
    call this%bc_y%init_from_components(this%coef)
    call this%bc_z%init_from_components(this%coef)
  end subroutine symmetry_aligned_init_from_components

  subroutine symmetry_aligned_finalize(this)
    class(symmetry_aligned_t), target, intent(inout) :: this
    integer :: i, facet, el
    type(tuple_i4_t), pointer :: bfp(:)
    real(kind=rp) :: sx, sy, sz
    real(kind=rp), parameter :: TOL = 1e-3_rp
    type(tuple_i4_t) :: bc_facet

    call this%finalize_base()

    associate(c => this%coef, nx => this%coef%nx, ny => this%coef%ny, &
         nz => this%coef%nz)
      bfp => this%marked_facet%array()
      do i = 1, this%marked_facet%size()
         bc_facet = bfp(i)
         facet = bc_facet%x(1)
         el = bc_facet%x(2)
         call this%get_normal_axis(sx, sy, sz, facet, el)

         if (sx .lt. TOL) call this%bc_x%mark_facet(facet, el)
         if (sy .lt. TOL) call this%bc_y%mark_facet(facet, el)
         if (sz .lt. TOL) call this%bc_z%mark_facet(facet, el)
      end do
    end associate
    call this%bc_x%finalize()
    call this%bc_y%finalize()
    call this%bc_z%finalize()
  end subroutine symmetry_aligned_finalize

  subroutine symmetry_aligned_get_normal_axis(this, sx, sy, sz, facet, el)
    class(symmetry_aligned_t), target, intent(inout) :: this
    real(kind=rp), intent(out) :: sx, sy, sz
    integer, intent(in) :: facet
    integer, intent(in) :: el
    integer :: j, l

    associate(c => this%coef, nx => this%coef%nx, ny => this%coef%ny, &
         nz => this%coef%nz)
      sx = 0.0_rp
      sy = 0.0_rp
      sz = 0.0_rp
      select case (facet)
      case (1, 2)
         do l = 2, c%Xh%lx - 1
            do j = 2, c%Xh%lx -1
               sx = sx + abs(abs(nx(l, j, facet, el)) - 1.0_rp)
               sy = sy + abs(abs(ny(l, j, facet, el)) - 1.0_rp)
               sz = sz + abs(abs(nz(l, j, facet, el)) - 1.0_rp)
            end do
         end do
      case (3, 4)
         do l = 2, c%Xh%lx - 1
            do j = 2, c%Xh%lx - 1
               sx = sx + abs(abs(nx(l, j, facet, el)) - 1.0_rp)
               sy = sy + abs(abs(ny(l, j, facet, el)) - 1.0_rp)
               sz = sz + abs(abs(nz(l, j, facet, el)) - 1.0_rp)
            end do
         end do
      case (5, 6)
         do l = 2, c%Xh%lx - 1
            do j = 2, c%Xh%lx - 1
               sx = sx + abs(abs(nx(l, j, facet, el)) - 1.0_rp)
               sy = sy + abs(abs(ny(l, j, facet, el)) - 1.0_rp)
               sz = sz + abs(abs(nz(l, j, facet, el)) - 1.0_rp)
            end do
         end do
      end select
      sx = sx / (c%Xh%lx - 2)**2
      sy = sy / (c%Xh%lx - 2)**2
      sz = sz / (c%Xh%lx - 2)**2
    end associate
  end subroutine symmetry_aligned_get_normal_axis

  subroutine symmetry_aligned_apply_scalar(this, x, n, time, strong)
    class(symmetry_aligned_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
  end subroutine symmetry_aligned_apply_scalar

  subroutine symmetry_aligned_apply_vector(this, x, y, z, n, time, strong)
    class(symmetry_aligned_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       call this%bc_x%apply_scalar(x, n, strong = .true.)
       call this%bc_y%apply_scalar(y, n, strong = .true.)
       call this%bc_z%apply_scalar(z, n, strong = .true.)
    end if
  end subroutine symmetry_aligned_apply_vector

  subroutine symmetry_aligned_apply_scalar_dev(this, x_d, time, strong, strm)
    class(symmetry_aligned_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
  end subroutine symmetry_aligned_apply_scalar_dev

  subroutine symmetry_aligned_apply_vector_dev(this, x_d, y_d, z_d, &
       time, strong, strm)
    class(symmetry_aligned_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_ .and. (this%msk(0) .gt. 0)) then
       call device_symmetry_aligned_apply_vector(this%bc_x%msk_d, &
            this%bc_y%msk_d, this%bc_z%msk_d, x_d, y_d, z_d, &
            this%bc_x%msk(0), this%bc_y%msk(0), this%bc_z%msk(0), strm)
    end if
  end subroutine symmetry_aligned_apply_vector_dev

  subroutine symmetry_aligned_free(this)
    class(symmetry_aligned_t), target, intent(inout) :: this

    call this%free_base()
    call this%bc_x%free()
    call this%bc_y%free()
    call this%bc_z%free()
  end subroutine symmetry_aligned_free

end module symmetry_aligned
