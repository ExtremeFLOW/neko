! Copyright (c) 2020-2021, The Neko Authors
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
!> Mixed Dirichlet-Neumann axis aligned symmetry plane
module symmetry
  use device_symmetry, only : device_symmetry_apply_vector
  use dirichlet, only : dirichlet_t
  use num_types, only : rp
  use bc, only : bc_t
  use tuple, only : tuple_i4_t
  use coefs, only : coef_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Mixed Dirichlet-Neumann symmetry plane condition
  type, public, extends(bc_t) :: symmetry_t
     type(dirichlet_t) :: bc_x
     type(dirichlet_t) :: bc_y
     type(dirichlet_t) :: bc_z
   contains
     procedure, pass(this) :: init => symmetry_init
     procedure, pass(this) :: apply_scalar => symmetry_apply_scalar
     procedure, pass(this) :: apply_vector => symmetry_apply_vector
     procedure, pass(this) :: apply_scalar_dev => symmetry_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => symmetry_apply_vector_dev
     !> Destructor.
     procedure, pass(this) :: free => symmetry_free
  end type symmetry_t

contains

  !> Initialize symmetry mask for each axis
  subroutine symmetry_init(this, coef)
    class(symmetry_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    integer :: i, j, l
    type(tuple_i4_t), pointer :: bfp(:)
    real(kind=rp) :: sx, sy, sz
    real(kind=rp), parameter :: TOL = 1d-3
    type(tuple_i4_t) :: bc_facet
    integer :: facet, el

    call this%bc_x%init_base(this%coef)
    call this%bc_y%init_base(this%coef)
    call this%bc_z%init_base(this%coef)

    associate(c => this%coef, nx => this%coef%nx, ny => this%coef%ny, &
              nz => this%coef%nz)
      bfp => this%marked_facet%array()
      do i = 1, this%marked_facet%size()
         bc_facet = bfp(i)
         facet = bc_facet%x(1)
         el = bc_facet%x(2)
         sx = 0d0
         sy = 0d0
         sz = 0d0
         select case (facet)
         case (1, 2)
            do l = 2, c%Xh%lx - 1
               do j = 2, c%Xh%lx -1
                  sx = sx + abs(abs(nx(l, j, facet, el)) - 1d0)
                  sy = sy + abs(abs(ny(l, j, facet, el)) - 1d0)
                  sz = sz + abs(abs(nz(l, j, facet, el)) - 1d0)
               end do
            end do
         case (3, 4)
            do l = 2, c%Xh%lx - 1
               do j = 2, c%Xh%lx - 1
                  sx = sx + abs(abs(nx(l, j, facet, el)) - 1d0)
                  sy = sy + abs(abs(ny(l, j, facet, el)) - 1d0)
                  sz = sz + abs(abs(nz(l, j, facet, el)) - 1d0)
               end do
            end do
         case (5, 6)
            do l = 2, c%Xh%lx - 1
               do j = 2, c%Xh%lx - 1
                  sx = sx + abs(abs(nx(l, j, facet, el)) - 1d0)
                  sy = sy + abs(abs(ny(l, j, facet, el)) - 1d0)
                  sz = sz + abs(abs(nz(l, j, facet, el)) - 1d0)
               end do
            end do
         end select
         sx = sx / (c%Xh%lx - 2)**2
         sy = sy / (c%Xh%lx - 2)**2
         sz = sz / (c%Xh%lx - 2)**2

         if (sx .lt. TOL) then
            call this%bc_x%mark_facet(facet, el)
         end if

         if (sy .lt. TOL) then
            call this%bc_y%mark_facet(facet, el)
         end if

         if (sz .lt. TOL) then
            call this%bc_z%mark_facet(facet, el)
         end if
      end do
    end associate
    call this%bc_x%finalize()
    call this%bc_x%set_g(0.0_rp)
    call this%bc_y%finalize()
    call this%bc_y%set_g(0.0_rp)
    call this%bc_z%finalize()
    call this%bc_z%set_g(0.0_rp)

  end subroutine symmetry_init

  !> No-op scalar apply
  subroutine symmetry_apply_scalar(this, x, n, t, tstep)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
  end subroutine symmetry_apply_scalar

  !> Apply symmetry conditions (axis aligned)
  subroutine symmetry_apply_vector(this, x, y, z, n, t, tstep)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call this%bc_x%apply_scalar(x, n)
    call this%bc_y%apply_scalar(y, n)
    call this%bc_z%apply_scalar(z, n)

  end subroutine symmetry_apply_vector

  !> No-op scalar apply (device version)
  subroutine symmetry_apply_scalar_dev(this, x_d, t, tstep)
    class(symmetry_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
  end subroutine symmetry_apply_scalar_dev

  !> Apply symmetry conditions (axis aligned) (device version)
  subroutine symmetry_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(symmetry_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call device_symmetry_apply_vector(this%bc_x%msk_d, this%bc_y%msk_d, &
                                      this%bc_z%msk_d, x_d, y_d, z_d, &
                                      this%bc_x%msk(0), &
                                      this%bc_y%msk(0), &
                                      this%bc_z%msk(0))

  end subroutine symmetry_apply_vector_dev

  !> Destructor
  subroutine symmetry_free(this)
    class(symmetry_t), target, intent(inout) :: this

    call this%free_base()
    call this%bc_x%free()
    call this%bc_y%free()
    call this%bc_z%free()

  end subroutine symmetry_free
end module symmetry
