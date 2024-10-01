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
!> Dirichlet condition on axis aligned plane in the non normal direction
module non_normal
  use symmetry, only : symmetry_t
  use num_types, only : rp
  use tuple, only : tuple_i4_t
  use coefs, only : coef_t
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Dirichlet condition in non normal direction of a plane
  type, public, extends(symmetry_t) :: non_normal_t
   contains
     !> Constructor.
     procedure, pass(this) :: init => non_normal_init
     !> Destructor.
     procedure, pass(this) :: free => non_normal_free
  end type non_normal_t

contains

  !> Constructor.
  subroutine non_normal_init(this, coef)
    class(non_normal_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    integer :: i, j, l
    type(tuple_i4_t), pointer :: bfp(:)
    real(kind=rp) :: sx, sy, sz
    real(kind=rp), parameter :: TOL = 1d-3
    type(tuple_i4_t) :: bc_facet
    integer :: facet, el

    call this%bc_x%free()
    call this%bc_y%free()
    call this%bc_z%free()
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
            call this%bc_y%mark_facet(facet, el)
            call this%bc_z%mark_facet(facet, el)
         end if

         if (sy .lt. TOL) then
            call this%bc_x%mark_facet(facet, el)
            call this%bc_z%mark_facet(facet, el)
         end if

         if (sz .lt. TOL) then
            call this%bc_y%mark_facet(facet, el)
            call this%bc_x%mark_facet(facet, el)
         end if
      end do
    end associate
    call this%bc_x%finalize()
    call this%bc_x%set_g(0.0_rp)
    call this%bc_y%finalize()
    call this%bc_y%set_g(0.0_rp)
    call this%bc_z%finalize()
    call this%bc_z%set_g(0.0_rp)
  end subroutine non_normal_init

  !> Destructor
  subroutine non_normal_free(this)
    class(non_normal_t), target, intent(inout) :: this

    call this%symmetry_t%free()
  end subroutine non_normal_free
end module non_normal
