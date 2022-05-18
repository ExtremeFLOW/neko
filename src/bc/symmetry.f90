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
  use device_symmetry
  use neko_config
  use num_types
  use dirichlet
  use bc
  use device
  use coefs
  use math
  use utils
  use stack
  use tuple
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Mixed Dirichlet-Neumann symmetry plane condition
  type, public, extends(bc_t) :: symmetry_t
     type(dirichlet_t) :: bc_x
     type(dirichlet_t) :: bc_y
     type(dirichlet_t) :: bc_z
   contains
     procedure, pass(this) :: init_msk => symmetry_init_msk
     procedure, pass(this) :: apply_scalar => symmetry_apply_scalar
     procedure, pass(this) :: apply_vector => symmetry_apply_vector
     procedure, pass(this) :: apply_scalar_dev => symmetry_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => symmetry_apply_vector_dev
     final :: symmetry_free
  end type symmetry_t

contains

  !> Initialize symmetry mask for each axis
  subroutine symmetry_init_msk(this, c)
    class(symmetry_t), intent(inout) :: this
    type(coef_t), intent(in) :: c
    integer :: i, m, j, k, l
    type(tuple_i4_t), pointer :: bfp(:)
    real(kind=rp) :: sx,sy,sz
    real(kind=rp), parameter :: TOL = 1d-3
    type(tuple_i4_t) :: bc_facet
    integer :: facet, el
    
    call symmetry_free(this)

    call this%bc_x%init(c%dof)
    call this%bc_y%init(c%dof)
    call this%bc_z%init(c%dof)
    
    associate(nx => c%nx, ny => c%ny, nz => c%nz)
      bfp => this%marked_facet%array()
      do i = 1, this%marked_facet%size()
         k = this%msk(i)
         bc_facet = bfp(i)
         facet = bc_facet%x(1)
         el = bc_facet%x(2)
         sx = 0d0
         sy = 0d0
         sz = 0d0
         select case (facet)               
         case(1,2)
            do l = 2, c%Xh%lx - 1
               do j = 2, c%Xh%lx -1
                  sx = sx + abs(abs(nx(l, j, facet, el)) - 1d0)
                  sy = sy + abs(abs(ny(l, j, facet, el)) - 1d0)
                  sz = sz + abs(abs(nz(l, j, facet, el)) - 1d0)
               end do
            end do
         case(3,4)
            do l = 2, c%Xh%lx - 1
               do j = 2, c%Xh%lx - 1
                  sx = sx + abs(abs(nx(l, j, facet, el)) - 1d0)
                  sy = sy + abs(abs(ny(l, j, facet, el)) - 1d0)
                  sz = sz + abs(abs(nz(l, j, facet, el)) - 1d0)
               end do
            end do
         case(5,6)
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

  end subroutine symmetry_init_msk
  
  subroutine symmetry_free(this)
    type(symmetry_t), intent(inout) :: this
    
    call this%bc_x%free()
    call this%bc_y%free()
    call this%bc_z%free()

  end subroutine symmetry_free
  
  !> No-op scalar apply
  subroutine symmetry_apply_scalar(this, x, n)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
  end subroutine symmetry_apply_scalar

  !> Apply symmetry conditions (axis aligned)
  subroutine symmetry_apply_vector(this, x, y, z, n)
    class(symmetry_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    call this%bc_x%apply_scalar(x,n)
    call this%bc_y%apply_scalar(y,n)
    call this%bc_z%apply_scalar(z,n)
    
  end subroutine symmetry_apply_vector

  !> No-op scalar apply (device version)
  subroutine symmetry_apply_scalar_dev(this, x_d)
    class(symmetry_t), intent(inout), target :: this
    type(c_ptr) :: x_d
  end subroutine symmetry_apply_scalar_dev

  !> Apply symmetry conditions (axis aligned) (device version)
  subroutine symmetry_apply_vector_dev(this, x_d, y_d, z_d)
    class(symmetry_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d

    call device_symmetry_apply_vector(this%bc_x%msk_d, this%bc_y%msk_d, &
                                      this%bc_z%msk_d, x_d, y_d, z_d, &
                                      this%bc_x%msk(0), &
                                      this%bc_y%msk(0), &
                                      this%bc_z%msk(0))

  end subroutine symmetry_apply_vector_dev
      
end module symmetry
