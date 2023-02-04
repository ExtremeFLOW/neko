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
!> Defines wall boundary conditions
module wall
  use device_wall
  use num_types
  use dirichlet
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> No-slip Wall boundary condition
  type, public, extends(dirichlet_t) :: no_slip_wall_t
   contains
     procedure, pass(this) :: apply_scalar => no_slip_wall_apply_scalar
     procedure, pass(this) :: apply_vector => no_slip_wall_apply_vector
     procedure, pass(this) :: apply_scalar_dev => no_slip_wall_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => no_slip_wall_apply_vector_dev
  end type no_slip_wall_t

contains

  !> Boundary condition apply for a no-slip wall condition
  !! to a vector @a x
  subroutine no_slip_wall_apply_scalar(this, x, n)
    class(no_slip_wall_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    integer :: i, m, k

    m = this%msk(0)
    !$omp do
    do i = 1, m
       k = this%msk(i)
       x(k) = 0d0
    end do
    !$omp end do
  end subroutine no_slip_wall_apply_scalar
  
  !> Boundary condition apply for a no-slip wall condition
  !! to vectors @a x, @a y and @a z
  subroutine no_slip_wall_apply_vector(this, x, y, z, n)
    class(no_slip_wall_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    integer :: i, m, k

    m = this%msk(0)
    !$omp do
    do i = 1, m
       k = this%msk(i)
       x(k) = 0d0
       y(k) = 0d0
       z(k) = 0d0
    end do
    !$omp end do
    
  end subroutine no_slip_wall_apply_vector

  !> Boundary condition apply for a no-slip wall condition
  !! to a vector @a x (device version)
  subroutine no_slip_wall_apply_scalar_dev(this, x_d)
    class(no_slip_wall_t), intent(inout), target :: this
    type(c_ptr) :: x_d

    call device_no_slip_wall_apply_scalar(this%msk_d, x_d, size(this%msk))
    
  end subroutine no_slip_wall_apply_scalar_dev
  
  !> Boundary condition apply for a no-slip wall condition
  !! to vectors @a x, @a y and @a z (device version)
  subroutine no_slip_wall_apply_vector_dev(this, x_d, y_d, z_d)
    class(no_slip_wall_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d

    call device_no_slip_wall_apply_vector(this%msk_d, x_d, y_d, z_d, &
                                          size(this%msk))
    
  end subroutine no_slip_wall_apply_vector_dev
  
end module wall
