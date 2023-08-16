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
!> Defines a dirichlet boundary condition
module dirichlet
  use device_dirichlet
  use num_types
  use bc
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Generic Dirichlet boundary condition
  !! \f$ x = g \f$ on \f$\partial \Omega\f$
  type, public, extends(bc_t) :: dirichlet_t
     real(kind=rp), private :: g
   contains
     procedure, pass(this) :: apply_scalar => dirichlet_apply_scalar
     procedure, pass(this) :: apply_vector => dirichlet_apply_vector
     procedure, pass(this) :: apply_scalar_dev => dirichlet_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => dirichlet_apply_vector_dev
     procedure, pass(this) :: set_g => dirichlet_set_g
  end type dirichlet_t

contains

  !> Boundary condition apply for a generic Dirichlet condition
  !! to a vector @a x
  subroutine dirichlet_apply_scalar(this, x, n, t, tstep)
    class(dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%g
    end do
  end subroutine dirichlet_apply_scalar

  !> Boundary condition apply for a generic Dirichlet condition
  !! to vectors @a x, @a y and @a z
  subroutine dirichlet_apply_vector(this, x, y, z, n, t, tstep)
    class(dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k

    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       x(k) = this%g
       y(k) = this%g
       z(k) = this%g
    end do
    
  end subroutine dirichlet_apply_vector

  !> Boundary condition apply for a generic Dirichlet condition
  !! to a vector @a x (device version)
  subroutine dirichlet_apply_scalar_dev(this, x_d, t, tstep)
    class(dirichlet_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call device_dirichlet_apply_scalar(this%msk_d, x_d, &
                                       this%g, size(this%msk))
    
  end subroutine dirichlet_apply_scalar_dev
  
  !> Boundary condition apply for a generic Dirichlet condition 
  !! to vectors @a x, @a y and @a z (device version)
  subroutine dirichlet_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(dirichlet_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    call device_dirichlet_apply_vector(this%msk_d, x_d, y_d, z_d, &
                                       this%g, size(this%msk))
    
  end subroutine dirichlet_apply_vector_dev

  !> Set value of \f$ g \f$
  subroutine dirichlet_set_g(this, g)
    class(dirichlet_t), intent(inout) :: this
    real(kind=rp), intent(in) :: g

    this%g = g
    
  end subroutine dirichlet_set_g
  
end module dirichlet
