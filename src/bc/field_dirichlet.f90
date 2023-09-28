! Copyright (c) 2020-2023, The Neko Authors
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
!> Defines inflow dirichlet conditions
module field_dirichlet
  use num_types
  use coefs
  use dirichlet
  use device
  use utils
  use field, only : field_t
  use field_list, only : field_list_t
  use math
  use device_math
  use dofmap, only : dofmap_t
  implicit none
  private
  
  !> User defined dirichlet condition for inlet (vector valued)
  !! Would be neat to add another class that contains all three 
  !! dirichlet bcs for the velocity, this bc would then implement
  !! apply_vector.
  type, public, extends(dirichlet_t) :: field_dirichlet_t
     type(field_t) :: field_bc
   contains
     procedure, pass(this) :: init_field => field_dirichlet_init
     procedure, pass(this) :: apply_scalar => field_dirichlet_apply_scalar
     procedure, pass(this) :: apply_vector => field_dirichlet_apply_vector
     procedure, pass(this) :: apply_vector_dev => field_dirichlet_apply_vector_dev
     procedure, pass(this) :: apply_scalar_dev => field_dirichlet_apply_scalar_dev
  end type field_dirichlet_t

  abstract interface
   
     !> Abstract interface defining a dirichlet condition on a list of fields
     !! @param field_bc_list list of fields that are used to extract values for field_dirichlet
     !! @param t Current time
     !! @param tstep Current time-step
     !! Should we pass down coef? and perhaps something more?
     subroutine field_dirichlet_update(field_bc_list, t, tstep)
       import rp
       import field_list_t
       type(field_list_t), intent(inout) :: field_bc_list
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
     end subroutine field_dirichlet_update
  end interface

  public :: field_dirichlet_update

contains
     
  !> Initialize symmetry mask for each axis
  subroutine field_dirichlet_init(this,bc_name)
    class(field_dirichlet_t), intent(inout) :: this
    character(len=*), intent(in) :: bc_name
    call this%field_bc%init(this%dof, bc_name)
  end subroutine field_dirichlet_init
  
  subroutine field_dirichlet_free(this)
    type(field_dirichlet_t), intent(inout) :: this

    call this%field_bc%free()
    
  end subroutine field_dirichlet_free
  
  !> No-op scalar apply 
  subroutine field_dirichlet_apply_scalar(this, x, n, t, tstep)
    class(field_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    if (this%msk(0) .gt. 0) then
       call masked_copy(x,this%field_bc%x,this%msk, n, this%msk(0))
    end if
  end subroutine field_dirichlet_apply_scalar
  
  !> No-op scalar apply (device version)
  subroutine field_dirichlet_apply_scalar_dev(this, x_d, t, tstep)
    class(field_dirichlet_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep

    if (this%msk(0) .gt. 0) then
       call device_masked_copy(x_d,this%field_bc%x_d,this%msk_d, &
            this%field_bc%dof%size(), this%msk(0))
    end if
  
  end subroutine field_dirichlet_apply_scalar_dev

  !> Apply field defined dirichlet conditions (vector valued)
  subroutine field_dirichlet_apply_vector(this, x, y, z, n, t, tstep)
    class(field_dirichlet_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k, idx(4), facet, tstep_
    real(kind=rp) :: t_

  end subroutine field_dirichlet_apply_vector

  subroutine field_dirichlet_apply_vector_dev(this, x_d, y_d, z_d, t, tstep)
    class(field_dirichlet_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    integer :: i, m, k, idx(4), facet, tstep_
    integer(c_size_t) :: s
    real(kind=rp) :: t_
    real(kind=rp), allocatable :: x(:)
    real(kind=rp), allocatable :: y(:)
    real(kind=rp), allocatable :: z(:)

   end subroutine field_dirichlet_apply_vector_dev

end module field_dirichlet
