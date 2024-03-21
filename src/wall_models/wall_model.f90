! Copyright (c) 2024, The Neko Authors
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
!
!> Implements `wall_model_t`.
module wall_model
  use num_types, only : rp
  use field, only : field_t, field_ptr_t
  use json_module, only : json_file
  use field_registry, only : neko_field_registry
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, HOST_TO_DEVICE
  use vector, only : vector_t
  use utils, only : neko_error, nonlinear_index

  implicit none
  private

  !> Base abstract type for wall-stress models for wall-modelled LES.
  type, abstract, public :: wall_model_t
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     type(dofmap_t), pointer :: dofmap => null()
     !> The boundary condition mask. First element, holds the array size!
     integer, pointer :: msk(:) => null()
     !> The boundary condition facet ids. First element holds the array size!
     integer, pointer :: facet(:) => null()
     !> The x component of the shear stress.
     real(kind=rp), allocatable :: tau_x(:)
     !> The y component of the shear stress.
     real(kind=rp), allocatable :: tau_y(:)
     !> The z component of the shear stress.
     real(kind=rp), allocatable :: tau_z(:)
     !> The r indices of the sampling points
     integer, allocatable :: ind_r(:)
     !> The s indices of the sampling points
     integer, allocatable :: ind_s(:)
     !> The t indices of the sampling points
     integer, allocatable :: ind_t(:)
     !> The element indices of the sampling points
     integer, allocatable :: ind_e(:)
     !> The sampling height
     type(vector_t) :: h
     !> Number of nodes in the boundary
     integer :: n_nodes = 0
   contains
     !> Constructor for the wall_model_t (base) class.
     procedure, pass(this) :: init_base => wall_model_init_base
     !> Destructor for the wall_model_t (base) class.
     procedure, pass(this) :: free_base => wall_model_free_base
     !> The common constructor.
     procedure(wall_model_init), pass(this), deferred :: init
     !> Destructor.
     procedure(wall_model_free), pass(this), deferred :: free
     !> Compute the wall shear stress.
     procedure(wall_model_compute), pass(this), deferred :: compute

     procedure, pass(this) :: find_points => wall_model_find_points
  end type wall_model_t

  abstract interface
     !> Compute wall shear stress.
     !! @param t The time value.
     !! @param tstep The current time-step.
     subroutine wall_model_compute(this, t, tstep)
       import wall_model_t, rp
       class(wall_model_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
     end subroutine wall_model_compute
  end interface

  abstract interface
     !> Common constructor.
     !! @param dofmap SEM map of degrees of freedom.
     !! @param coef SEM coefficients.
     !! @param json A dictionary with parameters.
     subroutine wall_model_init(this, dofmap, coef, msk, facet, json)
       import wall_model_t, json_file, dofmap_t, coef_t
       class(wall_model_t), intent(inout) :: this
       type(dofmap_t), intent(in) :: dofmap
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: msk(:)
       integer, intent(in) :: facet(:)
       type(json_file), intent(inout) :: json
     end subroutine wall_model_init
  end interface

  abstract interface
     !> Destructor.
     subroutine wall_model_free(this)
       import wall_model_t
       class(wall_model_t), intent(inout) :: this
     end subroutine wall_model_free
  end interface

contains
  !> Constructor for the wall_model_t (base) class.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param nu_name The name of the turbulent viscosity field.
  subroutine wall_model_init_base(this, dofmap, coef, msk, facet)
    class(wall_model_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), target, intent(in) :: coef
    integer, target, intent(in) :: msk(:)
    integer, target, intent(in) :: facet(:)

    this%coef => coef
    this%dofmap => coef%dof
    this%msk => msk
    this%facet => facet

    allocate(this%tau_x(this%msk(1)))
    allocate(this%tau_y(this%msk(1)))
    allocate(this%tau_z(this%msk(1)))

    allocate(this%ind_r(this%msk(1)))
    allocate(this%ind_s(this%msk(1)))
    allocate(this%ind_t(this%msk(1)))
    allocate(this%ind_e(this%msk(1)))

    call this%h%init(this%msk(1))

    call this%find_points

  end subroutine wall_model_init_base

  !> Destructor for the wall_model_t (base) class.
  subroutine wall_model_free_base(this)
    class(wall_model_t), intent(inout) :: this

    nullify(this%coef)
    nullify(this%msk)
    nullify(this%facet)

    if (allocated(this%tau_x)) then
      deallocate(this%tau_x)
    end if
    if (allocated(this%tau_y)) then
      deallocate(this%tau_y)
    end if
    if (allocated(this%tau_z)) then
      deallocate(this%tau_z)
    end if
    if (allocated(this%ind_r)) then
      deallocate(this%ind_r)
    end if

    call this%h%free()
  end subroutine wall_model_free_base

  subroutine wall_model_find_points(this)
    class(wall_model_t), intent(inout) :: this
    integer :: n_nodes, fid, idx(4), i, linear
    real(kind=rp) :: normal(3), p(3), x, y, z, xw, yw, zw, dot

    n_nodes = this%msk(1)
    this%n_nodes = n_nodes
    do i = 1, n_nodes
       ! Because these arrays in bc.f90 start from 0, we have to add 1 here.
       linear = this%msk(i + 1)
       fid = this%facet(i + 1)
       idx = nonlinear_index(linear, this%coef%Xh%lx, this%coef%Xh%lx,&
                             this%coef%Xh%lx)
       normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), fid)
       ! inward normal
       normal = -normal

       select case (fid)
       case(1)
         this%ind_r(i) = idx(1) + 1
         this%ind_s(i) = idx(2)
         this%ind_t(i) = idx(3)
       case(2)
         this%ind_r(i) = idx(1) - 1
         this%ind_s(i) = idx(2)
         this%ind_t(i) = idx(3)
       case(3)
         this%ind_r(i) = idx(1)
         this%ind_s(i) = idx(2) + 1
         this%ind_t(i) = idx(3)
       case(4)
         this%ind_r(i) = idx(1)
         this%ind_s(i) = idx(2) - 1
         this%ind_t(i) = idx(3)
       case(5)
         this%ind_r(i) = idx(1)
         this%ind_s(i) = idx(2)
         this%ind_t(i) = idx(3) + 1
       case(6)
         this%ind_r(i) = idx(1)
         this%ind_s(i) = idx(2)
         this%ind_t(i) = idx(3) - 1
       case default
         call neko_error("The face index is not correct")
       end select
       this%ind_e(i) = idx(4)

       ! Location of the wall node
       xw = this%dofmap%x(idx(1), idx(2), idx(3), idx(4))
       yw = this%dofmap%y(idx(1), idx(2), idx(3), idx(4))
       zw = this%dofmap%z(idx(1), idx(2), idx(3), idx(4))

       ! Location of the sampling point
!       write(*,*) "IND", this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i), fid
       x = this%dofmap%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
                         this%ind_e(i))
       y = this%dofmap%y(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
                         this%ind_e(i))
       z = this%dofmap%z(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
                         this%ind_e(i))


       ! vector from the sampling point to the wall
       p(1) = x - xw
       p(2) = y - yw
       p(3) = z - zw

       ! project on the normal direction to get h
       this%h%x(i) = p(1)*normal(1) + p(2)*normal(2) + p(3)*normal(3)
    end do

    !write(*,*) this%h%x(1:3), size(this%h%x)

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(this%h%x, this%h%x_d, size(this%h%x), HOST_TO_DEVICE,&
                         sync=.false.)
    end if
  end subroutine wall_model_find_points

end module wall_model