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
  use math, only : glmin, glmax
  use comm, only : pe_rank
  use logger, only : neko_log, NEKO_LOG_DEBUG

  implicit none
  private

  !> Base abstract type for wall-stress models for wall-modelled LES.
  type, abstract, public :: wall_model_t
     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Map of degrees of freedom.
     type(dofmap_t), pointer :: dof => null()
     !> The boundary condition mask. Stores the array size at index zero!
     integer, pointer :: msk(:) => null()
     !> The boundary condition facet ids. Stores the array size at index zero!
     integer, pointer :: facet(:) => null()
     !> The x component of the shear stress.
     real(kind=rp), allocatable :: tau_x(:)
     !> The y component of the shear stress.
     real(kind=rp), allocatable :: tau_y(:)
     !> The z component of the shear stress.
     real(kind=rp), allocatable :: tau_z(:)
     !> The x component of the normal.
     type(vector_t) :: n_x
     !> The y component of the normal.
     type(vector_t) :: n_y
     !> The z component of the normal.
     type(vector_t) :: n_z
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
     !> Sampling index
     integer :: h_index = 0
     !> Number of nodes in the boundary
     integer :: n_nodes = 0
     !> Kinematic viscosity value.
     real(kind=rp) :: nu = 0_rp
     !> The 3D field with the computed stress magnitude at the boundary.
     type(field_t), pointer :: tau_field => null()
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
     !> Find the sampling points based on the value of `h_index`.
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
     !! @param coef SEM coefficients.
     !! @param msk The boundary mask.
     !! @param facet The boundary facets.
     !! @param nu The molecular kinematic viscosity.
     !! @param h_index The off-wall index of the sampling cell.
     !! @param json A dictionary with parameters.
     subroutine wall_model_init(this, coef, msk, facet, nu, h_index, json)
       import wall_model_t, json_file, dofmap_t, coef_t, rp
       class(wall_model_t), intent(inout) :: this
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: msk(:)
       integer, intent(in) :: facet(:)
       real(kind=rp), intent(in) :: nu
       integer, intent(in) :: h_index
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

  interface
     !> Wall model factory. Both constructs and initializes the object.
     !! @param object The object to be allocated.
     !! @param coef SEM coefficients.
     !! @param msk The boundary mask.
     !! @param facet The boundary facets.
     !! @param nu The molecular kinematic viscosity.
     !! @param h_index The off-wall index of the sampling cell.
     !! @param json A dictionary with parameters.
     module subroutine wall_model_factory(object, coef, msk, facet, nu, &
          h_index, json)
       class(wall_model_t), allocatable, target, intent(inout) :: object
       type(coef_t), intent(in) :: coef
       integer, intent(in) :: msk(:)
       integer, intent(in) :: facet(:)
       real(kind=rp), intent(in) :: nu
       integer, intent(in) :: h_index
       type(json_file), intent(inout) :: json
     end subroutine wall_model_factory
  end interface

  public :: wall_model_factory

contains
  !> Constructor for the wall_model_t (base) class.
  !! @param dof SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param nu_name The name of the turbulent viscosity field.
  subroutine wall_model_init_base(this, coef, msk, facet, nu, index)
    class(wall_model_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    integer, target, intent(in) :: msk(0:)
    integer, target, intent(in) :: facet(:)
    real(kind=rp), intent(in) :: nu
    integer, intent(in) :: index

    call this%free_base

    this%coef => coef
    this%dof => coef%dof
    this%msk(0:msk(0)) => msk
    this%facet => facet
    this%nu = nu
    this%h_index = index

    call neko_field_registry%add_field(this%dof, "tau", &
                                       ignore_existing = .true.)

    this%tau_field => neko_field_registry%get_field("tau")

    allocate(this%tau_x(this%msk(0)))
    allocate(this%tau_y(this%msk(0)))
    allocate(this%tau_z(this%msk(0)))

    allocate(this%ind_r(this%msk(0)))
    allocate(this%ind_s(this%msk(0)))
    allocate(this%ind_t(this%msk(0)))
    allocate(this%ind_e(this%msk(0)))

    call this%h%init(this%msk(0))
    call this%n_x%init(this%msk(0))
    call this%n_y%init(this%msk(0))
    call this%n_z%init(this%msk(0))

    call this%find_points

  end subroutine wall_model_init_base

  !> Destructor for the wall_model_t (base) class.
  subroutine wall_model_free_base(this)
    class(wall_model_t), intent(inout) :: this

    nullify(this%coef)
    nullify(this%msk)
    nullify(this%facet)
    nullify(this%tau_field)

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
    call this%n_x%free()
    call this%n_y%free()
    call this%n_z%free()
  end subroutine wall_model_free_base

  subroutine wall_model_find_points(this)
    class(wall_model_t), intent(inout) :: this
    integer :: n_nodes, fid, idx(4), i, linear
    real(kind=rp) :: normal(3), p(3), x, y, z, xw, yw, zw, magp
    real(kind=rp) :: hmin, hmax

    n_nodes = this%msk(0)
    this%n_nodes = n_nodes
    do i = 1, n_nodes
       linear = this%msk(i)
       fid = this%facet(i)
       idx = nonlinear_index(linear, this%coef%Xh%lx, this%coef%Xh%lx,&
                             this%coef%Xh%lx)
       normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), fid)

       this%n_x%x(i) = normal(1)
       this%n_y%x(i) = normal(2)
       this%n_z%x(i) = normal(3)

       ! inward normal
       normal = -normal

       select case (fid)
       case (1)
         this%ind_r(i) = idx(1) + this%h_index
         this%ind_s(i) = idx(2)
         this%ind_t(i) = idx(3)
       case (2)
         this%ind_r(i) = idx(1) - this%h_index
         this%ind_s(i) = idx(2)
         this%ind_t(i) = idx(3)
       case (3)
         this%ind_r(i) = idx(1)
         this%ind_s(i) = idx(2) + this%h_index
         this%ind_t(i) = idx(3)
       case (4)
         this%ind_r(i) = idx(1)
         this%ind_s(i) = idx(2) - this%h_index
         this%ind_t(i) = idx(3)
       case (5)
         this%ind_r(i) = idx(1)
         this%ind_s(i) = idx(2)
         this%ind_t(i) = idx(3) + this%h_index
       case (6)
         this%ind_r(i) = idx(1)
         this%ind_s(i) = idx(2)
         this%ind_t(i) = idx(3) - this%h_index
       case default
         call neko_error("The face index is not correct")
       end select
       this%ind_e(i) = idx(4)

       ! Location of the wall node
       xw = this%dof%x(idx(1), idx(2), idx(3), idx(4))
       yw = this%dof%y(idx(1), idx(2), idx(3), idx(4))
       zw = this%dof%z(idx(1), idx(2), idx(3), idx(4))

       ! Location of the sampling point
!       write(*,*) "IND", this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i), fid
       x = this%dof%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
                      this%ind_e(i))
       y = this%dof%y(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
                      this%ind_e(i))
       z = this%dof%z(this%ind_r(i), this%ind_s(i), this%ind_t(i), &
                         this%ind_e(i))


       ! Vector from the sampling point to the wall
       p(1) = x - xw
       p(2) = y - yw
       p(3) = z - zw

       ! Total distance to the sampling point
       magp = sqrt(p(1)**2 + p(2)**2 + p(3)**2)

       ! Project on the normal direction to get h
       this%h%x(i) = p(1)*normal(1) + p(2)*normal(2) + p(3)*normal(3)

       ! Look at how much the total distance distance from the normal and warn
       ! if significant
       if ((this%h%x(i) - magp) / magp > 0.1 &
           .and. (neko_log%level_ .eq. NEKO_LOG_DEBUG)) then
          write(*,*) "Significant missalignment between wall normal and &
                   & sampling point direction at wall node", xw, yw, zw
       end if
    end do

    hmin = glmin(this%h%x, n_nodes)
    hmax = glmax(this%h%x, n_nodes)

    if (pe_rank .eq. 0) then
       write(*,*) "h min / max", hmin, hmax
    end if



    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(this%h%x, this%h%x_d, n_nodes, HOST_TO_DEVICE,&
                         sync = .false.)
      call device_memcpy(this%n_x%x, this%n_x%x_d, n_nodes, HOST_TO_DEVICE, &
                         sync = .false.)
      call device_memcpy(this%n_y%x, this%n_y%x_d, n_nodes, HOST_TO_DEVICE, &
                         sync = .false.)
      call device_memcpy(this%n_z%x, this%n_z%x_d, n_nodes, HOST_TO_DEVICE, &
                         sync = .true.)
    end if
  end subroutine wall_model_find_points

end module wall_model
