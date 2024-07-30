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
!> Implements `dynamic_smagorinsky_t`.
module dynamic_smagorinsky
  use num_types, only : rp
  use math
  use field_list, only : field_list_t
  use field, only : field_t
  use les_model, only : les_model_t
  use dofmap , only : dofmap_t
  use json_utils, only : json_get, json_get_or_default
  use json_module, only : json_file
  use utils, only : neko_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use coefs, only : coef_t
  use elementwise_filter, only : elementwise_filter_t
  use dynamic_smagorinsky_cpu, only : dynamic_smagorinsky_compute_cpu
  implicit none
  private

  !> Implements the dynamic Smagorinsky LES model.
  !! @note Reference DOI: 10.1063/1.857955
  type, public, extends(les_model_t) :: dynamic_smagorinsky_t
     !> Coefficient of the model.
     type(field_t) :: c_dyn
     !> Test filter.
     character(len=:), allocatable :: test_filter_type
     type(elementwise_filter_t) :: test_filter
     !> Mij
     !! index for tensor mij and lij:
     !! 1=>1,1; 2=>2,2; 3=>3,3; 4=>1,2; 5=>1,3; 6=>2,3;
     type(field_t) :: mij(6)
     !> Germano Identity
     type(field_t) :: lij(6)
     !> <M_ij L_ij>
     type(field_t) :: num
     !> <M_lm M_lm>
     type(field_t) :: den
  contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => dynamic_smagorinsky_init
     !> Destructor.
     procedure, pass(this) :: free => dynamic_smagorinsky_free
     !> Compute eddy viscosity.
     procedure, pass(this) :: compute => dynamic_smagorinsky_compute

  end type dynamic_smagorinsky_t

contains
  !> Constructor.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  subroutine dynamic_smagorinsky_init(this, dofmap, coef, json)
    class(dynamic_smagorinsky_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dofmap
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: nut_name !! The name of the SGS viscosity field.
    integer :: i
    character(len=:), allocatable :: delta_type

    call json_get_or_default(json, "nut_field", nut_name, "nut")
    call json_get_or_default(json, "delta_type", delta_type, "pointwise")

    call this%free()
    call this%init_base(dofmap, coef, nut_name, delta_type)
    this%test_filter_type = "nonBoyd"
    ! Filter assumes lx = ly = lz
    call this%test_filter%init(dofmap%xh%lx, this%test_filter_type)
    call set_ds_filt(this%test_filter)

    call this%c_dyn%init(dofmap, "ds_c_dyn")
    call this%num%init(dofmap, "ds_num")
    call this%den%init(dofmap, "ds_den")

    do i = 1, 6
       call this%mij(i)%init(dofmap)
       call this%lij(i)%init(dofmap)
    end do

  end subroutine dynamic_smagorinsky_init

  !> Destructor for the les_model_t (base) class.
  subroutine dynamic_smagorinsky_free(this)
    class(dynamic_smagorinsky_t), intent(inout) :: this
    integer :: i

    call this%c_dyn%free()
    do i = 1, 6
       call this%mij(i)%free()
       call this%lij(i)%free()
    end do
    call this%num%free()
    call this%den%free()
    call this%test_filter%free()
    call this%free_base()

  end subroutine dynamic_smagorinsky_free

  !> Compute eddy viscosity.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine dynamic_smagorinsky_compute(this, t, tstep)
    class(dynamic_smagorinsky_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    if (NEKO_BCKND_DEVICE .eq. 1) then
        call neko_error("Dynamic Smagorinsky model not implemented on &
             &accelarators.")
    else
        call dynamic_smagorinsky_compute_cpu(t, tstep, this%coef, this%nut, &
                                this%delta, this%c_dyn, this%test_filter, &
                                this%mij, this%lij, this%num, this%den)
    end if

  end subroutine dynamic_smagorinsky_compute

  !> Set up the test filter
  subroutine set_ds_filt(filter_1d)
    type(elementwise_filter_t), intent(inout) :: filter_1d
    integer :: i

    if (filter_1d%nx .le. 2) then
        call neko_error("Dynamic Smagorinsky model error: test filter is not &
             &defined for the current polynomial order")
    end if
    if (mod(filter_1d%nx,2) .eq. 0) then ! number of grid spacing is odd
       ! cutoff at polynomial order int((filter_1d%nx)/2)
       filter_1d%trnsfr(int((filter_1d%nx)/2)) = 0.95_rp
       filter_1d%trnsfr(int((filter_1d%nx)/2) + 1) = 0.50_rp
       filter_1d%trnsfr(int((filter_1d%nx)/2) + 2) = 0.05_rp
       if ((int((filter_1d%nx)/2) + 2) .lt. filter_1d%nx) then
          do i = int((filter_1d%nx)/2) + 3, filter_1d%nx
             filter_1d%trnsfr(i) = 0.0_rp
          end do
       end if
       ! make delta_ratio = (nx-1)/(nt-1) as close to 2
       filter_1d%nt = int(filter_1d%nx/2) + 1
    else ! number of grid spacing is even
       ! cutoff at polynomial order int((filter_1d%nx-1)/2)
       filter_1d%trnsfr(int((filter_1d%nx-1)/2)) = 0.95_rp
       filter_1d%trnsfr(int((filter_1d%nx-1)/2) + 1) = 0.50_rp
       filter_1d%trnsfr(int((filter_1d%nx-1)/2) + 2) = 0.05_rp
       if ((int((filter_1d%nx-1)/2) + 2) .lt. filter_1d%nx) then
          do i = int((filter_1d%nx-1)/2) + 3, filter_1d%nx
             filter_1d%trnsfr(i) = 0.0_rp
          end do
       end if
       ! make delta_ratio = (nx-1)/(nt-1) = 2
       filter_1d%nt = int((filter_1d%nx-1)/2) + 1
    end if

    call filter_1d%build_1d()

  end subroutine set_ds_filt

end module dynamic_smagorinsky
