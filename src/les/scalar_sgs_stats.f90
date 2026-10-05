! Copyright (c) 2026, The Neko Authors
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
!> Computes the subgrid-scale contributions for the scalar fluxes.
!! We use the Reynolds decomposition for a field u = <u> + u' = U + u'
module scalar_sgs_stats
  use mean_field, only : mean_field_t
  use num_types, only : rp
  use field_math, only : field_cmult2
  use operators, only : grad
  use coefs, only : coef_t
  use field, only : field_t
  use field_list, only : field_list_t
  use stats_quant, only : stats_quant_t
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private

  type, public, extends(stats_quant_t) :: scalar_sgs_stats_t
     !> Pointers to the instantenious quantities.
     !> Scalar diffusivity. Not associated when it is computed from nut, in
     !! which case it is only formed in a scratch field during update.
     type(field_t), pointer :: alphat => null()
     type(field_t), pointer :: nut => null() !< scalar diffusivity
     logical :: nut_dependency = .false.
     real(kind=rp) :: pr_turb !< turbulent Prandtl number
     type(field_t), pointer :: s => null() !< scalar

     type(mean_field_t) :: alphat_mean !< <alphat>

     type(mean_field_t) :: alphatdsdx !< <alphat*dsdx>
     type(mean_field_t) :: alphatdsdy !< <alphat*dsdy>
     type(mean_field_t) :: alphatdsdz !< <alphat*dsdz>

     !> SEM coefficients.
     type(coef_t), pointer :: coef => null()
     !> Number of statistical fields to be computed.
     integer :: n_stats = 4
     !> A list of size n_stats, whith entries pointing to the fields that will
     !! be output (the field components above.) Used to write the output.
     type(field_list_t) :: stat_fields
   contains
     generic :: init => init_alphat, init_nut
     !> Constructor.
     procedure, pass(this) :: init_alphat => scalar_sgs_stats_init_alphat
     procedure, pass(this) :: init_nut => scalar_sgs_stats_init_nut
     !> Destructor.
     procedure, pass(this) :: free => scalar_sgs_stats_free
     !> Update all the mean value fields with a new sample.
     procedure, pass(this) :: update => scalar_sgs_stats_update
     !> Reset all the computed means values and sampling times to zero.
     procedure, pass(this) :: reset => scalar_sgs_stats_reset
  end type scalar_sgs_stats_t

contains

  !> Constructor. Initialize the fields associated with scalar_sgs_stats.
  !> This version uses alphat directly.
  !! @param coef SEM coefficients. Optional.
  !! @param s_name The name of the scalar field.
  !! @param alphat_field Specifies the name of the alphat field.
  subroutine scalar_sgs_stats_init_alphat(this, coef, s_name, alphat_field)
    class(scalar_sgs_stats_t), intent(inout), target:: this
    type(coef_t), target, optional :: coef
    character(len=*), intent(in) :: s_name
    character(len=*), intent(in) :: alphat_field

    call this%free()
    this%coef => coef

    this%s => neko_registry%get_field(s_name)
    this%alphat => neko_registry%get_field(alphat_field)
    this%nut_dependency = .false.

    ! Initialize mean fields. The fluxes are means of products, updated with
    ! update_col3. The field passed to init only provides the dofmap.
    call this%alphat_mean%init(this%alphat)

    call this%alphatdsdx%init(this%s, trim(s_name) // 'alphatdsdx')
    call this%alphatdsdy%init(this%s, trim(s_name) // 'alphatdsdy')
    call this%alphatdsdz%init(this%s, trim(s_name) // 'alphatdsdz')

    allocate(this%stat_fields%items(this%n_stats))

    call this%stat_fields%assign_to_field(1, this%alphat_mean%mf)
    call this%stat_fields%assign_to_field(2, this%alphatdsdx%mf)
    call this%stat_fields%assign_to_field(3, this%alphatdsdy%mf)
    call this%stat_fields%assign_to_field(4, this%alphatdsdz%mf)

  end subroutine scalar_sgs_stats_init_alphat

  !> Constructor. Initialize the fields associated with scalar_sgs_stats.
  !> This version uses nut field and turbulent Prandtl number.
  !! @param coef SEM coefficients. Optional.
  !! @param s_name The name of the scalar field.
  !! @param nut_field Specifies the name of the nut field.
  !! @param pr_turb Turbulent Prandtl number.
  subroutine scalar_sgs_stats_init_nut(this, coef, s_name, nut_field, pr_turb)
    class(scalar_sgs_stats_t), intent(inout), target:: this
    type(coef_t), target, optional :: coef
    character(len=*), intent(in) :: s_name
    character(len=*), intent(in) :: nut_field
    real(kind=rp), intent(in) :: pr_turb

    call this%free()
    this%coef => coef

    this%s => neko_registry%get_field(s_name)
    this%nut => neko_registry%get_field(nut_field)
    this%pr_turb = pr_turb
    this%nut_dependency = .true.

    ! Initialize mean fields. alphat = nut / pr_turb only exists in a scratch
    ! field during update, so all means are updated with explicitly passed
    ! samples. The field passed to init only provides the dofmap.
    call this%alphat_mean%init(this%nut, &
         'mean_' // trim(s_name) // 'alphat_temp')

    call this%alphatdsdx%init(this%s, trim(s_name) // 'alphatdsdx')
    call this%alphatdsdy%init(this%s, trim(s_name) // 'alphatdsdy')
    call this%alphatdsdz%init(this%s, trim(s_name) // 'alphatdsdz')

    allocate(this%stat_fields%items(this%n_stats))

    call this%stat_fields%assign_to_field(1, this%alphat_mean%mf)
    call this%stat_fields%assign_to_field(2, this%alphatdsdx%mf)
    call this%stat_fields%assign_to_field(3, this%alphatdsdy%mf)
    call this%stat_fields%assign_to_field(4, this%alphatdsdz%mf)

  end subroutine scalar_sgs_stats_init_nut

  !> Updates all fields with a new sample.
  !! @param k Time elapsed since the last update.
  subroutine scalar_sgs_stats_update(this, k)
    class(scalar_sgs_stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k
    type(field_t), pointer :: alphat, dsdx, dsdy, dsdz
    integer :: temp_indices(4)

    call neko_scratch_registry%request_field(dsdx, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(dsdy, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(dsdz, temp_indices(3), .false.)

    if (this%nut_dependency) then
       call neko_scratch_registry%request_field(alphat, temp_indices(4), &
            .false.)
       call field_cmult2(alphat, this%nut, 1.0_rp / this%pr_turb)
    else
       alphat => this%alphat
    end if
    call this%alphat_mean%update_field(k, alphat)

    call grad(dsdx%x, dsdy%x, dsdz%x, this%s%x, this%coef)

    call this%alphatdsdx%update_col3(k, alphat, dsdx)
    call this%alphatdsdy%update_col3(k, alphat, dsdy)
    call this%alphatdsdz%update_col3(k, alphat, dsdz)

    if (this%nut_dependency) then
       call neko_scratch_registry%relinquish_field(temp_indices)
    else
       call neko_scratch_registry%relinquish_field(temp_indices(1:3))
    end if
  end subroutine scalar_sgs_stats_update


  !> Destructor.
  subroutine scalar_sgs_stats_free(this)
    class(scalar_sgs_stats_t), intent(inout) :: this

    call this%alphat_mean%free()

    call this%alphatdsdx%free()
    call this%alphatdsdy%free()
    call this%alphatdsdz%free()

    nullify(this%coef)
    nullify(this%s)
    nullify(this%alphat)
    nullify(this%nut)

    call this%stat_fields%free()

  end subroutine scalar_sgs_stats_free

  !> Resets all the computed means values and sampling times to zero.
  subroutine scalar_sgs_stats_reset(this)
    class(scalar_sgs_stats_t), intent(inout), target:: this

    call this%alphat_mean%reset()

    call this%alphatdsdx%reset()
    call this%alphatdsdy%reset()
    call this%alphatdsdz%reset()

  end subroutine scalar_sgs_stats_reset

end module scalar_sgs_stats
