! Copyright (c) 2025, The Neko Authors
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
!> Computes the subgrid-scale contributions for Reynolds stresses.
!! We use the Reynolds decomposition for a field u = <u> + u' = U + u'
!! Spatial derivatives i.e. du/dx we denote dudx
module fluid_sgs_stats
  use mean_field, only : mean_field_t
  use num_types, only : rp
  use field_math, only : field_cmult, field_col3
  use operators, only : strain_rate
  use coefs, only : coef_t
  use field, only : field_t
  use field_list, only : field_list_t
  use stats_quant, only : stats_quant_t
  use registry, only : neko_registry
  implicit none
  private

  type, public, extends(stats_quant_t) :: fluid_sgs_stats_t
     !> Work fields
     type(field_t) :: stats_work

     !> Pointers to the instantenious quantities.
     type(field_t), pointer :: nut !< nut
     type(field_t), pointer :: u !< u
     type(field_t), pointer :: v !< v
     type(field_t), pointer :: w !< w

     type(mean_field_t) :: nut_mean !< <nut>

     type(mean_field_t) :: uu_sgs !< <uu_sgs> = <2nut*S11>
     type(mean_field_t) :: vv_sgs !< <vv_sgs> = <2nut*S22>
     type(mean_field_t) :: ww_sgs !< <ww_sgs> = <2nut*S33>
     type(mean_field_t) :: uv_sgs !< <uv_sgs> = <2nut*S12>
     type(mean_field_t) :: uw_sgs !< <uw_sgs> = <2nut*S13>
     type(mean_field_t) :: vw_sgs !< <vw_sgs> = <2nut*S23>

     !> gradients
     type(field_t) :: double_s11
     type(field_t) :: double_s22
     type(field_t) :: double_s33
     type(field_t) :: double_s12
     type(field_t) :: double_s13
     type(field_t) :: double_s23

     !> SEM coefficients.
     type(coef_t), pointer :: coef

     !> Specify the name of the eddy viscosity field.
     character(10) :: nut_field

     !> Number of statistical fields to be computed.
     integer :: n_stats = 7

     !> A list of size n_stats, whith entries pointing to the fields that will
     !! be output (the field components above.) Used to write the output.
     type(field_list_t) :: stat_fields
   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_sgs_stats_init
     !> Destructor.
     procedure, pass(this) :: free => fluid_sgs_stats_free
     !> Update all the mean value fields with a new sample.
     procedure, pass(this) :: update => fluid_sgs_stats_update
     !> Reset all the computed means values and sampling times to zero.
     procedure, pass(this) :: reset => fluid_sgs_stats_reset
  end type fluid_sgs_stats_t

contains

  !> Constructor. Initialize the fields associated with fluid_sgs_stats.
  !! @param coef SEM coefficients.
  !! @param u The x component of velocity.
  !! @param v The y component of velocity.
  !! @param w The z component of velocity.
  !! @param nut_field Specifies the name of the nut field.
  !! Optional, defaults to `uut`.
  subroutine fluid_sgs_stats_init(this, coef, u, v, w, nut_field)
    class(fluid_sgs_stats_t), intent(inout), target:: this
    type(coef_t), target :: coef
    type(field_t), target, intent(in) :: u, v, w
    character(*), intent(in), optional :: nut_field

    call this%free()
    this%coef => coef

    this%u => u
    this%v => v
    this%w => w

    if (present(nut_field)) then
       this%nut_field = trim(nut_field)
    else
       this%nut_field = 'nut'
    end if
    this%nut => neko_registry%get_field_by_name(this%nut_field)

    ! Initialize work fields
    call this%stats_work%init(this%u%dof, 'stats')

    ! Initialize mean fields
    call this%nut_mean%init(this%stats_work)

    call this%uu_sgs%init(this%stats_work, 'uu_sgs')
    call this%vv_sgs%init(this%stats_work, 'vv_sgs')
    call this%ww_sgs%init(this%stats_work, 'ww_sgs')
    call this%uv_sgs%init(this%stats_work, 'uv_sgs')
    call this%uw_sgs%init(this%stats_work, 'uw_sgs')
    call this%vw_sgs%init(this%stats_work, 'vw_sgs')

    call this%double_s11%init(this%u%dof, 'double_s11')
    call this%double_s22%init(this%u%dof, 'double_s22')
    call this%double_s33%init(this%u%dof, 'double_s33')
    call this%double_s12%init(this%u%dof, 'double_s12')
    call this%double_s13%init(this%u%dof, 'double_s13')
    call this%double_s23%init(this%u%dof, 'double_s23')

    allocate(this%stat_fields%items(this%n_stats))

    call this%stat_fields%assign_to_field(1, this%nut_mean%mf)

    call this%stat_fields%assign_to_field(2, this%uu_sgs%mf)
    call this%stat_fields%assign_to_field(3, this%vv_sgs%mf)
    call this%stat_fields%assign_to_field(4, this%ww_sgs%mf)
    call this%stat_fields%assign_to_field(5, this%uv_sgs%mf)
    call this%stat_fields%assign_to_field(6, this%uw_sgs%mf)
    call this%stat_fields%assign_to_field(7, this%vw_sgs%mf)

  end subroutine fluid_sgs_stats_init

  !> Updates all fields with a new sample.
  !! @param k Time elapsed since the last update.
  subroutine fluid_sgs_stats_update(this, k)
    class(fluid_sgs_stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k
    integer :: n

    associate(stats_work => this%stats_work)
      n = stats_work%dof%size()

      call this%nut_mean%update(k)

      call strain_rate(this%double_s11%x, &
                       this%double_s22%x, &
                       this%double_s33%x, &
                       this%double_s12%x, &
                       this%double_s13%x, &
                       this%double_s23%x, this%u, this%v, this%w, this%coef)

      ! form the double sij tensor
      call field_cmult(this%double_s11, 2.0_rp)
      call field_cmult(this%double_s22, 2.0_rp)
      call field_cmult(this%double_s33, 2.0_rp)
      call field_cmult(this%double_s12, 2.0_rp)
      call field_cmult(this%double_s13, 2.0_rp)
      call field_cmult(this%double_s23, 2.0_rp)

      call field_col3(stats_work, this%nut, this%double_s11)
      call this%uu_sgs%update(k)
      call field_col3(stats_work, this%nut, this%double_s22)
      call this%vv_sgs%update(k)
      call field_col3(stats_work, this%nut, this%double_s33)
      call this%ww_sgs%update(k)
      call field_col3(stats_work, this%nut, this%double_s12)
      call this%uv_sgs%update(k)
      call field_col3(stats_work, this%nut, this%double_s13)
      call this%uw_sgs%update(k)
      call field_col3(stats_work, this%nut, this%double_s23)
      call this%vw_sgs%update(k)

    end associate

  end subroutine fluid_sgs_stats_update


  !> Destructor.
  subroutine fluid_sgs_stats_free(this)
    class(fluid_sgs_stats_t), intent(inout) :: this

    call this%stats_work%free()

    call this%nut_mean%free()

    call this%uu_sgs%free()
    call this%vv_sgs%free()
    call this%ww_sgs%free()
    call this%uv_sgs%free()
    call this%uw_sgs%free()
    call this%vw_sgs%free()

    call this%double_s11%free()
    call this%double_s22%free()
    call this%double_s33%free()
    call this%double_s12%free()
    call this%double_s13%free()
    call this%double_s23%free()

    nullify(this%coef)
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%nut)

  end subroutine fluid_sgs_stats_free

  !> Resets all the computed means values and sampling times to zero.
  subroutine fluid_sgs_stats_reset(this)
    class(fluid_sgs_stats_t), intent(inout), target:: this

    call this%nut_mean%reset()

    call this%uu_sgs%reset()
    call this%vv_sgs%reset()
    call this%ww_sgs%reset()
    call this%uv_sgs%reset()
    call this%uw_sgs%reset()
    call this%vw_sgs%reset()

  end subroutine fluid_sgs_stats_reset

end module fluid_sgs_stats
