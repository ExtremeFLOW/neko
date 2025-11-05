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
!> Computes various statistics for the scalar fields.
!! We use the Reynolds decomposition for a field u = <u> + u' = U + u'
!! Spatial derivatives i.e. du/dx we denote dudx
module scalar_stats
  use mean_field, only : mean_field_t
  use device_math, only : device_col3, device_col2, device_cfill, &
       device_invcol2, device_addcol3
  use num_types, only : rp
  use math, only : invers2, col2, addcol3, col3, copy, subcol3
  use operators, only : opgrad
  use coefs, only : coef_t
  use field, only : field_t
  use field_list, only : field_list_t
  use stats_quant, only : stats_quant_t
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_warning
  implicit none
  private

  type, public, extends(stats_quant_t) :: scalar_stats_t
     !> Work fields
     type(field_t) :: stats_ss
     type(field_t) :: stats_uiuj
     type(field_t) :: stats_work

     !> Pointers to the instantenious quantities.
     type(field_t), pointer :: s !< s
     type(field_t), pointer :: u !< u
     type(field_t), pointer :: v !< v
     type(field_t), pointer :: w !< w
     type(field_t), pointer :: p !< p

     type(mean_field_t) :: s_mean !< <s>

     type(mean_field_t) :: us !< <us>
     type(mean_field_t) :: vs !< <vs>
     type(mean_field_t) :: ws !< <ws>

     type(mean_field_t) :: ss !< <ss>
     type(mean_field_t) :: sss !< <sss>
     type(mean_field_t) :: ssss !< <ssss>

     type(mean_field_t) :: uss !< <uss>
     type(mean_field_t) :: vss !< <vss>
     type(mean_field_t) :: wss !< <wss>

     type(mean_field_t) :: uus !< <uus>
     type(mean_field_t) :: vvs !< <vvs>
     type(mean_field_t) :: wws !< <wws>
     type(mean_field_t) :: uvs !< <uvs>
     type(mean_field_t) :: uws !< <uws>
     type(mean_field_t) :: vws !< <vws>

     type(mean_field_t) :: ps !< <ps>

     type(mean_field_t) :: pdsdx
     type(mean_field_t) :: pdsdy
     type(mean_field_t) :: pdsdz

     type(mean_field_t) :: udsdx
     type(mean_field_t) :: udsdy
     type(mean_field_t) :: udsdz
     type(mean_field_t) :: vdsdx
     type(mean_field_t) :: vdsdy
     type(mean_field_t) :: vdsdz
     type(mean_field_t) :: wdsdx
     type(mean_field_t) :: wdsdy
     type(mean_field_t) :: wdsdz

     type(mean_field_t) :: sdudx
     type(mean_field_t) :: sdudy
     type(mean_field_t) :: sdudz
     type(mean_field_t) :: sdvdx
     type(mean_field_t) :: sdvdy
     type(mean_field_t) :: sdvdz
     type(mean_field_t) :: sdwdx
     type(mean_field_t) :: sdwdy
     type(mean_field_t) :: sdwdz

     type(mean_field_t) :: ess !< <dsdx**2 + dsdy**2 + dsdz**2>
     type(mean_field_t) :: eus !< <dudx*dsdx + dudy*dsdy + dudz*dsdz>
     type(mean_field_t) :: evs !< <dvdx*dsdx + dvdy*dsdy + dvdz*dsdz>
     type(mean_field_t) :: ews !< <dwdx*dsdx + dwdy*dsdy + dwdz*dsdz>

     !> gradients
     type(field_t) :: dsdx
     type(field_t) :: dsdy
     type(field_t) :: dsdz
     type(field_t) :: dudx
     type(field_t) :: dudy
     type(field_t) :: dudz
     type(field_t) :: dvdx
     type(field_t) :: dvdy
     type(field_t) :: dvdz
     type(field_t) :: dwdx
     type(field_t) :: dwdy
     type(field_t) :: dwdz

     !> SEM coefficients.
     type(coef_t), pointer :: coef
     !> Number of statistical fields to be computed.
     integer :: n_stats = 42
     !> Specifies a subset of the statistics to be collected. All 42 fields by
     !! default.
     character(5) :: stat_set
     !> A list of size n_stats, whith entries pointing to the fields that will
     !! be output (the field components above.) Used to write the output.
     type(field_list_t) :: stat_fields
   contains
     !> Constructor.
     procedure, pass(this) :: init => scalar_stats_init
     !> Destructor.
     procedure, pass(this) :: free => scalar_stats_free
     !> Update all the mean value fields with a new sample.
     procedure, pass(this) :: update => scalar_stats_update
     !> Reset all the computed means values and sampling times to zero.
     procedure, pass(this) :: reset => scalar_stats_reset
     ! Convert computed weak gradients to strong.
     procedure, pass(this) :: make_strong_grad => scalar_stats_make_strong_grad
  end type scalar_stats_t

contains

  !> Constructor. Initialize the fields associated with scalar_stats.
  !! @param coef SEM coefficients. Optional.
  !! @param s The scalar.
  !! @param u The x component of velocity.
  !! @param v The y component of velocity.
  !! @param w The z component of velocity.
  !! @param p The pressure.
  !! @param set Specifies the subset of the statistics to be collected.
  !! Optional. Either `basic` or `full`, defaults to `full`.
  subroutine scalar_stats_init(this, coef, s, u, v, w, p, set)
    class(scalar_stats_t), intent(inout), target:: this
    type(coef_t), target, optional :: coef
    type(field_t), target, intent(in) :: s, u, v, w, p
    character(*), intent(in), optional :: set

    call this%free()
    this%coef => coef

    this%s => s
    this%u => u
    this%v => v
    this%w => w
    this%p => p

    if (present(set)) then
       this%stat_set = trim(set)
       if (this%stat_set .eq. 'basic') then
          this%n_stats = 5
       end if
    else
       this%stat_set = 'full'
       this%n_stats = 42
    end if

    ! Initialize work fields
    call this%stats_work%init(this%u%dof, 'stats')
    call this%stats_ss%init(this%u%dof, 'ss temp')

    ! Initialize mean fields
    call this%s_mean%init(this%s)

    call this%us%init(this%stats_work, 'us')
    call this%vs%init(this%stats_work, 'vs')
    call this%ws%init(this%stats_work, 'ws')

    call this%ss%init(this%stats_ss, 'ss')

    if (this%n_stats .eq. 42) then
       ! Initialize req. work fields
       call this%stats_uiuj%init(this%u%dof, 'uiuj temp')

       ! Initialize req. gradient fields
       call this%dsdx%init(this%u%dof, 'dsdx')
       call this%dsdy%init(this%u%dof, 'dsdy')
       call this%dsdz%init(this%u%dof, 'dsdz')

       call this%dudx%init(this%u%dof, 'dudx')
       call this%dudy%init(this%u%dof, 'dudy')
       call this%dudz%init(this%u%dof, 'dudz')
       call this%dvdx%init(this%u%dof, 'dvdx')
       call this%dvdy%init(this%u%dof, 'dvdy')
       call this%dvdz%init(this%u%dof, 'dvdz')
       call this%dwdx%init(this%u%dof, 'dwdx')
       call this%dwdy%init(this%u%dof, 'dwdy')
       call this%dwdz%init(this%u%dof, 'dwdz')

       ! Initialize req. mean fields
       call this%sss%init(this%stats_work, 'sss')
       call this%ssss%init(this%stats_work, 'ssss')

       call this%uss%init(this%stats_work, 'uss')
       call this%vss%init(this%stats_work, 'vss')
       call this%wss%init(this%stats_work, 'wss')

       call this%uus%init(this%stats_work, 'uus')
       call this%vvs%init(this%stats_work, 'vvs')
       call this%wws%init(this%stats_work, 'wws')
       call this%uvs%init(this%stats_work, 'uvs')
       call this%uws%init(this%stats_work, 'uws')
       call this%vws%init(this%stats_work, 'vws')

       call this%ps%init(this%stats_work, 'ps')

       call this%pdsdx%init(this%stats_work, 'pdsdx')
       call this%pdsdy%init(this%stats_work, 'pdsdy')
       call this%pdsdz%init(this%stats_work, 'pdsdz')

       call this%udsdx%init(this%stats_work, 'udsdx')
       call this%udsdy%init(this%stats_work, 'udsdy')
       call this%udsdz%init(this%stats_work, 'udsdz')
       call this%vdsdx%init(this%stats_work, 'vdsdx')
       call this%vdsdy%init(this%stats_work, 'vdsdy')
       call this%vdsdz%init(this%stats_work, 'vdsdz')
       call this%wdsdx%init(this%stats_work, 'wdsdx')
       call this%wdsdy%init(this%stats_work, 'wdsdy')
       call this%wdsdz%init(this%stats_work, 'wdsdz')

       call this%sdudx%init(this%stats_work, 'sdudx')
       call this%sdudy%init(this%stats_work, 'sdudy')
       call this%sdudz%init(this%stats_work, 'sdudz')
       call this%sdvdx%init(this%stats_work, 'sdvdx')
       call this%sdvdy%init(this%stats_work, 'sdvdy')
       call this%sdvdz%init(this%stats_work, 'sdvdz')
       call this%sdwdx%init(this%stats_work, 'sdwdx')
       call this%sdwdy%init(this%stats_work, 'sdwdy')
       call this%sdwdz%init(this%stats_work, 'sdwdz')

       call this%ess%init(this%stats_work, 'ess')
       call this%eus%init(this%stats_work, 'eus')
       call this%evs%init(this%stats_work, 'evs')
       call this%ews%init(this%stats_work, 'ews')
    end if

    allocate(this%stat_fields%items(this%n_stats))

    call this%stat_fields%assign_to_field(1, this%s_mean%mf)

    call this%stat_fields%assign_to_field(2, this%us%mf)
    call this%stat_fields%assign_to_field(3, this%vs%mf)
    call this%stat_fields%assign_to_field(4, this%ws%mf)

    call this%stat_fields%assign_to_field(5, this%ss%mf)

    if (this%n_stats .eq. 42) then
       call this%stat_fields%assign_to_field(6, this%sss%mf)
       call this%stat_fields%assign_to_field(7, this%ssss%mf)

       call this%stat_fields%assign_to_field(8, this%uss%mf)
       call this%stat_fields%assign_to_field(9, this%vss%mf)
       call this%stat_fields%assign_to_field(10, this%wss%mf)
       call this%stat_fields%assign_to_field(11, this%uus%mf)
       call this%stat_fields%assign_to_field(12, this%vvs%mf)
       call this%stat_fields%assign_to_field(13, this%wws%mf)
       call this%stat_fields%assign_to_field(14, this%uvs%mf)
       call this%stat_fields%assign_to_field(15, this%uws%mf)
       call this%stat_fields%assign_to_field(16, this%vws%mf)

       call this%stat_fields%assign_to_field(17, this%ps%mf)

       call this%stat_fields%assign_to_field(18, this%pdsdx%mf)
       call this%stat_fields%assign_to_field(19, this%pdsdy%mf)
       call this%stat_fields%assign_to_field(20, this%pdsdz%mf)

       call this%stat_fields%assign_to_field(21, this%udsdx%mf)
       call this%stat_fields%assign_to_field(22, this%udsdy%mf)
       call this%stat_fields%assign_to_field(23, this%udsdz%mf)
       call this%stat_fields%assign_to_field(24, this%vdsdx%mf)
       call this%stat_fields%assign_to_field(25, this%vdsdy%mf)
       call this%stat_fields%assign_to_field(26, this%vdsdz%mf)
       call this%stat_fields%assign_to_field(27, this%wdsdx%mf)
       call this%stat_fields%assign_to_field(28, this%wdsdy%mf)
       call this%stat_fields%assign_to_field(29, this%wdsdz%mf)

       call this%stat_fields%assign_to_field(30, this%sdudx%mf)
       call this%stat_fields%assign_to_field(31, this%sdudy%mf)
       call this%stat_fields%assign_to_field(32, this%sdudz%mf)
       call this%stat_fields%assign_to_field(33, this%sdvdx%mf)
       call this%stat_fields%assign_to_field(34, this%sdvdy%mf)
       call this%stat_fields%assign_to_field(35, this%sdvdz%mf)
       call this%stat_fields%assign_to_field(36, this%sdwdx%mf)
       call this%stat_fields%assign_to_field(37, this%sdwdy%mf)
       call this%stat_fields%assign_to_field(38, this%sdwdz%mf)

       call this%stat_fields%assign_to_field(39, this%ess%mf)
       call this%stat_fields%assign_to_field(40, this%eus%mf)
       call this%stat_fields%assign_to_field(41, this%evs%mf)
       call this%stat_fields%assign_to_field(42, this%ews%mf)
    end if

  end subroutine scalar_stats_init

  !> Updates all fields with a new sample.
  !! @param k Time elapsed since the last update.
  subroutine scalar_stats_update(this, k)
    class(scalar_stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k
    integer :: n

    associate(stats_work => this%stats_work, stats_ss => this%stats_ss, &
         stats_uiuj => this%stats_uiuj)
      n = stats_work%dof%size()

      !> U%f is u and U%mf is <u>
      if (NEKO_BCKND_DEVICE .eq. 1) then

         call this%s_mean%update(k)

         call device_col3(stats_work%x_d, this%u%x_d, this%s%x_d, n)
         call this%us%update(k)
         call device_col3(stats_work%x_d, this%v%x_d, this%s%x_d, n)
         call this%vs%update(k)
         call device_col3(stats_work%x_d, this%w%x_d, this%s%x_d, n)
         call this%ws%update(k)

         call device_col3(stats_ss%x_d, this%s%x_d, this%s%x_d, n)
         call this%ss%update(k)

         if (this%n_stats .eq. 5) return

         call device_col3(stats_work%x_d, this%stats_ss%x_d, this%s%x_d, n)
         call this%sss%update(k)
         call device_col3(stats_work%x_d, this%stats_ss%x_d, this%stats_ss%x_d, n)
         call this%ssss%update(k)

         call device_col3(stats_work%x_d, this%u%x_d, this%stats_ss%x_d, n)
         call this%uss%update(k)
         call device_col3(stats_work%x_d, this%v%x_d, this%stats_ss%x_d, n)
         call this%vss%update(k)
         call device_col3(stats_work%x_d, this%w%x_d, this%stats_ss%x_d, n)
         call this%wss%update(k)

         call device_col3(stats_uiuj%x_d, this%u%x_d, this%u%x_d, n)
         call device_col3(stats_work%x_d, stats_uiuj%x_d, this%s%x_d, n)
         call this%uus%update(k)
         call device_col3(stats_uiuj%x_d, this%v%x_d, this%v%x_d, n)
         call device_col3(stats_work%x_d, stats_uiuj%x_d, this%s%x_d, n)
         call this%vvs%update(k)
         call device_col3(stats_uiuj%x_d, this%w%x_d, this%w%x_d, n)
         call device_col3(stats_work%x_d, stats_uiuj%x_d, this%s%x_d, n)
         call this%wws%update(k)
         call device_col3(stats_uiuj%x_d, this%u%x_d, this%v%x_d, n)
         call device_col3(stats_work%x_d, stats_uiuj%x_d, this%s%x_d, n)
         call this%uvs%update(k)
         call device_col3(stats_uiuj%x_d, this%u%x_d, this%w%x_d, n)
         call device_col3(stats_work%x_d, stats_uiuj%x_d, this%s%x_d, n)
         call this%uws%update(k)
         call device_col3(stats_uiuj%x_d, this%v%x_d, this%w%x_d, n)
         call device_col3(stats_work%x_d, stats_uiuj%x_d, this%s%x_d, n)
         call this%vws%update(k)

         call device_col3(stats_work%x_d, this%p%x_d, this%s%x_d, n)
         call this%ps%update(k)

      else

         call this%s_mean%update(k)

         call col3(stats_work%x, this%u%x, this%s%x, n)
         call this%us%update(k)
         call col3(stats_work%x, this%v%x, this%s%x, n)
         call this%vs%update(k)
         call col3(stats_work%x, this%w%x, this%s%x, n)
         call this%ws%update(k)

         call col3(stats_ss%x, this%s%x, this%s%x, n)
         call this%ss%update(k)

         if (this%n_stats .eq. 5) return

         call col3(stats_work%x, this%stats_ss%x, this%s%x, n)
         call this%sss%update(k)
         call col3(stats_work%x, this%stats_ss%x, this%stats_ss%x, n)
         call this%ssss%update(k)

         call col3(stats_work%x, this%u%x, this%stats_ss%x, n)
         call this%uss%update(k)
         call col3(stats_work%x, this%v%x, this%stats_ss%x, n)
         call this%vss%update(k)
         call col3(stats_work%x, this%w%x, this%stats_ss%x, n)
         call this%wss%update(k)

         call col3(stats_uiuj%x, this%u%x, this%u%x, n)
         call col3(stats_work%x, stats_uiuj%x, this%s%x, n)
         call this%uus%update(k)
         call col3(stats_uiuj%x, this%v%x, this%v%x, n)
         call col3(stats_work%x, stats_uiuj%x, this%s%x, n)
         call this%vvs%update(k)
         call col3(stats_uiuj%x, this%w%x, this%w%x, n)
         call col3(stats_work%x, stats_uiuj%x, this%s%x, n)
         call this%wws%update(k)
         call col3(stats_uiuj%x, this%u%x, this%v%x, n)
         call col3(stats_work%x, stats_uiuj%x, this%s%x, n)
         call this%uvs%update(k)
         call col3(stats_uiuj%x, this%u%x, this%w%x, n)
         call col3(stats_work%x, stats_uiuj%x, this%s%x, n)
         call this%uws%update(k)
         call col3(stats_uiuj%x, this%v%x, this%w%x, n)
         call col3(stats_work%x, stats_uiuj%x, this%s%x, n)
         call this%vws%update(k)

         call col3(stats_work%x, this%p%x, this%s%x, n)
         call this%ps%update(k)

      end if

      call opgrad(this%dsdx%x, this%dsdy%x, this%dsdz%x, this%s%x, this%coef)
      call opgrad(this%dudx%x, this%dudy%x, this%dudz%x, this%u%x, this%coef)
      call opgrad(this%dvdx%x, this%dvdy%x, this%dvdz%x, this%v%x, this%coef)
      call opgrad(this%dwdx%x, this%dwdy%x, this%dwdz%x, this%w%x, this%coef)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col3(stats_work%x_d, this%dsdx%x_d, this%p%x_d, n)
         call this%pdsdx%update(k)
         call device_col3(stats_work%x_d, this%dsdy%x_d, this%p%x_d, n)
         call this%pdsdy%update(k)
         call device_col3(stats_work%x_d, this%dsdz%x_d, this%p%x_d, n)
         call this%pdsdz%update(k)

         call device_col3(stats_work%x_d, this%dsdx%x_d, this%u%x_d, n)
         call this%udsdx%update(k)
         call device_col3(stats_work%x_d, this%dsdy%x_d, this%u%x_d, n)
         call this%udsdy%update(k)
         call device_col3(stats_work%x_d, this%dsdz%x_d, this%u%x_d, n)
         call this%udsdz%update(k)

         call device_col3(stats_work%x_d, this%dsdx%x_d, this%v%x_d, n)
         call this%vdsdx%update(k)
         call device_col3(stats_work%x_d, this%dsdy%x_d, this%v%x_d, n)
         call this%vdsdy%update(k)
         call device_col3(stats_work%x_d, this%dsdz%x_d, this%v%x_d, n)
         call this%vdsdz%update(k)

         call device_col3(stats_work%x_d, this%dsdx%x_d, this%w%x_d, n)
         call this%wdsdx%update(k)
         call device_col3(stats_work%x_d, this%dsdy%x_d, this%w%x_d, n)
         call this%wdsdy%update(k)
         call device_col3(stats_work%x_d, this%dsdz%x_d, this%w%x_d, n)
         call this%wdsdz%update(k)

         call device_col3(stats_work%x_d, this%dudx%x_d, this%s%x_d, n)
         call this%sdudx%update(k)
         call device_col3(stats_work%x_d, this%dudy%x_d, this%s%x_d, n)
         call this%sdudy%update(k)
         call device_col3(stats_work%x_d, this%dudz%x_d, this%s%x_d, n)
         call this%sdudz%update(k)

         call device_col3(stats_work%x_d, this%dvdx%x_d, this%s%x_d, n)
         call this%sdvdx%update(k)
         call device_col3(stats_work%x_d, this%dvdy%x_d, this%s%x_d, n)
         call this%sdvdy%update(k)
         call device_col3(stats_work%x_d, this%dvdz%x_d, this%s%x_d, n)
         call this%sdvdz%update(k)

         call device_col3(stats_work%x_d, this%dwdx%x_d, this%s%x_d, n)
         call this%sdwdx%update(k)
         call device_col3(stats_work%x_d, this%dwdy%x_d, this%s%x_d, n)
         call this%sdwdy%update(k)
         call device_col3(stats_work%x_d, this%dwdz%x_d, this%s%x_d, n)
         call this%sdwdz%update(k)

         call device_col3(stats_work%x_d, this%dsdx%x_d, &
              this%dsdx%x_d, n)
         call device_addcol3(stats_work%x_d, this%dsdy%x_d, &
              this%dsdy%x_d, n)
         call device_addcol3(stats_work%x_d, this%dsdz%x_d, &
              this%dsdz%x_d, n)
         call this%ess%update(k)

         call device_col3(stats_work%x_d, this%dudx%x_d, &
              this%dsdx%x_d, n)
         call device_addcol3(stats_work%x_d, this%dudy%x_d, &
              this%dsdy%x_d, n)
         call device_addcol3(stats_work%x_d, this%dudz%x_d, &
              this%dsdz%x_d, n)
         call this%eus%update(k)

         call device_col3(stats_work%x_d, this%dvdx%x_d, &
              this%dsdx%x_d, n)
         call device_addcol3(stats_work%x_d, this%dvdy%x_d, &
              this%dsdy%x_d, n)
         call device_addcol3(stats_work%x_d, this%dvdz%x_d, &
              this%dsdz%x_d, n)
         call this%evs%update(k)

         call device_col3(stats_work%x_d, this%dwdx%x_d, &
              this%dsdx%x_d, n)
         call device_addcol3(stats_work%x_d, this%dwdy%x_d, &
              this%dsdy%x_d, n)
         call device_addcol3(stats_work%x_d, this%dwdz%x_d, &
              this%dsdz%x_d, n)
         call this%ews%update(k)

      else

         call col3(stats_work%x, this%dsdx%x, this%p%x, n)
         call this%pdsdx%update(k)
         call col3(stats_work%x, this%dsdy%x, this%p%x, n)
         call this%pdsdy%update(k)
         call col3(stats_work%x, this%dsdz%x, this%p%x, n)
         call this%pdsdz%update(k)

         call col3(stats_work%x, this%dsdx%x, this%u%x, n)
         call this%udsdx%update(k)
         call col3(stats_work%x, this%dsdy%x, this%u%x, n)
         call this%udsdy%update(k)
         call col3(stats_work%x, this%dsdz%x, this%u%x, n)
         call this%udsdz%update(k)

         call col3(stats_work%x, this%dsdx%x, this%v%x, n)
         call this%vdsdx%update(k)
         call col3(stats_work%x, this%dsdy%x, this%v%x, n)
         call this%vdsdy%update(k)
         call col3(stats_work%x, this%dsdz%x, this%v%x, n)
         call this%vdsdz%update(k)

         call col3(stats_work%x, this%dsdx%x, this%w%x, n)
         call this%wdsdx%update(k)
         call col3(stats_work%x, this%dsdy%x, this%w%x, n)
         call this%wdsdy%update(k)
         call col3(stats_work%x, this%dsdz%x, this%w%x, n)
         call this%wdsdz%update(k)

         call col3(stats_work%x, this%dudx%x, this%s%x, n)
         call this%sdudx%update(k)
         call col3(stats_work%x, this%dudy%x, this%s%x, n)
         call this%sdudy%update(k)
         call col3(stats_work%x, this%dudz%x, this%s%x, n)
         call this%sdudz%update(k)

         call col3(stats_work%x, this%dvdx%x, this%s%x, n)
         call this%sdvdx%update(k)
         call col3(stats_work%x, this%dvdy%x, this%s%x, n)
         call this%sdvdy%update(k)
         call col3(stats_work%x, this%dvdz%x, this%s%x, n)
         call this%sdvdz%update(k)

         call col3(stats_work%x, this%dwdx%x, this%s%x, n)
         call this%sdwdx%update(k)
         call col3(stats_work%x, this%dwdy%x, this%s%x, n)
         call this%sdwdy%update(k)
         call col3(stats_work%x, this%dwdz%x, this%s%x, n)
         call this%sdwdz%update(k)

         call col3(stats_work%x, this%dsdx%x, &
              this%dsdx%x, n)
         call addcol3(stats_work%x, this%dsdy%x, &
              this%dsdy%x, n)
         call addcol3(stats_work%x, this%dsdz%x, &
              this%dsdz%x, n)
         call this%ess%update(k)

         call col3(stats_work%x, this%dudx%x, &
              this%dsdx%x, n)
         call addcol3(stats_work%x, this%dudy%x, &
              this%dsdy%x, n)
         call addcol3(stats_work%x, this%dudz%x, &
              this%dsdz%x, n)
         call this%eus%update(k)

         call col3(stats_work%x, this%dvdx%x, &
              this%dsdx%x, n)
         call addcol3(stats_work%x, this%dvdy%x, &
              this%dsdy%x, n)
         call addcol3(stats_work%x, this%dvdz%x, &
              this%dsdz%x, n)
         call this%evs%update(k)

         call col3(stats_work%x, this%dwdx%x, &
              this%dsdx%x, n)
         call addcol3(stats_work%x, this%dwdy%x, &
              this%dsdy%x, n)
         call addcol3(stats_work%x, this%dwdz%x, &
              this%dsdz%x, n)
         call this%ews%update(k)

      end if
    end associate

  end subroutine scalar_stats_update


  !> Destructor.
  subroutine scalar_stats_free(this)
    class(scalar_stats_t), intent(inout) :: this

    call this%stats_work%free()
    call this%stats_ss%free()
    call this%stats_uiuj%free()

    call this%s_mean%free()

    call this%us%free()
    call this%vs%free()
    call this%ws%free()

    call this%ss%free()

    if (this%n_stats .eq. 42) then
       call this%sss%free()
       call this%ssss%free()

       call this%uss%free()
       call this%vss%free()
       call this%wss%free()

       call this%uus%free()
       call this%vvs%free()
       call this%wws%free()
       call this%uvs%free()
       call this%uws%free()
       call this%vws%free()

       call this%ps%free()

       call this%pdsdx%free()
       call this%pdsdy%free()
       call this%pdsdz%free()

       call this%udsdx%free()
       call this%udsdy%free()
       call this%udsdz%free()
       call this%vdsdx%free()
       call this%vdsdy%free()
       call this%vdsdz%free()
       call this%wdsdx%free()
       call this%wdsdy%free()
       call this%wdsdz%free()

       call this%sdudx%free()
       call this%sdudy%free()
       call this%sdudz%free()
       call this%sdvdx%free()
       call this%sdvdy%free()
       call this%sdvdz%free()
       call this%sdwdx%free()
       call this%sdwdy%free()
       call this%sdwdz%free()

       call this%ess%free()
       call this%eus%free()
       call this%evs%free()
       call this%ews%free()

       call this%dsdx%free()
       call this%dsdy%free()
       call this%dsdz%free()
       call this%dudx%free()
       call this%dudy%free()
       call this%dudz%free()
       call this%dvdx%free()
       call this%dvdy%free()
       call this%dvdz%free()
       call this%dwdx%free()
       call this%dwdy%free()
       call this%dwdz%free()
    end if

    nullify(this%coef)
    nullify(this%s)
    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)

  end subroutine scalar_stats_free

  !> Resets all the computed means values and sampling times to zero.
  subroutine scalar_stats_reset(this)
    class(scalar_stats_t), intent(inout), target:: this

    call this%s_mean%reset()

    call this%us%reset()
    call this%vs%reset()
    call this%ws%reset()

    call this%ss%reset()

    if (this%n_stats .eq. 42) then
       call this%sss%reset()
       call this%ssss%reset()

       call this%uss%reset()
       call this%vss%reset()
       call this%wss%reset()

       call this%uus%reset()
       call this%vvs%reset()
       call this%wws%reset()
       call this%uvs%reset()
       call this%uws%reset()
       call this%vws%reset()

       call this%ps%reset()

       call this%pdsdx%reset()
       call this%pdsdy%reset()
       call this%pdsdz%reset()

       call this%udsdx%reset()
       call this%udsdy%reset()
       call this%udsdz%reset()
       call this%vdsdx%reset()
       call this%vdsdy%reset()
       call this%vdsdz%reset()
       call this%wdsdx%reset()
       call this%wdsdy%reset()
       call this%wdsdz%reset()

       call this%sdudx%reset()
       call this%sdudy%reset()
       call this%sdudz%reset()
       call this%sdvdx%reset()
       call this%sdvdy%reset()
       call this%sdvdz%reset()
       call this%sdwdx%reset()
       call this%sdwdy%reset()
       call this%sdwdz%reset()

       call this%ess%reset()
       call this%eus%reset()
       call this%evs%reset()
       call this%ews%reset()
    end if

  end subroutine scalar_stats_reset

  ! Convert computed weak gradients to strong.
  subroutine scalar_stats_make_strong_grad(this)
    class(scalar_stats_t) :: this
    integer :: n

    if (this%n_stats .eq. 5) return

    n = size(this%coef%B)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(this%stats_work%x_d, 1.0_rp, n)
       call device_invcol2(this%stats_work%x_d, this%coef%B_d, n)

       call device_col2(this%pdsdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdsdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdsdz%mf%x_d, this%stats_work%x_d, n)

       call device_col2(this%udsdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%udsdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%udsdz%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%vdsdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%vdsdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%vdsdz%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%wdsdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%wdsdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%wdsdz%mf%x_d, this%stats_work%x_d, n)

       call device_col2(this%sdudx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%sdudy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%sdudz%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%sdvdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%sdvdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%sdvdz%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%sdwdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%sdwdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%sdwdz%mf%x_d, this%stats_work%x_d, n)

       call device_col2(this%stats_work%x_d, this%stats_work%x_d, n)

       call device_col2(this%ess%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%eus%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%evs%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%ews%mf%x_d, this%stats_work%x_d, n)

    else

       call invers2(this%stats_work%x, this%coef%B, n)

       call col2(this%pdsdx%mf%x, this%stats_work%x, n)
       call col2(this%pdsdy%mf%x, this%stats_work%x, n)
       call col2(this%pdsdz%mf%x, this%stats_work%x, n)

       call col2(this%udsdx%mf%x, this%stats_work%x, n)
       call col2(this%udsdy%mf%x, this%stats_work%x, n)
       call col2(this%udsdz%mf%x, this%stats_work%x, n)
       call col2(this%vdsdx%mf%x, this%stats_work%x, n)
       call col2(this%vdsdy%mf%x, this%stats_work%x, n)
       call col2(this%vdsdz%mf%x, this%stats_work%x, n)
       call col2(this%wdsdx%mf%x, this%stats_work%x, n)
       call col2(this%wdsdy%mf%x, this%stats_work%x, n)
       call col2(this%wdsdz%mf%x, this%stats_work%x, n)

       call col2(this%sdudx%mf%x, this%stats_work%x, n)
       call col2(this%sdudy%mf%x, this%stats_work%x, n)
       call col2(this%sdudz%mf%x, this%stats_work%x, n)
       call col2(this%sdvdx%mf%x, this%stats_work%x, n)
       call col2(this%sdvdy%mf%x, this%stats_work%x, n)
       call col2(this%sdvdz%mf%x, this%stats_work%x, n)
       call col2(this%sdwdx%mf%x, this%stats_work%x, n)
       call col2(this%sdwdy%mf%x, this%stats_work%x, n)
       call col2(this%sdwdz%mf%x, this%stats_work%x, n)

       call col2(this%stats_work%x, this%stats_work%x, n)

       call col2(this%ess%mf%x, this%stats_work%x, n)
       call col2(this%eus%mf%x, this%stats_work%x, n)
       call col2(this%evs%mf%x, this%stats_work%x, n)
       call col2(this%ews%mf%x, this%stats_work%x, n)

    end if

  end subroutine scalar_stats_make_strong_grad

end module scalar_stats
