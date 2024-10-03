! Copyright (c) 2022-2024, The Neko Authors
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
!> Computes various statistics for the fluid fields.
!! We use the Reynolds decomposition for a field u = <u> + u' = U + u'
!! Spatial derivatives i.e. du/dx we denote dudx
module fluid_stats
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
  use device
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_warning
  implicit none
  private

  type, public, extends(stats_quant_t) :: fluid_stats_t
     !> Work fields
     type(field_t) :: stats_u
     type(field_t) :: stats_v
     type(field_t) :: stats_w
     type(field_t) :: stats_p
     type(field_t) :: stats_work

     !> Pointers to the instantenious quantities.
     type(field_t), pointer :: u !< u
     type(field_t), pointer :: v !< v
     type(field_t), pointer :: w !< w
     type(field_t), pointer :: p !< p

     type(mean_field_t) :: u_mean !< <u>
     type(mean_field_t) :: v_mean !< <v>
     type(mean_field_t) :: w_mean !< <w>
     type(mean_field_t) :: p_mean !< <p>
     !> Velocity squares
     type(mean_field_t) :: uu !< <uu>
     type(mean_field_t) :: vv !< <vv>
     type(mean_field_t) :: ww !< <ww>
     type(mean_field_t) :: uv !< <uv>
     type(mean_field_t) :: uw !< <uw>
     type(mean_field_t) :: vw !< <vw>
     !> Velocity cubes
     type(mean_field_t) :: uuu !< <uuu>
     type(mean_field_t) :: vvv !< <vvv>
     type(mean_field_t) :: www !< <www>
     type(mean_field_t) :: uuv !< <uuv>
     type(mean_field_t) :: uuw !< <uuw>
     type(mean_field_t) :: uvv !< <uvv>
     type(mean_field_t) :: uvw !< <uvv>
     type(mean_field_t) :: vvw !< <vvw>
     type(mean_field_t) :: uww !< <uww>
     type(mean_field_t) :: vww !< <vww>
     !> Velocity squares squared
     type(mean_field_t) :: uuuu !< <uuuu>
     type(mean_field_t) :: vvvv !< <vvvv>
     type(mean_field_t) :: wwww !< <wwww>
     !> Pressure
     type(mean_field_t) :: pp !< <pp>
     type(mean_field_t) :: ppp !< <ppp>
     type(mean_field_t) :: pppp !< <pppp>
     !> Pressure * velocity
     type(mean_field_t) :: pu !< <pu>
     type(mean_field_t) :: pv !< <pv>
     type(mean_field_t) :: pw !< <pw>

     !> Derivatives
     type(mean_field_t) :: pdudx
     type(mean_field_t) :: pdudy
     type(mean_field_t) :: pdudz
     type(mean_field_t) :: pdvdx
     type(mean_field_t) :: pdvdy
     type(mean_field_t) :: pdvdz
     type(mean_field_t) :: pdwdx
     type(mean_field_t) :: pdwdy
     type(mean_field_t) :: pdwdz

     !> Combinations of sums of duvwdxyz*duvwdxyz
     type(mean_field_t) :: e11
     type(mean_field_t) :: e22
     type(mean_field_t) :: e33
     type(mean_field_t) :: e12
     type(mean_field_t) :: e13
     type(mean_field_t) :: e23
     !> gradients
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
     integer :: n_stats = 44
     !> Specifies a subset of the statistics to be collected. All 44 fields by
     !! default.
     character(5) :: stat_set
     !> A list of size n_stats, whith entries pointing to the fields that will
     !! be output (the field components above.) Used to write the output.
     type(field_list_t)  :: stat_fields
   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_stats_init
     !> Destructor. 
     procedure, pass(this) :: free => fluid_stats_free
     !> Update all the mean value fields with a new sample.
     procedure, pass(this) :: update => fluid_stats_update
     !> Reset all the computed means values and sampling times to zero.
     procedure, pass(this) :: reset => fluid_stats_reset
     ! Convert computed weak gradients to strong.
     procedure, pass(this) :: make_strong_grad => fluid_stats_make_strong_grad
     !> Compute certain physical statistical quantities based on existing mean 
     !! fields.
     procedure, pass(this) :: post_process => fluid_stats_post_process
  end type fluid_stats_t

contains

  !> Constructor. Initialize the fields associated with fluid_stats.
  !! @param coef SEM coefficients. Optional.
  !! @param u The x component of velocity.
  !! @param v The y component of velocity.
  !! @param w The z component of velocity.
  !! @param p The pressure. 
  !! @param set Specifies the subset of the statistics to be collected. 
  !! Optional. Either `basic` or `full`, defaults to `full`.
  subroutine fluid_stats_init(this, coef, u, v, w, p, set)
    class(fluid_stats_t), intent(inout), target:: this
    type(coef_t), target, optional :: coef
    type(field_t), target, intent(inout) :: u, v, w, p
    character(*), intent(in), optional :: set

    call this%free()
    this%coef => coef

    this%u => u
    this%v => v
    this%w => w
    this%p => p

    if (present(set)) then 
       this%stat_set  = trim(set)
       if (this%stat_set .eq. 'basic') then
          this%n_stats = 11
       end if
    else
       this%stat_set = 'full'
       this%n_stats = 44
    end if
     
    call this%stats_work%init(this%u%dof, 'stats')
    call this%stats_u%init(this%u%dof, 'u temp')
    call this%stats_v%init(this%u%dof, 'v temp')
    call this%stats_w%init(this%u%dof, 'w temp')
    call this%stats_p%init(this%u%dof, 'p temp')
    call this%u_mean%init(this%u)
    call this%v_mean%init(this%v)
    call this%w_mean%init(this%w)
    call this%p_mean%init(this%p)
    call this%uu%init(this%stats_u, 'uu')
    call this%vv%init(this%stats_v, 'vv')
    call this%ww%init(this%stats_w, 'ww')
    call this%uv%init(this%stats_work, 'uv')
    call this%uw%init(this%stats_work, 'uw')
    call this%vw%init(this%stats_work, 'vw')
    call this%pp%init(this%stats_p, 'pp')
 
    if (this%n_stats .eq. 44) then
       call this%dudx%init(this%u%dof, 'dudx')
       call this%dudy%init(this%u%dof, 'dudy')
       call this%dudz%init(this%u%dof, 'dudz')
       call this%dvdx%init(this%u%dof, 'dvdx')
       call this%dvdy%init(this%u%dof, 'dvdy')
       call this%dvdz%init(this%u%dof, 'dvdz')
       call this%dwdx%init(this%u%dof, 'dwdx')
       call this%dwdy%init(this%u%dof, 'dwdy')
       call this%dwdz%init(this%u%dof, 'dwdz')

       call this%uuu%init(this%stats_work, 'uuu')
       call this%vvv%init(this%stats_work, 'vvv')
       call this%www%init(this%stats_work, 'www')
       call this%uuv%init(this%stats_work, 'uuv')
       call this%uuw%init(this%stats_work, 'uuw')
       call this%uvv%init(this%stats_work, 'uvv')
       call this%uvw%init(this%stats_work, 'uvw')
       call this%vvw%init(this%stats_work, 'vvw')
       call this%uww%init(this%stats_work, 'uww')
       call this%vww%init(this%stats_work, 'vww')
       call this%uuuu%init(this%stats_work, 'uuuu')
       call this%vvvv%init(this%stats_work, 'vvvv')
       call this%wwww%init(this%stats_work, 'wwww')
       !> Pressure
       call this%ppp%init(this%stats_work, 'ppp')
       call this%pppp%init(this%stats_work, 'pppp')
       !> Pressure * velocity
       call this%pu%init(this%stats_work, 'pu')
       call this%pv%init(this%stats_work, 'pv')
       call this%pw%init(this%stats_work, 'pw')

       call this%pdudx%init(this%stats_work, 'pdudx')
       call this%pdudy%init(this%stats_work, 'pdudy')
       call this%pdudz%init(this%stats_work, 'pdudz')
       call this%pdvdx%init(this%stats_work, 'pdvdx')
       call this%pdvdy%init(this%stats_work, 'pdvdy')
       call this%pdvdz%init(this%stats_work, 'pdvdz')
       call this%pdwdx%init(this%stats_work, 'pdwdx')
       call this%pdwdy%init(this%stats_work, 'pdwdy')
       call this%pdwdz%init(this%stats_work, 'pdwdz')

       call this%e11%init(this%stats_work, 'e11')
       call this%e22%init(this%stats_work, 'e22')
       call this%e33%init(this%stats_work, 'e33')
       call this%e12%init(this%stats_work, 'e12')
       call this%e13%init(this%stats_work, 'e13')
       call this%e23%init(this%stats_work, 'e23')
    end if

    allocate(this%stat_fields%items(this%n_stats))

    call this%stat_fields%assign_to_field(1, this%p_mean%mf)
    call this%stat_fields%assign_to_field(2, this%u_mean%mf)
    call this%stat_fields%assign_to_field(3, this%v_mean%mf)
    call this%stat_fields%assign_to_field(4, this%w_mean%mf)
    call this%stat_fields%assign_to_field(5, this%pp%mf)
    call this%stat_fields%assign_to_field(6, this%uu%mf)
    call this%stat_fields%assign_to_field(7, this%vv%mf)
    call this%stat_fields%assign_to_field(8, this%ww%mf)
    call this%stat_fields%assign_to_field(9, this%uv%mf)
    call this%stat_fields%assign_to_field(10, this%uw%mf)
    call this%stat_fields%assign_to_field(11, this%vw%mf)

    if (this%n_stats .eq. 44) then
       call this%stat_fields%assign_to_field(12, this%uuu%mf)
       call this%stat_fields%assign_to_field(13, this%vvv%mf)
       call this%stat_fields%assign_to_field(14, this%www%mf)
       call this%stat_fields%assign_to_field(15, this%uuv%mf)
       call this%stat_fields%assign_to_field(16, this%uuw%mf)
       call this%stat_fields%assign_to_field(17, this%uvv%mf)
       call this%stat_fields%assign_to_field(18, this%uvw%mf)
       call this%stat_fields%assign_to_field(19, this%vvw%mf)
       call this%stat_fields%assign_to_field(20, this%uww%mf)
       call this%stat_fields%assign_to_field(21, this%vww%mf)
       call this%stat_fields%assign_to_field(22, this%uuuu%mf)
       call this%stat_fields%assign_to_field(23, this%vvvv%mf)
       call this%stat_fields%assign_to_field(24, this%wwww%mf)
       call this%stat_fields%assign_to_field(25, this%ppp%mf)
       call this%stat_fields%assign_to_field(26, this%pppp%mf)
       call this%stat_fields%assign_to_field(27, this%pu%mf)
       call this%stat_fields%assign_to_field(28, this%pv%mf)
       call this%stat_fields%assign_to_field(29, this%pw%mf)

       call this%stat_fields%assign_to_field(30, this%pdudx%mf)
       call this%stat_fields%assign_to_field(31, this%pdudy%mf)
       call this%stat_fields%assign_to_field(32, this%pdudz%mf)
       call this%stat_fields%assign_to_field(33, this%pdvdx%mf)
       call this%stat_fields%assign_to_field(34, this%pdvdy%mf)
       call this%stat_fields%assign_to_field(35, this%pdvdz%mf)
       call this%stat_fields%assign_to_field(36, this%pdwdx%mf)
       call this%stat_fields%assign_to_field(37, this%pdwdy%mf)
       call this%stat_fields%assign_to_field(38, this%pdwdz%mf)
       call this%stat_fields%assign_to_field(39, this%e11%mf)
       call this%stat_fields%assign_to_field(40, this%e22%mf)
       call this%stat_fields%assign_to_field(41, this%e33%mf)
       call this%stat_fields%assign_to_field(42, this%e12%mf)
       call this%stat_fields%assign_to_field(43, this%e13%mf)
       call this%stat_fields%assign_to_field(44, this%e23%mf)
    end if

  end subroutine fluid_stats_init

  !> Updates all fields with a new sample.
  !! @param k Time elapsed since the last update.
  subroutine fluid_stats_update(this, k)
    class(fluid_stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k
    integer :: n

    associate(stats_work => this%stats_work, stats_u => this%stats_u, &
              stats_v => this%stats_v, stats_w => this%stats_w, &
              stats_p => this%stats_p)
      n = stats_work%dof%size()

      !> U%f is u and U%mf is <u>
      if (NEKO_BCKND_DEVICE .eq. 1) then

         call this%u_mean%update(k)
         call this%v_mean%update(k)
         call this%w_mean%update(k)
         call this%p_mean%update(k)

         call device_col3(stats_u%x_d, this%u%x_d, this%u%x_d, n)
         call device_col3(stats_v%x_d, this%v%x_d, this%v%x_d, n)
         call device_col3(stats_w%x_d, this%w%x_d, this%w%x_d, n)
         call device_col3(stats_p%x_d, this%p%x_d, this%p%x_d, n)

         call this%uu%update(k)
         call this%vv%update(k)
         call this%ww%update(k)
         call this%pp%update(k)

         call device_col3(stats_work%x_d, this%u%x_d, this%v%x_d, n)
         call this%uv%update(k)
         call device_col3(stats_work%x_d, this%u%x_d, this%w%x_d, n)
         call this%uw%update(k)
         call device_col3(stats_work%x_d, this%v%x_d, this%w%x_d, n)
         call this%vw%update(k)
         if (this%n_stats .eq. 11) return
         call device_col2(stats_work%x_d, this%u%x_d, n)
         call this%uvw%update(k)
         call device_col3(stats_work%x_d, this%stats_u%x_d, this%u%x_d, n)
         call this%uuu%update(k)
         call device_col3(stats_work%x_d, this%stats_v%x_d, this%v%x_d, n)
         call this%vvv%update(k)
         call device_col3(stats_work%x_d, this%stats_w%x_d, this%w%x_d, n)
         call this%www%update(k)
         call device_col3(stats_work%x_d, this%stats_u%x_d, this%v%x_d, n)
         call this%uuv%update(k)
         call device_col3(stats_work%x_d, this%stats_u%x_d, this%w%x_d, n)
         call this%uuw%update(k)
         call device_col3(stats_work%x_d, this%stats_v%x_d, this%u%x_d, n)
         call this%uvv%update(k)
         call device_col3(stats_work%x_d, this%stats_v%x_d, this%w%x_d, n)
         call this%vvw%update(k)
         call device_col3(stats_work%x_d, this%stats_w%x_d, this%u%x_d, n)
         call this%uww%update(k)
         call device_col3(stats_work%x_d, this%stats_w%x_d, this%v%x_d, n)
         call this%vww%update(k)

         call device_col3(stats_work%x_d, this%stats_u%x_d, this%stats_u%x_d, n)
         call this%uuuu%update(k)
         call device_col3(stats_work%x_d, this%stats_v%x_d, this%stats_v%x_d, n)
         call this%vvvv%update(k)
         call device_col3(stats_work%x_d, this%stats_w%x_d, this%stats_w%x_d, n)
         call this%wwww%update(k)

         call device_col3(stats_work%x_d, this%stats_p%x_d, this%p%x_d, n)
         call this%ppp%update(k)
         call device_col3(stats_work%x_d, this%stats_p%x_d, this%stats_p%x_d, n)
         call this%pppp%update(k)

         call device_col3(stats_work%x_d, this%p%x_d, this%u%x_d, n)
         call this%pu%update(k)
         call device_col3(stats_work%x_d, this%p%x_d, this%v%x_d, n)
         call this%pv%update(k)
         call device_col3(stats_work%x_d, this%p%x_d, this%w%x_d, n)
         call this%pw%update(k)

      else

         call this%u_mean%update(k)
         call this%v_mean%update(k)
         call this%w_mean%update(k)
         call this%p_mean%update(k)
         call col3(stats_u%x, this%u%x, this%u%x, n)
         call col3(stats_v%x, this%v%x, this%v%x, n)
         call col3(stats_w%x, this%w%x, this%w%x, n)
         call col3(stats_p%x, this%p%x, this%p%x, n)

         call this%uu%update(k)
         call this%vv%update(k)
         call this%ww%update(k)
         call this%pp%update(k)

         call col3(stats_work%x, this%u%x, this%v%x, n)
         call this%uv%update(k)
         call col3(stats_work%x, this%u%x, this%w%x, n)
         call this%uw%update(k)
         call col3(stats_work%x, this%v%x, this%w%x, n)
         call this%vw%update(k)

         if (this%n_stats .eq. 11) return

         call col2(stats_work%x, this%u%x, n)
         call this%uvw%update(k)
         call col3(stats_work%x, this%stats_u%x, this%u%x, n)
         call this%uuu%update(k)
         call col3(stats_work%x, this%stats_v%x, this%v%x, n)
         call this%vvv%update(k)
         call col3(stats_work%x, this%stats_w%x, this%w%x, n)
         call this%www%update(k)
         call col3(stats_work%x, this%stats_u%x, this%v%x, n)
         call this%uuv%update(k)
         call col3(stats_work%x, this%stats_u%x, this%w%x, n)
         call this%uuw%update(k)
         call col3(stats_work%x, this%stats_v%x, this%u%x, n)
         call this%uvv%update(k)
         call col3(stats_work%x, this%stats_v%x, this%w%x, n)
         call this%vvw%update(k)
         call col3(stats_work%x, this%stats_w%x, this%u%x, n)
         call this%uww%update(k)
         call col3(stats_work%x, this%stats_w%x, this%v%x, n)
         call this%vww%update(k)

         call col3(stats_work%x, this%stats_u%x, this%stats_u%x, n)
         call this%uuuu%update(k)
         call col3(stats_work%x, this%stats_v%x, this%stats_v%x, n)
         call this%vvvv%update(k)
         call col3(stats_work%x, this%stats_w%x, this%stats_w%x, n)
         call this%wwww%update(k)

         call col3(stats_work%x, this%stats_p%x, this%p%x, n)
         call this%ppp%update(k)
         call col3(stats_work%x, this%stats_p%x, this%stats_p%x, n)
         call this%pppp%update(k)

         call col3(stats_work%x, this%p%x, this%u%x,n)
         call this%pu%update(k)
         call col3(stats_work%x, this%p%x, this%v%x,n)
         call this%pv%update(k)
         call col3(stats_work%x, this%p%x, this%w%x,n)
         call this%pw%update(k)


      end if
      call opgrad(this%dudx%x, this%dudy%x, this%dudz%x, this%u%x, this%coef)
      call opgrad(this%dvdx%x, this%dvdy%x, this%dvdz%x, this%v%x, this%coef)
      call opgrad(this%dwdx%x, this%dwdy%x, this%dwdz%x, this%w%x, this%coef)

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col3(stats_work%x_d, this%dudx%x_d, this%p%x_d, n)
         call this%pdudx%update(k)
         call device_col3(stats_work%x_d, this%dudy%x_d, this%p%x_d, n)
         call this%pdudy%update(k)
         call device_col3(stats_work%x_d, this%dudz%x_d, this%p%x_d, n)
         call this%pdudz%update(k)

         call device_col3(stats_work%x_d, this%dvdx%x_d, this%p%x_d, n)
         call this%pdvdx%update(k)
         call device_col3(stats_work%x_d, this%dvdy%x_d, this%p%x_d, n)
         call this%pdvdy%update(k)
         call device_col3(stats_work%x_d, this%dvdz%x_d, this%p%x_d, n)
         call this%pdvdz%update(k)

         call device_col3(stats_work%x_d, this%dwdx%x_d, this%p%x_d, n)
         call this%pdwdx%update(k)
         call device_col3(stats_work%x_d, this%dwdy%x_d, this%p%x_d, n)
         call this%pdwdy%update(k)
         call device_col3(stats_work%x_d, this%dwdz%x_d, this%p%x_d, n)
         call this%pdwdz%update(k)

         call device_col3(this%stats_work%x_d, this%dudx%x_d, this%dudx%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dudy%x_d, &
              this%dudy%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dudz%x_d, &
              this%dudz%x_d, n)
         call this%e11%update(k)
         call device_col3(this%stats_work%x_d, this%dvdx%x_d, this%dvdx%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dvdy%x_d, & 
              this%dvdy%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dvdz%x_d, & 
              this%dvdz%x_d, n)
         call this%e22%update(k)
         call device_col3(this%stats_work%x_d, this%dwdx%x_d, this%dwdx%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dwdy%x_d, &
              this%dwdy%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dwdz%x_d, &
              this%dwdz%x_d, n)
         call this%e33%update(k)
         call device_col3(this%stats_work%x_d, this%dudx%x_d, &
              this%dvdx%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dudy%x_d, &
              this%dvdy%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dudz%x_d, &
              this%dvdz%x_d, n)
         call this%e12%update(k)
         call device_col3(this%stats_work%x_d, this%dudx%x_d, this%dwdx%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dudy%x_d, &
              this%dwdy%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dudz%x_d, &
              this%dwdz%x_d, n)
         call this%e13%update(k)
         call device_col3(this%stats_work%x_d, this%dvdx%x_d, this%dwdx%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dvdy%x_d, &
              this%dwdy%x_d, n)
         call device_addcol3(this%stats_work%x_d, this%dvdz%x_d, &
              this%dwdz%x_d, n)
         call this%e23%update(k)
      else
         call col3(stats_work%x, this%dudx%x, this%p%x, n)
         call this%pdudx%update(k)
         call col3(stats_work%x, this%dudy%x, this%p%x, n)
         call this%pdudy%update(k)
         call col3(stats_work%x, this%dudz%x, this%p%x, n)
         call this%pdudz%update(k)

         call col3(stats_work%x, this%dvdx%x, this%p%x, n)
         call this%pdvdx%update(k)
         call col3(stats_work%x, this%dvdy%x, this%p%x, n)
         call this%pdvdy%update(k)
         call col3(stats_work%x, this%dvdz%x, this%p%x, n)
         call this%pdvdz%update(k)

         call col3(stats_work%x, this%dwdx%x, this%p%x, n)
         call this%pdwdx%update(k)
         call col3(stats_work%x, this%dwdy%x, this%p%x, n)
         call this%pdwdy%update(k)
         call col3(stats_work%x, this%dwdz%x, this%p%x, n)
         call this%pdwdz%update(k)

         call col3(this%stats_work%x, this%dudx%x, this%dudx%x, n)
         call addcol3(this%stats_work%x, this%dudy%x, this%dudy%x, n)
         call addcol3(this%stats_work%x, this%dudz%x, this%dudz%x, n)
         call this%e11%update(k)
         call col3(this%stats_work%x, this%dvdx%x, this%dvdx%x, n)
         call addcol3(this%stats_work%x, this%dvdy%x, this%dvdy%x, n)
         call addcol3(this%stats_work%x, this%dvdz%x, this%dvdz%x, n)
         call this%e22%update(k)
         call col3(this%stats_work%x, this%dwdx%x, this%dwdx%x, n)
         call addcol3(this%stats_work%x, this%dwdy%x, this%dwdy%x, n)
         call addcol3(this%stats_work%x, this%dwdz%x, this%dwdz%x, n)
         call this%e33%update(k)
         call col3(this%stats_work%x, this%dudx%x, this%dvdx%x, n)
         call addcol3(this%stats_work%x, this%dudy%x, this%dvdy%x, n)
         call addcol3(this%stats_work%x, this%dudz%x, this%dvdz%x, n)
         call this%e12%update(k)
         call col3(this%stats_work%x,this%dudx%x, this%dwdx%x,n)
         call addcol3(this%stats_work%x,this%dudy%x, this%dwdy%x,n)
         call addcol3(this%stats_work%x,this%dudz%x, this%dwdz%x,n)
         call this%e13%update(k)
         call col3(this%stats_work%x, this%dvdx%x, this%dwdx%x, n)
         call addcol3(this%stats_work%x, this%dvdy%x, this%dwdy%x, n)
         call addcol3(this%stats_work%x, this%dvdz%x, this%dwdz%x, n)
         call this%e23%update(k)

      end if
    end associate

  end subroutine fluid_stats_update


  !> Destructor.
  subroutine fluid_stats_free(this)
    class(fluid_stats_t), intent(inout) :: this

    call this%stats_work%free()
    call this%stats_u%free()
    call this%stats_v%free()
    call this%stats_w%free()

    call this%u_mean%free()
    call this%v_mean%free()
    call this%w_mean%free()
    call this%p_mean%free()

    call this%uu%free()
    call this%vv%free()
    call this%ww%free()
    call this%uv%free()
    call this%uw%free()
    call this%vw%free()
    call this%pp%free()

    call this%dUdx%free()
    call this%dUdy%free()
    call this%dUdz%free()
    call this%dVdx%free()
    call this%dVdy%free()
    call this%dVdz%free()
    call this%dWdx%free()
    call this%dWdy%free()
    call this%dWdz%free()

  end subroutine fluid_stats_free

  !> Resets all the computed means values and sampling times to zero.
  subroutine fluid_stats_reset(this)
    class(fluid_stats_t), intent(inout), target:: this

    call this%p_mean%reset()
    call this%u_mean%reset()
    call this%v_mean%reset()
    call this%w_mean%reset()

    call this%uu%reset()
    call this%vv%reset()
    call this%ww%reset()
    call this%uv%reset()
    call this%uw%reset()
    call this%vw%reset()
    call this%pp%reset()
    if (this%n_stats .eq. 44) then
       call this%uuu%reset()
       call this%vvv%reset()
       call this%www%reset()
       call this%uuv%reset()
       call this%uuw%reset()
       call this%uvv%reset()
       call this%uvw%reset()
       call this%vvw%reset()
       call this%uww%reset()
       call this%vww%reset()
       call this%uuuu%reset()
       call this%vvvv%reset()
       call this%wwww%reset()
       call this%ppp%reset()
       call this%pppp%reset()
       call this%pu%reset()
       call this%pv%reset()
       call this%pw%reset()

       call this%pdudx%reset()
       call this%pdudy%reset()
       call this%pdudz%reset()
       call this%pdvdx%reset()
       call this%pdvdy%reset()
       call this%pdvdz%reset()
       call this%pdwdx%reset()
       call this%pdwdy%reset()
       call this%pdwdz%reset()

       call this%e11%reset()
       call this%e22%reset()
       call this%e33%reset()
       call this%e12%reset()
       call this%e13%reset()
       call this%e23%reset()
    end if

  end subroutine fluid_stats_reset

  ! Convert computed weak gradients to strong.
  subroutine fluid_stats_make_strong_grad(this)
    class(fluid_stats_t) :: this
    integer :: n
 
    if (this%n_stats .eq. 11) return
 
    n = size(this%coef%B)
 
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(this%stats_work%x_d, 1.0_rp, n)
       call device_invcol2(this%stats_work%x_d, this%coef%B_d, n)
       call device_col2(this%pdudx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdudy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdudz%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdvdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdvdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdvdz%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdwdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdwdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdwdz%mf%x_d, this%stats_work%x_d, n)

       call device_col2(this%stats_work%x_d, this%stats_work%x_d, n)
       call device_col2(this%e11%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%e22%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%e33%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%e12%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%e13%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%e23%mf%x_d, this%stats_work%x_d, n)


    else
       call invers2(this%stats_work%x, this%coef%B, n)
       call col2(this%pdudx%mf%x, this%stats_work%x, n)
       call col2(this%pdudy%mf%x, this%stats_work%x, n)
       call col2(this%pdudz%mf%x, this%stats_work%x, n)
       call col2(this%pdvdx%mf%x, this%stats_work%x, n)
       call col2(this%pdvdy%mf%x, this%stats_work%x, n)
       call col2(this%pdvdz%mf%x, this%stats_work%x, n)
       call col2(this%pdwdx%mf%x, this%stats_work%x, n)
       call col2(this%pdwdy%mf%x, this%stats_work%x, n)
       call col2(this%pdwdz%mf%x, this%stats_work%x, n)

       call col2(this%stats_work%x, this%stats_work%x, n)
       call col2(this%e11%mf%x, this%stats_work%x, n)
       call col2(this%e22%mf%x, this%stats_work%x, n)
       call col2(this%e33%mf%x, this%stats_work%x, n)
       call col2(this%e12%mf%x, this%stats_work%x, n)
       call col2(this%e13%mf%x, this%stats_work%x, n)
       call col2(this%e23%mf%x, this%stats_work%x, n)

    end if

  end subroutine fluid_stats_make_strong_grad

  !> Compute certain physical statistical quantities based on existing mean 
  !! fields.
  subroutine fluid_stats_post_process(this, mean, reynolds, pressure_flatness,&
       pressure_skewness, skewness_tensor, mean_vel_grad, dissipation_tensor)
    class(fluid_stats_t) :: this
    type(field_list_t), intent(inout), optional :: mean
    type(field_list_t), intent(inout), optional :: reynolds
    type(field_list_t), intent(inout), optional :: pressure_skewness
    type(field_list_t), intent(inout), optional :: pressure_flatness
    type(field_list_t), intent(inout), optional :: skewness_tensor
    type(field_list_t), intent(inout), optional :: mean_vel_grad
    type(field_list_t), intent(inout), optional :: dissipation_tensor
    integer :: n

    if (present(mean)) then
       n = mean%item_size(1)
       call copy(mean%items(1)%ptr%x, this%u_mean%mf%x, n)
       call copy(mean%items(2)%ptr%x, this%v_mean%mf%x, n)
       call copy(mean%items(3)%ptr%x, this%w_mean%mf%x, n)
       call copy(mean%items(4)%ptr%x, this%p_mean%mf%x, n)
    end if

    if (present(reynolds)) then
       n = reynolds%item_size(1)
       call copy(reynolds%items(1)%ptr%x, this%pp%mf%x, n)
       call subcol3(reynolds%items(1)%ptr%x, this%p_mean%mf%x, &
            this%p_mean%mf%x, n)

       call copy(reynolds%items(2)%ptr%x, this%uu%mf%x, n)
       call subcol3(reynolds%items(2)%ptr%x, this%u_mean%mf%x, &
            this%u_mean%mf%x, n)

       call copy(reynolds%items(3)%ptr%x, this%vv%mf%x, n)
       call subcol3(reynolds%items(3)%ptr%x, this%v_mean%mf%x, &
            this%v_mean%mf%x,n)

       call copy(reynolds%items(4)%ptr%x, this%ww%mf%x, n)
       call subcol3(reynolds%items(4)%ptr%x, this%w_mean%mf%x, &
            this%w_mean%mf%x,n)

       call copy(reynolds%items(5)%ptr%x, this%uv%mf%x, n)
       call subcol3(reynolds%items(5)%ptr%x, this%u_mean%mf%x, &
            this%v_mean%mf%x, n)

       call copy(reynolds%items(6)%ptr%x, this%uw%mf%x, n)
       call subcol3(reynolds%items(6)%ptr%x, this%u_mean%mf%x, &
            this%w_mean%mf%x, n)

       call copy(reynolds%items(7)%ptr%x, this%vw%mf%x, n)
       call subcol3(reynolds%items(7)%ptr%x, this%v_mean%mf%x, &
            this%w_mean%mf%x, n)
    end if
    if (present(pressure_skewness)) then

       call neko_warning('Presssure skewness stat not implemented'// &
                         ' in fluid_stats, process stats in python instead')

    end if

    if (present(pressure_flatness)) then
       call neko_warning('Presssure flatness stat not implemented'// &
                         ' in fluid_stats, process stats in python instead')

    end if

    if (present(skewness_tensor)) then
       call neko_warning('Skewness tensor stat not implemented'// &
                         ' in fluid_stats, process stats in python instead')
    end if

    if (present(mean_vel_grad)) then
       !Compute gradient of mean flow
       n = mean_vel_grad%item_size(1)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%u_mean%mf%x, this%u_mean%mf%x_d, n, &
                             HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%v_mean%mf%x, this%v_mean%mf%x_d, n, &
                             HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%w_mean%mf%x, this%w_mean%mf%x_d, n, &
                             HOST_TO_DEVICE, sync = .false.)
          call opgrad(this%dudx%x, this%dudy%x, this%dudz%x, &
                      this%u_mean%mf%x, this%coef)
          call opgrad(this%dvdx%x, this%dvdy%x, this%dvdz%x, &
                      this%v_mean%mf%x, this%coef)
          call opgrad(this%dwdx%x, this%dwdy%x, this%dwdz%x, &
                      this%w_mean%mf%x, this%coef)
          call device_memcpy(this%dudx%x, this%dudx%x_d, n, &
                             DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(this%dvdx%x, this%dvdx%x_d, n, &
                             DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(this%dwdx%x, this%dwdx%x_d, n, &
                             DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(this%dudy%x, this%dudy%x_d, n, &
                             DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(this%dvdy%x, this%dvdy%x_d, n, &
                             DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(this%dwdy%x, this%dwdy%x_d, n, &
                             DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(this%dudz%x, this%dudz%x_d, n, &
                             DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(this%dvdz%x, this%dvdz%x_d, n, &
                             DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(this%dwdz%x, this%dwdz%x_d, n, &
                             DEVICE_TO_HOST, sync = .true.)
       else
          call opgrad(this%dudx%x, this%dudy%x, this%dudz%x, &
                      this%u_mean%mf%x, this%coef)
          call opgrad(this%dvdx%x, this%dvdy%x, this%dvdz%x, &
                      this%v_mean%mf%x, this%coef)
          call opgrad(this%dwdx%x, this%dwdy%x, this%dwdz%x, & 
                      this%w_mean%mf%x, this%coef)
       end if
       call invers2(this%stats_work%x, this%coef%B,n)
       call col3(mean_vel_grad%items(1)%ptr%x, this%dudx%x, & 
                 this%stats_work%x, n)
       call col3(mean_vel_grad%items(2)%ptr%x, this%dudy%x, &
                 this%stats_work%x, n)
       call col3(mean_vel_grad%items(3)%ptr%x, this%dudz%x, &
                 this%stats_work%x, n)
       call col3(mean_vel_grad%items(4)%ptr%x, this%dvdx%x, &
                 this%stats_work%x, n)
       call col3(mean_vel_grad%items(5)%ptr%x, this%dvdy%x, &
                 this%stats_work%x, n)
       call col3(mean_vel_grad%items(6)%ptr%x, this%dvdz%x, &
                 this%stats_work%x, n)
       call col3(mean_vel_grad%items(7)%ptr%x, this%dwdx%x, &
                 this%stats_work%x, n)
       call col3(mean_vel_grad%items(8)%ptr%x, this%dwdy%x, &
                 this%stats_work%x, n)
       call col3(mean_vel_grad%items(9)%ptr%x, this%dwdz%x, &
                 this%stats_work%x, n)

    end if

    if (present(dissipation_tensor)) then
       call neko_warning('Dissipation tensor stat not implemented'// &
                         ' in fluid_stats, process stats in python instead')
    end if

  end subroutine fluid_stats_post_process

end module fluid_stats
