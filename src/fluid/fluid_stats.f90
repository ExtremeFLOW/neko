! Copyright (c) 2022-2026, The Neko Authors
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
  use device_math, only : device_col2, device_cfill, device_invcol2
  use num_types, only : rp
  use math, only : copy, subcol3
  use operators, only : opgrad
  use coefs, only : coef_t
  use field, only : field_t
  use field_list, only : field_list_t
  use field_math, only : field_col3
  use stats_quant, only : stats_quant_t
  use scratch_registry, only : neko_scratch_registry
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_warning
  implicit none
  private

  type, public, extends(stats_quant_t) :: fluid_stats_t
     !> Pointers to the instantaneous quantities.
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
     !> SEM coefficients.
     type(coef_t), pointer :: coef
     !> Number of statistical fields to be computed.
     integer :: n_stats = 44
     !> Specifies a subset of the statistics to be collected. All 44 fields by
     !! default.
     character(5) :: stat_set
     !> A list of size n_stats, with entries pointing to the fields that will
     !! be output (the field components above.) Used to write the output.
     type(field_list_t) :: stat_fields
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
  subroutine fluid_stats_init(this, coef, u, v, w, p, set, name)
    class(fluid_stats_t), intent(inout), target:: this
    type(coef_t), target, optional :: coef
    type(field_t), target, intent(in) :: u, v, w, p
    character(*), intent(in), optional :: set
    character(*), intent(in), optional :: name

    character(len=1024) :: unique_name
    unique_name = ""

    call this%free()
    this%coef => coef

    this%u => u
    this%v => v
    this%w => w
    this%p => p

    if (present(set)) then
       this%stat_set = trim(set)
       if (this%stat_set .eq. 'basic') then
          this%n_stats = 11
       end if
    else
       this%stat_set = 'full'
       this%n_stats = 44
    end if

    if (present(name)) then
       unique_name = name // "/"
    else
       unique_name = "fluid_stats/"
    end if

    call this%u_mean%init(this%u, trim(unique_name) // 'mean_u')
    call this%v_mean%init(this%v, trim(unique_name) // 'mean_v')
    call this%w_mean%init(this%w, trim(unique_name) // 'mean_w')
    call this%p_mean%init(this%p, trim(unique_name) // 'mean_p')
    ! The remaining means are of products of the instantaneous fields, and are
    ! updated with update_col3/update_dot3. The field passed to init only
    ! provides the dofmap.
    call this%uu%init(this%u, trim(unique_name) // 'mean_uu')
    call this%vv%init(this%u, trim(unique_name) // 'mean_vv')
    call this%ww%init(this%u, trim(unique_name) // 'mean_ww')
    call this%uv%init(this%u, trim(unique_name) // 'mean_uv')
    call this%uw%init(this%u, trim(unique_name) // 'mean_uw')
    call this%vw%init(this%u, trim(unique_name) // 'mean_vw')
    call this%pp%init(this%u, trim(unique_name) // 'mean_pp')

    if (this%n_stats .eq. 44) then
       call this%uuu%init(this%u, trim(unique_name) // 'mean_uuu')
       call this%vvv%init(this%u, trim(unique_name) // 'mean_vvv')
       call this%www%init(this%u, trim(unique_name) // 'mean_www')
       call this%uuv%init(this%u, trim(unique_name) // 'mean_uuv')
       call this%uuw%init(this%u, trim(unique_name) // 'mean_uuw')
       call this%uvv%init(this%u, trim(unique_name) // 'mean_uvv')
       call this%uvw%init(this%u, trim(unique_name) // 'mean_uvw')
       call this%vvw%init(this%u, trim(unique_name) // 'mean_vvw')
       call this%uww%init(this%u, trim(unique_name) // 'mean_uww')
       call this%vww%init(this%u, trim(unique_name) // 'mean_vww')
       call this%uuuu%init(this%u, trim(unique_name) // 'mean_uuuu')
       call this%vvvv%init(this%u, trim(unique_name) // 'mean_vvvv')
       call this%wwww%init(this%u, trim(unique_name) // 'mean_wwww')
       !> Pressure
       call this%ppp%init(this%u, trim(unique_name) // 'mean_ppp')
       call this%pppp%init(this%u, trim(unique_name) // 'mean_pppp')
       !> Pressure * velocity
       call this%pu%init(this%u, trim(unique_name) // 'mean_pu')
       call this%pv%init(this%u, trim(unique_name) // 'mean_pv')
       call this%pw%init(this%u, trim(unique_name) // 'mean_pw')

       call this%pdudx%init(this%u, trim(unique_name) // 'mean_pdudx')
       call this%pdudy%init(this%u, trim(unique_name) // 'mean_pdudy')
       call this%pdudz%init(this%u, trim(unique_name) // 'mean_pdudz')
       call this%pdvdx%init(this%u, trim(unique_name) // 'mean_pdvdx')
       call this%pdvdy%init(this%u, trim(unique_name) // 'mean_pdvdy')
       call this%pdvdz%init(this%u, trim(unique_name) // 'mean_pdvdz')
       call this%pdwdx%init(this%u, trim(unique_name) // 'mean_pdwdx')
       call this%pdwdy%init(this%u, trim(unique_name) // 'mean_pdwdy')
       call this%pdwdz%init(this%u, trim(unique_name) // 'mean_pdwdz')

       call this%e11%init(this%u, trim(unique_name) // 'mean_e11')
       call this%e22%init(this%u, trim(unique_name) // 'mean_e22')
       call this%e33%init(this%u, trim(unique_name) // 'mean_e33')
       call this%e12%init(this%u, trim(unique_name) // 'mean_e12')
       call this%e13%init(this%u, trim(unique_name) // 'mean_e13')
       call this%e23%init(this%u, trim(unique_name) // 'mean_e23')
    end if

    call this%stat_fields%init(this%n_stats)

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
    type(field_t), pointer :: work
    type(field_t), pointer :: dudx, dudy, dudz
    type(field_t), pointer :: dvdx, dvdy, dvdz
    type(field_t), pointer :: dwdx, dwdy, dwdz
    integer :: idx(10)

    associate(u => this%u, v => this%v, w => this%w, p => this%p)

      call this%u_mean%update(k)
      call this%v_mean%update(k)
      call this%w_mean%update(k)
      call this%p_mean%update(k)

      call this%uu%update_col3(k, u, u)
      call this%vv%update_col3(k, v, v)
      call this%ww%update_col3(k, w, w)
      call this%pp%update_col3(k, p, p)

      call this%uv%update_col3(k, u, v)
      call this%uw%update_col3(k, u, w)
      call this%vw%update_col3(k, v, w)

      if (this%n_stats .eq. 11) return

      call neko_scratch_registry%request_field(work, idx(1), &
           this%u%dof, .false.)

      call field_col3(work, u, u)
      call this%uuu%update_col3(k, work, u)
      call this%uuv%update_col3(k, work, v)
      call this%uuw%update_col3(k, work, w)
      call this%uuuu%update_col3(k, work, work)

      call field_col3(work, v, v)
      call this%vvv%update_col3(k, work, v)
      call this%uvv%update_col3(k, work, u)
      call this%vvw%update_col3(k, work, w)
      call this%vvvv%update_col3(k, work, work)

      call field_col3(work, w, w)
      call this%www%update_col3(k, work, w)
      call this%uww%update_col3(k, work, u)
      call this%vww%update_col3(k, work, v)
      call this%wwww%update_col3(k, work, work)

      call field_col3(work, v, w)
      call this%uvw%update_col3(k, work, u)

      call field_col3(work, p, p)
      call this%ppp%update_col3(k, work, p)
      call this%pppp%update_col3(k, work, work)

      call neko_scratch_registry%relinquish_field(idx(1))

      call this%pu%update_col3(k, p, u)
      call this%pv%update_col3(k, p, v)
      call this%pw%update_col3(k, p, w)

      call neko_scratch_registry%request_field(dudx, idx(2), &
           this%u%dof, .false.)
      call neko_scratch_registry%request_field(dudy, idx(3), &
           this%u%dof, .false.)
      call neko_scratch_registry%request_field(dudz, idx(4), &
           this%u%dof, .false.)
      call neko_scratch_registry%request_field(dvdx, idx(5), &
           this%u%dof, .false.)
      call neko_scratch_registry%request_field(dvdy, idx(6), &
           this%u%dof, .false.)
      call neko_scratch_registry%request_field(dvdz, idx(7), &
           this%u%dof, .false.)
      call neko_scratch_registry%request_field(dwdx, idx(8), &
           this%u%dof, .false.)
      call neko_scratch_registry%request_field(dwdy, idx(9), &
           this%u%dof, .false.)
      call neko_scratch_registry%request_field(dwdz, idx(10), &
           this%u%dof, .false.)

      call opgrad(dudx%x, dudy%x, dudz%x, u%x, this%coef)
      call opgrad(dvdx%x, dvdy%x, dvdz%x, v%x, this%coef)
      call opgrad(dwdx%x, dwdy%x, dwdz%x, w%x, this%coef)

      call this%pdudx%update_col3(k, dudx, p)
      call this%pdudy%update_col3(k, dudy, p)
      call this%pdudz%update_col3(k, dudz, p)
      call this%pdvdx%update_col3(k, dvdx, p)
      call this%pdvdy%update_col3(k, dvdy, p)
      call this%pdvdz%update_col3(k, dvdz, p)
      call this%pdwdx%update_col3(k, dwdx, p)
      call this%pdwdy%update_col3(k, dwdy, p)
      call this%pdwdz%update_col3(k, dwdz, p)

      call this%e11%update_dot3(k, dudx, dudy, dudz, dudx, dudy, dudz)
      call this%e22%update_dot3(k, dvdx, dvdy, dvdz, dvdx, dvdy, dvdz)
      call this%e33%update_dot3(k, dwdx, dwdy, dwdz, dwdx, dwdy, dwdz)
      call this%e12%update_dot3(k, dudx, dudy, dudz, dvdx, dvdy, dvdz)
      call this%e13%update_dot3(k, dudx, dudy, dudz, dwdx, dwdy, dwdz)
      call this%e23%update_dot3(k, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)

      call neko_scratch_registry%relinquish_field(idx(2:10))
    end associate

  end subroutine fluid_stats_update


  !> Destructor.
  subroutine fluid_stats_free(this)
    class(fluid_stats_t), intent(inout) :: this

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

    call this%uuu%free()
    call this%vvv%free()
    call this%www%free()
    call this%uuv%free()
    call this%uuw%free()
    call this%uvv%free()
    call this%uvw%free()
    call this%vvw%free()
    call this%uww%free()
    call this%vww%free()

    call this%uuuu%free()
    call this%vvvv%free()
    call this%wwww%free()

    call this%pp%free()
    call this%ppp%free()
    call this%pppp%free()

    call this%pu%free()
    call this%pv%free()
    call this%pw%free()

    call this%pdudx%free()
    call this%pdudy%free()
    call this%pdudz%free()
    call this%pdvdx%free()
    call this%pdvdy%free()
    call this%pdvdz%free()
    call this%pdwdx%free()
    call this%pdwdy%free()
    call this%pdwdz%free()

    call this%e11%free()
    call this%e22%free()
    call this%e33%free()
    call this%e12%free()
    call this%e13%free()
    call this%e23%free()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)
    nullify(this%coef)

    call this%stat_fields%free()

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
    integer :: n, i
    real(kind=rp) :: wrk, wrk_sqr
    type(field_t), pointer :: work
    integer :: idx

    if (this%n_stats .eq. 11) return

    n = size(this%coef%B)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_scratch_registry%request_field(work, idx, &
            this%u%dof, .false.)
       call device_cfill(work%x_d, 1.0_rp, n)
       call device_invcol2(work%x_d, this%coef%B_d, n)
       call device_col2(this%pdudx%mf%x_d, work%x_d, n)
       call device_col2(this%pdudy%mf%x_d, work%x_d, n)
       call device_col2(this%pdudz%mf%x_d, work%x_d, n)
       call device_col2(this%pdvdx%mf%x_d, work%x_d, n)
       call device_col2(this%pdvdy%mf%x_d, work%x_d, n)
       call device_col2(this%pdvdz%mf%x_d, work%x_d, n)
       call device_col2(this%pdwdx%mf%x_d, work%x_d, n)
       call device_col2(this%pdwdy%mf%x_d, work%x_d, n)
       call device_col2(this%pdwdz%mf%x_d, work%x_d, n)

       call device_col2(work%x_d, work%x_d, n)
       call device_col2(this%e11%mf%x_d, work%x_d, n)
       call device_col2(this%e22%mf%x_d, work%x_d, n)
       call device_col2(this%e33%mf%x_d, work%x_d, n)
       call device_col2(this%e12%mf%x_d, work%x_d, n)
       call device_col2(this%e13%mf%x_d, work%x_d, n)
       call device_col2(this%e23%mf%x_d, work%x_d, n)
       call neko_scratch_registry%relinquish_field(idx)

    else
       !$omp parallel do private(i, wrk, wrk_sqr)
       do i = 1, n
          wrk = 1.0_rp / this%coef%B(i,1,1,1)
          this%pdudx%mf%x(i,1,1,1) = this%pdudx%mf%x(i,1,1,1) * wrk
          this%pdudy%mf%x(i,1,1,1) = this%pdudy%mf%x(i,1,1,1) * wrk
          this%pdudz%mf%x(i,1,1,1) = this%pdudz%mf%x(i,1,1,1) * wrk

          this%pdvdx%mf%x(i,1,1,1) = this%pdvdx%mf%x(i,1,1,1) * wrk
          this%pdvdy%mf%x(i,1,1,1) = this%pdvdy%mf%x(i,1,1,1) * wrk
          this%pdvdz%mf%x(i,1,1,1) = this%pdvdz%mf%x(i,1,1,1) * wrk

          this%pdwdx%mf%x(i,1,1,1) = this%pdwdx%mf%x(i,1,1,1) * wrk
          this%pdwdy%mf%x(i,1,1,1) = this%pdwdy%mf%x(i,1,1,1) * wrk
          this%pdwdz%mf%x(i,1,1,1) = this%pdwdz%mf%x(i,1,1,1) * wrk

          wrk_sqr = wrk * wrk
          this%e11%mf%x(i,1,1,1) = this%e11%mf%x(i,1,1,1) * wrk_sqr
          this%e22%mf%x(i,1,1,1) = this%e22%mf%x(i,1,1,1) * wrk_sqr
          this%e33%mf%x(i,1,1,1) = this%e33%mf%x(i,1,1,1) * wrk_sqr

          this%e12%mf%x(i,1,1,1) = this%e12%mf%x(i,1,1,1) * wrk_sqr
          this%e13%mf%x(i,1,1,1) = this%e13%mf%x(i,1,1,1) * wrk_sqr
          this%e23%mf%x(i,1,1,1) = this%e23%mf%x(i,1,1,1) * wrk_sqr
       end do
       !$omp end parallel do
    end if

  end subroutine fluid_stats_make_strong_grad

  !> Compute certain physical statistical quantities based on existing mean
  !! fields.
  subroutine fluid_stats_post_process(this, mean, reynolds, pressure_flatness,&
       pressure_skewness, skewness_tensor, mean_vel_grad, dissipation_tensor)
    class(fluid_stats_t) :: this
    type(field_list_t), intent(inout), optional :: mean
    type(field_list_t), intent(inout), optional :: reynolds
    type(field_list_t), intent(in), optional :: pressure_skewness
    type(field_list_t), intent(in), optional :: pressure_flatness
    type(field_list_t), intent(in), optional :: skewness_tensor
    type(field_list_t), intent(inout), optional :: mean_vel_grad
    type(field_list_t), intent(in), optional :: dissipation_tensor
    integer :: n, i
    real(kind=rp) :: wrk
    type(field_t), pointer :: dudx, dudy, dudz
    type(field_t), pointer :: dvdx, dvdy, dvdz
    type(field_t), pointer :: dwdx, dwdy, dwdz
    integer :: idx(9)

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
       call neko_scratch_registry%request_field(dudx, idx(1), &
            this%u%dof, .false.)
       call neko_scratch_registry%request_field(dudy, idx(2), &
            this%u%dof, .false.)
       call neko_scratch_registry%request_field(dudz, idx(3), &
            this%u%dof, .false.)
       call neko_scratch_registry%request_field(dvdx, idx(4), &
            this%u%dof, .false.)
       call neko_scratch_registry%request_field(dvdy, idx(5), &
            this%u%dof, .false.)
       call neko_scratch_registry%request_field(dvdz, idx(6), &
            this%u%dof, .false.)
       call neko_scratch_registry%request_field(dwdx, idx(7), &
            this%u%dof, .false.)
       call neko_scratch_registry%request_field(dwdy, idx(8), &
            this%u%dof, .false.)
       call neko_scratch_registry%request_field(dwdz, idx(9), &
            this%u%dof, .false.)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%u_mean%mf%x, this%u_mean%mf%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%v_mean%mf%x, this%v_mean%mf%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%w_mean%mf%x, this%w_mean%mf%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
          call opgrad(dudx%x, dudy%x, dudz%x, &
               this%u_mean%mf%x, this%coef)
          call opgrad(dvdx%x, dvdy%x, dvdz%x, &
               this%v_mean%mf%x, this%coef)
          call opgrad(dwdx%x, dwdy%x, dwdz%x, &
               this%w_mean%mf%x, this%coef)
          call device_memcpy(dudx%x, dudx%x_d, n, &
               DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(dvdx%x, dvdx%x_d, n, &
               DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(dwdx%x, dwdx%x_d, n, &
               DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(dudy%x, dudy%x_d, n, &
               DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(dvdy%x, dvdy%x_d, n, &
               DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(dwdy%x, dwdy%x_d, n, &
               DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(dudz%x, dudz%x_d, n, &
               DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(dvdz%x, dvdz%x_d, n, &
               DEVICE_TO_HOST, sync = .false.)
          call device_memcpy(dwdz%x, dwdz%x_d, n, &
               DEVICE_TO_HOST, sync = .true.)
       else
          call opgrad(dudx%x, dudy%x, dudz%x, &
               this%u_mean%mf%x, this%coef)
          call opgrad(dvdx%x, dvdy%x, dvdz%x, &
               this%v_mean%mf%x, this%coef)
          call opgrad(dwdx%x, dwdy%x, dwdz%x, &
               this%w_mean%mf%x, this%coef)
       end if

       !$omp parallel do private(i, wrk)
       do i = 1, n
          wrk = 1.0_rp / this%coef%B(i,1,1,1)
          mean_vel_grad%items(1)%ptr%x(i,1,1,1) = dudx%x(i,1,1,1) * wrk
          mean_vel_grad%items(2)%ptr%x(i,1,1,1) = dudy%x(i,1,1,1) * wrk
          mean_vel_grad%items(3)%ptr%x(i,1,1,1) = dudz%x(i,1,1,1) * wrk

          mean_vel_grad%items(4)%ptr%x(i,1,1,1) = dvdx%x(i,1,1,1) * wrk
          mean_vel_grad%items(5)%ptr%x(i,1,1,1) = dvdy%x(i,1,1,1) * wrk
          mean_vel_grad%items(6)%ptr%x(i,1,1,1) = dvdz%x(i,1,1,1) * wrk

          mean_vel_grad%items(7)%ptr%x(i,1,1,1) = dwdx%x(i,1,1,1) * wrk
          mean_vel_grad%items(8)%ptr%x(i,1,1,1) = dwdy%x(i,1,1,1) * wrk
          mean_vel_grad%items(9)%ptr%x(i,1,1,1) = dwdz%x(i,1,1,1) * wrk
       end do
       !$omp end parallel do

       call neko_scratch_registry%relinquish_field(idx)
    end if

    if (present(dissipation_tensor)) then
       call neko_warning('Dissipation tensor stat not implemented'// &
            ' in fluid_stats, process stats in python instead')
    end if

  end subroutine fluid_stats_post_process

end module fluid_stats
