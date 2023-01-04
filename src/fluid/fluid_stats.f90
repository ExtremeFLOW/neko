! Copyright (c) 2022, The Neko Authors
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
!! 

module fluid_stats
  use mean_field
  use device_math
  use device_mathops
  use mathops
  use math
  use operators
  use coefs
  use field_registry
  use gather_scatter
  implicit none

  type, extends(stats_quant_t) :: fluid_stats_t
     type(field_t) :: stats_u !< Not reasonable to allocate 20 something fields
     type(field_t) :: stats_v !< Not reasonable to allocate 20 something fields
     type(field_t) :: stats_w !< Not reasonable to allocate 20 something fields
     type(field_t) :: stats_p !< Not reasonable to allocate 20 something fields
     type(field_t) :: stats_work !< Not reasonable to allocate 20 something fields
     type(field_t), pointer :: u !< u
     type(field_t), pointer :: v !< v
     type(field_t), pointer :: w !< w
     type(field_t), pointer :: p !< p

     type(field_t), pointer :: u_mean !< <u>
     type(field_t), pointer :: v_mean !< <v>
     type(field_t), pointer :: w_mean !< <w>
     type(field_t), pointer :: p_mean !< <p>
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



     type(coef_t), pointer :: coef
     integer :: n_stats = 40
     type(field_list_t)  :: stat_fields
   contains
     procedure, pass(this) :: init => fluid_stats_init
     procedure, pass(this) :: free => fluid_stats_free
     procedure, pass(this) :: update => fluid_stats_update
     procedure, pass(this) :: reset => fluid_stats_reset
     procedure, pass(this) :: make_strong_grad => fluid_stats_make_strong_grad
  end type fluid_stats_t

contains

  !> Initialize a mean flow field
  subroutine fluid_stats_init(this, coef)
    class(fluid_stats_t), intent(inout), target:: this
    type(coef_t), target :: coef
    this%coef => coef
  
    this%u_mean => neko_field_registry%get_field('mean_u')
    this%v_mean => neko_field_registry%get_field('mean_v')
    this%w_mean => neko_field_registry%get_field('mean_w')
    this%p_mean => neko_field_registry%get_field('mean_p')
    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')
    this%p => neko_field_registry%get_field('p')

    call this%stats_work%init(this%u%dof, 'stats')
    call this%stats_u%init(this%u%dof, 'u temp')
    call this%stats_v%init(this%u%dof, 'v temp')
    call this%stats_w%init(this%u%dof, 'w temp')
    call this%stats_p%init(this%u%dof, 'p temp')
 
    call this%dudx%init(this%u%dof, 'dudx')
    call this%dudy%init(this%u%dof, 'dudy')
    call this%dudz%init(this%u%dof, 'dudz')
    call this%dvdx%init(this%u%dof, 'dvdx')
    call this%dvdy%init(this%u%dof, 'dvdy')
    call this%dvdz%init(this%u%dof, 'dvdz')
    call this%dwdx%init(this%u%dof, 'dwdx')
    call this%dwdy%init(this%u%dof, 'dwdy')
    call this%dwdz%init(this%u%dof, 'dwdz')
    
    call this%uu%init(this%stats_u, 'uu')
    call this%vv%init(this%stats_v, 'vv')
    call this%ww%init(this%stats_w, 'ww')
    call this%uv%init(this%stats_work, 'uv')
    call this%uw%init(this%stats_work, 'uw')
    call this%vw%init(this%stats_work, 'vw')
    call this%uuu%init(this%stats_work, 'uuu')!< <uuu>
    call this%vvv%init(this%stats_work, 'vvv')!< <vvv>
    call this%www%init(this%stats_work, 'www')!< <www>
    call this%uuv%init(this%stats_work, 'uuv')!< <uuv>
    call this%uuw%init(this%stats_work, 'uuw')!< <uuw>
    call this%uvv%init(this%stats_work, 'uvv')!< <uvv>
    call this%uvw%init(this%stats_work, 'uvw')!< <uvv>
    call this%vvw%init(this%stats_work, 'vvw')!< <vvw>
    call this%uww%init(this%stats_work, 'uww')!< <uww>
    call this%vww%init(this%stats_work, 'vww')!< <vww>
    call this%uuuu%init(this%stats_work, 'uuuu') !< <uuuu>
    call this%vvvv%init(this%stats_work, 'vvvv') !< <vvvv>
    call this%wwww%init(this%stats_work, 'wwww') !< <wwww>
    !> Pressure
    call this%pp%init(this%stats_p, 'pp')
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

    allocate(this%stat_fields%fields(this%n_stats))

    this%stat_fields%fields(1)%f => this%pp%mf
    this%stat_fields%fields(2)%f => this%uu%mf
    this%stat_fields%fields(3)%f => this%vv%mf
    this%stat_fields%fields(4)%f => this%ww%mf
    this%stat_fields%fields(5)%f => this%uv%mf
    this%stat_fields%fields(6)%f => this%uw%mf
    this%stat_fields%fields(7)%f => this%vw%mf
    this%stat_fields%fields(8)%f => this%uuu%mf !< <uuu>
    this%stat_fields%fields(9)%f => this%vvv%mf !< <vvv>
    this%stat_fields%fields(10)%f => this%www%mf !< <www>
    this%stat_fields%fields(11)%f => this%uuv%mf !< <uuv>
    this%stat_fields%fields(12)%f => this%uuw%mf !< <uuw>
    this%stat_fields%fields(13)%f => this%uvv%mf !< <uvv>
    this%stat_fields%fields(14)%f => this%uvw%mf !< <uvv>
    this%stat_fields%fields(15)%f => this%vvw%mf !< <vvw>
    this%stat_fields%fields(16)%f => this%uww%mf !< <uww>
    this%stat_fields%fields(17)%f => this%vww%mf !< <vww>
    this%stat_fields%fields(18)%f => this%uuuu%mf !< <uuuu>
    this%stat_fields%fields(19)%f => this%vvvv%mf !< <vvvv>
    this%stat_fields%fields(20)%f => this%wwww%mf !< <wwww>
    this%stat_fields%fields(21)%f => this%ppp%mf
    this%stat_fields%fields(22)%f => this%pppp%mf
    this%stat_fields%fields(23)%f => this%pu%mf
    this%stat_fields%fields(24)%f => this%pv%mf
    this%stat_fields%fields(25)%f => this%pw%mf

    this%stat_fields%fields(26)%f => this%pdudx%mf
    this%stat_fields%fields(27)%f => this%pdudy%mf
    this%stat_fields%fields(28)%f => this%pdudz%mf
    this%stat_fields%fields(29)%f => this%pdvdx%mf
    this%stat_fields%fields(30)%f => this%pdvdy%mf
    this%stat_fields%fields(31)%f => this%pdvdz%mf
    this%stat_fields%fields(32)%f => this%pdwdx%mf
    this%stat_fields%fields(33)%f => this%pdwdy%mf
    this%stat_fields%fields(34)%f => this%pdwdz%mf
    this%stat_fields%fields(35)%f => this%e11%mf
    this%stat_fields%fields(36)%f => this%e22%mf
    this%stat_fields%fields(37)%f => this%e33%mf
    this%stat_fields%fields(38)%f => this%e12%mf
    this%stat_fields%fields(39)%f => this%e13%mf
    this%stat_fields%fields(40)%f => this%e23%mf


  end subroutine fluid_stats_init
  !> Updates all fields
  subroutine fluid_stats_update(this, k)
    class(fluid_stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k
    integer :: n



    associate(stats_work => this%stats_work, stats_u => this%stats_u,&
              stats_v => this%stats_v, stats_w => this%stats_w, stats_p => this%stats_p)
    n = stats_work%dof%size()

    !> U%f is u and U%mf is <u>
    if (NEKO_BCKND_DEVICE .eq. 1) then

       call device_col3(stats_u%x_d,this%u%x_d, this%u%x_d,n)
       call device_col3(stats_v%x_d,this%v%x_d, this%v%x_d,n)
       call device_col3(stats_w%x_d,this%w%x_d, this%w%x_d,n)
       call device_col3(stats_p%x_d,this%p%x_d, this%p%x_d,n)

       call this%uu%update(k)
       call this%vv%update(k)
       call this%ww%update(k)
       call this%pp%update(k)

       call device_col3(stats_work%x_d,this%u%x_d, this%v%x_d,n)
       call this%uv%update(k)
       call device_col3(stats_work%x_d,this%u%x_d, this%w%x_d,n)
       call this%uw%update(k)
       call device_col3(stats_work%x_d,this%v%x_d, this%w%x_d,n)
       call this%vw%update(k)

       call device_col2(stats_work%x_d, this%u%x_d,n)
       call this%uvw%update(k)
       call device_col3(stats_work%x_d,this%stats_u%x_d, this%u%x_d,n)
       call this%uuu%update(k)
       call device_col3(stats_work%x_d,this%stats_v%x_d, this%v%x_d,n)
       call this%vvv%update(k)
       call device_col3(stats_work%x_d,this%stats_w%x_d, this%w%x_d,n)
       call this%www%update(k)
       call device_col3(stats_work%x_d,this%stats_u%x_d, this%v%x_d,n)
       call this%uuv%update(k)
       call device_col3(stats_work%x_d,this%stats_u%x_d, this%w%x_d,n)
       call this%uuw%update(k)
       call device_col3(stats_work%x_d,this%stats_v%x_d, this%u%x_d,n)
       call this%uvv%update(k)
       call device_col3(stats_work%x_d,this%stats_v%x_d, this%w%x_d,n)
       call this%vvw%update(k)
       call device_col3(stats_work%x_d,this%stats_w%x_d, this%u%x_d,n)
       call this%uww%update(k)
       call device_col3(stats_work%x_d,this%stats_w%x_d, this%v%x_d,n)
       call this%vww%update(k)

       call device_col3(stats_work%x_d,this%stats_u%x_d, this%stats_u%x_d,n)
       call this%uuuu%update(k)
       call device_col3(stats_work%x_d,this%stats_v%x_d, this%stats_v%x_d,n)
       call this%vvvv%update(k)
       call device_col3(stats_work%x_d,this%stats_w%x_d, this%stats_w%x_d,n)
       call this%wwww%update(k)

       call device_col3(stats_work%x_d,this%stats_p%x_d, this%p%x_d,n)
       call this%ppp%update(k)
       call device_col3(stats_work%x_d,this%stats_p%x_d, this%stats_p%x_d,n)
       call this%pppp%update(k)

       call device_col3(stats_work%x_d,this%p%x_d, this%u%x_d,n)
       call this%pu%update(k)
       call device_col3(stats_work%x_d,this%p%x_d, this%v%x_d,n)
       call this%pv%update(k)
       call device_col3(stats_work%x_d,this%p%x_d, this%w%x_d,n)
       call this%pw%update(k)

    else

       call col3(stats_u%x,this%u%x, this%u%x,n)
       call col3(stats_v%x,this%v%x, this%v%x,n)
       call col3(stats_w%x,this%w%x, this%w%x,n)
       call col3(stats_p%x,this%p%x, this%p%x,n)

       call this%uu%update(k)
       call this%vv%update(k)
       call this%ww%update(k)
       call this%pp%update(k)

       call col3(stats_work%x,this%u%x, this%v%x,n)
       call this%uv%update(k)
       call col3(stats_work%x,this%u%x, this%w%x,n)
       call this%uw%update(k)
       call col3(stats_work%x,this%v%x, this%w%x,n)
       call this%vw%update(k)

       call col2(stats_work%x, this%u%x,n)
       call this%uvw%update(k)
       call col3(stats_work%x,this%stats_u%x, this%u%x,n)
       call this%uuu%update(k)
       call col3(stats_work%x,this%stats_v%x, this%v%x,n)
       call this%vvv%update(k)
       call col3(stats_work%x,this%stats_w%x, this%w%x,n)
       call this%www%update(k)
       call col3(stats_work%x,this%stats_u%x, this%v%x,n)
       call this%uuv%update(k)
       call col3(stats_work%x,this%stats_u%x, this%w%x,n)
       call this%uuw%update(k)
       call col3(stats_work%x,this%stats_v%x, this%u%x,n)
       call this%uvv%update(k)
       call col3(stats_work%x,this%stats_v%x, this%w%x,n)
       call this%vvw%update(k)
       call col3(stats_work%x,this%stats_w%x, this%u%x,n)
       call this%uww%update(k)
       call col3(stats_work%x,this%stats_w%x, this%v%x,n)
       call this%vww%update(k)

       call col3(stats_work%x,this%stats_u%x, this%stats_u%x,n)
       call this%uuuu%update(k)
       call col3(stats_work%x,this%stats_v%x, this%stats_v%x,n)
       call this%vvvv%update(k)
       call col3(stats_work%x,this%stats_w%x, this%stats_w%x,n)
       call this%wwww%update(k)

       call col3(stats_work%x,this%stats_p%x, this%p%x,n)
       call this%ppp%update(k)
       call col3(stats_work%x,this%stats_p%x, this%stats_p%x,n)
       call this%pppp%update(k)

       call col3(stats_work%x,this%p%x, this%u%x,n)
       call this%pu%update(k)
       call col3(stats_work%x,this%p%x, this%v%x,n)
       call this%pv%update(k)
       call col3(stats_work%x,this%p%x, this%w%x,n)
       call this%pw%update(k)


    end if

    call opgrad(this%dudx%x,this%dudy%x, this%dudz%x,this%u%x,this%coef)
    call opgrad(this%dvdx%x,this%dvdy%x, this%dvdz%x,this%v%x,this%coef)
    call opgrad(this%dwdx%x,this%dwdy%x, this%dwdz%x,this%w%x,this%coef)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3(stats_work%x_d,this%dudx%x_d, this%p%x_d,n)
       call this%pdudx%update(k)
       call device_col3(stats_work%x_d,this%dudy%x_d, this%p%x_d,n)
       call this%pdudy%update(k)
       call device_col3(stats_work%x_d,this%dudz%x_d, this%p%x_d,n)
       call this%pdudz%update(k)

       call device_col3(stats_work%x_d,this%dvdx%x_d, this%p%x_d,n)
       call this%pdvdx%update(k)
       call device_col3(stats_work%x_d,this%dvdy%x_d, this%p%x_d,n)
       call this%pdvdy%update(k)
       call device_col3(stats_work%x_d,this%dvdz%x_d, this%p%x_d,n)
       call this%pdvdz%update(k)

       call device_col3(stats_work%x_d,this%dwdx%x_d, this%p%x_d,n)
       call this%pdwdx%update(k)
       call device_col3(stats_work%x_d,this%dwdy%x_d, this%p%x_d,n)
       call this%pdwdy%update(k)
       call device_col3(stats_work%x_d,this%dwdz%x_d, this%p%x_d,n)
       call this%pdwdz%update(k)

       call device_col3(this%stats_work%x_d,this%dudx%x_d, this%dudx%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dudy%x_d, this%dudy%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dudz%x_d, this%dudz%x_d,n)
       call this%e11%update(k)
       call device_col3(this%stats_work%x_d,this%dvdx%x_d, this%dvdx%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dvdy%x_d, this%dvdy%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dvdz%x_d, this%dvdz%x_d,n)
       call this%e22%update(k)
       call device_col3(this%stats_work%x_d,this%dwdx%x_d, this%dwdx%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dwdy%x_d, this%dwdy%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dwdz%x_d, this%dwdz%x_d,n)
       call this%e33%update(k)
       call device_col3(this%stats_work%x_d,this%dudx%x_d, this%dvdx%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dudy%x_d, this%dvdy%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dudz%x_d, this%dvdz%x_d,n)
       call this%e12%update(k)
       call device_col3(this%stats_work%x_d,this%dvdx%x_d, this%dwdx%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dvdy%x_d, this%dwdy%x_d,n)
       call device_addcol3(this%stats_work%x_d,this%dvdz%x_d, this%dwdz%x_d,n)
       call this%e23%update(k)


    else
       call col3(stats_work%x,this%dudx%x, this%p%x,n)
       call this%pdudx%update(k)
       call col3(stats_work%x,this%dudy%x, this%p%x,n)
       call this%pdudy%update(k)
       call col3(stats_work%x,this%dudz%x, this%p%x,n)
       call this%pdudz%update(k)

       call col3(stats_work%x,this%dvdx%x, this%p%x,n)
       call this%pdvdx%update(k)
       call col3(stats_work%x,this%dvdy%x, this%p%x,n)
       call this%pdvdy%update(k)
       call col3(stats_work%x,this%dvdz%x, this%p%x,n)
       call this%pdvdz%update(k)

       call col3(stats_work%x,this%dwdx%x, this%p%x,n)
       call this%pdwdx%update(k)
       call col3(stats_work%x,this%dwdy%x, this%p%x,n)
       call this%pdwdy%update(k)
       call col3(stats_work%x,this%dwdz%x, this%p%x,n)
       call this%pdwdz%update(k)

       call col3(this%stats_work%x,this%dudx%x, this%dudx%x,n)
       call addcol3(this%stats_work%x,this%dudy%x, this%dudy%x,n)
       call addcol3(this%stats_work%x,this%dudz%x, this%dudz%x,n)
       call this%e11%update(k)
       call col3(this%stats_work%x,this%dvdx%x, this%dvdx%x,n)
       call addcol3(this%stats_work%x,this%dvdy%x, this%dvdy%x,n)
       call addcol3(this%stats_work%x,this%dvdz%x, this%dvdz%x,n)
       call this%e22%update(k)
       call col3(this%stats_work%x,this%dwdx%x, this%dwdx%x,n)
       call addcol3(this%stats_work%x,this%dwdy%x, this%dwdy%x,n)
       call addcol3(this%stats_work%x,this%dwdz%x, this%dwdz%x,n)
       call this%e33%update(k)
       call col3(this%stats_work%x,this%dudx%x, this%dvdx%x,n)
       call addcol3(this%stats_work%x,this%dudy%x, this%dvdy%x,n)
       call addcol3(this%stats_work%x,this%dudz%x, this%dvdz%x,n)
       call this%e12%update(k)
       call col3(this%stats_work%x,this%dvdx%x, this%dwdx%x,n)
       call addcol3(this%stats_work%x,this%dvdy%x, this%dwdy%x,n)
       call addcol3(this%stats_work%x,this%dvdz%x, this%dwdz%x,n)
       call this%e23%update(k)

    end if

  end associate

  end subroutine fluid_stats_update


  !> Deallocates a mean flow field
  subroutine fluid_stats_free(this)
    class(fluid_stats_t), intent(inout) :: this

    call this%stats_work%free()
    call this%stats_u%free()
    call this%stats_v%free()
    call this%stats_w%free()
    
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
 
  !> Initialize a mean flow field
  subroutine fluid_stats_reset(this)
    class(fluid_stats_t), intent(inout), target:: this
    
    call this%uu%reset()
    call this%vv%reset()
    call this%ww%reset()
    call this%uv%reset()
    call this%uw%reset()
    call this%vw%reset()
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
    !> Pressure
    call this%pp%reset()
    call this%ppp%reset()
    call this%pppp%reset()
     !> Pressure * velocity
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

  end subroutine fluid_stats_reset

  subroutine fluid_stats_make_strong_grad(this)
    class(fluid_stats_t) :: this
    integer :: n
    n = size(this%coef%B)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(this%stats_work%x_d, 1.0_rp,n)
       call device_invcol2(this%stats_work%x_d, this%coef%B_d,n)
       call device_col2(this%pdudx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdudy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdudz%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdvdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdvdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdvdz%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdwdx%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdwdy%mf%x_d, this%stats_work%x_d, n)
       call device_col2(this%pdwdz%mf%x_d, this%stats_work%x_d, n)

       call device_col2(this%stats_work%x_d, this%stats_work%x_d,n)
       call device_col2(this%e11%mf%x_d,this%stats_work%x_d, n)
       call device_col2(this%e22%mf%x_d,this%stats_work%x_d, n)
       call device_col2(this%e33%mf%x_d,this%stats_work%x_d, n)
       call device_col2(this%e12%mf%x_d,this%stats_work%x_d, n)
       call device_col2(this%e13%mf%x_d,this%stats_work%x_d, n)
       call device_col2(this%e23%mf%x_d,this%stats_work%x_d, n)


    else 
       call invers2(this%stats_work%x, this%coef%B,n)
       call col2(this%pdudx%mf%x, this%stats_work%x, n)
       call col2(this%pdudy%mf%x, this%stats_work%x, n)
       call col2(this%pdudz%mf%x, this%stats_work%x, n)
       call col2(this%pdvdx%mf%x, this%stats_work%x, n)
       call col2(this%pdvdy%mf%x, this%stats_work%x, n)
       call col2(this%pdvdz%mf%x, this%stats_work%x, n)
       call col2(this%pdwdx%mf%x, this%stats_work%x, n)
       call col2(this%pdwdy%mf%x, this%stats_work%x, n)
       call col2(this%pdwdz%mf%x, this%stats_work%x, n)

       call col2(this%stats_work%x, this%stats_work%x,n)
       call col2(this%e11%mf%x,this%stats_work%x, n)
       call col2(this%e22%mf%x,this%stats_work%x, n)
       call col2(this%e33%mf%x,this%stats_work%x, n)
       call col2(this%e12%mf%x,this%stats_work%x, n)
       call col2(this%e13%mf%x,this%stats_work%x, n)
       call col2(this%e23%mf%x,this%stats_work%x, n)

    end if

  end subroutine fluid_stats_make_strong_grad
 
end module fluid_stats
