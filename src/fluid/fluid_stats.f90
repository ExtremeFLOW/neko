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
!!
!! Without an averaging direction, every statistic is a `mean_field_t`,
!! a 3D field in the registry. With an averaging direction, the products
!! sampled at every step are summed over the homogeneous direction(s) on
!! the device and accumulated in a 2D (`map_2d_t`) or 1D (`map_1d_t`)
!! accumulator instead, so that no 3D statistics fields are kept, unless
!! the 3D fields are explicitly kept and averaged when written.
module fluid_stats
  use mean_field, only : mean_field_t
  use device_math, only : device_col2, device_cfill, device_invcol2, &
       device_glsc2
  use field_math, only : field_col3, field_col2, field_addcol3, &
       field_copy, field_cadd
  use num_types, only : rp
  use math, only : copy, subcol3, glsc2, invcol2
  use operators, only : opgrad
  use coefs, only : coef_t
  use field, only : field_t
  use field_list, only : field_list_t
  use stats_quant, only : stats_quant_t
  use scratch_registry, only : neko_scratch_registry
  use map_1d, only : map_1d_t
  use map_2d, only : map_2d_t
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_warning, neko_error
  implicit none
  private

  !> Slots of the statistics, which is also their order in the output.
  integer, parameter :: S_P = 1, S_U = 2, S_V = 3, S_W = 4, S_PP = 5, &
       S_UU = 6, S_VV = 7, S_WW = 8, S_UV = 9, S_UW = 10, S_VW = 11, &
       S_UUU = 12, S_VVV = 13, S_WWW = 14, S_UUV = 15, S_UUW = 16, &
       S_UVV = 17, S_UVW = 18, S_VVW = 19, S_UWW = 20, S_VWW = 21, &
       S_UUUU = 22, S_VVVV = 23, S_WWWW = 24, S_PPP = 25, S_PPPP = 26, &
       S_PU = 27, S_PV = 28, S_PW = 29, S_PDUDX = 30, S_PDUDY = 31, &
       S_PDUDZ = 32, S_PDVDX = 33, S_PDVDY = 34, S_PDVDZ = 35, S_PDWDX = 36, &
       S_PDWDY = 37, S_PDWDZ = 38, S_E11 = 39, S_E22 = 40, S_E33 = 41, &
       S_E12 = 42, S_E13 = 43, S_E23 = 44

  !> Pointer to a mean field.
  type :: mean_field_ptr_t
     type(mean_field_t), pointer :: ptr => null()
  end type mean_field_ptr_t

  type, public, extends(stats_quant_t) :: fluid_stats_t
     !> Work fields borrowed from the scratch registry while sampling: the
     !! sampled product and the square it is built from.
     type(field_t), pointer :: prod => null()
     type(field_t), pointer :: sqr => null()
     !> Pressure in the requested gauge, borrowed while sampling.
     type(field_t), pointer :: p_gauged => null()

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
     !> The mean fields in the order of the output, used without averaging.
     type(mean_field_ptr_t), allocatable :: means(:)
     !> Gradient fields borrowed while sampling, holding the gradients of
     !! two velocity components (a and b) at a time, see `update`.
     type(field_t), pointer :: dadx => null()
     type(field_t), pointer :: dady => null()
     type(field_t), pointer :: dadz => null()
     type(field_t), pointer :: dbdx => null()
     type(field_t), pointer :: dbdy => null()
     type(field_t), pointer :: dbdz => null()

     !> SEM coefficients.
     type(coef_t), pointer :: coef
     !> Number of statistical fields to be computed.
     integer :: n_stats = 44
     !> Specifies a subset of the statistics to be collected. All 44 fields by
     !! default.
     character(5) :: stat_set
     !> A list of size n_stats, with entries pointing to the fields that will
     !! be output (the field components above.) Used to write the output
     !! without averaging.
     type(field_list_t) :: stat_fields
     !> Subtract the volume-weighted mean of the pressure before sampling.
     logical :: volume_mean_gauge = .false.
     !> Scratch registry indices of the fields borrowed while sampling.
     integer, allocatable :: work_idx(:)
     !> Dimension of the stored statistics: 3 for 3D mean fields, 2 or 1
     !! when accumulated in the averaged space.
     integer :: avg_dim = 3
     !> Dimension of the written statistics: 3 without averaging, 2 when
     !! averaged in one direction and 1 when averaged in two directions.
     integer :: output_dim = 3
     !> Map of the averaging in one direction, accumulating the statistics
     !! when they are not kept as 3D fields.
     type(map_2d_t) :: map_2d
     !> Map of the averaging in two directions, accumulating the statistics
     !! when they are not kept as 3D fields.
     type(map_1d_t) :: map_1d
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
     !> Borrow the work fields from the scratch registry.
     procedure, private, pass(this) :: acquire_work => fluid_stats_acquire_work
     !> Return the work fields to the scratch registry.
     procedure, private, pass(this) :: release_work => fluid_stats_release_work
     !> Add a sample of one statistic.
     procedure, private, pass(this) :: sample => fluid_stats_sample
     !> Compute the gradient of a field for the statistics.
     procedure, private, pass(this) :: gradient => fluid_stats_gradient
     !> Sample the products of the pressure with a gradient.
     procedure, private, pass(this) :: sample_products => &
          fluid_stats_sample_products
     !> Sample the dot product of two gradients.
     procedure, private, pass(this) :: sample_dot => fluid_stats_sample_dot
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
  !! @param name Name of the statistics, used to prefix the mean fields in
  !! the registry. Optional.
  !! @param pressure_gauge Gauge of the pressure entering the statistics.
  !! Optional. Either `solver`, the pressure as computed by the solver, or
  !! `volume_mean`, the pressure shifted to have a zero volume-weighted mean
  !! at every sample. Defaults to `solver`.
  !! @param avg_direction Direction(s) to average the statistics in,
  !! `x`, `y`, `z`, `xy`, `xz`, `yz` or `none`. Optional, defaults to `none`.
  !! With an averaging direction the statistics are accumulated directly in
  !! the averaged space and no 3D mean fields are created, unless
  !! `keep_3d_fields` is true: the 3D mean fields are then kept, as without
  !! an averaging direction, and averaged when written.
  !! @param keep_3d_fields Keep the statistics as 3D mean fields with an
  !! averaging direction. Optional, false by default.
  subroutine fluid_stats_init(this, coef, u, v, w, p, set, name, &
       pressure_gauge, avg_direction, keep_3d_fields)
    class(fluid_stats_t), intent(inout), target:: this
    type(coef_t), target, optional :: coef
    type(field_t), target, intent(in) :: u, v, w, p
    character(*), intent(in), optional :: set
    character(*), intent(in), optional :: name
    character(*), intent(in), optional :: pressure_gauge
    character(*), intent(in), optional :: avg_direction
    logical, intent(in), optional :: keep_3d_fields

    character(len=1024) :: unique_name
    logical :: keep_3d
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

    if (present(pressure_gauge)) then
       select case (trim(pressure_gauge))
       case ('solver')
          this%volume_mean_gauge = .false.
       case ('volume_mean')
          this%volume_mean_gauge = .true.
       case default
          call neko_error("fluid_stats: unknown pressure_gauge '" // &
               trim(pressure_gauge) // "', use 'solver' or 'volume_mean'")
       end select
    end if

    keep_3d = .false.
    if (present(keep_3d_fields)) keep_3d = keep_3d_fields

    this%output_dim = 3
    if (present(avg_direction)) then
       select case (trim(avg_direction))
       case ('x', 'y', 'z')
          this%output_dim = 2
          call this%map_2d%init_char(coef, avg_direction, 1e-7_rp)
          if (.not. keep_3d) call this%map_2d%accumulate_init(this%n_stats)
       case ('xy', 'yx', 'xz', 'zx', 'yz', 'zy')
          this%output_dim = 1
          call this%map_1d%init_char(coef, avg_direction, 1e-7_rp)
          if (.not. keep_3d) call this%map_1d%accumulate_init(this%n_stats)
       case ('none', '')
          this%output_dim = 3
       case default
          call neko_error("fluid_stats: unknown avg_direction '" // &
               trim(avg_direction) // "'")
       end select
    end if
    this%avg_dim = 3
    if (.not. keep_3d) this%avg_dim = this%output_dim

    if (this%avg_dim .ne. 3) return

    ! The field sampled by each product statistic is a work field that is
    ! borrowed from the scratch registry while sampling, see acquire_work.
    ! Here the means are only created, with the velocity as a placeholder.
    call this%u_mean%init(this%u, trim(unique_name) // 'mean_u')
    call this%v_mean%init(this%v, trim(unique_name) // 'mean_v')
    call this%w_mean%init(this%w, trim(unique_name) // 'mean_w')
    call this%p_mean%init(this%p, trim(unique_name) // 'mean_p')
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

    allocate(this%means(this%n_stats))
    this%means(S_P)%ptr => this%p_mean
    this%means(S_U)%ptr => this%u_mean
    this%means(S_V)%ptr => this%v_mean
    this%means(S_W)%ptr => this%w_mean
    this%means(S_PP)%ptr => this%pp
    this%means(S_UU)%ptr => this%uu
    this%means(S_VV)%ptr => this%vv
    this%means(S_WW)%ptr => this%ww
    this%means(S_UV)%ptr => this%uv
    this%means(S_UW)%ptr => this%uw
    this%means(S_VW)%ptr => this%vw

    if (this%n_stats .eq. 44) then
       this%means(S_UUU)%ptr => this%uuu
       this%means(S_VVV)%ptr => this%vvv
       this%means(S_WWW)%ptr => this%www
       this%means(S_UUV)%ptr => this%uuv
       this%means(S_UUW)%ptr => this%uuw
       this%means(S_UVV)%ptr => this%uvv
       this%means(S_UVW)%ptr => this%uvw
       this%means(S_VVW)%ptr => this%vvw
       this%means(S_UWW)%ptr => this%uww
       this%means(S_VWW)%ptr => this%vww
       this%means(S_UUUU)%ptr => this%uuuu
       this%means(S_VVVV)%ptr => this%vvvv
       this%means(S_WWWW)%ptr => this%wwww
       this%means(S_PPP)%ptr => this%ppp
       this%means(S_PPPP)%ptr => this%pppp
       this%means(S_PU)%ptr => this%pu
       this%means(S_PV)%ptr => this%pv
       this%means(S_PW)%ptr => this%pw
       this%means(S_PDUDX)%ptr => this%pdudx
       this%means(S_PDUDY)%ptr => this%pdudy
       this%means(S_PDUDZ)%ptr => this%pdudz
       this%means(S_PDVDX)%ptr => this%pdvdx
       this%means(S_PDVDY)%ptr => this%pdvdy
       this%means(S_PDVDZ)%ptr => this%pdvdz
       this%means(S_PDWDX)%ptr => this%pdwdx
       this%means(S_PDWDY)%ptr => this%pdwdy
       this%means(S_PDWDZ)%ptr => this%pdwdz
       this%means(S_E11)%ptr => this%e11
       this%means(S_E22)%ptr => this%e22
       this%means(S_E33)%ptr => this%e33
       this%means(S_E12)%ptr => this%e12
       this%means(S_E13)%ptr => this%e13
       this%means(S_E23)%ptr => this%e23
    end if

    call this%stat_fields%init(this%n_stats)
    block
      integer :: i
      do i = 1, this%n_stats
         call this%stat_fields%assign_to_field(i, this%means(i)%ptr%mf)
      end do
    end block

  end subroutine fluid_stats_init

  !> Borrows the work fields from the scratch registry: two for the
  !! products, six for the gradients with the full set, and one for the
  !! gauged pressure.
  subroutine fluid_stats_acquire_work(this)
    class(fluid_stats_t), intent(inout) :: this
    integer :: n_work, i

    n_work = 2
    if (this%n_stats .eq. 44) n_work = n_work + 6
    if (this%volume_mean_gauge) n_work = n_work + 1
    allocate(this%work_idx(n_work))

    i = 1
    call neko_scratch_registry%request_field(this%prod, this%work_idx(i), &
         .false.)
    i = i + 1
    call neko_scratch_registry%request_field(this%sqr, this%work_idx(i), &
         .false.)

    if (this%n_stats .eq. 44) then
       i = i + 1
       call neko_scratch_registry%request_field(this%dadx, &
            this%work_idx(i), .false.)
       i = i + 1
       call neko_scratch_registry%request_field(this%dady, &
            this%work_idx(i), .false.)
       i = i + 1
       call neko_scratch_registry%request_field(this%dadz, &
            this%work_idx(i), .false.)
       i = i + 1
       call neko_scratch_registry%request_field(this%dbdx, &
            this%work_idx(i), .false.)
       i = i + 1
       call neko_scratch_registry%request_field(this%dbdy, &
            this%work_idx(i), .false.)
       i = i + 1
       call neko_scratch_registry%request_field(this%dbdz, &
            this%work_idx(i), .false.)
    end if

    if (this%volume_mean_gauge) then
       i = i + 1
       call neko_scratch_registry%request_field(this%p_gauged, &
            this%work_idx(i), .false.)
    end if

  end subroutine fluid_stats_acquire_work

  !> Returns the work fields to the scratch registry.
  subroutine fluid_stats_release_work(this)
    class(fluid_stats_t), intent(inout) :: this

    if (allocated(this%work_idx)) then
       call neko_scratch_registry%relinquish_field(this%work_idx)
       deallocate(this%work_idx)
    end if

    nullify(this%prod, this%sqr, this%p_gauged)
    nullify(this%dadx, this%dady, this%dadz)
    nullify(this%dbdx, this%dbdy, this%dbdz)

  end subroutine fluid_stats_release_work

  !> Adds a sample of the statistic in slot `slot`, either to its mean
  !! field or, with an averaging direction, to the accumulated averages.
  !! @param slot Slot of the statistic.
  !! @param f Sampled field.
  !! @param k Time elapsed since the last sample.
  subroutine fluid_stats_sample(this, slot, f, k)
    class(fluid_stats_t), intent(inout) :: this
    integer, intent(in) :: slot
    type(field_t), intent(in), target :: f
    real(kind=rp), intent(in) :: k

    select case (this%avg_dim)
    case (3)
       this%means(slot)%ptr%f => f
       call this%means(slot)%ptr%update(k)
    case (2)
       call this%map_2d%accumulate(f, slot, k)
    case (1)
       call this%map_1d%accumulate(f, slot, k)
    end select

  end subroutine fluid_stats_sample

  !> Updates all fields with a new sample.
  !! @param k Time elapsed since the last update.
  !!
  !! The products are formed in two borrowed work fields, and the gradient
  !! statistics in six borrowed gradient fields that hold the gradients of
  !! two velocity components at a time: u and v first, then w in place of
  !! v, and finally v in place of u. The gradient of v is thus computed
  !! twice, but eight fields are borrowed at once instead of fourteen.
  subroutine fluid_stats_update(this, k)
    class(fluid_stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k
    type(field_t), pointer :: p
    real(kind=rp) :: p_mean_vol
    integer :: n

    call this%acquire_work()

    associate(prod => this%prod, sqr => this%sqr, u => this%u, v => this%v, &
         w => this%w, dadx => this%dadx, dady => this%dady, &
         dadz => this%dadz, dbdx => this%dbdx, dbdy => this%dbdy, &
         dbdz => this%dbdz)
      n = prod%dof%size()

      ! The averages are normalised by the accumulated volumes.
      if (this%avg_dim .eq. 2) call this%map_2d%accumulate_volume(k)
      if (this%avg_dim .eq. 1) call this%map_1d%accumulate_volume(k)

      ! Shift the pressure to the requested gauge.
      if (this%volume_mean_gauge) then
         p => this%p_gauged
         if (NEKO_BCKND_DEVICE .eq. 1) then
            p_mean_vol = device_glsc2(this%p%x_d, this%coef%B_d, n) / &
                 this%coef%volume
         else
            p_mean_vol = glsc2(this%p%x, this%coef%B, n) / this%coef%volume
         end if
         call field_copy(p, this%p, n)
         call field_cadd(p, -p_mean_vol, n)
      else
         p => this%p
      end if

      call this%sample(S_U, u, k)
      call this%sample(S_V, v, k)
      call this%sample(S_W, w, k)
      call this%sample(S_P, p, k)

      ! Second moments, with the third and fourth moments formed from the
      ! square while it is available.
      call field_col3(sqr, u, u, n)
      call this%sample(S_UU, sqr, k)
      if (this%n_stats .eq. 44) then
         call field_col3(prod, sqr, u, n)
         call this%sample(S_UUU, prod, k)
         call field_col3(prod, sqr, v, n)
         call this%sample(S_UUV, prod, k)
         call field_col3(prod, sqr, w, n)
         call this%sample(S_UUW, prod, k)
         call field_col3(prod, sqr, sqr, n)
         call this%sample(S_UUUU, prod, k)
      end if

      call field_col3(sqr, v, v, n)
      call this%sample(S_VV, sqr, k)
      if (this%n_stats .eq. 44) then
         call field_col3(prod, sqr, v, n)
         call this%sample(S_VVV, prod, k)
         call field_col3(prod, sqr, u, n)
         call this%sample(S_UVV, prod, k)
         call field_col3(prod, sqr, w, n)
         call this%sample(S_VVW, prod, k)
         call field_col3(prod, sqr, sqr, n)
         call this%sample(S_VVVV, prod, k)
      end if

      call field_col3(sqr, w, w, n)
      call this%sample(S_WW, sqr, k)
      if (this%n_stats .eq. 44) then
         call field_col3(prod, sqr, w, n)
         call this%sample(S_WWW, prod, k)
         call field_col3(prod, sqr, u, n)
         call this%sample(S_UWW, prod, k)
         call field_col3(prod, sqr, v, n)
         call this%sample(S_VWW, prod, k)
         call field_col3(prod, sqr, sqr, n)
         call this%sample(S_WWWW, prod, k)
      end if

      call field_col3(sqr, p, p, n)
      call this%sample(S_PP, sqr, k)
      if (this%n_stats .eq. 44) then
         call field_col3(prod, sqr, p, n)
         call this%sample(S_PPP, prod, k)
         call field_col3(prod, sqr, sqr, n)
         call this%sample(S_PPPP, prod, k)
      end if

      call field_col3(prod, u, v, n)
      call this%sample(S_UV, prod, k)
      if (this%n_stats .eq. 44) then
         call field_col2(prod, w, n)
         call this%sample(S_UVW, prod, k)
      end if
      call field_col3(prod, u, w, n)
      call this%sample(S_UW, prod, k)
      call field_col3(prod, v, w, n)
      call this%sample(S_VW, prod, k)

      if (this%n_stats .eq. 44) then
         call field_col3(prod, p, u, n)
         call this%sample(S_PU, prod, k)
         call field_col3(prod, p, v, n)
         call this%sample(S_PV, prod, k)
         call field_col3(prod, p, w, n)
         call this%sample(S_PW, prod, k)

         ! The gradients of u (a) and v (b).
         call this%gradient(u, dadx, dady, dadz)
         call this%gradient(v, dbdx, dbdy, dbdz)
         call this%sample_products(p, dadx, dady, dadz, &
              S_PDUDX, S_PDUDY, S_PDUDZ, k)
         call this%sample_products(p, dbdx, dbdy, dbdz, &
              S_PDVDX, S_PDVDY, S_PDVDZ, k)
         call this%sample_dot(dadx, dady, dadz, dadx, dady, dadz, S_E11, k)
         call this%sample_dot(dbdx, dbdy, dbdz, dbdx, dbdy, dbdz, S_E22, k)
         call this%sample_dot(dadx, dady, dadz, dbdx, dbdy, dbdz, S_E12, k)

         ! The gradient of w (b) replaces that of v.
         call this%gradient(w, dbdx, dbdy, dbdz)
         call this%sample_products(p, dbdx, dbdy, dbdz, &
              S_PDWDX, S_PDWDY, S_PDWDZ, k)
         call this%sample_dot(dbdx, dbdy, dbdz, dbdx, dbdy, dbdz, S_E33, k)
         call this%sample_dot(dadx, dady, dadz, dbdx, dbdy, dbdz, S_E13, k)

         ! The gradient of v (a) replaces that of u.
         call this%gradient(v, dadx, dady, dadz)
         call this%sample_dot(dadx, dady, dadz, dbdx, dbdy, dbdz, S_E23, k)
      end if
    end associate

    call this%release_work()

  end subroutine fluid_stats_update

  !> Computes the gradient of @a f into @a dfdx, @a dfdy and @a dfdz.
  !! The weak gradient of `opgrad` is converted to the strong one here
  !! when the statistics are accumulated in the averaged space, and at
  !! the output otherwise, see `make_strong_grad`.
  subroutine fluid_stats_gradient(this, f, dfdx, dfdy, dfdz)
    class(fluid_stats_t), intent(inout) :: this
    type(field_t), intent(in) :: f
    type(field_t), intent(inout) :: dfdx, dfdy, dfdz
    integer :: n

    call opgrad(dfdx%x, dfdy%x, dfdz%x, f%x, this%coef)

    if (this%avg_dim .ne. 3) then
       n = f%dof%size()
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_invcol2(dfdx%x_d, this%coef%B_d, n)
          call device_invcol2(dfdy%x_d, this%coef%B_d, n)
          call device_invcol2(dfdz%x_d, this%coef%B_d, n)
       else
          call invcol2(dfdx%x, this%coef%B, n)
          call invcol2(dfdy%x, this%coef%B, n)
          call invcol2(dfdz%x, this%coef%B, n)
       end if
    end if

  end subroutine fluid_stats_gradient

  !> Samples the products of @a p with the three gradient components into
  !! the slots @a sx, @a sy and @a sz.
  subroutine fluid_stats_sample_products(this, p, gx, gy, gz, sx, sy, sz, k)
    class(fluid_stats_t), intent(inout) :: this
    type(field_t), intent(in) :: p, gx, gy, gz
    integer, intent(in) :: sx, sy, sz
    real(kind=rp), intent(in) :: k
    integer :: n

    n = p%dof%size()
    call field_col3(this%prod, p, gx, n)
    call this%sample(sx, this%prod, k)
    call field_col3(this%prod, p, gy, n)
    call this%sample(sy, this%prod, k)
    call field_col3(this%prod, p, gz, n)
    call this%sample(sz, this%prod, k)

  end subroutine fluid_stats_sample_products

  !> Samples the dot product of two gradients into @a slot.
  subroutine fluid_stats_sample_dot(this, ax, ay, az, bx, by, bz, slot, k)
    class(fluid_stats_t), intent(inout) :: this
    type(field_t), intent(in) :: ax, ay, az, bx, by, bz
    integer, intent(in) :: slot
    real(kind=rp), intent(in) :: k
    integer :: n

    n = ax%dof%size()
    call field_col3(this%prod, ax, bx, n)
    call field_addcol3(this%prod, ay, by, n)
    call field_addcol3(this%prod, az, bz, n)
    call this%sample(slot, this%prod, k)

  end subroutine fluid_stats_sample_dot


  !> Destructor.
  subroutine fluid_stats_free(this)
    class(fluid_stats_t), intent(inout) :: this

    call this%release_work()

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

    if (allocated(this%means)) deallocate(this%means)

    call this%map_2d%free()
    call this%map_1d%free()
    this%avg_dim = 3
    this%output_dim = 3

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
    integer :: i

    select case (this%avg_dim)
    case (3)
       do i = 1, this%n_stats
          call this%means(i)%ptr%reset()
       end do
    case (2)
       call this%map_2d%accumulate_reset()
    case (1)
       call this%map_1d%accumulate_reset()
    end select

  end subroutine fluid_stats_reset

  ! Convert computed weak gradients to strong. Only without averaging,
  ! otherwise the gradients are made strong when sampling.
  subroutine fluid_stats_make_strong_grad(this)
    class(fluid_stats_t) :: this
    type(field_t), pointer :: work
    integer :: n, i, idx
    real(kind=rp) :: wrk, wrk_sqr

    if (this%n_stats .eq. 11 .or. this%avg_dim .ne. 3) return

    n = size(this%coef%B)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_scratch_registry%request_field(work, idx, .false.)
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
  !! fields. Only available without an averaging direction.
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
    type(field_t), pointer :: dudx, dudy, dudz, dvdx, dvdy, dvdz
    type(field_t), pointer :: dwdx, dwdy, dwdz
    integer :: grad_idx(9)
    integer :: n, i
    real(kind=rp) :: wrk

    if (this%avg_dim .ne. 3) then
       call neko_error('fluid_stats: post_process requires statistics ' // &
            'without an averaging direction')
    end if

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
       call neko_scratch_registry%request_field(dudx, grad_idx(1), .false.)
       call neko_scratch_registry%request_field(dudy, grad_idx(2), .false.)
       call neko_scratch_registry%request_field(dudz, grad_idx(3), .false.)
       call neko_scratch_registry%request_field(dvdx, grad_idx(4), .false.)
       call neko_scratch_registry%request_field(dvdy, grad_idx(5), .false.)
       call neko_scratch_registry%request_field(dvdz, grad_idx(6), .false.)
       call neko_scratch_registry%request_field(dwdx, grad_idx(7), .false.)
       call neko_scratch_registry%request_field(dwdy, grad_idx(8), .false.)
       call neko_scratch_registry%request_field(dwdz, grad_idx(9), .false.)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%u_mean%mf%x, this%u_mean%mf%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%v_mean%mf%x, this%v_mean%mf%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%w_mean%mf%x, this%w_mean%mf%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
          call opgrad(dudx%x, dudy%x, dudz%x, this%u_mean%mf%x, this%coef)
          call opgrad(dvdx%x, dvdy%x, dvdz%x, this%v_mean%mf%x, this%coef)
          call opgrad(dwdx%x, dwdy%x, dwdz%x, this%w_mean%mf%x, this%coef)
          call device_memcpy(dudx%x, dudx%x_d, n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dvdx%x, dvdx%x_d, n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dwdx%x, dwdx%x_d, n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dudy%x, dudy%x_d, n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dvdy%x, dvdy%x_d, n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dwdy%x, dwdy%x_d, n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dudz%x, dudz%x_d, n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dvdz%x, dvdz%x_d, n, DEVICE_TO_HOST, &
               sync = .false.)
          call device_memcpy(dwdz%x, dwdz%x_d, n, DEVICE_TO_HOST, &
               sync = .true.)
       else
          call opgrad(dudx%x, dudy%x, dudz%x, this%u_mean%mf%x, this%coef)
          call opgrad(dvdx%x, dvdy%x, dvdz%x, this%v_mean%mf%x, this%coef)
          call opgrad(dwdx%x, dwdy%x, dwdz%x, this%w_mean%mf%x, this%coef)
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
       call neko_scratch_registry%relinquish_field(grad_idx)

    end if

    if (present(dissipation_tensor)) then
       call neko_warning('Dissipation tensor stat not implemented'// &
            ' in fluid_stats, process stats in python instead')
    end if

  end subroutine fluid_stats_post_process

end module fluid_stats
