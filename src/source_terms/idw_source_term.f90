! Copyright (c) 2023-2026, The Neko Authors
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
!> Implements an inverse distance weighting based source term
module idw_source_term
  use num_types, only : rp, dp
  use field_list, only : field_list_t
  use json_module, only : json_file, json_value, json_core
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use source_term, only : source_term_t
  use field_math, only : field_col2, field_col3, field_copy, field_rzero, &
       field_add2
  use coefs, only : coef_t
  use field, only : field_t
  use utils, only : neko_error
  use tri_mesh, only : tri_mesh_t
  use registry, only : neko_registry
  use field, only : field_t
  use file, only : file_t
  use comm, only : NEKO_COMM, pe_rank, MPI_REAL_PRECISION
  use intersection_detector, only : intersect_detector_t
  use mesh, only : mesh_t
  use stack, only : stack_i4_t, stack_pt_t
  use global_interpolation, only : global_interpolation_t
  use point, only : point_t
  use math, only : NEKO_EPS, glsc2, col2
  use aabb, only : aabb_t, get_aabb
  use time_state, only : time_state_t
  use time_based_controller, only : time_based_controller_t
  use vector, only : vector_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use elementwise_filter, only : elementwise_filter_t
  use PDE_filter, only : PDE_filter_t
  use filter, only : filter_t
  use device_math, only : device_col2, device_glsc2, device_col3, &
       device_cmult, device_add2, device_rzero
  use device_mathops, only : device_opcolv
  use mathops, only : opcolv
  use profiler, only : profiler_start_region, profiler_end_region
  use device, only : device_free, device_map, device_memcpy, &
       device_event_sync, glb_cmd_event, HOST_TO_DEVICE, DEVICE_TO_HOST
  use device_idw_source_term, only : device_idw_gather, &
       device_idw_interp_partials
  use gather_scatter, only : gs_t, GS_OP_ADD
  use mpi_f08, only : MPI_Allreduce, MPI_IN_PLACE, MPI_MIN, &
       MPI_MAX, MPI_INTEGER, MPI_SUM, MPI_DOUBLE_PRECISION
  use logger, only : neko_log, LOG_SIZE
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Inverse distance weighting source term.
  type, public, extends(source_term_t) :: idw_source_term_t
     !> Smallest distance between between points and dofs
     real(kind=dp) :: ds_min
     real(kind=dp) :: ds_max
     type(intersect_detector_t) :: intersect
     type(global_interpolation_t) :: global_interp
     type(point_t), allocatable :: lag_pts(:)
     type(point_t), allocatable :: lag_nrm(:)
     type(stack_i4_t), allocatable :: lag_el(:)
     real(kind=rp), allocatable :: xyz(:,:)
     real(kind=rp), allocatable :: fu_ib(:)
     real(kind=rp), allocatable :: fv_ib(:)
     real(kind=rp), allocatable :: fw_ib(:)
     real(kind=rp), allocatable :: fum_ib(:)
     real(kind=rp), allocatable :: fvm_ib(:)
     real(kind=rp), allocatable :: fwm_ib(:)
     type(c_ptr) :: fu_ib_d = C_NULL_PTR
     type(c_ptr) :: fv_ib_d = C_NULL_PTR
     type(c_ptr) :: fw_ib_d = C_NULL_PTR
     type(c_ptr) :: fum_ib_d = C_NULL_PTR
     type(c_ptr) :: fvm_ib_d = C_NULL_PTR
     type(c_ptr) :: fwm_ib_d = C_NULL_PTR
     !> CSR transpose of lag_el (per element -> lag points) for the gather kernel
     integer, allocatable :: el_off(:)
     integer, allocatable :: el_lag(:)
     integer, allocatable :: active_el(:)
     !> Lagrangian point coordinates, unpacked for the device kernel
     real(kind=rp), allocatable :: lpx(:), lpy(:), lpz(:)
     type(c_ptr) :: el_off_d = C_NULL_PTR
     type(c_ptr) :: el_lag_d = C_NULL_PTR
     type(c_ptr) :: active_el_d = C_NULL_PTR
     type(c_ptr) :: lpx_d = C_NULL_PTR
     type(c_ptr) :: lpy_d = C_NULL_PTR
     type(c_ptr) :: lpz_d = C_NULL_PTR
     !> CSR of lag_el (per lag point -> elements, 0-based) walked by the
     !! device interpolation kernel
     integer, allocatable :: lag_off(:)
     integer, allocatable :: lag_els(:)
     !> Shepard partial sums, 8 per lag point, as the device kernel
     !! returns them (see idw_interp_shepard_partials for the layout)
     real(kind=rp), allocatable :: part(:,:)
     type(c_ptr) :: lag_off_d = C_NULL_PTR
     type(c_ptr) :: lag_els_d = C_NULL_PTR
     type(c_ptr) :: part_d = C_NULL_PTR
     type(c_ptr) :: Bm_d = C_NULL_PTR
     integer :: n_active = 0
     integer :: lx3 = 0
     real(kind=rp) :: pwr_param
     real(kind=rp) :: rmax
     type(field_t) :: w
     type(field_t) :: wm
     type(field_t) :: ds
     type(field_t) :: mmsk
     type(field_t) :: pmsk
     type(field_t) :: tmp
     !> Immersed boundary forcing, spread, assembled and filtered on its
     !! own before being added to the right-hand side fields
     type(field_t) :: ib_fx
     type(field_t) :: ib_fy
     type(field_t) :: ib_fz
     !> Write the force on the immersed objects to a csv file
     logical :: force_output = .false.
     !> Force on the immersed objects, -rho times the integral of the
     !! IB forcing, times force_scale
     real(kind=rp) :: force(3) = 0.0_rp
     real(kind=rp) :: force_scale = 1.0_rp
     !> Registry name of the fluid density
     character(len=:), allocatable :: rho_name
     type(file_t) :: force_file
     type(vector_t) :: force_row
     type(time_based_controller_t) :: force_controller
     type(gs_t) :: gs
     logical :: one_sided
     !> Use Shepard IDW interpolation onto the markers instead of the
     !! global (rst-based) interpolation
     logical :: idw_interp = .false.
     !> Shepard weights K B mult / w, which make the interpolation the
     !! mass-weighted adjoint of the spread so that the forcing never adds
     !! kinetic energy (implies `idw_interp`, and `interp_rmax` = `rmax`)
     logical :: adjoint_interp = .false.
     !> Spread the marker forcing with the mass-weighted adjoint of the
     !! spectral interpolation, f = -(1/dt) mask Binv gs_add(I^T D u_m),
     !! instead of the IDW kernel
     logical :: adjoint_spread = .false.
     !> Take the adjoint in the pointwise GLL mass inner product (L2
     !! projection of the point forces, response ~ 1/B at every node)
     !! instead of the element-lumped one (bounded node response)
     logical :: adjoint_mass = .false.
     !> Inverse of the element-lumped mass (element volume over node
     !! count, averaged over the copies of a shared node), the weight of
     !! the inner product the adjoint spread is taken in
     type(field_t) :: winv
     !> Lumped marker weights D of the adjoint spread, + and - side
     real(kind=rp), allocatable :: dp_mark(:)
     real(kind=rp), allocatable :: dm_mark(:)
     !> Constant-mode gain per marker, (G D 1)_m, + and - side (diagnostic)
     real(kind=rp), allocatable :: gp_mark(:)
     real(kind=rp), allocatable :: gm_mark(:)
     !> One-sided interpolation weight sums c_m = I(swp 1), I(swm 1) and the
     !! per-marker normalisation s_m = 1 / max(c_m, cmin) (0 for c_m <= 0)
     real(kind=rp), allocatable :: cp_mark(:), cm_mark(:)
     real(kind=rp), allocatable :: sp_mark(:), sm_mark(:)
     !> Deposit weights D_m s_m handed to the transpose
     real(kind=rp), allocatable :: dsp_mark(:), dsm_mark(:)
     type(c_ptr) :: sp_mark_d = C_NULL_PTR
     type(c_ptr) :: sm_mark_d = C_NULL_PTR
     type(c_ptr) :: dsp_mark_d = C_NULL_PTR
     type(c_ptr) :: dsm_mark_d = C_NULL_PTR
     !> Side weights of the adjoint spread, pmsk/(pmsk+mmsk) and
     !! mmsk/(pmsk+mmsk): a surface node shared by both sides counts half
     !! on each
     type(field_t) :: swp
     type(field_t) :: swm
     !> Half-width of the band of nodes on the surface that belong to both
     !! sides, in local grid spacings
     real(kind=rp) :: mask_band = 0.0_rp
     !> Floor of the one-sided weight sum below which the normalisation is
     !! clamped
     real(kind=rp) :: cmin = 0.25_rp
     !> Gain factor on the lumped weights of the adjoint spread
     real(kind=rp) :: spread_gain = 1.0_rp
     !> Scaled marker values -D u_m / dt handed to the transpose
     real(kind=rp), allocatable :: vals_adj(:)
     type(c_ptr) :: dp_mark_d = C_NULL_PTR
     type(c_ptr) :: dm_mark_d = C_NULL_PTR
     type(c_ptr) :: vals_adj_d = C_NULL_PTR
     !> Assembled mass matrix (1/Binv): the weight of the mass inner product
     !! in which the adjoint interpolation is the transpose of the spread
     real(kind=rp), allocatable :: Bm(:,:,:,:)
     !> Cutoff radius (in units of ds) for the Shepard interpolation
     real(kind=rp) :: interp_rmax
     !> Reduction slot per marker for the Shepard interpolation: markers
     !! held by several ranks share a slot, 0 marks a single-holder marker
     integer, allocatable :: shared_slot(:)
     !> Global number of markers held by more than one rank
     integer :: n_shared_glb = 0
     class(filter_t), allocatable :: fltr
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => idw_source_term_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => idw_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => idw_source_term_compute
     !> Initialise lagrangian from a boundary mesh
     procedure, pass(this) :: init_boundary_mesh => idw_init_boundary_mesh
     !> Initialise the csv output of the force on the immersed objects
     procedure, pass(this) :: init_force_output => idw_init_force_output
     !> Compute and write the force on the immersed objects
     procedure, pass(this) :: write_force => idw_write_force
  end type idw_source_term_t

  public :: idw_interp_shepard, idw_interp_shepard_partials, &
       idw_interp_shepard_normalize, idw_build_shared_slots, inv_dist_weight

contains

  subroutine idw_source_term_init_from_json(this, json, fields, &
       coef, variable_name)
    class(idw_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), target, intent(in) :: coef
    character(len=*), intent(in) :: variable_name
    type(json_value), pointer :: json_object_list
    type(json_core) :: core
    type(json_file) :: object_settings
    character(len=:), allocatable :: object_type
    integer :: n_regions, i, j, k, e, n_lags
    integer :: msk_zeros(2)
    logical :: marker_output
    integer :: np_min, nm_min, np_sum, nm_sum, np_empty, nm_empty, n_lag_glb
    integer :: gid_offset
    integer, allocatable :: np_i(:), nm_i(:), cnt_shr(:,:)
    character(len=LOG_SIZE) :: log_buf
    real(kind=dp) :: aabb_padding,dx_max, dy_max, dz_max, ds_max, ds_min
    real(kind=dp) :: diam_min, dxe, dye, dze
    type(stack_i4_t) :: overlaps
    type(stack_pt_t) :: lagrangian_points
    type(stack_pt_t) :: lagrangian_normals
    type(stack_i4_t) :: lagrangian_gids
    real(kind=rp) :: start_time, end_time
    character(len=:), allocatable :: filter_type
    character(len=:), allocatable :: interp_scheme
    character(len=:), allocatable :: spread_scheme
    type(json_file) :: filter_subdict

    ! Mandatory fields for the general source term
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    ! The forcing feeds the velocity back with gain 1/dt. Extrapolating it
    ! in time with the other explicit terms is unstable (BDF3/EXT3 allows a
    ! gain of at most 20/21), so the scheme applies it as computed.
    this%extrapolate = .false.

    call neko_log%section('Inverse distance weighting')

    call json_get_or_default(json, "rmax", this%rmax, 1.0_rp)
    write(log_buf, '(A,f5.2)') 'Rmax       : ', this%rmax
    call neko_log%message(log_buf)

    call json_get_or_default(json, "padding", aabb_padding, 0.125_dp)
    write(log_buf, '(A,f5.2)') 'Padding    : ', aabb_padding
    call neko_log%message(log_buf)
    
    call json_get_or_default(json, "power_parameter", this%pwr_param, 0.5_rp)
    write(log_buf, '(A,f5.2)') 'IDW Power  : ', this%pwr_param
    call neko_log%message(log_buf)
    call json_get_or_default(json, "one_sided", this%one_sided, .true.)
    write(log_buf, '(A,L1)') 'One sided  : ', this%one_sided
    call neko_log%message(log_buf)

    call json_get_or_default(json, "interpolation", interp_scheme, &
         "spectral")
    select case (trim(interp_scheme))
    case ("spectral", "barycentric")
       ! Point evaluation of the element's polynomial interpolant at the
       ! marker; "barycentric" is the historical name, after the formula
       ! used to evaluate the Lagrange basis
       this%idw_interp = .false.
    case ("idw")
       this%idw_interp = .true.
    case ("adjoint")
       ! Shepard interpolation with the spread's own stencil and weights
       ! K B mult / w: the interpolation is then the mass-weighted adjoint
       ! of the spread, interp o spread is symmetric positive semi-definite
       ! and the forcing can only remove kinetic energy
       this%idw_interp = .true.
       this%adjoint_interp = .true.
    case default
       call neko_error('IDW source term unknown interpolation scheme: ' &
            // trim(interp_scheme))
    end select
    call neko_log%message('Interp     : '// trim(interp_scheme))

    call json_get_or_default(json, "spread", spread_scheme, "idw")
    select case (trim(spread_scheme))
    case ("idw")
       this%adjoint_spread = .false.
    case ("adjoint", "adjoint_mass")
       ! The spread is the adjoint of the spectral interpolation, so that
       ! interp o spread is symmetric positive semi-definite and the
       ! forcing dissipates the energy of the inner product it is taken
       ! in. "adjoint" uses the element-lumped mass (element volume over
       ! node count): constant within an element, so the node response
       ! has no 1/B amplification at the low-mass GLL corner nodes, and
       ! it follows the physical mass across elements, so the dissipated
       ! quantity is the kinetic energy up to the GLL weight ratios.
       ! "adjoint_mass" uses the pointwise GLL mass, the L2 projection of
       ! the point forces: variationally exact, but a corner node of a
       ! cut element, far out in the fluid, is kicked by 1/B.
       if (this%idw_interp) then
          call neko_error('IDW source term: the adjoint spread pairs with &
               &the spectral interpolation')
       end if
       this%adjoint_spread = .true.
       this%adjoint_mass = trim(spread_scheme) .eq. 'adjoint_mass'
    case default
       call neko_error('IDW source term unknown spread scheme: ' &
            // trim(spread_scheme))
    end select
    call neko_log%message('Spread     : '// trim(spread_scheme))
    if (this%adjoint_spread) then
       call json_get_or_default(json, "spread_gain", this%spread_gain, 1.0_rp)
       write(log_buf, '(A,f5.2)') 'Spread gain: ', this%spread_gain
       call neko_log%message(log_buf)
       call json_get_or_default(json, "one_sided_min_weight", this%cmin, &
            0.25_rp)
       write(log_buf, '(A,f5.2)') 'Min weight : ', this%cmin
       call neko_log%message(log_buf)
    end if
    ! Nodes within mask_band grid spacings of the surface belong to both
    ! sides. At an element face the Lagrange basis is supported on the
    ! face alone, so a face-aligned wall seen one-sidedly has no footprint
    ! on the fluid side; sharing the surface nodes restores it. Default on
    ! for the adjoint spread only, the kernel spread keeps its masks.
    if (this%adjoint_spread) then
       call json_get_or_default(json, "mask_band", this%mask_band, 0.1_rp)
    else
       call json_get_or_default(json, "mask_band", this%mask_band, 0.0_rp)
    end if
    if (this%one_sided) then
       write(log_buf, '(A,f5.2)') 'Mask band  : ', this%mask_band
       call neko_log%message(log_buf)
    end if

    call json_get_or_default(json, "interpolation_rmax", this%interp_rmax, &
         2.0_rp)
    ! The adjoint pairing needs the stencil of the spread
    if (this%adjoint_interp) this%interp_rmax = this%rmax
    if (this%idw_interp) then
       write(log_buf, '(A,f5.2)') 'Interp rmax: ', this%interp_rmax
       call neko_log%message(log_buf)
    end if

    ! The mass inner product sums the local B over the copies of a dof;
    ! Binv inverts exactly that sum, so 1/Binv is the assembled mass on
    ! every copy
    allocate(this%Bm(coef%Xh%lx, coef%Xh%ly, coef%Xh%lz, coef%msh%nelv))
    this%Bm = 1.0_rp / coef%Binv

    ! Element-lumped mass for the adjoint spread: the local mass of an
    ! element spread evenly over its nodes, then averaged over the copies
    ! of a shared node so that the weight is the same on every copy
    if (this%adjoint_spread) then
       call this%winv%init(coef%dof, "ib_winv")
       do e = 1, coef%msh%nelv
          this%winv%x(:,:,:,e) = sum(coef%B(:,:,:,e)) / real(coef%Xh%lxyz, rp)
       end do
    end if


    call json_get_or_default(json, 'filter.type', filter_type, 'none')
    select case (filter_type)
    case ('PDE')
       allocate(PDE_filter_t::this%fltr)       
    case ('elementwise')
       allocate(elementwise_filter_t::this%fltr)
    case ('none')       
    case default
       call neko_error('IDW source term unknown filter type')
    end select
    
    if (allocated(this%fltr)) then
       call json_get(json, 'filter', filter_subdict)
        call this%fltr%init(filter_subdict, coef)
     end if

    if (json%valid_path('force_output')) then
       call this%init_force_output(json, variable_name)
    end if


    ! Naive apporach to find the smallest distance between two dofs in the mesh

    ds_min = huge(0.0_rp)
    ds_max = -huge(0.0_rp)

    call this%ds%init(coef%dof)

    associate (x => coef%dof%x%x, y => coef%dof%y%x, z => coef%dof%z%x, &
         lx => coef%Xh%lx, ds => this%ds%x)

      !$omp parallel do private(i, j, k, dx_max, dy_max, dz_max) &
      !$omp reduction(max:ds_max) reduction(min:ds_min)
      do e = 1, coef%msh%nelv
         do k = 2, lx-1
            do j = 2, lx-1
               do i = 2, lx-1
                  dx_max = max(abs(x(i,j,k,e) - x(i+1,j,k,e)), &
                       abs(x(i-1,j,k,e) - x(i,j,k,e)), &
                       abs(x(i,j,k,e) - x(i,j+1,k,e)), &
                       abs(x(i,j-1,k,e) - x(i,j,k,e)), &
                       abs(x(i,j,k,e) - x(i,j,k+1,e)), &
                       abs(x(i,j,k-1,e) - x(i,j,k,e)), &
                       abs(x(i-1, j-1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j-1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j+1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j+1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j-1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j-1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j+1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j+1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i, j+1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i, j+1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i, j-1, k+1, e) - x(i,j,k,e)), &
                       abs(x(i, j-1, k-1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j, k+1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j, k-1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j, k+1, e) - x(i,j,k,e)), &
                       abs(x(i+1, j, k-1, e) - x(i,j,k,e)), &
                       abs(x(i-1, j-1, k, e) - x(i,j,k,e)), &
                       abs(x(i-1, j+1, k, e) - x(i,j,k,e)), &
                       abs(x(i+1, j-1, k, e) - x(i,j,k,e)), &
                       abs(x(i+1, j+1, k, e) - x(i,j,k,e)))

                  dy_max = max(abs(y(i,j,k,e) - y(i+1,j,k,e)), &
                       abs(y(i-1,j,k,e) - y(i,j,k,e)), &
                       abs(y(i,j,k,e) - y(i,j+1,k,e)), &
                       abs(y(i,j-1,k,e) - y(i,j,k,e)), &
                       abs(y(i,j,k,e) - y(i,j,k+1,e)), &
                       abs(y(i,j,k-1,e) - y(i,j,k,e)), &
                       abs(y(i-1, j-1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j-1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j+1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j+1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j-1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j-1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j+1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j+1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i, j+1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i, j+1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i, j-1, k+1, e) - y(i,j,k,e)), &
                       abs(y(i, j-1, k-1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j, k+1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j, k-1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j, k+1, e) - y(i,j,k,e)), &
                       abs(y(i+1, j, k-1, e) - y(i,j,k,e)), &
                       abs(y(i-1, j-1, k, e) - y(i,j,k,e)), &
                       abs(y(i-1, j+1, k, e) - y(i,j,k,e)), &
                       abs(y(i+1, j-1, k, e) - y(i,j,k,e)), &
                       abs(y(i+1, j+1, k, e) - y(i,j,k,e)))


                  dz_max = max(abs(z(i,j,k,e) - z(i+1,j,k,e)), &
                       abs(z(i-1,j,k,e) - z(i,j,k,e)), &
                       abs(z(i,j,k,e) - z(i,j+1,k,e)), &
                       abs(z(i,j-1,k,e) - z(i,j,k,e)), &
                       abs(z(i,j,k,e) - z(i,j,k+1,e)), &
                       abs(z(i,j,k-1,e) - z(i,j,k,e)), &
                       abs(z(i-1, j-1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j-1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j+1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j+1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j-1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j-1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j+1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j+1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i, j+1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i, j+1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i, j-1, k+1, e) - z(i,j,k,e)), &
                       abs(z(i, j-1, k-1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j, k+1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j, k-1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j, k+1, e) - z(i,j,k,e)), &
                       abs(z(i+1, j, k-1, e) - z(i,j,k,e)), &
                       abs(z(i-1, j-1, k, e) - z(i,j,k,e)), &
                       abs(z(i-1, j+1, k, e) - z(i,j,k,e)), &
                       abs(z(i+1, j-1, k, e) - z(i,j,k,e)), &
                       abs(z(i+1, j+1, k, e) - z(i,j,k,e)))
                  ds(i,j,k,e) = (dx_max + dy_max + dz_max) / 3.0_rp
               end do
            end do
         end do


         ds(1,:,:,e) = ds(2,:,:,e)
         ds(lx,:,:,e) = ds(lx-1,:,:,e)

         ds(:,1,:,e) = ds(:,2,:,e)
         ds(:,lx,:,e) = ds(:,lx-1,:,e)

         ds(:,:,1,e) = ds(:,:,2,e)
         ds(:,:,lx,e) = ds(:,:,lx-1,e)


         ds_max = max(ds_max, maxval(ds(:,:,:,e)))
         ds_min = min(ds_min, minval(ds(:,:,:,e)))

      end do
      !$omp end parallel do
    end associate


    this%ds_min = ds_min
    this%ds_max = ds_max

    call MPI_Allreduce(MPI_IN_PLACE, this%ds_min, 1, &
         MPI_DOUBLE_PRECISION, MPI_MIN, NEKO_COMM)
    write(log_buf, '(A,ES13.6)') 'Minimum ds :', this%ds_min
    call neko_log%message(log_buf)

    call MPI_Allreduce(MPI_IN_PLACE, this%ds_max, 1, &
         MPI_DOUBLE_PRECISION, MPI_MAX, NEKO_COMM)
    write(log_buf, '(A,ES13.6)') 'Maximum ds :', this%ds_max
    call neko_log%message(log_buf)


    call this%intersect%init(coef%msh, aabb_padding)
    call lagrangian_points%init()
    call lagrangian_normals%init()
    call lagrangian_gids%init()
    gid_offset = 0

    call json%get('objects', json_object_list)
    call json%info('objects', n_children=n_regions)
    call json%get_core(core)

    if (n_regions .lt. 10) then
       write(log_buf, '(A, I1)') 'Objects    : ', n_regions
    else if (n_regions .ge. 100) then
       write(log_buf, '(A, I2)') 'Objects    : ', n_regions
    else
       write(log_buf, '(A, I3)') 'Objects    : ', n_regions
    end if
    call neko_log%message(log_buf)

    call neko_log%begin()
    do i = 1, n_regions
       call neko_log%begin()
       call json_extract_item(core, json_object_list, i , object_settings)
       call json_get_or_default(object_settings, 'type', object_type, 'none')

       call neko_log%message('Type       : '// trim(object_type))
       select case (object_type)
       case ('boundary_mesh')
          call this%init_boundary_mesh(lagrangian_points, lagrangian_normals, &
               lagrangian_gids, gid_offset, object_settings)
       case ('none')
          call neko_error('IDW source term objects require a region type')
       case default
          call neko_error('IDW source term unkown region type')
       end select
       call neko_log%end()
    end do
    call neko_log%end()

    ! Report total number of lagrangian points generated, this differs
    ! from the number of triangles due to refinement and clipping against
    ! the fluid mesh
    n_lags = lagrangian_points%size()
    call MPI_Allreduce(MPI_IN_PLACE, n_lags, 1, MPI_INTEGER, MPI_SUM, NEKO_COMM)
    if (n_lags .lt. 1e1) then
       call neko_log%message('Type       : '// trim(object_type))
       write(log_buf, '(A, I1)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e2) then
       write(log_buf, '(A, I2)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e3) then
       write(log_buf, '(A, I3)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e4) then
       write(log_buf, '(A, I4)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e5) then
       write(log_buf, '(A, I5)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e6) then
       write(log_buf, '(A, I6)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e7) then
       write(log_buf, '(A, I7)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e8) then
       write(log_buf, '(A, I8)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e9) then
       write(log_buf, '(A, I9)') 'Tot lagpts : ', n_lags
    else if (n_lags .lt. 1e10) then
       write(log_buf, '(A, I10)') 'Tot lagpts : ', n_lags
    end if
    call neko_log%message(log_buf)

    allocate(this%xyz(3, lagrangian_points%size()))
    allocate(this%fu_ib(lagrangian_points%size()))
    allocate(this%fv_ib(lagrangian_points%size()))
    allocate(this%fw_ib(lagrangian_points%size()))
    allocate(this%fum_ib(lagrangian_points%size()))
    allocate(this%fvm_ib(lagrangian_points%size()))
    allocate(this%fwm_ib(lagrangian_points%size()))
    allocate(this%lag_pts(lagrangian_points%size()))
    allocate(this%lag_nrm(lagrangian_normals%size()))

    select type(pt => lagrangian_points%data)
      type is (point_t)
       do i = 1, lagrangian_points%size()
          this%xyz(1, i) = pt(i)%x(1)
          this%xyz(2, i) = pt(i)%x(2)
          this%xyz(3, i) = pt(i)%x(3)
          this%lag_pts(i) = pt(i)
       end do
    end select

    select type(pt => lagrangian_normals%data)
      type is (point_t)
       do i = 1, lagrangian_normals%size()
          this%lag_nrm(i) = pt(i)
       end do
    end select

    n_lags = lagrangian_points%size()

    ! The Shepard (idw) interpolation never evaluates through the global
    ! (rst-based) interpolation, so the global search is only set up for
    ! the barycentric scheme. The Shepard scheme instead needs a reduction
    ! slot for every marker held by more than one rank, so that the partial
    ! sums of stencils straddling an MPI boundary can be summed across the
    ! holder ranks.
    if (.not. this%idw_interp) then
       call this%global_interp%init(coef%dof, NEKO_COMM)
       call this%global_interp%find_points_xyz(this%xyz, n_lags)
    else
       allocate(this%shared_slot(n_lags))
       this%shared_slot = 0
       select type (gd => lagrangian_gids%data)
         type is (integer)
          call idw_build_shared_slots(gd(1:n_lags), gid_offset, &
               this%shared_slot, this%n_shared_glb)
       end select
       write(log_buf, '(A,I9)') 'Shared mrks: ', this%n_shared_glb
       call neko_log%message(log_buf)
    end if

    ! Construct list of overlapping elements for each lagrangian particle
    allocate(this%lag_el(lagrangian_points%size()))
    call overlaps%init()
    do i = 1, size(this%lag_pts)
       call this%lag_el(i)%init()
       call this%intersect%overlap(this%lag_pts(i), overlaps)
       do while(.not. overlaps%is_empty())
          e = overlaps%pop()
          call this%lag_el(i)%push(e)
       end do
    end do
    call overlaps%free()

    call idw_build_el_csr(this, coef%msh%nelv)

    ! Construct weight field
    call this%w%init(coef%dof, "ib_weight")
    call this%wm%init(coef%dof, "ib_mweight")

    call this%gs%init(coef%dof)

    call idw_assemble(this%gs, this%ds, coef%mult)
    call idw_assemble(this%gs, this%ds, coef%mult)

    call this%mmsk%init(coef%dof, "ib_mmask")
    call this%pmsk%init(coef%dof, "ib_pmask")
    call this%tmp%init(coef%dof, "ib_tmp")
    call this%ib_fx%init(coef%dof, "ib_fx")
    call this%ib_fy%init(coef%dof, "ib_fy")
    call this%ib_fz%init(coef%dof, "ib_fz")

    if (this%one_sided) then
       call idw_compute_mask(this%mmsk, this%pmsk, this%lag_pts, &
            this%lag_nrm, this%active_el, this%el_off, this%el_lag, &
            coef%dof%x%x, coef%dof%y%x, coef%dof%z%x, this%ds%x, &
            this%mask_band, coef%Xh%lx, coef%msh%nelv)
    else
       ! Assign the host arrays directly: the field_t defined assignment only
       ! fills the device buffer on a device build, and idw_assemble below and
       ! idw_compute_weight both read the host side.
       this%mmsk%x = 0.0_rp
       this%pmsk%x = 1.0_rp
    end if

    call idw_assemble(this%gs, this%pmsk, coef%mult)
    call idw_assemble(this%gs, this%mmsk, coef%mult)

    if (this%adjoint_spread) then
       call idw_assemble(this%gs, this%winv, coef%mult)
       this%winv%x = 1.0_rp / this%winv%x
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%winv%x, this%winv%x_d, this%winv%size(), &
               HOST_TO_DEVICE, sync = .true.)
       end if
    end if

    ! Masked dofs per side (reduced over ranks, copies counted once per
    ! rank): both counts zero means the side split is not in effect
    if (this%one_sided) then
       msk_zeros(1) = count(this%pmsk%x .lt. 0.5_rp)
       msk_zeros(2) = count(this%mmsk%x .lt. 0.5_rp)
       call MPI_Allreduce(MPI_IN_PLACE, msk_zeros, 2, MPI_INTEGER, MPI_SUM, &
            NEKO_COMM)
       write(log_buf, '(A,I0,A,I0)') 'Mask zeros : + side ', msk_zeros(1), &
            ', - side ', msk_zeros(2)
       call neko_log%message(log_buf)
    end if

    call idw_compute_weight(this%w, this%wm, this%pmsk, this%lag_pts, &
         this%active_el, this%el_off, this%el_lag, &
         coef%dof%x%x, coef%dof%y%x, coef%dof%z%x, this%ds%x, this%rmax, &
         this%pwr_param, coef%Xh%lx,coef%msh%nelv)

    call idw_assemble(this%gs, this%w, coef%mult)
    call idw_assemble(this%gs, this%wm, coef%mult)

    if (this%adjoint_spread) call idw_adjoint_spread_weights(this)

    ! Optional dump of the markers this rank holds, for inspection next to
    ! the solution: position, normal and the lumped weights of the adjoint
    ! spread (zero for the kernel spread). One csv per rank.
    call json_get_or_default(json, 'marker_output', marker_output, .false.)
    if (marker_output) call idw_write_markers(this)

    if (this%idw_interp) then
       ! Per-marker local stencil counts; the counts of markers held by
       ! several ranks are then summed over their holders so the diagnostic
       ! reflects the stencil the reduced interpolation actually sees.
       allocate(np_i(n_lags), nm_i(n_lags))
       call idw_interp_stencil_counts(this%lag_pts, this%lag_el, &
            this%pmsk%x, coef%dof%x%x, coef%dof%y%x, coef%dof%z%x, this%ds%x, &
            this%interp_rmax, coef%Xh%lx, coef%msh%nelv, np_i, nm_i)

       if (this%n_shared_glb .gt. 0) then
          allocate(cnt_shr(2, this%n_shared_glb))
          cnt_shr = 0
          do i = 1, n_lags
             if (this%shared_slot(i) .gt. 0) then
                cnt_shr(1, this%shared_slot(i)) = np_i(i)
                cnt_shr(2, this%shared_slot(i)) = nm_i(i)
             end if
          end do
          call MPI_Allreduce(MPI_IN_PLACE, cnt_shr, 2 * this%n_shared_glb, &
               MPI_INTEGER, MPI_SUM, NEKO_COMM)
       end if

       ! Stats over single-holder markers, reduced over ranks
       np_min = huge(0)
       nm_min = huge(0)
       np_sum = 0
       nm_sum = 0
       np_empty = 0
       nm_empty = 0
       n_lag_glb = 0
       do i = 1, n_lags
          if (this%shared_slot(i) .eq. 0) then
             np_min = min(np_min, np_i(i))
             nm_min = min(nm_min, nm_i(i))
             np_sum = np_sum + np_i(i)
             nm_sum = nm_sum + nm_i(i)
             if (np_i(i) .eq. 0) np_empty = np_empty + 1
             if (nm_i(i) .eq. 0) nm_empty = nm_empty + 1
             n_lag_glb = n_lag_glb + 1
          end if
       end do

       call MPI_Allreduce(MPI_IN_PLACE, np_min, 1, MPI_INTEGER, &
            MPI_MIN, NEKO_COMM)
       call MPI_Allreduce(MPI_IN_PLACE, nm_min, 1, MPI_INTEGER, &
            MPI_MIN, NEKO_COMM)
       call MPI_Allreduce(MPI_IN_PLACE, np_sum, 1, MPI_INTEGER, &
            MPI_SUM, NEKO_COMM)
       call MPI_Allreduce(MPI_IN_PLACE, nm_sum, 1, MPI_INTEGER, &
            MPI_SUM, NEKO_COMM)
       call MPI_Allreduce(MPI_IN_PLACE, np_empty, 1, MPI_INTEGER, &
            MPI_SUM, NEKO_COMM)
       call MPI_Allreduce(MPI_IN_PLACE, nm_empty, 1, MPI_INTEGER, &
            MPI_SUM, NEKO_COMM)
       call MPI_Allreduce(MPI_IN_PLACE, n_lag_glb, 1, MPI_INTEGER, &
            MPI_SUM, NEKO_COMM)

       ! Fold in each shared marker once; the reduced buffer is identical
       ! on every rank, so the stats stay rank independent
       do i = 1, this%n_shared_glb
          np_min = min(np_min, cnt_shr(1, i))
          nm_min = min(nm_min, cnt_shr(2, i))
          np_sum = np_sum + cnt_shr(1, i)
          nm_sum = nm_sum + cnt_shr(2, i)
          if (cnt_shr(1, i) .eq. 0) np_empty = np_empty + 1
          if (cnt_shr(2, i) .eq. 0) nm_empty = nm_empty + 1
       end do
       n_lag_glb = n_lag_glb + this%n_shared_glb

       if (allocated(cnt_shr)) deallocate(cnt_shr)
       deallocate(np_i, nm_i)

       write(log_buf, '(A,I6,A,F7.1)') 'Interp + side: min nodes ', &
            np_min, ', mean ', real(np_sum, rp) / max(n_lag_glb, 1)
       call neko_log%message(log_buf)
       write(log_buf, '(A,I6)') 'Interp + side: empty markers ', np_empty
       call neko_log%message(log_buf)
       if (this%one_sided) then
          write(log_buf, '(A,I6,A,F7.1)') 'Interp - side: min nodes ', &
               nm_min, ', mean ', real(nm_sum, rp) / max(n_lag_glb, 1)
          call neko_log%message(log_buf)
          write(log_buf, '(A,I6)') 'Interp - side: empty markers ', nm_empty
          call neko_log%message(log_buf)
       end if

       diam_min = huge(0.0_dp)
       do e = 1, coef%msh%nelv
          dxe = maxval(coef%dof%x%x(:,:,:,e)) - minval(coef%dof%x%x(:,:,:,e))
          dye = maxval(coef%dof%y%x(:,:,:,e)) - minval(coef%dof%y%x(:,:,:,e))
          dze = maxval(coef%dof%z%x(:,:,:,e)) - minval(coef%dof%z%x(:,:,:,e))
          diam_min = min(diam_min, sqrt(dxe**2 + dye**2 + dze**2))
       end do
       call MPI_Allreduce(MPI_IN_PLACE, diam_min, 1, &
            MPI_DOUBLE_PRECISION, MPI_MIN, NEKO_COMM)

       if (this%interp_rmax * this%ds_max .gt. aabb_padding * diam_min) then
          call neko_log%warning('interpolation_rmax*ds may exceed the &
               &element search reach; consider increasing padding')
       end if
    end if

    call lagrangian_points%free()
    call lagrangian_normals%free()
    call lagrangian_gids%free()

    if (NEKO_BCKND_DEVICE .eq. 1) call idw_build_device_maps(this)

    call neko_log%end_section()

  end subroutine idw_source_term_init_from_json

  !> Assemble a host-built field to a continuous representation,
  subroutine idw_assemble(gs, fld, mult)
    type(gs_t), intent(inout) :: gs
    type(field_t), intent(inout) :: fld
    real(kind=rp), dimension(:,:,:,:), intent(in) :: mult
    integer :: n

    n = fld%size()
    fld%x = fld%x * mult

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(fld%x, fld%x_d, n, HOST_TO_DEVICE, sync = .true.)
       call gs%op(fld, GS_OP_ADD)
       call device_memcpy(fld%x, fld%x_d, n, DEVICE_TO_HOST, sync = .true.)
    else
       call gs%op(fld, GS_OP_ADD)
    end if

  end subroutine idw_assemble

  !> Write the markers held by this rank to `idw_markers_<rank>.csv`:
  !! x, y, z, nx, ny, nz, d_plus, d_minus, gain_plus, gain_minus (the
  !! lumped weights of the adjoint spread and its constant-mode gain per
  !! marker, zero for the kernel spread). Markers held by several ranks
  !! appear in every holder's file.
  subroutine idw_write_markers(this)
    class(idw_source_term_t), intent(inout) :: this
    character(len=64) :: fname
    real(kind=rp) :: dp_m, dm_m, gp_m, gm_m, cp_m, cm_m
    integer :: i, unit, n_lag

    n_lag = size(this%lag_pts)
    write(fname, '(A,I0,A)') 'idw_markers_', pe_rank, '.csv'
    open(newunit = unit, file = trim(fname), status = 'replace', &
         action = 'write')
    write(unit, '(A)') 'x,y,z,nx,ny,nz,d_plus,d_minus,gain_plus,' // &
         'gain_minus,c_plus,c_minus'
    do i = 1, n_lag
       dp_m = 0.0_rp
       dm_m = 0.0_rp
       gp_m = 0.0_rp
       gm_m = 0.0_rp
       cp_m = 0.0_rp
       cm_m = 0.0_rp
       if (allocated(this%dp_mark)) dp_m = this%dp_mark(i)
       if (allocated(this%dm_mark)) dm_m = this%dm_mark(i)
       if (allocated(this%gp_mark)) gp_m = this%gp_mark(i)
       if (allocated(this%gm_mark)) gm_m = this%gm_mark(i)
       if (allocated(this%cp_mark)) cp_m = this%cp_mark(i)
       if (allocated(this%cm_mark)) cm_m = this%cm_mark(i)
       write(unit, '(11(ES16.8,","),ES16.8)') this%lag_pts(i)%x(1), &
            this%lag_pts(i)%x(2), this%lag_pts(i)%x(3), &
            this%lag_nrm(i)%x(1), this%lag_nrm(i)%x(2), &
            this%lag_nrm(i)%x(3), dp_m, dm_m, gp_m, gm_m, cp_m, cm_m
    end do
    close(unit)

  end subroutine idw_write_markers

  !> Gather-scatter sum of a field without the multiplicity scaling, on
  !! the host array (the device mirror is used for the exchange on device
  !! builds)
  subroutine idw_assemble_sum(gs, fld)
    type(gs_t), intent(inout) :: gs
    type(field_t), intent(inout) :: fld
    integer :: n

    n = fld%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(fld%x, fld%x_d, n, HOST_TO_DEVICE, sync = .true.)
       call gs%op(fld, GS_OP_ADD)
       call device_memcpy(fld%x, fld%x_d, n, DEVICE_TO_HOST, sync = .true.)
    else
       call gs%op(fld, GS_OP_ADD)
    end if

  end subroutine idw_assemble_sum

  !> One side of the adjoint spectral spread: `out += msk * W^-1 *
  !! gs_add(I^T vals)`, the adjoint of the masked interpolation
  !! `I (msk * u)` in the inner product with weight W, so that for every
  !! field u: sum_m vals(m) * I(msk u)(m) = sum_i W(i) u(i) out(i) (out
  !! zero on entry, dofs counted once). W is the element-lumped mass
  !! (`winv`), or the assembled GLL mass (`Binv`) with `adjoint_mass`.
  !! Without `msk` the unmasked operator.
  !! With `on_host` false on a device build everything runs on the device
  !! (`vals` device-mapped); the point values still travel through the
  !! host inside the transpose.
  subroutine idw_adjoint_apply(this, vals, out, on_host, msk)
    class(idw_source_term_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vals(:)
    type(field_t), intent(inout) :: out
    logical, intent(in) :: on_host
    type(field_t), intent(in), optional :: msk
    integer :: n

    n = this%tmp%size()
    if (NEKO_BCKND_DEVICE .eq. 1 .and. .not. on_host) then
       call device_rzero(this%tmp%x_d, n)
       call this%global_interp%evaluate_transpose(vals, this%tmp%x, .false.)
       call this%gs%op(this%tmp, GS_OP_ADD)
       if (this%adjoint_mass) then
          call device_col2(this%tmp%x_d, this%coef%Binv_d, n)
       else
          call device_col2(this%tmp%x_d, this%winv%x_d, n)
       end if
       if (present(msk)) call device_col2(this%tmp%x_d, msk%x_d, n)
       call device_add2(out%x_d, this%tmp%x_d, n)
    else
       this%tmp%x = 0.0_rp
       call this%global_interp%evaluate_transpose(vals, this%tmp%x, .true.)
       call idw_assemble_sum(this%gs, this%tmp)
       if (this%adjoint_mass) then
          this%tmp%x = this%coef%Binv * this%tmp%x
       else
          this%tmp%x = this%winv%x * this%tmp%x
       end if
       if (present(msk)) then
          out%x = out%x + msk%x * this%tmp%x
       else
          out%x = out%x + this%tmp%x
       end if
    end if

  end subroutine idw_adjoint_apply

  !> Lumped marker weights of the adjoint spectral spread. With
  !! G = I Binv I^T the marker Gram matrix of one side (I the masked
  !! interpolation), the weights are D = gain / (|G| 1): by Gershgorin
  !! every eigenvalue of the gain operator interp o spread = G D then lies
  !! in [0, gain], whatever the marker density. Lumping by the signed row
  !! sums (D = 1/(G 1)) would give gain exactly one on a constant marker
  !! velocity, but G is the Christoffel-Darboux kernel of the polynomial
  !! space and oscillates: markers one node spacing apart couple with
  !! near-zero or negative weight, the signed row sum under-counts the
  !! diagonal and the oscillatory marker modes get gains of 2-4 and more
  !! (measured 1.4-4.2 on a flat sheet in one degree-7 element), past the
  !! stability limit of the time integration. The price of the absolute
  !! lumping is a constant-mode gain below one, reported in the log; the
  !! `spread_gain` factor scales it back up as long as it stays under the
  !! limit of the scheme.
  !! Markers held by more than one rank are deposited once per holder; the
  !! row sum sees the duplicates too, so the lumping absorbs them exactly.
  !! Markers whose absolute row sum is zero (no unmasked node reaches them)
  !! get zero weight.
  subroutine idw_adjoint_spread_weights(this)
    class(idw_source_term_t), intent(inout) :: this
    real(kind=rp), allocatable :: ones(:), rowsum(:)
    integer :: n_lag, nv, n

    n_lag = size(this%lag_pts)
    nv = max(n_lag, 1)
    n = this%tmp%size()
    allocate(this%dp_mark(n_lag), this%dm_mark(n_lag))
    allocate(this%gp_mark(n_lag), this%gm_mark(n_lag))
    allocate(this%cp_mark(n_lag), this%cm_mark(n_lag))
    allocate(this%sp_mark(nv), this%sm_mark(nv))
    allocate(this%dsp_mark(nv), this%dsm_mark(nv))
    allocate(ones(n_lag), rowsum(n_lag))
    ones = 1.0_rp
    this%dp_mark = 0.0_rp
    this%dm_mark = 0.0_rp
    this%gp_mark = 0.0_rp
    this%gm_mark = 0.0_rp
    this%cp_mark = 0.0_rp
    this%cm_mark = 0.0_rp
    this%sp_mark = 0.0_rp
    this%sm_mark = 0.0_rp
    this%dsp_mark = 0.0_rp
    this%dsm_mark = 0.0_rp

    ! Side weights: pmsk and mmsk are complementary 0/1 except on the
    ! surface band, where both are 1 and each side counts half. Host
    ! arithmetic, then the device mirror.
    call this%swp%init(this%coef%dof, "ib_swp")
    call this%swm%init(this%coef%dof, "ib_swm")
    if (this%one_sided) then
       this%swp%x = this%pmsk%x / (this%pmsk%x + this%mmsk%x)
       this%swm%x = this%mmsk%x / (this%pmsk%x + this%mmsk%x)
    else
       this%swp%x = 1.0_rp
       this%swm%x = 0.0_rp
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%swp%x, this%swp%x_d, n, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%swm%x, this%swm%x_d, n, HOST_TO_DEVICE, &
            sync = .true.)
    end if

    ! One-sided weight sums c_m = I(sw 1) and the normalisation s_m that
    ! makes the one-sided interpolation reproduce a constant on its side:
    ! s = 1/c for c >= cmin, clamped to 1/cmin for small positive c, and 0
    ! for c <= 0 (a marker whose side sum is negative would push a uniform
    ! flow the wrong way). The scaling is a per-marker factor on both the
    ! interpolation and the deposit, so the pair stays adjoint.
    this%tmp%x = this%swp%x
    call this%global_interp%evaluate(this%cp_mark, this%tmp%x, .true.)
    call idw_adjoint_normalisation(this%cp_mark, this%sp_mark, this%cmin, &
         'Adj. spread + side')
    if (this%one_sided) then
       this%tmp%x = this%swm%x
       call this%global_interp%evaluate(this%cm_mark, this%tmp%x, .true.)
       call idw_adjoint_normalisation(this%cm_mark, this%sm_mark, &
            this%cmin, 'Adj. spread - side')
    end if

    if (this%one_sided) then
       call idw_adjoint_lump(this, ones, rowsum, this%dp_mark, &
            'Adj. spread + side', this%sp_mark, this%swp)
       this%gp_mark = rowsum * this%dp_mark
       call idw_adjoint_lump(this, ones, rowsum, this%dm_mark, &
            'Adj. spread - side', this%sm_mark, this%swm)
       this%gm_mark = rowsum * this%dm_mark
    else
       call idw_adjoint_lump(this, ones, rowsum, this%dp_mark, &
            'Adj. spread', this%sp_mark)
       this%gp_mark = rowsum * this%dp_mark
    end if
    this%dsp_mark(1:n_lag) = this%dp_mark * this%sp_mark(1:n_lag)
    this%dsm_mark(1:n_lag) = this%dm_mark * this%sm_mark(1:n_lag)

    deallocate(ones, rowsum)

    ! Device mirrors: the interpolation scaling and the deposit weights
    ! are read by the per-step path, the scaled values are handed to the
    ! transpose as a mapped array
    allocate(this%vals_adj(nv))
    this%vals_adj = 0.0_rp
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%vals_adj, this%vals_adj_d, nv)
       call device_memcpy(this%vals_adj, this%vals_adj_d, nv, &
            HOST_TO_DEVICE, sync = .false.)
       call device_map(this%sp_mark, this%sp_mark_d, nv)
       call device_memcpy(this%sp_mark, this%sp_mark_d, nv, &
            HOST_TO_DEVICE, sync = .false.)
       call device_map(this%sm_mark, this%sm_mark_d, nv)
       call device_memcpy(this%sm_mark, this%sm_mark_d, nv, &
            HOST_TO_DEVICE, sync = .false.)
       call device_map(this%dsp_mark, this%dsp_mark_d, nv)
       call device_memcpy(this%dsp_mark, this%dsp_mark_d, nv, &
            HOST_TO_DEVICE, sync = .false.)
       call device_map(this%dsm_mark, this%dsm_mark_d, nv)
       call device_memcpy(this%dsm_mark, this%dsm_mark_d, nv, &
            HOST_TO_DEVICE, sync = .true.)
    end if

    if (this%one_sided) then
       call idw_adjoint_check_response(this, this%dp_mark, this%sp_mark, &
            'Adj. spread + side', this%swp)
       call idw_adjoint_check_response(this, this%dm_mark, this%sm_mark, &
            'Adj. spread - side', this%swm)
    else
       call idw_adjoint_check_response(this, this%dp_mark, this%sp_mark, &
            'Adj. spread')
    end if

  end subroutine idw_adjoint_spread_weights

  !> One forcing step applied to a unit velocity field, du = S D I(msk 1),
  !! on the host and, on a device build, through the per-step device path
  !! as well. The gain bound holds in the mass norm only, so on a graded
  !! mesh low-mass nodes can still overshoot: max|du| of order one is
  !! healthy, orders of magnitude above one blows a run up in its first
  !! step. The device result must match the host one to round-off; a
  !! difference points at the device data flow (deposit kernel, CSR
  !! mirrors, marker exchange), not at the method.
  subroutine idw_adjoint_check_response(this, d, scal, label, msk)
    class(idw_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: d(:), scal(:)
    character(len=*), intent(in) :: label
    type(field_t), intent(in), optional :: msk
    real(kind=rp), allocatable :: vals(:), du_host(:,:,:,:)
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: m(3)
    integer :: n_lag, n, nv

    n_lag = size(d)
    nv = max(n_lag, 1)
    n = this%tmp%size()
    allocate(vals(nv))

    ! Marker velocities of the (masked) unit field, host side
    if (present(msk)) then
       this%tmp%x = msk%x
    else
       this%tmp%x = 1.0_rp
    end if
    vals = 0.0_rp
    call this%global_interp%evaluate(vals, this%tmp%x, .true.)
    ! Normalised marker velocity s I(msk 1), then the deposit weight d s
    vals(1:n_lag) = d * scal(1:n_lag)**2 * vals(1:n_lag)

    this%ib_fx%x = 0.0_rp
    call idw_adjoint_apply(this, vals, this%ib_fx, .true., msk)
    m(1) = maxval(abs(this%ib_fx%x))
    m(2) = m(1)
    m(3) = 0.0_rp

    if (NEKO_BCKND_DEVICE .eq. 1) then
       du_host = this%ib_fx%x
       this%vals_adj(1:nv) = vals
       call device_memcpy(this%vals_adj, this%vals_adj_d, nv, &
            HOST_TO_DEVICE, sync = .true.)
       call device_rzero(this%ib_fx%x_d, n)
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fx, .false., msk)
       call device_memcpy(this%ib_fx%x, this%ib_fx%x_d, n, &
            DEVICE_TO_HOST, sync = .true.)
       m(2) = maxval(abs(this%ib_fx%x))
       m(3) = maxval(abs(this%ib_fx%x - du_host))
       call device_rzero(this%ib_fx%x_d, n)
       deallocate(du_host)
    end if
    this%ib_fx%x = 0.0_rp

    call MPI_Allreduce(MPI_IN_PLACE, m, 3, MPI_REAL_PRECISION, MPI_MAX, &
         NEKO_COMM)
    ! Two lines: the log buffer is LOG_SIZE characters
    write(log_buf, '(A,A,ES10.3)') trim(label), &
         ': unit response max|du| host ', m(1)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A,ES10.3,A,ES10.3)') trim(label), &
         ': unit response device ', m(2), ', max diff ', m(3)
    call neko_log%message(log_buf)

    ! Pointwise stability: a node whose velocity changes by more than
    ! twice its value in one step flips sign and grows, whatever the
    ! mass-norm spectrum says. The unit response scales with the gain, so
    ! the gain must stay below 2 * spread_gain / max|du|.
    if (m(1) .gt. 0.0_rp) then
       write(log_buf, '(A,A,F6.2)') trim(label), &
            ': node-stable gain bound ', 2.0_rp * this%spread_gain / m(1)
       call neko_log%message(log_buf)
       if (m(1) .gt. 2.0_rp) then
          call neko_log%warning('spread_gain exceeds the pointwise &
               &stability bound of the adjoint spread')
       end if
    end if

    deallocate(vals)

  end subroutine idw_adjoint_check_response

  !> Normalisation of the one-sided interpolation from its weight sums.
  !! Logs how many markers are clamped and how many switched off.
  subroutine idw_adjoint_normalisation(c, sc, cmin, label)
    real(kind=rp), intent(in) :: c(:)
    real(kind=rp), intent(inout) :: sc(:)
    real(kind=rp), intent(in) :: cmin
    character(len=*), intent(in) :: label
    character(len=LOG_SIZE) :: log_buf
    integer :: i, cnt(2)

    cnt = 0
    sc = 0.0_rp
    do i = 1, size(c)
       if (c(i) .ge. cmin) then
          sc(i) = 1.0_rp / c(i)
       else if (c(i) .gt. 0.0_rp) then
          sc(i) = 1.0_rp / cmin
          cnt(1) = cnt(1) + 1
       else
          sc(i) = 0.0_rp
          cnt(2) = cnt(2) + 1
       end if
    end do
    call MPI_Allreduce(MPI_IN_PLACE, cnt, 2, MPI_INTEGER, MPI_SUM, NEKO_COMM)
    write(log_buf, '(A,A,I0,A,I0)') trim(label), &
         ': weight sum clamped ', cnt(1), ', switched off ', cnt(2)
    call neko_log%message(log_buf)

  end subroutine idw_adjoint_normalisation

  !> Lumped weights of one side, `d = gain / (|G| 1)`. The absolute row
  !! sums are bounded from above by applying the interpolation and its
  !! transpose with the absolute values of the Lagrange weights,
  !! `|I| msk Binv gs_add(|I|^T 1)`, which needs no assembled G; the signed
  !! row sums `I msk Binv gs_add(I^T 1)` give the resulting gain on a
  !! constant marker velocity, `gain * (G 1) / (|G| 1)`, logged as min and
  !! mean over all markers together with the zero-weight count.
  subroutine idw_adjoint_lump(this, ones, rowsum, d, label, scal, msk)
    class(idw_source_term_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: ones(:), rowsum(:), d(:)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: scal(:)
    type(field_t), intent(in), optional :: msk
    real(kind=rp), allocatable :: wr(:,:), ws(:,:), wt(:,:), rowabs(:)
    real(kind=rp), allocatable :: v(:), gv(:)
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: rmax, gmin, gsum, nrm(2), lam
    integer :: i, n_lag, n_zero, n_glb, it
    integer, parameter :: n_power = 20

    n_lag = size(d)
    allocate(rowabs(n_lag))

    ! |G| 1: the same operators with the Lagrange weights replaced by
    ! their absolute values. The lumping runs on the host at start-up, so
    ! only the host weight arrays are swapped and the device copies stay
    associate (li => this%global_interp%local_interp)
      wr = li%weights_r
      ws = li%weights_s
      wt = li%weights_t
      li%weights_r(:,:) = abs(wr)
      li%weights_s(:,:) = abs(ws)
      li%weights_t(:,:) = abs(wt)
      call idw_adjoint_rowsum(this, ones, rowabs, scal, msk)
      li%weights_r(:,:) = wr
      li%weights_s(:,:) = ws
      li%weights_t(:,:) = wt
    end associate
    deallocate(wr, ws, wt)

    ! G 1: the signed row sums, for the constant-mode gain
    call idw_adjoint_rowsum(this, ones, rowsum, scal, msk)

    rmax = 0.0_rp
    do i = 1, n_lag
       rmax = max(rmax, rowabs(i))
    end do
    call MPI_Allreduce(MPI_IN_PLACE, rmax, 1, MPI_REAL_PRECISION, MPI_MAX, &
         NEKO_COMM)

    n_zero = 0
    gmin = huge(0.0_rp)
    gsum = 0.0_rp
    do i = 1, n_lag
       if (rowabs(i) .gt. 1.0e-6_rp * rmax) then
          d(i) = this%spread_gain / rowabs(i)
          gmin = min(gmin, rowsum(i) * d(i))
          gsum = gsum + rowsum(i) * d(i)
       else
          d(i) = 0.0_rp
          n_zero = n_zero + 1
       end if
    end do
    n_glb = n_lag - n_zero
    call MPI_Allreduce(MPI_IN_PLACE, n_zero, 1, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM)
    call MPI_Allreduce(MPI_IN_PLACE, n_glb, 1, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM)
    call MPI_Allreduce(MPI_IN_PLACE, gmin, 1, MPI_REAL_PRECISION, MPI_MIN, &
         NEKO_COMM)
    call MPI_Allreduce(MPI_IN_PLACE, gsum, 1, MPI_REAL_PRECISION, MPI_SUM, &
         NEKO_COMM)

    write(log_buf, '(A,A,F6.3,A,F6.3)') trim(label), &
         ': const. gain min ', gmin, ', mean ', gsum / max(n_glb, 1)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A,I6)') trim(label), ': zero-weight markers ', n_zero
    call neko_log%message(log_buf)

    ! Largest gain of interp o spread, G D, by power iteration: the
    ! Gershgorin bound above is pessimistic where the rows of G alternate
    ! in sign, and this is the number to hold under the stability limit
    ! of the time integration (4.5 for the un-extrapolated BDF3 forcing)
    ! when raising `spread_gain`. Markers held by several ranks enter the
    ! norms once per holder, which does not change the estimate's limit.
    allocate(v(n_lag), gv(n_lag))
    v = 1.0_rp
    lam = 0.0_rp
    do it = 1, n_power
       gv = d * v
       call idw_adjoint_rowsum(this, gv, rowsum, scal, msk)
       nrm(1) = 0.0_rp
       nrm(2) = 0.0_rp
       do i = 1, n_lag
          nrm(1) = nrm(1) + rowsum(i)**2
          nrm(2) = nrm(2) + v(i)**2
       end do
       call MPI_Allreduce(MPI_IN_PLACE, nrm, 2, MPI_REAL_PRECISION, MPI_SUM, &
            NEKO_COMM)
       if (nrm(1) .le. 0.0_rp .or. nrm(2) .le. 0.0_rp) exit
       lam = sqrt(nrm(1) / nrm(2))
       v = rowsum / sqrt(nrm(1))
    end do
    write(log_buf, '(A,A,F6.3,A,I0,A)') trim(label), &
         ': largest gain ', lam, ' (', n_power, ' power iterations)'
    call neko_log%message(log_buf)
    deallocate(v, gv)

    ! Leave the signed row sums in `rowsum` for the caller (the power
    ! iteration used it as scratch)
    call idw_adjoint_rowsum(this, ones, rowsum, scal, msk)

    deallocate(rowabs)

  end subroutine idw_adjoint_lump

  !> Row sums of one side's marker Gram matrix with the current weights:
  !! `rowsum = I(msk S(ones))`, S applied through the scratch forcing field.
  !! `rowsum = s I(msk S(s vals))` with the per-marker scaling s on both
  !! sides of the operator pair.
  subroutine idw_adjoint_rowsum(this, vals, rowsum, scal, msk)
    class(idw_source_term_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: vals(:), rowsum(:)
    real(kind=rp), intent(in) :: scal(:)
    type(field_t), intent(in), optional :: msk
    real(kind=rp), allocatable :: sv(:)
    integer :: n_lag

    n_lag = size(vals)
    allocate(sv(n_lag))
    sv = scal(1:n_lag) * vals

    this%ib_fx%x = 0.0_rp
    call idw_adjoint_apply(this, sv, this%ib_fx, .true., msk)
    rowsum = 0.0_rp
    ! Host arrays throughout: field_col3 would act on the device mirrors on
    ! a device build and leave the host side untouched
    if (present(msk)) then
       this%tmp%x = this%ib_fx%x * msk%x
       call this%global_interp%evaluate(rowsum, this%tmp%x, .true.)
    else
       call this%global_interp%evaluate(rowsum, this%ib_fx%x, .true.)
    end if
    rowsum = scal(1:n_lag) * rowsum
    this%ib_fx%x = 0.0_rp
    deallocate(sv)

  end subroutine idw_adjoint_rowsum

  !> Adjoint spectral spread of the direct forcing: per side and
  !! component, `ib += -(1/dt) msk Binv gs_add(I^T (D u_m))` with the
  !! marker velocities `fu_ib` (+ side) and `fum_ib` (- side) interpolated
  !! before. The result is continuous (assembled inside the operator).
  !! With `on_host` false the marker velocities are read from their device
  !! buffers and the deposit runs on the device.
  subroutine idw_spread_adjoint(this, dt, on_host)
    class(idw_source_term_t), intent(inout) :: this
    real(kind=dp), intent(in) :: dt
    logical, intent(in) :: on_host
    real(kind=rp) :: idt
    integer :: n_lag

    n_lag = size(this%fu_ib)
    idt = -1.0_rp / real(dt, kind=rp)

    ! The marker velocities fu_ib etc. are already normalised (s I(sw u));
    ! the deposit carries d s, so the pair interp/spread stays adjoint
    if (NEKO_BCKND_DEVICE .eq. 1 .and. .not. on_host) then
       call idw_spread_adjoint_side(this, this%fu_ib_d, this%dsp_mark_d, &
            this%ib_fx, idt, n_lag, this%one_sided, this%swp)
       call idw_spread_adjoint_side(this, this%fv_ib_d, this%dsp_mark_d, &
            this%ib_fy, idt, n_lag, this%one_sided, this%swp)
       call idw_spread_adjoint_side(this, this%fw_ib_d, this%dsp_mark_d, &
            this%ib_fz, idt, n_lag, this%one_sided, this%swp)
       if (this%one_sided) then
          call idw_spread_adjoint_side(this, this%fum_ib_d, &
               this%dsm_mark_d, this%ib_fx, idt, n_lag, .true., this%swm)
          call idw_spread_adjoint_side(this, this%fvm_ib_d, &
               this%dsm_mark_d, this%ib_fy, idt, n_lag, .true., this%swm)
          call idw_spread_adjoint_side(this, this%fwm_ib_d, &
               this%dsm_mark_d, this%ib_fz, idt, n_lag, .true., this%swm)
       end if
    else if (this%one_sided) then
       this%vals_adj(1:n_lag) = idt * this%dsp_mark(1:n_lag) * this%fu_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fx, .true., &
            this%swp)
       this%vals_adj(1:n_lag) = idt * this%dsp_mark(1:n_lag) * this%fv_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fy, .true., &
            this%swp)
       this%vals_adj(1:n_lag) = idt * this%dsp_mark(1:n_lag) * this%fw_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fz, .true., &
            this%swp)

       this%vals_adj(1:n_lag) = idt * this%dsm_mark(1:n_lag) * this%fum_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fx, .true., &
            this%swm)
       this%vals_adj(1:n_lag) = idt * this%dsm_mark(1:n_lag) * this%fvm_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fy, .true., &
            this%swm)
       this%vals_adj(1:n_lag) = idt * this%dsm_mark(1:n_lag) * this%fwm_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fz, .true., &
            this%swm)
    else
       this%vals_adj(1:n_lag) = idt * this%dsp_mark(1:n_lag) * this%fu_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fx, .true.)
       this%vals_adj(1:n_lag) = idt * this%dsp_mark(1:n_lag) * this%fv_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fy, .true.)
       this%vals_adj(1:n_lag) = idt * this%dsp_mark(1:n_lag) * this%fw_ib
       call idw_adjoint_apply(this, this%vals_adj, this%ib_fz, .true.)
    end if

  end subroutine idw_spread_adjoint

  !> One component and side of the device adjoint spread: scale the
  !! marker velocities `u_d` by `idt * d_d` into the mapped scratch and
  !! apply the (masked) adjoint.
  subroutine idw_spread_adjoint_side(this, u_d, d_d, out, idt, n_lag, &
       masked, msk)
    class(idw_source_term_t), intent(inout) :: this
    type(c_ptr), intent(in) :: u_d, d_d
    type(field_t), intent(inout) :: out
    real(kind=rp), intent(in) :: idt
    integer, intent(in) :: n_lag
    logical, intent(in) :: masked
    type(field_t), intent(in) :: msk

    if (n_lag .gt. 0) then
       call device_col3(this%vals_adj_d, d_d, u_d, n_lag)
       call device_cmult(this%vals_adj_d, idt, n_lag)
    end if
    if (masked) then
       call idw_adjoint_apply(this, this%vals_adj, out, .false., msk)
    else
       call idw_adjoint_apply(this, this%vals_adj, out, .false.)
    end if

  end subroutine idw_spread_adjoint_side

  !> Host interpolation of the adjoint spread: the side-weighted spectral
  !! evaluate normalised per marker, `s I(sw u)`, both sides under
  !! `one_sided`. Host arrays throughout (CPU builds only).
  subroutine idw_interp_adjoint_host(this, u, v, w)
    class(idw_source_term_t), intent(inout) :: this
    type(field_t), intent(inout) :: u, v, w
    integer :: n_lag

    n_lag = size(this%fu_ib)
    if (this%one_sided) then
       this%tmp%x = u%x * this%swp%x
       call this%global_interp%evaluate(this%fu_ib, this%tmp%x, .true.)
       this%tmp%x = v%x * this%swp%x
       call this%global_interp%evaluate(this%fv_ib, this%tmp%x, .true.)
       this%tmp%x = w%x * this%swp%x
       call this%global_interp%evaluate(this%fw_ib, this%tmp%x, .true.)
       this%tmp%x = u%x * this%swm%x
       call this%global_interp%evaluate(this%fum_ib, this%tmp%x, .true.)
       this%tmp%x = v%x * this%swm%x
       call this%global_interp%evaluate(this%fvm_ib, this%tmp%x, .true.)
       this%tmp%x = w%x * this%swm%x
       call this%global_interp%evaluate(this%fwm_ib, this%tmp%x, .true.)
       this%fu_ib = this%sp_mark(1:n_lag) * this%fu_ib
       this%fv_ib = this%sp_mark(1:n_lag) * this%fv_ib
       this%fw_ib = this%sp_mark(1:n_lag) * this%fw_ib
       this%fum_ib = this%sm_mark(1:n_lag) * this%fum_ib
       this%fvm_ib = this%sm_mark(1:n_lag) * this%fvm_ib
       this%fwm_ib = this%sm_mark(1:n_lag) * this%fwm_ib
    else
       call this%global_interp%evaluate(this%fu_ib, u%x, .true.)
       call this%global_interp%evaluate(this%fv_ib, v%x, .true.)
       call this%global_interp%evaluate(this%fw_ib, w%x, .true.)
       this%fu_ib = this%sp_mark(1:n_lag) * this%fu_ib
       this%fv_ib = this%sp_mark(1:n_lag) * this%fv_ib
       this%fw_ib = this%sp_mark(1:n_lag) * this%fw_ib
    end if

  end subroutine idw_interp_adjoint_host

  !> Device path of the adjoint spectral spread: the spectral (masked)
  !! interpolation evaluated straight into the device marker buffers, as
  !! in `idw_compute_device`, then the device spread.
  subroutine idw_compute_device_adjoint(this, u, v, w, time)
    class(idw_source_term_t), intent(inout) :: this
    type(field_t), intent(inout) :: u, v, w
    type(time_state_t), intent(in) :: time
    integer :: nt, n_lag

    nt = this%tmp%size()
    n_lag = size(this%fu_ib)

    ! Side-weighted interpolation, then the per-marker normalisation s
    if (this%one_sided) then
       call idw_interp_masked(this%global_interp, this%tmp, u, this%swp, &
            this%fu_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, v, this%swp, &
            this%fv_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, w, this%swp, &
            this%fw_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, u, this%swm, &
            this%fum_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, v, this%swm, &
            this%fvm_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, w, this%swm, &
            this%fwm_ib, nt)
       if (n_lag .gt. 0) then
          call device_col2(this%fu_ib_d, this%sp_mark_d, n_lag)
          call device_col2(this%fv_ib_d, this%sp_mark_d, n_lag)
          call device_col2(this%fw_ib_d, this%sp_mark_d, n_lag)
          call device_col2(this%fum_ib_d, this%sm_mark_d, n_lag)
          call device_col2(this%fvm_ib_d, this%sm_mark_d, n_lag)
          call device_col2(this%fwm_ib_d, this%sm_mark_d, n_lag)
       end if
    else
       call this%global_interp%evaluate(this%fu_ib, u%x, .false.)
       call this%global_interp%evaluate(this%fv_ib, v%x, .false.)
       call this%global_interp%evaluate(this%fw_ib, w%x, .false.)
       if (n_lag .gt. 0) then
          call device_col2(this%fu_ib_d, this%sp_mark_d, n_lag)
          call device_col2(this%fv_ib_d, this%sp_mark_d, n_lag)
          call device_col2(this%fw_ib_d, this%sp_mark_d, n_lag)
       end if
    end if

    call idw_spread_adjoint(this, time%dt, .false.)

  end subroutine idw_compute_device_adjoint

  !> Build the CSR transpose of `lag_el` (per element -> lag points). The
  !! lag points of an element are listed in increasing order, so a gather
  !! over an element's list adds the contributions to each dof in the same
  !! order as a scatter over the lag points. This lets the host loops run
  !! element-parallel without races, and gives the device gather kernel
  !! its connectivity.
  subroutine idw_build_el_csr(this, nelv)
    class(idw_source_term_t), intent(inout) :: this
    integer, intent(in) :: nelv
    integer :: n_lag, i, ee, e, k, n_csr
    integer, allocatable :: cursor(:)

    n_lag = size(this%lag_pts)

    ! Histogram of contributions per element (stored shifted by one slot)
    allocate(this%el_off(nelv + 1))
    this%el_off = 0
    do i = 1, n_lag
       select type (el => this%lag_el(i)%data)
         type is (integer)
          do ee = 1, this%lag_el(i)%size()
             e = el(ee)
             this%el_off(e + 1) = this%el_off(e + 1) + 1
          end do
       end select
    end do

    ! Prefix sum -> 0-based CSR offsets; element e occupies [el_off(e), el_off(e+1))
    do e = 1, nelv
       this%el_off(e + 1) = this%el_off(e + 1) + this%el_off(e)
    end do
    n_csr = this%el_off(nelv + 1)

    allocate(this%el_lag(n_csr))
    allocate(cursor(nelv))
    cursor(1:nelv) = this%el_off(1:nelv)

    ! Fill the CSR buffer (0-based lag indices for the kernel)
    do i = 1, n_lag
       select type (el => this%lag_el(i)%data)
         type is (integer)
          do ee = 1, this%lag_el(i)%size()
             e = el(ee)
             this%el_lag(cursor(e) + 1) = i - 1
             cursor(e) = cursor(e) + 1
          end do
       end select
    end do
    deallocate(cursor)

    ! Compact list of elements that actually receive contributions
    this%n_active = count(this%el_off(2:nelv + 1) - this%el_off(1:nelv) > 0)
    allocate(this%active_el(this%n_active))
    k = 0
    do e = 1, nelv
       if (this%el_off(e + 1) > this%el_off(e)) then
          k = k + 1
          this%active_el(k) = e - 1
       end if
    end do

    ! The lists themselves as a CSR (per lag point -> elements, 0-based),
    ! which the device interpolation kernel walks marker by marker
    allocate(this%lag_off(n_lag + 1))
    this%lag_off(1) = 0
    do i = 1, n_lag
       this%lag_off(i + 1) = this%lag_off(i) + this%lag_el(i)%size()
    end do
    allocate(this%lag_els(max(this%lag_off(n_lag + 1), 1)))
    do i = 1, n_lag
       select type (el => this%lag_el(i)%data)
         type is (integer)
          do ee = 1, this%lag_el(i)%size()
             this%lag_els(this%lag_off(i) + ee) = el(ee) - 1
          end do
       end select
    end do

  end subroutine idw_build_el_csr

  !> Build the device data structures for the IDW source term: upload the
  !! CSR built by idw_build_el_csr, the unpacked Lagrangian coordinates, and
  !! the mask/weight fields that are assembled on the host during init.
  subroutine idw_build_device_maps(this)
    class(idw_source_term_t), intent(inout) :: this
    integer :: nelv, n_lag, n, i, n_csr

    nelv = size(this%el_off) - 1
    n_lag = size(this%lag_pts)
    n_csr = this%el_off(nelv + 1)
    this%lx3 = this%w%Xh%lx**3

    ! Unpack Lagrangian coordinates
    allocate(this%lpx(n_lag), this%lpy(n_lag), this%lpz(n_lag))
    do i = 1, n_lag
       this%lpx(i) = this%lag_pts(i)%x(1)
       this%lpy(i) = this%lag_pts(i)%x(2)
       this%lpz(i) = this%lag_pts(i)%x(3)
    end do

    ! Map and upload the connectivity / per-point arrays (static after init)
    call device_map(this%el_off, this%el_off_d, nelv + 1)
    call device_memcpy(this%el_off, this%el_off_d, nelv + 1, &
         HOST_TO_DEVICE, sync = .false.)

    if (this%n_active > 0) then
       call device_map(this%active_el, this%active_el_d, this%n_active)
       call device_memcpy(this%active_el, this%active_el_d, this%n_active, &
            HOST_TO_DEVICE, sync = .false.)
    end if

    if (n_csr > 0) then
       call device_map(this%el_lag, this%el_lag_d, n_csr)
       call device_memcpy(this%el_lag, this%el_lag_d, n_csr, &
            HOST_TO_DEVICE, sync = .false.)
    end if

    ! Partial sums of the device interpolation kernel; allocated for every
    ! rank, since the finalisation that reads them is collective
    allocate(this%part(8, n_lag))

    if (n_lag > 0) then
       call device_map(this%lpx, this%lpx_d, n_lag)
       call device_map(this%lpy, this%lpy_d, n_lag)
       call device_map(this%lpz, this%lpz_d, n_lag)
       call device_memcpy(this%lpx, this%lpx_d, n_lag, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%lpy, this%lpy_d, n_lag, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%lpz, this%lpz_d, n_lag, HOST_TO_DEVICE, &
            sync = .false.)

       ! Device buffers for the interpolated values (refilled every step)
       call device_map(this%fu_ib, this%fu_ib_d, n_lag)
       call device_map(this%fv_ib, this%fv_ib_d, n_lag)
       call device_map(this%fw_ib, this%fw_ib_d, n_lag)
       call device_map(this%fum_ib, this%fum_ib_d, n_lag)
       call device_map(this%fvm_ib, this%fvm_ib_d, n_lag)
       call device_map(this%fwm_ib, this%fwm_ib_d, n_lag)

       ! Marker -> element CSR and partial-sum buffer of the interpolation
       ! kernel
       call device_map(this%lag_off, this%lag_off_d, n_lag + 1)
       call device_memcpy(this%lag_off, this%lag_off_d, n_lag + 1, &
            HOST_TO_DEVICE, sync = .false.)
       call device_map(this%lag_els, this%lag_els_d, size(this%lag_els))
       call device_memcpy(this%lag_els, this%lag_els_d, size(this%lag_els), &
            HOST_TO_DEVICE, sync = .false.)
       call device_map(this%part, this%part_d, 8 * n_lag)
    end if

    ! Upload the mask / weight fields assembled on the host during init
    n = this%w%size()
    call device_map(this%Bm, this%Bm_d, n)
    call device_memcpy(this%Bm, this%Bm_d, n, HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%w%x, this%w%x_d, n, HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%wm%x, this%wm%x_d, n, HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%ds%x, this%ds%x_d, n, HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%pmsk%x, this%pmsk%x_d, n, HOST_TO_DEVICE, &
         sync = .false.)
    call device_memcpy(this%mmsk%x, this%mmsk%x_d, n, HOST_TO_DEVICE, &
         sync = .true.)

  end subroutine idw_build_device_maps

  !> Interpolate a masked field at the Lagrangian points
  subroutine idw_interp_masked(global_interp, tmp, fld, msk, ib, nt)
    type(global_interpolation_t), intent(inout) :: global_interp
    type(field_t), intent(inout) :: tmp
    type(field_t), intent(in) :: fld, msk
    real(kind=rp), intent(inout) :: ib(:)
    integer, intent(in) :: nt

    call field_col3(tmp, fld, msk, nt)
    call global_interp%evaluate(ib, tmp%x, NEKO_BCKND_DEVICE .ne. 1)
  end subroutine idw_interp_masked

  !> Device path of the source term: interpolate on the device (the Shepard
  !! and adjoint partial sums with a marker-parallel kernel, the barycentric
  !! values with the global interpolation's device path), finalise the
  !! marker values on the host, then assemble the contributions with the
  !! atomic-free gather kernel. Mirrors the host branches of
  !! idw_source_term_compute. The only per-step host traffic is marker
  !! sized: the 8 partial sums per marker down and the 6 values up.
  subroutine idw_compute_device(this, fu, fv, fw, u, v, w, time)
    class(idw_source_term_t), intent(inout) :: this
    type(field_t), intent(inout) :: fu, fv, fw
    type(field_t), intent(inout) :: u, v, w
    type(time_state_t), intent(in) :: time
    integer :: n_lag, nt
    real(kind=rp) :: wtol

    n_lag = size(this%fu_ib)
    nt = this%tmp%size()

    this%fu_ib = 0.0_rp
    this%fv_ib = 0.0_rp
    this%fw_ib = 0.0_rp

    if (this%idw_interp) then
       ! Shepard / adjoint partial sums on the device, one thread block per
       ! marker; the 8 sums per marker come back for the shared-marker
       ! reduction and the normalisation, which stay on the host
       if (n_lag > 0) then
          call device_idw_interp_partials(this%part_d, u%x_d, v%x_d, w%x_d, &
               this%w%dof%x%x_d, this%w%dof%y%x_d, this%w%dof%z%x_d, &
               this%ds%x_d, this%pmsk%x_d, this%coef%mult_d, this%Bm_d, &
               this%w%x_d, this%wm%x_d, this%lpx_d, this%lpy_d, this%lpz_d, &
               this%lag_off_d, this%lag_els_d, n_lag, this%lx3, &
               this%interp_rmax, this%pwr_param, NEKO_EPS, 1.0e-10_rp, &
               this%adjoint_interp)
          call device_memcpy(this%part, this%part_d, 8 * n_lag, &
               DEVICE_TO_HOST, sync = .true.)
       end if
       call idw_interp_shepard_finalize(this%fu_ib, this%fv_ib, this%fw_ib, &
            this%fum_ib, this%fvm_ib, this%fwm_ib, this%part, &
            this%shared_slot, this%n_shared_glb)
    else if (this%one_sided) then
       this%fum_ib = 0.0_rp
       this%fvm_ib = 0.0_rp
       this%fwm_ib = 0.0_rp

       call idw_interp_masked(this%global_interp, this%tmp, u, this%pmsk, &
            this%fu_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, v, this%pmsk, &
            this%fv_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, w, this%pmsk, &
            this%fw_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, u, this%mmsk, &
            this%fum_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, v, this%mmsk, &
            this%fvm_ib, nt)
       call idw_interp_masked(this%global_interp, this%tmp, w, this%mmsk, &
            this%fwm_ib, nt)
    else
       ! Unmasked barycentric interpolation, evaluated on the device straight
       ! into the device buffers of the marker arrays
       call this%global_interp%evaluate(this%fu_ib, u%x, .false.)
       call this%global_interp%evaluate(this%fv_ib, v%x, .false.)
       call this%global_interp%evaluate(this%fw_ib, w%x, .false.)
    end if

    ! Stage the per-point values the gather kernel reads. The Shepard and
    ! adjoint values are finalised on the host and uploaded; the barycentric
    ! ones were evaluated straight into the device buffers. Both the uploads
    ! and the evaluation sit on the command queue the gather kernel runs on,
    ! so stream order is enough and no host synchronisation is needed here;
    ! the host arrays are not touched again before the next step's
    ! synchronous download of the partial sums.
    if (n_lag > 0 .and. this%idw_interp) then
       call device_memcpy(this%fu_ib, this%fu_ib_d, n_lag, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%fv_ib, this%fv_ib_d, n_lag, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%fw_ib, this%fw_ib_d, n_lag, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%fum_ib, this%fum_ib_d, n_lag, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%fvm_ib, this%fvm_ib_d, n_lag, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%fwm_ib, this%fwm_ib_d, n_lag, &
            HOST_TO_DEVICE, sync = .false.)
    end if

    ! one_sided uses the pmsk split with tol 1e-10; the unmasked path has
    ! pmsk == 1 everywhere and matches the host tol of 1e-12.
    wtol = merge(1.0e-10_rp, 1.0e-12_rp, this%one_sided)

    call device_idw_gather(fu%x_d, fv%x_d, fw%x_d, &
         this%fu_ib_d, this%fv_ib_d, this%fw_ib_d, &
         this%fum_ib_d, this%fvm_ib_d, this%fwm_ib_d, &
         this%w%dof%x%x_d, this%w%dof%y%x_d, this%w%dof%z%x_d, this%ds%x_d, &
         this%pmsk%x_d, this%w%x_d, this%wm%x_d, &
         this%lpx_d, this%lpy_d, this%lpz_d, &
         this%active_el_d, this%el_off_d, this%el_lag_d, &
         this%n_active, this%lx3, real(time%dt, kind=rp), this%rmax, &
         this%pwr_param, &
         NEKO_EPS, wtol)

  end subroutine idw_compute_device

  subroutine idw_source_term_free(this)
    class(idw_source_term_t), intent(inout) :: this
    integer :: i

    call this%free_base()

    call this%ds%free()

    if (allocated(this%xyz)) then
       deallocate(this%xyz)
    end if

    if (allocated(this%dp_mark)) deallocate(this%dp_mark)
    if (allocated(this%dm_mark)) deallocate(this%dm_mark)
    if (allocated(this%gp_mark)) deallocate(this%gp_mark)
    if (allocated(this%gm_mark)) deallocate(this%gm_mark)
    if (allocated(this%cp_mark)) deallocate(this%cp_mark)
    if (allocated(this%cm_mark)) deallocate(this%cm_mark)
    if (allocated(this%sp_mark)) deallocate(this%sp_mark)
    if (allocated(this%sm_mark)) deallocate(this%sm_mark)
    if (allocated(this%dsp_mark)) deallocate(this%dsp_mark)
    if (allocated(this%dsm_mark)) deallocate(this%dsm_mark)
    if (c_associated(this%sp_mark_d)) call device_free(this%sp_mark_d)
    if (c_associated(this%sm_mark_d)) call device_free(this%sm_mark_d)
    if (c_associated(this%dsp_mark_d)) call device_free(this%dsp_mark_d)
    if (c_associated(this%dsm_mark_d)) call device_free(this%dsm_mark_d)
    call this%swp%free()
    call this%swm%free()
    if (allocated(this%vals_adj)) deallocate(this%vals_adj)
    call this%winv%free()
    if (c_associated(this%dp_mark_d)) call device_free(this%dp_mark_d)
    if (c_associated(this%dm_mark_d)) call device_free(this%dm_mark_d)
    if (c_associated(this%vals_adj_d)) call device_free(this%vals_adj_d)

    if (allocated(this%fu_ib)) then
       deallocate(this%fu_ib)
    end if

    if (allocated(this%fv_ib)) then
       deallocate(this%fv_ib)
    end if

    if (allocated(this%fw_ib)) then
       deallocate(this%fw_ib)
    end if

    if (allocated(this%fum_ib)) then
       deallocate(this%fum_ib)
    end if

    if (allocated(this%fvm_ib)) then
       deallocate(this%fvm_ib)
    end if

    if (allocated(this%fwm_ib)) then
       deallocate(this%fwm_ib)
    end if

    if (c_associated(this%fu_ib_d)) then
       call device_free(this%fu_ib_d)
    end if
    
    if (c_associated(this%fv_ib_d)) then
       call device_free(this%fv_ib_d)
    end if
    
    if (c_associated(this%fw_ib_d)) then
       call device_free(this%fw_ib_d)
    end if
    
    if (c_associated(this%fum_ib_d)) then
       call device_free(this%fum_ib_d)
    end if
    
    if (c_associated(this%fvm_ib_d)) then
       call device_free(this%fvm_ib_d)
    end if
    
    if (c_associated(this%fwm_ib_d)) then
       call device_free(this%fwm_ib_d)
    end if
    
    if (allocated(this%lag_pts)) then
       deallocate(this%lag_pts)
    end if
    
    if (allocated(this%lag_el)) then
       do i = 1, size(this%lag_el)
          call this%lag_el(i)%free()
       end do
       deallocate(this%lag_el)
    end if

    if (allocated(this%fltr)) then
       call this%fltr%free()
       deallocate(this%fltr)
    end if

    if (c_associated(this%el_off_d)) call device_free(this%el_off_d)
    if (c_associated(this%el_lag_d)) call device_free(this%el_lag_d)
    if (c_associated(this%active_el_d)) call device_free(this%active_el_d)
    if (c_associated(this%lpx_d)) call device_free(this%lpx_d)
    if (c_associated(this%lpy_d)) call device_free(this%lpy_d)
    if (c_associated(this%lpz_d)) call device_free(this%lpz_d)

    if (allocated(this%el_off)) deallocate(this%el_off)
    if (allocated(this%el_lag)) deallocate(this%el_lag)
    if (allocated(this%active_el)) deallocate(this%active_el)
    if (allocated(this%lpx)) deallocate(this%lpx)
    if (allocated(this%lpy)) deallocate(this%lpy)
    if (allocated(this%lpz)) deallocate(this%lpz)

    if (allocated(this%shared_slot)) deallocate(this%shared_slot)
    this%n_shared_glb = 0
    if (allocated(this%Bm)) deallocate(this%Bm)
    if (allocated(this%lag_off)) deallocate(this%lag_off)
    if (allocated(this%lag_els)) deallocate(this%lag_els)
    if (allocated(this%part)) deallocate(this%part)
    if (c_associated(this%lag_off_d)) call device_free(this%lag_off_d)
    if (c_associated(this%lag_els_d)) call device_free(this%lag_els_d)
    if (c_associated(this%part_d)) call device_free(this%part_d)
    if (c_associated(this%Bm_d)) call device_free(this%Bm_d)

    call this%gs%free()

    call this%ib_fx%free()
    call this%ib_fy%free()
    call this%ib_fz%free()

    call this%force_file%free()
    call this%force_row%free()
    call this%force_controller%free()
    if (allocated(this%rho_name)) deallocate(this%rho_name)
    this%force_output = .false.
    this%force = 0.0_rp
    this%force_scale = 1.0_rp

  end subroutine idw_source_term_free

  subroutine idw_source_term_compute(this, time)
    class(idw_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, fu, fv, fw
    integer :: n

    n = this%fields%item_size(1)

    u => neko_registry%get_field('u')
    v => neko_registry%get_field('v')
    w => neko_registry%get_field('w')

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)

    ! The IB forcing is spread into its own fields, so that the assembly
    ! and the filter below only act on it and its integral can be taken
    call field_rzero(this%ib_fx)
    call field_rzero(this%ib_fy)
    call field_rzero(this%ib_fz)

    associate(global_interp => this%global_interp, &
         fu_ib => this%fu_ib, fv_ib => this%fv_ib, fw_ib => this%fw_ib, &
         fum_ib => this%fum_ib, fvm_ib => this%fvm_ib, fwm_ib => this%fwm_ib, &
         lag_pts => this%lag_pts, tmp => this%tmp, &
         ds => this%ds%x, x => this%w%dof%x%x, &
         y => this%w%dof%y%x, z => this%w%dof%z%x, lx => this%w%Xh%lx, &
         ibx => this%ib_fx%x, iby => this%ib_fy%x, ibz => this%ib_fz%x)


      fu_ib = 0.0_rp
      fv_ib = 0.0_rp
      fw_ib = 0.0_rp

      call profiler_start_region('IDW interpolation and spread')
      if (NEKO_BCKND_DEVICE .eq. 1) then
         if (this%adjoint_spread) then
            call idw_compute_device_adjoint(this, u, v, w, time)
         else
            call idw_compute_device(this, this%ib_fx, this%ib_fy, &
                 this%ib_fz, u, v, w, time)
         end if
      else
         ! Interpolation stage: velocity at the Lagrangian points
         if (this%idw_interp) then
            call idw_interp_shepard(fu_ib, fv_ib, fw_ib, &
                 fum_ib, fvm_ib, fwm_ib, lag_pts, this%lag_el, &
                 u%x, v%x, w%x, this%pmsk%x, this%coef%mult, x, y, z, ds, &
                 this%interp_rmax, this%pwr_param, lx, this%coef%msh%nelv, &
                 this%shared_slot, this%n_shared_glb, &
                 adjoint = this%adjoint_interp, B = this%Bm, &
                 sw_p = this%w%x, sw_m = this%wm%x)
         else if (this%adjoint_spread) then
            ! Side-weighted interpolation normalised per marker, s I(sw u)
            call idw_interp_adjoint_host(this, u, v, w)
         else if (this%one_sided) then
            fum_ib = 0.0_rp
            fvm_ib = 0.0_rp
            fwm_ib = 0.0_rp

            call field_col3(tmp, u, this%pmsk, tmp%size())
            call global_interp%evaluate(fu_ib, tmp%x, .true.)

            call field_col3(tmp, v, this%pmsk, tmp%size())
            call global_interp%evaluate(fv_ib, tmp%x, .true.)

            call field_col3(tmp, w, this%pmsk, tmp%size())
            call global_interp%evaluate(fw_ib, tmp%x, .true.)

            call field_col3(tmp, u, this%mmsk, tmp%size())
            call global_interp%evaluate(fum_ib, tmp%x, .true.)

            call field_col3(tmp, v, this%mmsk, tmp%size())
            call global_interp%evaluate(fvm_ib, tmp%x, .true.)

            call field_col3(tmp, w, this%mmsk, tmp%size())
            call global_interp%evaluate(fwm_ib, tmp%x, .true.)
         else
            call global_interp%evaluate(fu_ib, u%x, .true.)
            call global_interp%evaluate(fv_ib, v%x, .true.)
            call global_interp%evaluate(fw_ib, w%x, .true.)
         end if

         ! Spread stage: accumulate into the IB forcing fields
         if (this%adjoint_spread) then
            call idw_spread_adjoint(this, time%dt, .true.)
         else
            call idw_spread(ibx, iby, ibz, fu_ib, fv_ib, fw_ib, &
                 fum_ib, fvm_ib, fwm_ib, lag_pts, this%active_el, &
                 this%el_off, this%el_lag, x, y, z, ds, this%pmsk%x, &
                 this%w%x, this%wm%x, this%rmax, this%pwr_param, time%dt, &
                 this%one_sided, lx, this%coef%msh%nelv)
         end if
      end if

    end associate
    call profiler_end_region('IDW interpolation and spread')

    ! Assemble the IB forcing to a continuous representation: the three
    ! components scaled in one kernel and exchanged in one halo round.
    ! The adjoint spread is assembled (Binv gs_add) inside the operator
    ! and is already continuous.
    if (.not. this%adjoint_spread) then
       call profiler_start_region('IDW assembly')
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_opcolv(this%ib_fx%x_d, this%ib_fy%x_d, &
               this%ib_fz%x_d, this%coef%mult_d, 3, n)
       else
          call opcolv(this%ib_fx%x, this%ib_fy%x, this%ib_fz%x, &
               this%coef%mult, 3, n)
       end if

       call this%gs%op(this%ib_fx%x, this%ib_fy%x, this%ib_fz%x, n, &
            GS_OP_ADD, glb_cmd_event)
       call device_event_sync(glb_cmd_event)
       call profiler_end_region('IDW assembly')
    end if

    if (allocated(this%fltr)) then
       call profiler_start_region('IDW filter')
       call field_copy(this%tmp, this%ib_fx)
       call this%fltr%apply(this%ib_fx, this%tmp)

       call field_copy(this%tmp, this%ib_fy)
       call this%fltr%apply(this%ib_fy, this%tmp)

       call field_copy(this%tmp, this%ib_fz)
       call this%fltr%apply(this%ib_fz, this%tmp)
       call profiler_end_region('IDW filter')
    end if

    call field_add2(fu, this%ib_fx)
    call field_add2(fv, this%ib_fy)
    call field_add2(fw, this%ib_fz)

    if (this%force_output) then
       call profiler_start_region('IDW force')
       call this%write_force(time)
       call profiler_end_region('IDW force')
    end if

  end subroutine idw_source_term_compute

  !> Parse the optional `force_output` object and set up the csv file.
  !! @param json The JSON object of the source term.
  !! @param variable_name Name of the fluid scheme, used to look up the
  !! density `<variable_name>_rho` in the registry.
  subroutine idw_init_force_output(this, json, variable_name)
    class(idw_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: variable_name
    character(len=:), allocatable :: fname, control
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp) :: value
    integer :: nsteps
    logical :: overwrite

    call json_get_or_default(json, 'force_output.output_file', fname, &
         'idw_force.csv')
    call json_get_or_default(json, 'force_output.output_control', control, &
         'tsteps')
    call json_get_or_default(json, 'force_output.scale', this%force_scale, &
         1.0_rp)
    ! Append by default, so that a restart keeps the earlier history
    call json_get_or_default(json, 'force_output.overwrite', overwrite, &
         .false.)

    ! A source term does not know the end time of the case, so 'nsamples'
    ! is not supported
    select case (trim(control))
    case ('tsteps')
       call json_get_or_default(json, 'force_output.output_value', nsteps, 1)
       if (nsteps .lt. 1) then
          call neko_error('IDW force_output: output_value must be positive')
       end if
       value = real(nsteps, rp)
    case ('simulationtime')
       call json_get_or_default(json, 'force_output.output_value', value, &
            1.0_rp)
       if (value .le. 0.0_rp) then
          call neko_error('IDW force_output: output_value must be positive')
       end if
    case default
       call neko_error('IDW force_output: output_control must be tsteps &
            &or simulationtime')
    end select

    call this%force_controller%init(real(this%start_time, dp), &
         real(this%end_time, dp), control, real(value, dp))

    this%rho_name = trim(variable_name) // '_rho'

    call this%force_file%init(trim(fname), header = 'tstep,time,Fx,Fy,Fz', &
         overwrite = overwrite)
    call this%force_row%init(5)
    this%force_output = .true.

    call neko_log%message('Force file : ' // trim(fname))
    if (trim(control) .eq. 'tsteps') then
       write(log_buf, '(A,I0,A)') 'Force out  : every ', nsteps, ' tsteps'
    else
       write(log_buf, '(A,ES13.6)') 'Force out  : every ', value
    end if
    call neko_log%message(log_buf)

  end subroutine idw_init_force_output

  !> Compute the force on the immersed objects and write it to the csv
  !! file, if the output controller says so. The force on the objects is
  !! minus the momentum the IB forcing adds to the fluid,
  !! \f$ F = -\rho \int f_{IB} \, dV \f$ (times `scale`), integrated with
  !! the local mass matrix over the assembled and filtered forcing, i.e.
  !! the forcing that is added to the right-hand side. The forcing is built
  !! from the velocity at the start of the step and the scheme applies it as
  !! computed, without time extrapolation, so this is the force applied in
  !! the step. The rate of change of the
  !! fluid momentum inside the objects is not included, and forcing on
  !! nodes with strong velocity boundary conditions is included although
  !! the velocity solve discards it, which over-reports the force on
  !! objects touching such boundaries. Collective, all ranks must call it.
  subroutine idw_write_force(this, time)
    class(idw_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rho
    integer :: n

    if (.not. this%force_controller%check(time)) return

    n = this%ib_fx%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       this%force(1) = device_glsc2(this%ib_fx%x_d, this%coef%B_d, n)
       this%force(2) = device_glsc2(this%ib_fy%x_d, this%coef%B_d, n)
       this%force(3) = device_glsc2(this%ib_fz%x_d, this%coef%B_d, n)
    else
       this%force(1) = glsc2(this%ib_fx%x, this%coef%B, n)
       this%force(2) = glsc2(this%ib_fy%x, this%coef%B, n)
       this%force(3) = glsc2(this%ib_fz%x, this%coef%B, n)
    end if

    rho => neko_registry%get_field(this%rho_name)
    this%force = -this%force_scale * rho%x(1,1,1,1) * this%force

    this%force_row%x(1) = real(time%tstep, rp)
    this%force_row%x(2) = real(time%t, rp)
    this%force_row%x(3:5) = this%force
    call this%force_file%write(this%force_row)
    call this%force_controller%register_execution(time)

  end subroutine idw_write_force

  !> Shepard (inverse-distance-weighted) interpolation of the velocity onto
  !! the Lagrangian points, mirroring the spread operator: same kernel and
  !! stencil, hard pmsk side test, row-normalized per marker and side.
  !! Contributions carry the dof multiplicity weight so that shared nodes,
  !! visited once per adjoining stencil element, count once. Sides whose
  !! weight sum is below tolerance return zero.
  !!
  !! A marker whose stencil straddles an MPI boundary is held by every rank
  !! whose padded elements it overlaps, each copy seeing only the local part
  !! of the stencil. When `shared_slot` and `n_shared` are passed, the
  !! partial sums of such markers are summed over their holder ranks before
  !! normalization, so every copy computes the full-stencil value. The
  !! reduction is collective: all ranks in NEKO_COMM must call this routine.
  !!
  !! With `adjoint` set (and `B`, `sw_p`, `sw_m` passed) every weight is
  !! multiplied by the mass matrix and divided by the assembled spread
  !! weight of the node's side, so that with `rmax_i` equal to the spread
  !! radius the interpolation is the mass-weighted adjoint of the spread.
  subroutine idw_interp_shepard(fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, fwm_ib, &
       lag_pts, lag_el, u, v, w, pmsk, mult, x, y, z, ds, rmax_i, p, lx, ne, &
       shared_slot, n_shared, adjoint, B, sw_p, sw_m)
    real(kind=rp), intent(inout) :: fu_ib(:), fv_ib(:), fw_ib(:)
    real(kind=rp), intent(inout) :: fum_ib(:), fvm_ib(:), fwm_ib(:)
    type(point_t), intent(in) :: lag_pts(:)
    type(stack_i4_t), intent(inout) :: lag_el(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: u, v, w
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: pmsk, mult, x, y, z, ds
    real(kind=rp), intent(in) :: rmax_i, p
    integer, intent(in), optional :: shared_slot(:)
    integer, intent(in), optional :: n_shared
    logical, intent(in), optional :: adjoint
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in), optional :: B, sw_p, sw_m
    real(kind=rp), allocatable :: part(:,:)

    allocate(part(8, size(lag_pts)))
    call idw_interp_shepard_partials(part, lag_pts, lag_el, u, v, w, pmsk, &
         mult, x, y, z, ds, rmax_i, p, lx, ne, adjoint, B, sw_p, sw_m)
    call idw_interp_shepard_finalize(fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, &
         fwm_ib, part, shared_slot, n_shared)
    deallocate(part)

  end subroutine idw_interp_shepard

  !> Turn the local partial sums into the marker values: sum the partials
  !! of markers held by several ranks over their holders (collective when
  !! `shared_slot` and `n_shared` are passed), then normalise. Shared by the
  !! host and the device interpolation paths.
  subroutine idw_interp_shepard_finalize(fu_ib, fv_ib, fw_ib, fum_ib, &
       fvm_ib, fwm_ib, part, shared_slot, n_shared)
    real(kind=rp), intent(inout) :: fu_ib(:), fv_ib(:), fw_ib(:)
    real(kind=rp), intent(inout) :: fum_ib(:), fvm_ib(:), fwm_ib(:)
    real(kind=rp), intent(inout) :: part(:,:)
    integer, intent(in), optional :: shared_slot(:)
    integer, intent(in), optional :: n_shared
    real(kind=rp), allocatable :: shr(:,:)
    integer :: i, n_lag

    n_lag = size(part, 2)

    if (present(shared_slot) .and. present(n_shared)) then
       if (n_shared .gt. 0) then
          allocate(shr(8, n_shared))
          shr = 0.0_rp
          do i = 1, n_lag
             if (shared_slot(i) .gt. 0) then
                shr(:, shared_slot(i)) = part(:, i)
             end if
          end do
          call MPI_Allreduce(MPI_IN_PLACE, shr, 8 * n_shared, &
               MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
          do i = 1, n_lag
             if (shared_slot(i) .gt. 0) then
                part(:, i) = shr(:, shared_slot(i))
             end if
          end do
          deallocate(shr)
       end if
    end if

    call idw_interp_shepard_normalize(fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, &
         fwm_ib, part)

  end subroutine idw_interp_shepard_finalize

  !> Accumulate the per-marker Shepard partial sums over the local stencil
  !! elements. Row layout of `part`: the +side numerators and weight sum
  !! (rows 1-4: u, v, w, weight), then the -side (rows 5-8). Because every
  !! contribution carries the dof multiplicity weight, the partials of one
  !! marker are exactly summable over the ranks holding a copy of it.
  !!
  !! With `adjoint` set the weight of a node is K B mult / sw, where B is
  !! the assembled mass matrix (1/Binv, identical on every copy of a dof)
  !! and sw the assembled spread weight of the node's side (`sw_p` where
  !! pmsk > 0, `sw_m` elsewhere); nodes the spread leaves untouched (sw
  !! below the spread's tolerance) get weight zero. Summed over the stencil
  !! elements this is exactly the transpose of the spread in the mass inner
  !! product, including the markers listed by only some of a dof's elements.
  subroutine idw_interp_shepard_partials(part, lag_pts, lag_el, u, v, w, &
       pmsk, mult, x, y, z, ds, rmax_i, p, lx, ne, adjoint, B, sw_p, sw_m)
    real(kind=rp), intent(inout) :: part(:,:)
    type(point_t), intent(in) :: lag_pts(:)
    type(stack_i4_t), intent(inout) :: lag_el(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: u, v, w
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: pmsk, mult, x, y, z, ds
    real(kind=rp), intent(in) :: rmax_i, p
    logical, intent(in), optional :: adjoint
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in), optional :: B, sw_p, sw_m
    real(kind=rp), parameter :: wtol = 1e-10_rp
    integer :: i, j, k, l, e, ee
    real(kind=rp) :: r, wgt, acc(8)
    real(kind=dp) :: xp, yp, zp
    logical :: use_adjoint

    use_adjoint = .false.
    if (present(adjoint)) use_adjoint = adjoint
    if (use_adjoint) then
       if (.not. (present(B) .and. present(sw_p) .and. present(sw_m))) then
          call neko_error('idw_interp_shepard_partials: adjoint weights &
               &need B, sw_p and sw_m')
       end if
    end if

    ! Every marker owns its column of part, accumulated in registers and
    ! stored once; the stencil sizes differ between markers, hence the
    ! dynamic schedule
    !$omp parallel do private(i, j, k, l, e, ee, r, wgt, acc, xp, yp, zp) &
    !$omp schedule(dynamic)
    do i = 1, size(lag_pts)
       acc = 0.0_rp
       xp = lag_pts(i)%x(1)
       yp = lag_pts(i)%x(2)
       zp = lag_pts(i)%x(3)
       select type (el => lag_el(i)%data)
         type is (integer)
          do ee = 1, lag_el(i)%size()
             e = el(ee)
             do l = 1, lx
                do k = 1, lx
                   do j = 1, lx
                      r = sqrt((x(j,k,l,e) - xp)**2 &
                           + (y(j,k,l,e) - yp)**2 &
                           + (z(j,k,l,e) - zp)**2)
                      r = r / ds(j,k,l,e)
                      wgt = inv_dist_weight(r, rmax_i, p) * mult(j,k,l,e)
                      if (pmsk(j,k,l,e) .gt. 0.0_rp) then
                         if (use_adjoint) then
                            if (abs(sw_p(j,k,l,e)) .gt. wtol) then
                               wgt = wgt * B(j,k,l,e) / sw_p(j,k,l,e)
                            else
                               wgt = 0.0_rp
                            end if
                         end if
                         acc(1) = acc(1) + wgt * u(j,k,l,e)
                         acc(2) = acc(2) + wgt * v(j,k,l,e)
                         acc(3) = acc(3) + wgt * w(j,k,l,e)
                         acc(4) = acc(4) + wgt
                      else
                         if (use_adjoint) then
                            if (abs(sw_m(j,k,l,e)) .gt. wtol) then
                               wgt = wgt * B(j,k,l,e) / sw_m(j,k,l,e)
                            else
                               wgt = 0.0_rp
                            end if
                         end if
                         acc(5) = acc(5) + wgt * u(j,k,l,e)
                         acc(6) = acc(6) + wgt * v(j,k,l,e)
                         acc(7) = acc(7) + wgt * w(j,k,l,e)
                         acc(8) = acc(8) + wgt
                      end if
                   end do
                end do
             end do
          end do
       end select
       part(:, i) = acc
    end do
    !$omp end parallel do

  end subroutine idw_interp_shepard_partials

  !> Turn the (possibly rank-reduced) Shepard partial sums into the
  !! interpolated values; a side whose weight sum is below tolerance
  !! returns zero.
  subroutine idw_interp_shepard_normalize(fu_ib, fv_ib, fw_ib, fum_ib, &
       fvm_ib, fwm_ib, part)
    real(kind=rp), intent(inout) :: fu_ib(:), fv_ib(:), fw_ib(:)
    real(kind=rp), intent(inout) :: fum_ib(:), fvm_ib(:), fwm_ib(:)
    real(kind=rp), intent(in) :: part(:,:)
    integer :: i

    !$omp parallel do private(i)
    do i = 1, size(part, 2)
       if (part(4, i) .gt. 1e-12_rp) then
          fu_ib(i) = part(1, i) / part(4, i)
          fv_ib(i) = part(2, i) / part(4, i)
          fw_ib(i) = part(3, i) / part(4, i)
       else
          fu_ib(i) = 0.0_rp
          fv_ib(i) = 0.0_rp
          fw_ib(i) = 0.0_rp
       end if

       if (part(8, i) .gt. 1e-12_rp) then
          fum_ib(i) = part(5, i) / part(8, i)
          fvm_ib(i) = part(6, i) / part(8, i)
          fwm_ib(i) = part(7, i) / part(8, i)
       else
          fum_ib(i) = 0.0_rp
          fvm_ib(i) = 0.0_rp
          fwm_ib(i) = 0.0_rp
       end if
    end do
    !$omp end parallel do

  end subroutine idw_interp_shepard_normalize

  !> Assign a reduction-buffer slot to every marker held by more than one
  !! rank. `lag_gid` holds the global ids of this rank's markers (their
  !! position in the source boundary meshes); the ids are consistent across
  !! ranks because every rank scans the same boundary meshes. Slots are a
  !! prefix count over the reduced holder array, hence identical on every
  !! rank. Collective over NEKO_COMM.
  subroutine idw_build_shared_slots(lag_gid, n_gid_glb, shared_slot, n_shared)
    integer, intent(in) :: lag_gid(:)
    integer, intent(in) :: n_gid_glb
    integer, intent(out) :: shared_slot(:)
    integer, intent(out) :: n_shared
    integer, allocatable :: holders(:)
    integer :: i, g

    allocate(holders(n_gid_glb))
    holders = 0
    do i = 1, size(lag_gid)
       holders(lag_gid(i)) = 1
    end do
    call MPI_Allreduce(MPI_IN_PLACE, holders, n_gid_glb, MPI_INTEGER, &
         MPI_SUM, NEKO_COMM)

    ! Prefix count: holders becomes the slot per gid (0 = single holder)
    n_shared = 0
    do g = 1, n_gid_glb
       if (holders(g) .gt. 1) then
          n_shared = n_shared + 1
          holders(g) = n_shared
       else
          holders(g) = 0
       end if
    end do

    do i = 1, size(lag_gid)
       shared_slot(i) = holders(lag_gid(i))
    end do
    deallocate(holders)

  end subroutine idw_build_shared_slots

  subroutine idw_init_boundary_mesh(this, lag_pts, lag_nrm, lag_gid, &
       gid_offset, json)
    class(idw_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(stack_pt_t), intent(inout) :: lag_pts
    type(stack_pt_t), intent(inout) :: lag_nrm
    type(stack_i4_t), intent(inout) :: lag_gid
    integer, intent(inout) :: gid_offset
    type(file_t) :: mesh_file
    type(tri_mesh_t) :: boundary_mesh
    character(len=:), allocatable :: mesh_file_name
    character(len=LOG_SIZE) :: log_buf
    integer :: i, j, k, e, n, a, b, m, el_idx, n_sub_max
    integer(kind=8) :: n_markers
    type(stack_i4_t) :: overlaps
    type(point_t) :: tri_nrm, tri_cntr, tri_pt
    type(point_t), pointer :: p1, p2, p3
    real(kind=dp) :: spacing_fac, sb, tb, h_glb
    real(kind=dp) :: box_lo(3), box_hi(3)
    real(kind=dp), allocatable :: h_tri(:)
    integer, allocatable :: n_sub(:)
    logical :: uniform
    character(len=:), allocatable :: refine_mode

    real(kind=dp), dimension(:), allocatable :: box_min, box_max
    real(kind=rp), dimension(3) :: scaling, translation
    type(aabb_t) :: mesh_box, target_box
    character(len=:), allocatable :: mesh_transform
    integer :: idx_p
    logical :: keep_aspect_ratio


    call json_get(json, 'name', mesh_file_name)
    call mesh_file%init(mesh_file_name)
    call neko_log%message('Filename   : '// trim(mesh_file_name))
    call mesh_file%read(boundary_mesh)
    if (boundary_mesh%mpts .lt. 1e1) then
       write(log_buf, '(A, I1)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e2) then
       write(log_buf, '(A, I2)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e3) then
       write(log_buf, '(A, I3)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e4) then
       write(log_buf, '(A, I4)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e5) then
       write(log_buf, '(A, I5)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e6) then
       write(log_buf, '(A, I6)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e7) then
       write(log_buf, '(A, I7)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e8) then
       write(log_buf, '(A, I8)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e9) then
       write(log_buf, '(A, I9)') ' `-Points  : ', boundary_mesh%mpts
    else if (boundary_mesh%mpts .lt. 1e10) then
       write(log_buf, '(A, I10)') ' `-Points  : ', boundary_mesh%mpts
    end if
    call neko_log%message(log_buf)

    if (boundary_mesh%nelv .eq. 0) then
       call neko_error('No elements in the boundary mesh')
    end if

    call json_get_or_default(json, 'mesh_transform.type', &
         mesh_transform, 'none')
    mesh_box = get_aabb(boundary_mesh)

    select case (mesh_transform)
      case ('none')
       ! Do nothing
      case ('bounding_box')
       call json_get(json, 'mesh_transform.box_min', box_min)
       call json_get(json, 'mesh_transform.box_max', box_max)
       call json_get_or_default(json, 'mesh_transform.keep_aspect_ratio', &
            keep_aspect_ratio, .true.)

       if (size(box_min) .ne. 3 .or. size(box_max) .ne. 3) then
          call neko_error('Case file: mesh_transform. &
               &box_min and box_max must be 3 element arrays of reals')
       end if

       call target_box%init(box_min, box_max)


       scaling = target_box%get_diagonal() / mesh_box%get_diagonal()
       if (keep_aspect_ratio) then
          scaling = minval(scaling)
       end if

       translation = - scaling * mesh_box%get_min() + target_box%get_min()

       do idx_p = 1, boundary_mesh%mpts
          boundary_mesh%points(idx_p)%x = &
               scaling * boundary_mesh%points(idx_p)%x + translation
       end do

      case default
       call neko_error('Unknown mesh transform')
    end select

    ! Marker spacing relative to the local grid spacing ds. A triangle
    ! larger than the elements it crosses is refined into sub-triangles
    ! of edge length <= marker_spacing * ds before seeding, one marker
    ! per sub-triangle, so that the spread (reach rmax * ds) leaves no
    ! gaps in the body. A large value reproduces one marker per triangle.
    call json_get_or_default(json, 'marker_spacing', spacing_fac, 1.0_dp)
    if (spacing_fac .le. 0.0_dp) then
       call neko_error('IDW source term: marker_spacing must be positive')
    end if
    ! 'uniform' (default): one target spacing for the whole surface, the
    ! smallest grid spacing under any of its triangles, so the marker
    ! density is the same everywhere. A patch of denser markers acts as a
    ! stiffer piece of wall under the adjoint spread (clustered markers
    ! lump with less cancellation), and the step in wall stiffness sheds
    ! spurious vorticity. 'local': per-triangle spacing, fewer markers.
    call json_get_or_default(json, 'refinement', refine_mode, 'uniform')
    select case (trim(refine_mode))
    case ('uniform')
       uniform = .true.
    case ('local')
       uniform = .false.
    case default
       call neko_error('IDW source term: refinement must be uniform or local')
    end select

    call overlaps%init()

    ! Local grid spacing under each triangle: the smallest ds over the
    ! elements overlapping its centroid or vertices, reduced over all
    ! ranks so that every rank refines identically and the global marker
    ! ids stay aligned
    allocate(h_tri(boundary_mesh%nelv), n_sub(boundary_mesh%nelv))
    h_tri = huge(0.0_dp)
    do i = 1, boundary_mesh%nelv
       tri_cntr = boundary_mesh%el(i)%centroid()
       call this%intersect%overlap(tri_cntr, overlaps)
       do k = 1, 3
          call this%intersect%overlap(boundary_mesh%el(i)%p(k), overlaps)
       end do

       select type (el => overlaps%data)
       type is (integer)
          do j = 1, overlaps%size()
             e = el(j)
             h_tri(i) = min(h_tri(i), &
                  real(minval(this%ds%x(:,:,:,e)), kind=dp))
          end do
       end select
       call overlaps%clear()
    end do

    call MPI_Allreduce(MPI_IN_PLACE, h_tri, boundary_mesh%nelv, &
         MPI_DOUBLE_PRECISION, MPI_MIN, NEKO_COMM)

    if (uniform) then
       h_glb = minval(h_tri)
       where (h_tri .lt. huge(0.0_dp)) h_tri = h_glb
    end if

    n_sub_max = 1
    n_markers = 0
    box_lo = huge(0.0_dp)
    box_hi = -huge(0.0_dp)
    do i = 1, boundary_mesh%nelv
       if (h_tri(i) .lt. huge(0.0_dp)) then
          n_sub(i) = max(1, ceiling(boundary_mesh%el(i)%diameter() &
               / (spacing_fac * h_tri(i))))
       else
          ! Not overlapping the fluid mesh on any rank, dropped below
          n_sub(i) = 1
       end if
       if (n_sub(i) .gt. 1) then
          tri_cntr = boundary_mesh%el(i)%centroid()
          box_lo = min(box_lo, tri_cntr%x)
          box_hi = max(box_hi, tri_cntr%x)
       end if
       n_sub_max = max(n_sub_max, n_sub(i))
       n_markers = n_markers + int(n_sub(i), 8)**2
    end do
    deallocate(h_tri)

    write(log_buf, '(A,I0,A,I0,A,A,A)') ' `-Markers : ', n_markers, &
         ' (max ', n_sub_max, ' per edge, ', trim(refine_mode), ')'
    call neko_log%message(log_buf)
    if (n_sub_max .gt. 1) then
       write(log_buf, '(A,3ES11.3)') ' `-Refined : from ', box_lo
       call neko_log%message(log_buf)
       write(log_buf, '(A,3ES11.3)') ' `-          to   ', box_hi
       call neko_log%message(log_buf)
    end if

    ! Seed one marker per sub-triangle on a barycentric lattice with n
    ! subdivisions per edge: n(n+1)/2 upward and n(n-1)/2 downward
    ! sub-triangles, n**2 in total. n = 1 gives the centroid.
    do i = 1, boundary_mesh%nelv

       n = n_sub(i)
       tri_nrm = boundary_mesh%el(i)%normal()
       p1 => boundary_mesh%el(i)%p(1)
       p2 => boundary_mesh%el(i)%p(2)
       p3 => boundary_mesh%el(i)%p(3)

       m = 0
       do a = 0, n - 1
          do b = 0, n - 1 - a
             ! Upward sub-triangle (a, b), (a+1, b), (a, b+1)
             m = m + 1
             sb = (real(a, dp) + 1.0_dp / 3.0_dp) / real(n, dp)
             tb = (real(b, dp) + 1.0_dp / 3.0_dp) / real(n, dp)
             tri_pt%x = p1%x + sb * (p2%x - p1%x) + tb * (p3%x - p1%x)

             call this%intersect%overlap(tri_pt, overlaps)
             if (overlaps%size() .gt. 0) then
                call lag_pts%push(tri_pt)
                call lag_nrm%push(tri_nrm)
                el_idx = gid_offset + m
                call lag_gid%push(el_idx)
             end if
             call overlaps%clear()

             ! Downward sub-triangle (a+1, b), (a+1, b+1), (a, b+1)
             if (a + b .le. n - 2) then
                m = m + 1
                sb = (real(a, dp) + 2.0_dp / 3.0_dp) / real(n, dp)
                tb = (real(b, dp) + 2.0_dp / 3.0_dp) / real(n, dp)
                tri_pt%x = p1%x + sb * (p2%x - p1%x) + tb * (p3%x - p1%x)

                call this%intersect%overlap(tri_pt, overlaps)
                if (overlaps%size() .gt. 0) then
                   call lag_pts%push(tri_pt)
                   call lag_nrm%push(tri_nrm)
                   el_idx = gid_offset + m
                   call lag_gid%push(el_idx)
                end if
                call overlaps%clear()
             end if
          end do
       end do

       ! Every rank reads the same boundary mesh and refines it the same
       ! way, so the global marker ids stay aligned across ranks and
       ! consecutive objects
       gid_offset = gid_offset + n * n

    end do
    deallocate(n_sub)
    call overlaps%free()



!    call boundary_mesh%free()

  end subroutine idw_init_boundary_mesh

  !> Compute IB weight field. Gathers per element from the CSR transpose
  !! of `lag_el`, so every dof has one writer and the lag points are added
  !! in increasing order, independent of the number of threads.
  subroutine idw_compute_weight(w, wm, msk, lag_pts, active_el, el_off, &
       el_lag, x, y, z, ds, rmax, p, lx, ne)
    type(field_t), intent(inout) :: w, wm, msk
    type(point_t), intent(in) :: lag_pts(:)
    integer, intent(in) :: active_el(:), el_off(:), el_lag(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: x, y, z, ds
    real(kind=rp), intent(inout) :: p
    real(kind=rp), intent(inout) :: rmax
    integer :: a, i, ii, j, k, l, e
    real(kind=rp) :: r
    real(kind=dp) :: xp, yp, zp

    w%x = 0.0_rp
    wm%x = 0.0_rp

    !$omp parallel do private(a, i, ii, j, k, l, e, r, xp, yp, zp) &
    !$omp schedule(dynamic)
    do a = 1, size(active_el)
       e = active_el(a) + 1
       do ii = el_off(e) + 1, el_off(e + 1)
          i = el_lag(ii) + 1
          xp = lag_pts(i)%x(1)
          yp = lag_pts(i)%x(2)
          zp = lag_pts(i)%x(3)
          do l = 1, lx
             do k = 1, lx
                do j = 1, lx
                   r = sqrt((x(j,k,l,e) - xp)**2 &
                        + (y(j,k,l,e) - yp)**2 &
                        + (z(j,k,l,e) - zp)**2)
                   r = r / ds(j,k,l,e)
                   if (msk%x(j,k,l,e) .gt. 0) then
                      w%x(j, k, l, e) = w%x(j, k, l, e) &
                           + inv_dist_weight(r, rmax, p)
                   else
                      wm%x(j, k, l, e) = wm%x(j, k, l, e) &
                           + inv_dist_weight(r, rmax, p)
                   end if
                end do
             end do
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine idw_compute_weight

  !> Spread the per-marker velocities into the IB forcing. Gathers per
  !! element from the CSR transpose of `lag_el`, so every dof has one
  !! writer and the lag points are added in increasing order, i.e. the
  !! result is that of a scatter over the lag points, independent of the
  !! number of threads.
  subroutine idw_spread(ibx, iby, ibz, fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, &
       fwm_ib, lag_pts, active_el, el_off, el_lag, x, y, z, ds, pmsk, w, wm, &
       rmax, p, dt, one_sided, lx, ne)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(inout) :: ibx, iby, ibz
    real(kind=rp), intent(in) :: fu_ib(:), fv_ib(:), fw_ib(:)
    real(kind=rp), intent(in) :: fum_ib(:), fvm_ib(:), fwm_ib(:)
    type(point_t), intent(in) :: lag_pts(:)
    integer, intent(in) :: active_el(:), el_off(:), el_lag(:)
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: x, y, z, ds
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: pmsk, w, wm
    real(kind=rp), intent(in) :: rmax, p
    real(kind=dp), intent(in) :: dt
    logical, intent(in) :: one_sided
    integer :: a, i, ii, j, k, l, e
    real(kind=rp) :: r, idw
    real(kind=dp) :: xp, yp, zp

    if (one_sided) then
       !$omp parallel do private(a, i, ii, j, k, l, e, r, idw, xp, yp, zp) &
       !$omp schedule(dynamic)
       do a = 1, size(active_el)
          e = active_el(a) + 1
          do ii = el_off(e) + 1, el_off(e + 1)
             i = el_lag(ii) + 1
             xp = lag_pts(i)%x(1)
             yp = lag_pts(i)%x(2)
             zp = lag_pts(i)%x(3)
             do l = 1, lx
                do k = 1, lx
                   do j = 1, lx
                      r = sqrt((x(j,k,l,e) - xp)**2 &
                           + (y(j,k,l,e) - yp)**2 &
                           + (z(j,k,l,e) - zp)**2)
                      r = r / ds(j,k,l,e)
                      idw = inv_dist_weight(r, rmax, p)

                      if (pmsk(j,k,l,e) .gt. 0) then
                         if (abs(w(j,k,l,e)) .gt. 1e-10_rp) then
                            ibx(j,k,l,e) = ibx(j,k,l,e) &
                                 + (-fu_ib(i) * idw) / (w(j,k,l,e) * dt)

                            iby(j,k,l,e) = iby(j,k,l,e) &
                                 + (-fv_ib(i) * idw) / (w(j,k,l,e) * dt)

                            ibz(j,k,l,e) = ibz(j,k,l,e) &
                                 + (-fw_ib(i) * idw) / (w(j,k,l,e) * dt)
                         end if
                      else
                         if (abs(wm(j,k,l,e)) .gt. 1e-10_rp) then
                            ibx(j,k,l,e) = ibx(j,k,l,e) &
                                 + (-fum_ib(i) * idw) / (wm(j,k,l,e) * dt)

                            iby(j,k,l,e) = iby(j,k,l,e) &
                                 + (-fvm_ib(i) * idw) / (wm(j,k,l,e) * dt)

                            ibz(j,k,l,e) = ibz(j,k,l,e) &
                                 + (-fwm_ib(i) * idw) / (wm(j,k,l,e) * dt)
                         end if
                      end if
                   end do
                end do
             end do
          end do
       end do
       !$omp end parallel do
    else
       !$omp parallel do private(a, i, ii, j, k, l, e, r, idw, xp, yp, zp) &
       !$omp schedule(dynamic)
       do a = 1, size(active_el)
          e = active_el(a) + 1
          do ii = el_off(e) + 1, el_off(e + 1)
             i = el_lag(ii) + 1
             xp = lag_pts(i)%x(1)
             yp = lag_pts(i)%x(2)
             zp = lag_pts(i)%x(3)
             do l = 1, lx
                do k = 1, lx
                   do j = 1, lx
                      if (w(j,k,l,e) .gt. 1e-12_rp) then
                         r = sqrt((x(j,k,l,e) - xp)**2 &
                              + (y(j,k,l,e) - yp)**2 &
                              + (z(j,k,l,e) - zp)**2)
                         r = r / ds(j,k,l,e)
                         idw = inv_dist_weight(r, rmax, p)

                         ibx(j,k,l,e) = ibx(j,k,l,e) &
                              + (-fu_ib(i) * idw) / (w(j,k,l,e) * dt)

                         iby(j,k,l,e) = iby(j,k,l,e) &
                              + (-fv_ib(i) * idw) / (w(j,k,l,e) * dt)

                         ibz(j,k,l,e) = ibz(j,k,l,e) &
                              + (-fw_ib(i) * idw) / (w(j,k,l,e) * dt)
                      end if
                   end do
                end do
             end do
          end do
       end do
       !$omp end parallel do
    end if

  end subroutine idw_spread

  !> Count, for every marker and side, the local stencil nodes inside the
  !! interpolation cutoff. Used for the init-time diagnostic of the
  !! Shepard interpolation (spec: sparse read stencils do not
  !! self-attenuate, so surface them loudly). The counts of markers held
  !! by several ranks are summed over the holders by the caller.
  subroutine idw_interp_stencil_counts(lag_pts, lag_el, pmsk, x, y, z, ds, &
       rmax_i, lx, ne, np_i, nm_i)
    type(point_t), intent(in) :: lag_pts(:)
    type(stack_i4_t), intent(inout) :: lag_el(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: pmsk, x, y, z, ds
    real(kind=rp), intent(in) :: rmax_i
    integer, intent(out) :: np_i(:), nm_i(:)
    integer :: i, j, k, l, e, ee, np, nm
    real(kind=rp) :: r

    !$omp parallel do private(i, j, k, l, e, ee, np, nm, r) schedule(dynamic)
    do i = 1, size(lag_pts)
       np = 0
       nm = 0
       select type (el => lag_el(i)%data)
         type is (integer)
          do ee = 1, lag_el(i)%size()
             e = el(ee)
             do l = 1, lx
                do k = 1, lx
                   do j = 1, lx
                      r = sqrt((x(j,k,l,e) - lag_pts(i)%x(1))**2 &
                           + (y(j,k,l,e) - lag_pts(i)%x(2))**2 &
                           + (z(j,k,l,e) - lag_pts(i)%x(3))**2)
                      if (r / ds(j,k,l,e) .lt. rmax_i) then
                         if (pmsk(j,k,l,e) .gt. 0.0_rp) then
                            np = np + 1
                         else
                            nm = nm + 1
                         end if
                      end if
                   end do
                end do
             end do
          end do
       end select
       np_i(i) = np
       nm_i(i) = nm
    end do
    !$omp end parallel do

  end subroutine idw_interp_stencil_counts

  !> Compute IB mask fields. Gathers per element from the CSR transpose
  !! of `lag_el`: the nearest-marker distance of an element's dofs only
  !! depends on the markers overlapping it, so each thread keeps it in a
  !! private element-sized buffer.
  subroutine idw_compute_mask(mmsk, pmsk, lag_pts, lag_nrm, active_el, &
       el_off, el_lag, x, y, z, ds, band, lx, ne)
    type(field_t), intent(inout) :: mmsk, pmsk
    type(point_t), intent(in) :: lag_pts(:)
    type(point_t), intent(in) :: lag_nrm(:)
    integer, intent(in) :: active_el(:), el_off(:), el_lag(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: x, y, z, ds
    real(kind=rp), intent(in) :: band
    real(kind=rp) :: dist(lx,lx,lx)
    real(kind=rp) :: euler_pt(3)
    integer :: a, i, ii, j, k, l, e
    real(kind=rp) :: r, dn, nn
    real(kind=dp) :: xp, yp, zp

    mmsk%x = 1.0_rp
    pmsk%x = 1.0_rp

    !$omp parallel do private(a, i, ii, j, k, l, e, r, xp, yp, zp, &
    !$omp dist, euler_pt) schedule(dynamic)
    do a = 1, size(active_el)
       e = active_el(a) + 1
       dist = huge(0.0_rp)

       do ii = el_off(e) + 1, el_off(e + 1)
          i = el_lag(ii) + 1
          xp = lag_pts(i)%x(1)
          yp = lag_pts(i)%x(2)
          zp = lag_pts(i)%x(3)
          do l = 1, lx
             do k = 1, lx
                do j = 1, lx
                   r = sqrt((x(j,k,l,e) - xp)**2 &
                        + (y(j,k,l,e) - yp)**2 &
                        + (z(j,k,l,e) - zp)**2)
                   dist(j,k,l) = min(dist(j,k,l), r)
                end do
             end do
          end do
       end do

       do ii = el_off(e) + 1, el_off(e + 1)
          i = el_lag(ii) + 1
          xp = lag_pts(i)%x(1)
          yp = lag_pts(i)%x(2)
          zp = lag_pts(i)%x(3)
          do l = 1, lx
             do k = 1, lx
                do j = 1, lx
                   r = sqrt((x(j,k,l,e) - xp)**2 &
                        + (y(j,k,l,e) - yp)**2 &
                        + (z(j,k,l,e) - zp)**2)

                   if (r .le. dist(j,k,l)) then

                      euler_pt(1) = (x(j,k,l,e) - xp)
                      euler_pt(2) = (y(j,k,l,e) - yp)
                      euler_pt(3) = (z(j,k,l,e) - zp)

                      ! Signed distance to the marker's tangent plane; a
                      ! node within the band is on the surface and stays
                      ! on both sides
                      nn = sqrt(sum(lag_nrm(i)%x**2))
                      dn = sum(euler_pt * lag_nrm(i)%x) / max(nn, tiny(nn))
                      if (abs(dn) .le. band * ds(j,k,l,e)) then
                         continue
                      else if (dn .gt. 0) then
                         mmsk%x(j,k,l,e) = 0.0_rp
                      else
                         pmsk%x(j,k,l,e) = 0.0_rp
                      end if
                   end if
                end do
             end do
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine idw_compute_mask

  !> Inverse distance weighting coefficient
  !! @param r Radial distance to Lagrangian point.
  !! @param rmax Radial distance for sphere of influence.
  pure function inv_dist_weight(r, rmax, p) result(idw)
    real(kind=rp), intent(in) :: r
    real(kind=rp), intent(in) :: rmax
    real(kind=rp), intent(in) :: p
    real(kind=rp) :: idw

    if(r .ge. rmax) then
       idw = 0.0_rp
    else
       idw = ((rmax-r)/(rmax * r + NEKO_EPS))**p
    end if

  end function inv_dist_weight

end module idw_source_term
