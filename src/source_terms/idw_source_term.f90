! Copyright (c) 2023-2024, The Neko Authors
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
  use field_math, only : field_col2, field_col3, field_copy
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
  use math, only : NEKO_EPS
  use aabb, only : aabb_t, get_aabb
  use time_state, only : time_state_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use elementwise_filter, only : elementwise_filter_t
  use PDE_filter, only : PDE_filter_t
  use filter, only : filter_t
  use device_math, only : device_col2
  use device, only : device_free, device_map, device_memcpy, &
       HOST_TO_DEVICE, DEVICE_TO_HOST
  use device_idw_source_term, only : device_idw_gather
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
     type(gs_t) :: gs
     logical :: one_sided
     !> Use Shepard IDW interpolation onto the markers instead of the
     !! global (rst-based) interpolation
     logical :: idw_interp = .false.
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
  end type idw_source_term_t

  public :: idw_interp_shepard, idw_interp_shepard_partials, &
       idw_interp_shepard_normalize, idw_build_shared_slots

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

    ! Mandatory fields for the general source term
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)
    
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
         "barycentric")
    select case (trim(interp_scheme))
    case ("barycentric")
       this%idw_interp = .false.
    case ("idw")
       this%idw_interp = .true.
    case default
       call neko_error('IDW source term unknown interpolation scheme: ' &
            // trim(interp_scheme))
    end select
    call neko_log%message('Interp     : '// trim(interp_scheme))

    call json_get_or_default(json, "interpolation_rmax", this%interp_rmax, &
         2.0_rp)
    if (this%idw_interp) then
       write(log_buf, '(A,f5.2)') 'Interp rmax: ', this%interp_rmax
       call neko_log%message(log_buf)
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
        call this%fltr%init(json, coef)
     end if


    ! Naive apporach to find the smallest distance between two dofs in the mesh

    ds_min = huge(0.0_rp)
    ds_max = -huge(0.0_rp)

    call this%ds%init(coef%dof)

    associate (x => coef%dof%x, y => coef%dof%y, z => coef%dof%z, &
         lx => coef%Xh%lx, ds => this%ds%x)

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

    ! Report total number of lagrangian points generated, this might differ
    ! from the STL sources due to refinement (not implemented yet!)
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

    ! Construct weight field
    call this%w%init(coef%dof, "ib_weight")
    call this%wm%init(coef%dof, "ib_mweight")

    call this%gs%init(coef%dof)

    this%ds%x = this%ds%x * coef%mult

    call this%gs%op(this%ds, GS_OP_ADD)


    this%ds%x = this%ds%x * coef%mult
    call this%gs%op(this%ds, GS_OP_ADD)

    call this%mmsk%init(coef%dof, "ib_mmask")
    call this%pmsk%init(coef%dof, "ib_pmask")
    call this%tmp%init(coef%dof, "ib_tmp")

    if (this%one_sided) then
       call idw_compute_mask(this%mmsk, this%pmsk, this%lag_pts, this%lag_el, &
            this%lag_nrm, coef%dof%x, coef%dof%y, coef%dof%z, &
            coef%Xh%lx, coef%msh%nelv)
    else
       this%mmsk = 0.0_rp
       this%pmsk = 1.0_rp
    end if


    this%pmsk%x = this%pmsk%x * coef%mult

    call this%gs%op(this%pmsk, GS_OP_ADD)

    this%mmsk%x = this%mmsk%x * coef%mult

    call this%gs%op(this%mmsk, GS_OP_ADD)


    call idw_compute_weight(this%w, this%wm, this%pmsk, this%lag_pts, this%lag_el, &
         coef%dof%x, coef%dof%y, coef%dof%z, this%ds%x, this%rmax, &
         this%pwr_param, coef%Xh%lx,coef%msh%nelv)

    this%w%x = this%w%x * coef%mult

    call this%gs%op(this%w, GS_OP_ADD)

    this%wm%x = this%wm%x * coef%mult

    call this%gs%op(this%wm, GS_OP_ADD)

    if (this%idw_interp) then
       ! Per-marker local stencil counts; the counts of markers held by
       ! several ranks are then summed over their holders so the diagnostic
       ! reflects the stencil the reduced interpolation actually sees.
       allocate(np_i(n_lags), nm_i(n_lags))
       call idw_interp_stencil_counts(this%lag_pts, this%lag_el, &
            this%pmsk%x, coef%dof%x, coef%dof%y, coef%dof%z, this%ds%x, &
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
          dxe = maxval(coef%dof%x(:,:,:,e)) - minval(coef%dof%x(:,:,:,e))
          dye = maxval(coef%dof%y(:,:,:,e)) - minval(coef%dof%y(:,:,:,e))
          dze = maxval(coef%dof%z(:,:,:,e)) - minval(coef%dof%z(:,:,:,e))
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

  !> Build the device data structures for the IDW source term: the CSR
  !! transpose of `lag_el` (per element -> lag points) used by the gather
  !! kernel, the unpacked Lagrangian coordinates, and upload the mask/weight
  !! fields that are assembled on the host during init.
  subroutine idw_build_device_maps(this)
    class(idw_source_term_t), intent(inout) :: this
    integer :: nelv, n_lag, n, i, ee, e, k, n_csr
    integer, allocatable :: cursor(:)

    nelv = size(this%w%dof%x, 4)
    n_lag = size(this%lag_pts)
    this%lx3 = this%w%Xh%lx**3

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
    end if

    ! Upload the mask / weight fields assembled on the host during init
    n = this%w%size()
    call device_memcpy(this%w%x, this%w%x_d, n, HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%wm%x, this%wm%x_d, n, HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%ds%x, this%ds%x_d, n, HOST_TO_DEVICE, sync = .false.)
    call device_memcpy(this%pmsk%x, this%pmsk%x_d, n, HOST_TO_DEVICE, &
         sync = .false.)
    call device_memcpy(this%mmsk%x, this%mmsk%x_d, n, HOST_TO_DEVICE, &
         sync = .true.)

  end subroutine idw_build_device_maps

  !> Interpolate a masked field at the Lagrangian points on the host. The
  !! masked product `field_col3` runs on the device, so the result is synced
  !! back before the (host) gslib evaluation.
  subroutine idw_interp_masked(global_interp, tmp, fld, msk, ib, nt)
    type(global_interpolation_t), intent(inout) :: global_interp
    type(field_t), intent(inout) :: tmp
    type(field_t), intent(in) :: fld, msk
    real(kind=rp), intent(inout) :: ib(:)
    integer, intent(in) :: nt

    call field_col3(tmp, fld, msk, nt)
    if (NEKO_BCKND_DEVICE .eq. 1) &
         call device_memcpy(tmp%x, tmp%x_d, nt, DEVICE_TO_HOST, sync = .true.)
    call global_interp%evaluate(ib, tmp%x, .true.)
  end subroutine idw_interp_masked

  !> Device path of the source term: interpolate on the host (gslib), upload
  !! the per-point values, then assemble the contributions with the atomic-free
  !! gather kernel. Mirrors the host branches of idw_source_term_compute.
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
       ! Shepard interpolation runs on the host: refresh u,v,w host mirrors
       call device_memcpy(u%x, u%x_d, u%size(), DEVICE_TO_HOST, &
            sync = .false.)
       call device_memcpy(v%x, v%x_d, v%size(), DEVICE_TO_HOST, &
            sync = .false.)
       call device_memcpy(w%x, w%x_d, w%size(), DEVICE_TO_HOST, &
            sync = .true.)
       call idw_interp_shepard(this%fu_ib, this%fv_ib, this%fw_ib, &
            this%fum_ib, this%fvm_ib, this%fwm_ib, this%lag_pts, &
            this%lag_el, u%x, v%x, w%x, this%pmsk%x, this%coef%mult, &
            this%w%dof%x, this%w%dof%y, this%w%dof%z, this%ds%x, &
            this%interp_rmax, this%pwr_param, this%w%Xh%lx, &
            size(this%w%dof%x, 4), this%shared_slot, this%n_shared_glb)
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
       ! Unmasked interpolation reads the host velocity directly; refresh it.
       call device_memcpy(u%x, u%x_d, u%size(), DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(v%x, v%x_d, v%size(), DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(w%x, w%x_d, w%size(), DEVICE_TO_HOST, sync = .true.)
       call this%global_interp%evaluate(this%fu_ib, u%x, .true.)
       call this%global_interp%evaluate(this%fv_ib, v%x, .true.)
       call this%global_interp%evaluate(this%fw_ib, w%x, .true.)
    end if

    ! Upload the interpolated per-point values (last copy is synchronous so the
    ! kernel, enqueued on the same queue afterwards, sees the data).
    if (n_lag > 0) then
       call device_memcpy(this%fu_ib, this%fu_ib_d, n_lag, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%fv_ib, this%fv_ib_d, n_lag, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%fw_ib, this%fw_ib_d, n_lag, HOST_TO_DEVICE, &
            sync = .not. (this%one_sided .or. this%idw_interp))
       if (this%one_sided .or. this%idw_interp) then
          call device_memcpy(this%fum_ib, this%fum_ib_d, n_lag, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%fvm_ib, this%fvm_ib_d, n_lag, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%fwm_ib, this%fwm_ib_d, n_lag, &
               HOST_TO_DEVICE, sync = .true.)
       end if
    end if

    ! one_sided uses the pmsk split with tol 1e-10; the unmasked path has
    ! pmsk == 1 everywhere and matches the host tol of 1e-12.
    wtol = merge(1.0e-10_rp, 1.0e-12_rp, this%one_sided)

    call device_idw_gather(fu%x_d, fv%x_d, fw%x_d, &
         this%fu_ib_d, this%fv_ib_d, this%fw_ib_d, &
         this%fum_ib_d, this%fvm_ib_d, this%fwm_ib_d, &
         this%w%dof%x_d, this%w%dof%y_d, this%w%dof%z_d, this%ds%x_d, &
         this%pmsk%x_d, this%w%x_d, this%wm%x_d, &
         this%lpx_d, this%lpy_d, this%lpz_d, &
         this%active_el_d, this%el_off_d, this%el_lag_d, &
         this%n_active, this%lx3, time%dt, this%rmax, this%pwr_param, &
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

    call this%gs%free()

  end subroutine idw_source_term_free

  subroutine idw_source_term_compute(this, time)
    class(idw_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, fu, fv, fw
    integer :: i, j, k, l, e, ee, n
    real(kind=rp) :: r, idw

    n = this%fields%item_size(1)

    u => neko_registry%get_field('u')
    v => neko_registry%get_field('v')
    w => neko_registry%get_field('w')

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)

    associate(global_interp => this%global_interp, &
         fu_ib => this%fu_ib, fv_ib => this%fv_ib, fw_ib => this%fw_ib, &
         fum_ib => this%fum_ib, fvm_ib => this%fvm_ib, fwm_ib => this%fwm_ib, &
         lag_pts => this%lag_pts, tmp => this%tmp, &
         ds => this%ds%x, x => this%w%dof%x, &
         y => this%w%dof%y, z => this%w%dof%z, lx => this%w%Xh%lx)


      fu_ib = 0.0_rp
      fv_ib = 0.0_rp
      fw_ib = 0.0_rp

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call idw_compute_device(this, fu, fv, fw, u, v, w, time)
      else
         ! Interpolation stage: velocity at the Lagrangian points
         if (this%idw_interp) then
            call idw_interp_shepard(fu_ib, fv_ib, fw_ib, &
                 fum_ib, fvm_ib, fwm_ib, lag_pts, this%lag_el, &
                 u%x, v%x, w%x, this%pmsk%x, this%coef%mult, x, y, z, ds, &
                 this%interp_rmax, this%pwr_param, lx, this%coef%msh%nelv, &
                 this%shared_slot, this%n_shared_glb)
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

         ! Spread stage (unchanged)
         if (this%one_sided) then
         do i = 1, size(this%lag_pts)
            select type (el => this%lag_el(i)%data)
              type is (integer)
               do ee = 1, this%lag_el(i)%size()
                  e = el(ee)
                  do l = 1, lx
                     do k = 1, lx
                        do j = 1, lx
                           r = sqrt((x(j,k,l,e) - lag_pts(i)%x(1))**2 &
                                + (y(j,k,l,e) - lag_pts(i)%x(2))**2 &
                                + (z(j,k,l,e) - lag_pts(i)%x(3))**2)
                           r = r / ds(j,k,l,e)
                           idw = inv_dist_weight(r, this%rmax, this%pwr_param)

                           if (this%pmsk%x(j,k,l,e) .gt. 0) then
                              if (abs(this%w%x(j,k,l,e)) .gt. 1e-10_rp) then
                                 fu%x(j,k,l,e) = fu%x(j,k,l,e) &
                                      + (-fu_ib(i) * idw) / (this%w%x(j,k,l,e) &
                                      * time%dt)

                                 fv%x(j,k,l,e) = fv%x(j,k,l,e) &
                                      + (-fv_ib(i) * idw) / (this%w%x(j,k,l,e) &
                                      * time%dt)

                                 fw%x(j,k,l,e) = fw%x(j,k,l,e) &
                                      + (-fw_ib(i) * idw) / (this%w%x(j,k,l,e) &
                                      * time%dt)
                              end if
                           else
                              if (abs(this%wm%x(j,k,l,e)) .gt. 1e-10_rp) then
                                 fu%x(j,k,l,e) = fu%x(j,k,l,e) &
                                      + (-fum_ib(i) * idw) / (this%wm%x(j,k,l,e) &
                                      * time%dt)

                                 fv%x(j,k,l,e) = fv%x(j,k,l,e) &
                                      + (-fvm_ib(i) * idw) / (this%wm%x(j,k,l,e) &
                                      * time%dt)

                                 fw%x(j,k,l,e) = fw%x(j,k,l,e) &
                                      + (-fwm_ib(i) * idw) / (this%wm%x(j,k,l,e) &
                                      * time%dt)
                              end if
                           end if
                        end do
                     end do
                  end do
               end do
            end select
         end do
         else
         do i = 1, size(this%lag_pts)
            select type (el => this%lag_el(i)%data)
              type is (integer)
               do ee = 1, this%lag_el(i)%size()
                  e = el(ee)
                  do l = 1, lx
                     do k = 1, lx
                        do j = 1, lx
                           if (this%w%x(j,k,l,e) .gt. 1e-12_rp) then
                              r = sqrt((x(j,k,l,e) - lag_pts(i)%x(1))**2 &
                                   + (y(j,k,l,e) - lag_pts(i)%x(2))**2 &
                                   + (z(j,k,l,e) - lag_pts(i)%x(3))**2)
                              r = r / ds(j,k,l,e)
                              idw = inv_dist_weight(r, this%rmax, this%pwr_param)
                              
                              fu%x(j,k,l,e) = fu%x(j,k,l,e) &
                                   + (-fu_ib(i) * idw) / (this%w%x(j,k,l,e) &
                                   * time%dt)

                              fv%x(j,k,l,e) = fv%x(j,k,l,e) &
                                   + (-fv_ib(i) * idw) / (this%w%x(j,k,l,e) &
                                   * time%dt)

                              fw%x(j,k,l,e) = fw%x(j,k,l,e) &
                                   + (-fw_ib(i) * idw) / (this%w%x(j,k,l,e) &
                                   * time%dt)
                           end if
                        end do
                     end do
                  end do
               end do
            end select
         end do
         end if
      end if

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(fu%x_d, this%coef%mult_d, fu%size())
         call device_col2(fv%x_d, this%coef%mult_d, fu%size())
         call device_col2(fw%x_d, this%coef%mult_d, fu%size())         
      else
         do i = 1, n
            fu%x(i,1,1,1) = fu%x(i,1,1,1) * this%coef%mult(i,1,1,1)
            fv%x(i,1,1,1) = fv%x(i,1,1,1) * this%coef%mult(i,1,1,1)
            fw%x(i,1,1,1) = fw%x(i,1,1,1) * this%coef%mult(i,1,1,1)
         end do
      end if
      
      call this%gs%op(fu, GS_OP_ADD)
      call this%gs%op(fv, GS_OP_ADD)
      call this%gs%op(fw, GS_OP_ADD)

      if (allocated(this%fltr)) then
         call field_copy(tmp, fu)
         call this%fltr%apply(fu, tmp)
         
         call field_copy(tmp, fv)
         call this%fltr%apply(fv, tmp)
         
         call field_copy(tmp, fw)
         call this%fltr%apply(fw, tmp)
      end if

    end associate

  end subroutine idw_source_term_compute

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
  subroutine idw_interp_shepard(fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, fwm_ib, &
       lag_pts, lag_el, u, v, w, pmsk, mult, x, y, z, ds, rmax_i, p, lx, ne, &
       shared_slot, n_shared)
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
    real(kind=rp), allocatable :: part(:,:), shr(:,:)
    integer :: i

    allocate(part(8, size(lag_pts)))
    call idw_interp_shepard_partials(part, lag_pts, lag_el, u, v, w, pmsk, &
         mult, x, y, z, ds, rmax_i, p, lx, ne)

    if (present(shared_slot) .and. present(n_shared)) then
       if (n_shared .gt. 0) then
          allocate(shr(8, n_shared))
          shr = 0.0_rp
          do i = 1, size(lag_pts)
             if (shared_slot(i) .gt. 0) then
                shr(:, shared_slot(i)) = part(:, i)
             end if
          end do
          call MPI_Allreduce(MPI_IN_PLACE, shr, 8 * n_shared, &
               MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
          do i = 1, size(lag_pts)
             if (shared_slot(i) .gt. 0) then
                part(:, i) = shr(:, shared_slot(i))
             end if
          end do
          deallocate(shr)
       end if
    end if

    call idw_interp_shepard_normalize(fu_ib, fv_ib, fw_ib, fum_ib, fvm_ib, &
         fwm_ib, part)
    deallocate(part)

  end subroutine idw_interp_shepard

  !> Accumulate the per-marker Shepard partial sums over the local stencil
  !! elements. Row layout of `part`: the +side numerators and weight sum
  !! (rows 1-4: u, v, w, weight), then the -side (rows 5-8). Because every
  !! contribution carries the dof multiplicity weight, the partials of one
  !! marker are exactly summable over the ranks holding a copy of it.
  subroutine idw_interp_shepard_partials(part, lag_pts, lag_el, u, v, w, &
       pmsk, mult, x, y, z, ds, rmax_i, p, lx, ne)
    real(kind=rp), intent(inout) :: part(:,:)
    type(point_t), intent(in) :: lag_pts(:)
    type(stack_i4_t), intent(inout) :: lag_el(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: u, v, w
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: pmsk, mult, x, y, z, ds
    real(kind=rp), intent(in) :: rmax_i, p
    integer :: i, j, k, l, e, ee
    real(kind=rp) :: r, wgt

    do i = 1, size(lag_pts)
       part(:, i) = 0.0_rp
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
                      r = r / ds(j,k,l,e)
                      wgt = inv_dist_weight(r, rmax_i, p) * mult(j,k,l,e)
                      if (pmsk(j,k,l,e) .gt. 0.0_rp) then
                         part(1, i) = part(1, i) + wgt * u(j,k,l,e)
                         part(2, i) = part(2, i) + wgt * v(j,k,l,e)
                         part(3, i) = part(3, i) + wgt * w(j,k,l,e)
                         part(4, i) = part(4, i) + wgt
                      else
                         part(5, i) = part(5, i) + wgt * u(j,k,l,e)
                         part(6, i) = part(6, i) + wgt * v(j,k,l,e)
                         part(7, i) = part(7, i) + wgt * w(j,k,l,e)
                         part(8, i) = part(8, i) + wgt
                      end if
                   end do
                end do
             end do
          end do
       end select
    end do

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
    integer :: i, j, el_idx
    type(stack_i4_t) :: overlaps
    type(point_t) :: tri_nrm, tri_cntr

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

    call overlaps%init()

    do i = 1, boundary_mesh%nelv

       tri_cntr = boundary_mesh%el(i)%centroid()
       tri_nrm = boundary_mesh%el(i)%normal()


       call this%intersect%overlap(tri_cntr, overlaps)

       if (overlaps%size() .gt. 0) then

          call lag_pts%push(tri_cntr)
          call lag_nrm%push(tri_nrm)
          el_idx = gid_offset + i
          call lag_gid%push(el_idx)

!          do while (.not. overlaps%is_empty())
!             el_idx = overlaps%pop()
!          end do
       end if
       call overlaps%clear()

    end do

    ! Every rank reads the same boundary mesh, so the global marker ids
    ! stay aligned across ranks and consecutive objects
    gid_offset = gid_offset + boundary_mesh%nelv


!    call boundary_mesh%free()

  end subroutine idw_init_boundary_mesh

  !> Compute IB weight field
  subroutine idw_compute_weight(w, wm, msk, lag_pts, lag_el, x, y, z, &
       ds, rmax, p, lx, ne)
    type(field_t), intent(inout) :: w, wm, msk
    type(point_t), allocatable, intent(inout) :: lag_pts(:)
    type(stack_i4_t), allocatable, intent(inout) :: lag_el(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: x, y, z, ds
    real(kind=rp), intent(inout) :: p
    real(kind=rp), intent(inout) :: rmax
    integer :: i, j, k, l, e, ee
    real(kind=rp) :: r

    w%x = 0.0_rp
    wm%x = 0.0_rp

    do i = 1, size(lag_pts)
       select type(el => lag_el(i)%data)
         type is (integer)
          do ee = 1, lag_el(i)%size()
             e = el(ee)
             do l = 1, lx
                do k = 1, lx
                   do j = 1, lx
                      r = sqrt((x(j,k,l,e) - lag_pts(i)%x(1))**2 &
                           + (y(j,k,l,e) - lag_pts(i)%x(2))**2 &
                           + (z(j,k,l,e) - lag_pts(i)%x(3))**2)
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
       end select
    end do

  end subroutine idw_compute_weight

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

  end subroutine idw_interp_stencil_counts

  !> Compute IB mask fields
  subroutine idw_compute_mask(mmsk, pmsk, lag_pts, lag_el, lag_nrm, x, y, z, lx, ne)
    type(field_t), intent(inout) :: mmsk, pmsk
    type(point_t), allocatable, intent(inout) :: lag_pts(:)
    type(stack_i4_t), allocatable, intent(inout) :: lag_el(:)
    type(point_t), allocatable, intent(inout) :: lag_nrm(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne), intent(in) :: x, y, z
    real(kind=rp), allocatable :: dist(:,:,:,:)
    real(kind=rp) :: euler_pt(3)
    integer :: i, j, k, l, e, ee
    real(kind=rp) :: r

    mmsk%x = 1.0_rp
    pmsk%x = 1.0_rp

    allocate(dist(lx,lx,lx,ne))
    dist = huge(0.0_rp)

    do i = 1, size(lag_pts)
       select type(el => lag_el(i)%data)
         type is (integer)
          do ee = 1, lag_el(i)%size()
             e = el(ee)
             do l = 1, lx
                do k = 1, lx
                   do j = 1, lx
                      r = sqrt((x(j,k,l,e) - lag_pts(i)%x(1))**2 &
                           + (y(j,k,l,e) - lag_pts(i)%x(2))**2 &
                           + (z(j,k,l,e) - lag_pts(i)%x(3))**2)
                      dist(j,k,l,e) = min(dist(j,k,l,e), r)
                   end do
                end do
             end do
          end do
       end select
    end do

    do i = 1, size(lag_pts)
       select type(el => lag_el(i)%data)
         type is (integer)
          do ee = 1, lag_el(i)%size()
             e = el(ee)
             do l = 1, lx
                do k = 1, lx
                   do j = 1, lx
                      r = sqrt((x(j,k,l,e) - lag_pts(i)%x(1))**2 &
                           + (y(j,k,l,e) - lag_pts(i)%x(2))**2 &
                           + (z(j,k,l,e) - lag_pts(i)%x(3))**2)

                      if (r .le. dist(j,k,l,e)) then

                         euler_pt(1) = (x(j,k,l,e) - lag_pts(i)%x(1))
                         euler_pt(2) = (y(j,k,l,e) - lag_pts(i)%x(2))
                         euler_pt(3) = (z(j,k,l,e) - lag_pts(i)%x(3))

                         if (sum(euler_pt * lag_nrm(i)%x) .gt. 0) then
                            mmsk%x(j,k,l,e) = 0.0_rp
                         else
                            pmsk%x(j,k,l,e) = 0.0_rp
                         end if
                      end if
                   end do
                end do
             end do
          end do
       end select
    end do

    deallocate(dist)

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
