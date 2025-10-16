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
  use coefs, only : coef_t
  use field, only : field_t
  use utils, only : neko_error
  use tri_mesh, only : tri_mesh_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use file, only : file_t
  use comm, only : MPI_REAL_PRECISION, NEKO_COMM
  use intersection_detector, only : intersect_detector_t
  use mesh, only : mesh_t
  use stack, only : stack_i4_t, stack_pt_t
  use global_interpolation, only : global_interpolation_t
  use point, only : point_t
  use math, only : NEKO_EPS
  use aabb, only : aabb_t, get_aabb
  use time_state, only : time_state_t
  use gather_scatter
  use mpi_f08
  use logger
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
     real(kind=dp), allocatable :: xyz(:,:)
     real(kind=rp), allocatable :: fu_ib(:)
     real(kind=rp), allocatable :: fv_ib(:)
     real(kind=rp), allocatable :: fw_ib(:)
     real(kind=rp), allocatable :: fum_ib(:)
     real(kind=rp), allocatable :: fvm_ib(:)
     real(kind=rp), allocatable :: fwm_ib(:)
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
    character(len=LOG_SIZE) :: log_buf
    real(kind=dp) :: aabb_padding,dx_max, dy_max, dz_max, ds_max, ds_min
    type(stack_i4_t) :: overlaps
    type(stack_pt_t) :: lagrangian_points
    type(stack_pt_t) :: lagrangian_normals
    real(kind=rp) :: start_time, end_time

    ! Mandatory fields for the general source term
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    call neko_log%section('Inverse distance weighting')

    call json_get_or_default(json, "rmax", this%rmax, 1.0_rp)
    write(log_buf, '(A,f5.2)') 'Rmax       : ', this%rmax
    call neko_log%message(log_buf)

    call json_get_or_default(json, "padding", aabb_padding, 0.125_rp)
    write(log_buf, '(A,f5.2)') 'Padding    : ', aabb_padding
    call neko_log%message(log_buf)
    
    call json_get_or_default(json, "power_parameter", this%pwr_param, 0.5_rp)
    write(log_buf, '(A,f5.2)') 'IDW Power  : ', this%pwr_param
    call neko_log%message(log_buf)
    call json_get_or_default(json, "one_sided", this%one_sided, .false.)
    write(log_buf, '(A,L1)') 'One sided  : ', this%one_sided
    call neko_log%message(log_buf)

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
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM)
    write(log_buf, '(A,ES13.6)') 'Minimum ds :', this%ds_min
    call neko_log%message(log_buf)

    call MPI_Allreduce(MPI_IN_PLACE, this%ds_max, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
    write(log_buf, '(A,ES13.6)') 'Maximum ds :', this%ds_max
    call neko_log%message(log_buf)


    call this%intersect%init(coef%msh, aabb_padding)
    call lagrangian_points%init()
    call lagrangian_normals%init()

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
               object_settings)
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

    call this%global_interp%init(coef%dof, NEKO_COMM)

    n_lags = lagrangian_points%size()
    call this%global_interp%find_points_and_redist(this%xyz, n_lags)

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
       this%mmsk%x = 0.0_rp
       this%pmsk%x = 1.0_rp
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

    call lagrangian_points%free()
    call lagrangian_normals%free()

  end subroutine idw_source_term_init_from_json

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

    if (allocated(this%lag_pts)) then
       deallocate(this%lag_pts)
    end if

    if (allocated(this%lag_el)) then
       do i = 1, size(this%lag_el)
          call this%lag_el(i)%free()
       end do
       deallocate(this%lag_el)
    end if

    call this%gs%free()

  end subroutine idw_source_term_free

  subroutine idw_source_term_compute(this, time)
    class(idw_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, fu, fv, fw
    integer :: i, j, k, l, e, ee, n
    real(kind=rp) :: r, idw

    n = this%fields%item_size(1)

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

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

      if (this%one_sided) then

         fum_ib = 0.0_rp
         fvm_ib = 0.0_rp
         fwm_ib = 0.0_rp

         do i = 1, n
            tmp%x(i,1,1,1) = u%x(i,1,1,1) * this%pmsk%x(i,1,1,1)
         end do
         call global_interp%evaluate(fu_ib, tmp%x, .false.)

         do i = 1, n
            tmp%x(i,1,1,1) = v%x(i,1,1,1) * this%pmsk%x(i,1,1,1)
         end do
         call global_interp%evaluate(fv_ib, tmp%x, .false.)

         do i = 1, n
            tmp%x(i,1,1,1) = w%x(i,1,1,1) * this%pmsk%x(i,1,1,1)
         end do
         call global_interp%evaluate(fw_ib, tmp%x, .false.)

         do i = 1, n
            tmp%x(i,1,1,1) = u%x(i,1,1,1) * this%mmsk%x(i,1,1,1)
         end do
         call global_interp%evaluate(fum_ib, tmp%x, .false.)

         do i = 1, n
            tmp%x(i,1,1,1) = v%x(i,1,1,1) * this%mmsk%x(i,1,1,1)
         end do
         call global_interp%evaluate(fvm_ib, tmp%x, .false.)

         do i = 1, n
            tmp%x(i,1,1,1) = w%x(i,1,1,1) * this%mmsk%x(i,1,1,1)
         end do
         call global_interp%evaluate(fwm_ib, tmp%x, .false.)

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
                              if (abs(this%w%x(j,k,l,e)) .gt. 1e-8_rp) then
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
                              if (abs(this%wm%x(j,k,l,e)) .gt. 1e-8_rp) then
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

         call global_interp%evaluate(fu_ib, u%x, .false.)
         call global_interp%evaluate(fv_ib, v%x, .false.)
         call global_interp%evaluate(fw_ib, w%x, .false.)


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

      do i = 1, n
         fu%x(i,1,1,1) = fu%x(i,1,1,1) * this%coef%mult(i,1,1,1)
         fv%x(i,1,1,1) = fv%x(i,1,1,1) * this%coef%mult(i,1,1,1)
         fw%x(i,1,1,1) = fw%x(i,1,1,1) * this%coef%mult(i,1,1,1)
      end do

      call this%gs%op(fu, GS_OP_ADD)
      call this%gs%op(fv, GS_OP_ADD)
      call this%gs%op(fw, GS_OP_ADD)

    end associate

  end subroutine idw_source_term_compute

  subroutine idw_init_boundary_mesh(this, lag_pts, lag_nrm, json)
    class(idw_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(stack_pt_t), intent(inout) :: lag_pts
    type(stack_pt_t), intent(inout) :: lag_nrm
    type(file_t) :: mesh_file
    type(tri_mesh_t) :: boundary_mesh
    character(len=:), allocatable :: mesh_file_name
    character(len=LOG_SIZE) :: log_buf
    integer :: i, j, el_idx
    type(stack_i4_t) :: overlaps
    type(point_t) :: tri_nrm, tri_cntr

    real(kind=rp), dimension(:), allocatable :: box_min, box_max
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

!          do while (.not. overlaps%is_empty())
!             el_idx = overlaps%pop()
!          end do
       end if
       call overlaps%clear()

    end do



!    call boundary_mesh%free()

  end subroutine idw_init_boundary_mesh

  !> Compute IB weight field
  subroutine idw_compute_weight(w, wm, msk, lag_pts, lag_el, x, y, z, &
       ds, rmax, p, lx, ne)
    type(field_t), intent(inout) :: w, wm, msk
    type(point_t), allocatable, intent(inout) :: lag_pts(:)
    type(stack_i4_t), allocatable, intent(inout) :: lag_el(:)
    integer, intent(in) :: lx, ne
    real(kind=rp), dimension(lx,lx,lx,ne) :: x, y, z, ds
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
