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
!
!> Creates a 1d GLL point map along a specified direction based on the
!! connectivity in the mesh.
module map_1d
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use space, only : space_t
  use dofmap, only : dofmap_t
  use gather_scatter, only : GS_OP_ADD
  use mesh, only : mesh_t
  use device, only : device_memcpy, device_map, device_unmap, &
       device_sync, HOST_TO_DEVICE, DEVICE_TO_HOST
  use comm, only : pe_size, pe_rank, NEKO_COMM, MPI_REAL_PRECISION
  use coefs, only : coef_t
  use field_list, only : field_list_t
  use field, only : field_t
  use device_math, only : device_slab_sum, device_gather_add, device_rzero
  use matrix, only : matrix_t
  use vector, only : vector_ptr_t
  use utils, only : neko_error, neko_warning
  use math, only : glmax, glmin, glimax, cmult, add2s1, col2
  use math, only : slab_sum, gather_add
  use mpi_f08, only : MPI_Allreduce, MPI_SUM, MPI_Barrier, MPI_IN_PLACE
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Map every GLL point in the mesh to a level in one physical direction.
  !! Can be used to average across the two remaining directions.
  !! @details
  !! This type is used to average fields over planes normal to one physical
  !! direction. The constructor first determines which local tensor-product
  !! direction (`r`, `s`, or `t`) is most closely aligned with the requested
  !! physical direction (`x`, `y`, or `z`) for each element. It then assigns a
  !! one-dimensional element level to each element and expands that element
  !! level into GLL-point levels along the aligned local direction.
  !!
  !! The element-level assignment is based on propagation of the global minimum
  !! coordinate through the gather-scatter connectivity. Elements whose minimum
  !! coordinate matches the global minimum are assigned to the first level. On
  !! each following pass, the current lower face of every element is masked to
  !! the global maximum, a gather-scatter add exchanges boundary values, and the
  !! result is converted into an effective boundary minimum. Elements that have
  !! newly received the global minimum are assigned to the next level. This
  !! repeats until the propagated minimum reaches the global maximum level.
  !!
  !! Once element levels are known, every GLL point receives a global 1D level
  !! (`pt_lvl`) ordered consistently with the local element orientation. The
  !! constructor also accumulates `volume_per_gll_lvl`, the quadrature weight
  !! volume associated with each 1D level. Averaging then sums each field into
  !! these levels and divides by the corresponding level volume.
  !!
  !! The propagation algorithm assumes that the connectivity and coordinate
  !! layout allow a unique stack of element levels in the requested direction.
  !! If an element level remains unassigned, the resulting point levels are not
  !! valid for the volume accumulation or averaging steps.
  !! The map is built once at initialisation from the mesh coordinates
  !! and assumes that the mesh does not move.
  !! @remark Could also be rather easily extended to say polar coordinates
  !! as well (I think). Martin Karp
  type, public :: map_1d_t
     !> Local tensor-product direction aligned with the requested physical
     !! direction for each element. Values 1, 2, and 3 correspond to `r`, `s`,
     !! and `t`, respectively.
     integer, allocatable :: dir_el(:)
     !> One-dimensional element level assigned by propagating the global
     !! minimum coordinate through the mesh connectivity.
     integer, allocatable :: el_lvl(:)
     !> One-dimensional GLL level for each local point in each element.
     !! Levels run from 1 to `n_gll_lvls` and are used as row indices in the
     !! averaged output matrix.
     integer, allocatable :: pt_lvl(:, :, :, :)
     !> Number of element levels in the requested physical direction.
     integer :: n_el_lvls
     !> Number of total GLL levels in the requested physical direction.
     !! For a polynomial space with `lx` points per element direction, this is
     !! `n_el_lvls * lx`.
     integer :: n_gll_lvls
     !> Dofmap that owns the physical coordinates used to build the map.
     type(dofmap_t), pointer :: dof => null()
     !> SEM coefficients that provide quadrature weights and connectivity.
     type(coef_t), pointer :: coef => null()
     !> Mesh associated with `dof` and `coef`.
     type(mesh_t), pointer :: msh => null()
     !> Requested physical direction of the 1D mapping.
     !! Values 1, 2, and 3 correspond to `x`, `y`, and `z`, respectively.
     integer :: dir
     !> Coordinate comparison tolerance used when identifying propagated levels.
     real(kind=rp) :: tol = 1e-7
     !> Integrated quadrature weight volume associated with each GLL level.
     !! Used as the denominator when computing plane averages.
     real(kind=rp), allocatable :: volume_per_gll_lvl(:)
     !> Volume-averaged coordinate, in the requested direction, of each GLL
     !! level. Written as the first column of the averaged output.
     real(kind=rp), allocatable :: coord_per_gll_lvl(:)
     !> Number of fields accumulated per sample, see `accumulate`.
     integer :: n_acc = 0
     !> 0-based strides of the first reduction of each element, over the
     !! first homogeneous direction: the two kept indices and the summed one.
     integer, allocatable :: el_sa(:), el_sb(:), el_sh(:)
     !> 0-based strides of the second reduction, over the second homogeneous
     !! direction of the `(lx, lx)` slab sums, keeping the level index.
     integer, allocatable :: el_sa2(:), el_sb2(:), el_sh2(:)
     !> 0-based permutation codes, all zero, and identity tables of the two
     !! reductions.
     integer, allocatable :: el_code0(:), tbl0_lxy(:), tbl0_lx(:)
     !> Rows (levels) of the level sums: `csr_ptr` holds `(n_gll_lvls + 1)`
     !! 0-based offsets into `csr_list`, the 0-based entries of `tmp1`
     !! summed by each level.
     integer, allocatable :: csr_ptr(:), csr_list(:)
     !> Slab sums of the elements after the first, `(lxy * nelv)`, and the
     !! second, `(lx * nelv)`, reduction.
     real(kind=rp), allocatable :: tmp(:), tmp1(:)
     !> Accumulated level sums, `(n_gll_lvls, 0:n_acc)`. Slot 0 holds the
     !! accumulated level volumes.
     real(kind=rp), allocatable :: acc(:,:)
     type(c_ptr) :: el_sa_d = C_NULL_PTR
     type(c_ptr) :: el_sb_d = C_NULL_PTR
     type(c_ptr) :: el_sh_d = C_NULL_PTR
     type(c_ptr) :: el_sa2_d = C_NULL_PTR
     type(c_ptr) :: el_sb2_d = C_NULL_PTR
     type(c_ptr) :: el_sh2_d = C_NULL_PTR
     type(c_ptr) :: el_code0_d = C_NULL_PTR
     type(c_ptr) :: tbl0_lxy_d = C_NULL_PTR
     type(c_ptr) :: tbl0_lx_d = C_NULL_PTR
     type(c_ptr) :: csr_ptr_d = C_NULL_PTR
     type(c_ptr) :: csr_list_d = C_NULL_PTR
     type(c_ptr) :: tmp_d = C_NULL_PTR
     type(c_ptr) :: tmp1_d = C_NULL_PTR
     type(c_ptr) :: acc_d = C_NULL_PTR
   contains
     !> Constructor
     procedure, pass(this) :: init_int => map_1d_init
     procedure, pass(this) :: init_char => map_1d_init_char
     generic :: init => init_int, init_char
     !> Destructor
     procedure, pass(this) :: free => map_1d_free
     !> Average field list along planes
     procedure, pass(this) :: average_planes_fld_lst => &
          map_1d_average_field_list
     procedure, pass(this) :: average_planes_vec_ptr => &
          map_1d_average_vector_ptr
     generic :: average_planes => average_planes_fld_lst, average_planes_vec_ptr
     !> Prepare the accumulation of fields averaged over the planes.
     procedure, pass(this) :: accumulate_init => map_1d_accumulate_init
     !> Add the plane sums of a field, scaled, to the accumulated sums.
     procedure, pass(this) :: accumulate => map_1d_accumulate
     !> Add the plane volumes, scaled, to the accumulated sums.
     procedure, pass(this) :: accumulate_volume => map_1d_accumulate_volume
     !> Reset the accumulated sums to zero.
     procedure, pass(this) :: accumulate_reset => map_1d_accumulate_reset
     !> Output the accumulated averages as a matrix.
     procedure, pass(this) :: accumulated_average => &
          map_1d_accumulated_average
     !> Free the accumulation data.
     procedure, pass(this) :: accumulate_free => map_1d_accumulate_free
  end type map_1d_t


contains

  subroutine map_1d_init(this, coef, dir, tol)
    class(map_1d_t) :: this
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    integer :: nelv, lx, n, i, e, lvl, ierr
    real(kind=rp), contiguous, pointer :: line(:, :, :, :)
    real(kind=rp), allocatable :: min_vals(:, :, :, :)
    real(kind=rp), allocatable :: min_temp(:, :, :, :)
    type(c_ptr) :: min_vals_d = c_null_ptr
    real(kind=rp) :: el_dim(3, 3), glb_min, glb_max, el_min, atol

    call this%free()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       if (pe_rank .eq. 0) then
          call neko_warning('map_1d does not copy indices to device, ' // &
               ' but ok if used on cpu and for io')
       end if
    end if

    this%dir = dir
    this%tol = tol
    this%dof => coef%dof
    this%coef => coef
    this%msh => coef%msh
    nelv = this%msh%nelv
    lx = this%dof%Xh%lx
    n = this%dof%size()

    if (dir .eq. 1) then
       line => this%dof%x%x
    else if (dir .eq. 2) then
       line => this%dof%y%x
    else if (dir .eq. 3) then
       line => this%dof%z%x
    else
       call neko_error('Invalid dir for geopmetric comm')
    end if
    allocate(this%dir_el(nelv))
    allocate(this%el_lvl(nelv))
    allocate(this%pt_lvl(lx, lx, lx, nelv))
    allocate(min_vals(lx, lx, lx, nelv))
    allocate(min_temp(lx, lx, lx, nelv))
    call MPI_BARRIER(NEKO_COMM)
    if (NEKO_BCKND_DEVICE .eq. 1) then

       call device_map(min_vals, min_vals_d, n)
    end if
    call MPI_BARRIER(NEKO_COMM)

    do i = 1, nelv
       ! Store which direction r,s,t corresponds to specified direction, x,y,z
       ! we assume elements are stacked on each other...
       ! Check which one of the normalized vectors are closest to dir
       ! If we want to incorporate other directions, we should look here

       ! This follows the point ordering defined in hex.f90
       el_dim(1, :) = abs(this%msh%elements(i)%e%pts(1)%p%x - &
            this%msh%elements(i)%e%pts(2)%p%x)
       el_dim(1, :) = el_dim(1, :)/norm2(el_dim(1, :))
       el_dim(2, :) = abs(this%msh%elements(i)%e%pts(1)%p%x - &
            this%msh%elements(i)%e%pts(3)%p%x)
       el_dim(2, :) = el_dim(2, :)/norm2(el_dim(2, :))
       el_dim(3, :) = abs(this%msh%elements(i)%e%pts(1)%p%x - &
            this%msh%elements(i)%e%pts(5)%p%x)
       el_dim(3, :) = el_dim(3, :)/norm2(el_dim(3, :))
       ! Checks which directions in rst the xyz corresponds to
       ! 1 corresponds to r, 2 to s, 3 to t and are stored in dir_el
       this%dir_el(i) = maxloc(el_dim(:, this%dir), dim = 1)
    end do

    glb_min = glmin(line, n)
    glb_max = glmax(line, n)
    ! Tolerance for a coordinate to count as the minimum, relative to the
    ! extent of the domain rather than to the minimum itself, which may be 0.
    atol = this%tol * (glb_max - glb_min)

    i = 1
    this%el_lvl = -1
    ! Check what the minimum value in each element and put in min_vals
    do e = 1, nelv
       el_min = minval(line(:, :, :, e))
       min_vals(:, :, :, e) = el_min
       ! Check if this element is on the bottom,
       ! in this case assign el_lvl = i = 1
       if (abs(el_min - glb_min) .le. atol) then
          if (this%el_lvl(e) .eq. -1) this%el_lvl(e) = i
       end if
    end do

    ! While loop where at each iteration the global maximum value
    ! propagates down one level.
    ! When the minimum value has propagated to the highest level this stops.
    ! Only works when the bottom plate of the domain is flat.
    do while (abs(glmax(min_vals, n) - glb_min) .gt. atol)

       ! The propagation passes the minimum between layers through the
       ! face-interior nodes, so it never finishes without them (polynomial
       ! order 1) or when the mesh is not stacked in the requested direction.
       ! There cannot be more levels than elements.
       if (i .gt. this%msh%glb_nelv) then
          call neko_error('map_1d: the element levels could not be ' // &
               'determined, the mesh must be stacked in the requested ' // &
               'direction and the polynomial order at least 2')
       end if

       ! This is the assigned level
       i = i + 1

       do e = 1, nelv
          !Sets the value at the bottom of each element to glb_max
          if (this%dir_el(e) .eq. 1) then
             if (line(1, 1, 1, e) .gt. line(lx, 1, 1, e)) then
                min_vals(lx, :, :, e) = glb_max
             else
                min_vals(1, :, :, e) = glb_max
             end if
          end if
          if (this%dir_el(e) .eq. 2) then
             if (line(1, 1, 1, e) .gt. line(1, lx, 1, e)) then
                min_vals(:, lx, :, e) = glb_max
             else
                min_vals(:, 1, :, e) = glb_max
             end if
          end if
          if (this%dir_el(e) .eq. 3) then
             if (line(1, 1, 1, e) .gt. line(1, 1, lx, e)) then
                min_vals(:, :, lx, e) = glb_max
             else
                min_vals(:, :, 1, e) = glb_max
             end if
          end if
       end do

       !Make sketchy min as GS_OP_MIN is not supported with device mpi
       min_temp = min_vals

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(min_vals, min_vals_d, n, HOST_TO_DEVICE, &
               sync = .false.)
       end if

       !Propagates the minimum value along the element boundary.
       call coef%gs_h%op(min_vals, n, GS_OP_ADD)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(min_vals, min_vals_d, n, DEVICE_TO_HOST, &
               sync = .true.)
       end if

       !Obtain average along boundary
       call col2(min_vals, coef%mult, n)
       call cmult(min_temp, -1.0_rp, n)
       call add2s1(min_vals, min_temp, 2.0_rp, n)

       !Checks the new minimum value on each element
       !Assign this value to all points in this element in min_val
       !If the element has not already been assigned a level,
       !and it has obtained the minval, set el_lvl = i
       do e = 1, nelv
          el_min = minval(min_vals(:, :, :, e))
          min_vals(:, :, :, e) = el_min
          if (abs(el_min - glb_min) .le. atol) then
             if (this%el_lvl(e) .eq. -1) this%el_lvl(e) = i
          end if
       end do
    end do
    this%n_el_lvls = glimax(this%el_lvl, nelv)
    this%n_gll_lvls = this%n_el_lvls*lx

    ! Every element must have received a level, otherwise the point levels
    ! would index outside the level arrays.
    if (glimax(abs(min(this%el_lvl, 0)), nelv) .gt. 0) then
       call neko_error('map_1d: an element was not assigned a level, the ' // &
            'mesh must be stacked in the requested direction')
    end if

    !Numbers the points in each element based on the element level
    !and its orientation
    do e = 1, nelv
       do i = 1, lx
          lvl = lx * (this%el_lvl(e) - 1) + i
          if (this%dir_el(e) .eq. 1) then
             if (line(1, 1, 1, e) .gt. line(lx, 1, 1, e)) then
                this%pt_lvl(lx-i+1, :, :, e) = lvl
             else
                this%pt_lvl(i, :, :, e) = lvl
             end if
          end if
          if (this%dir_el(e) .eq. 2) then
             if (line(1, 1, 1, e) .gt. line(1, lx, 1, e)) then
                this%pt_lvl(:, lx-i+1, :, e) = lvl
             else
                this%pt_lvl(:, i, :, e) = lvl
             end if
          end if
          if (this%dir_el(e) .eq. 3) then
             if (line(1, 1, 1, e) .gt. line(1, 1, lx, e)) then
                this%pt_lvl(:, :, lx-i+1, e) = lvl
             else
                this%pt_lvl(:, :, i, e) = lvl
             end if
          end if
       end do
    end do

    if (allocated(min_vals)) then
       if (c_associated(min_vals_d)) call device_unmap(min_vals, min_vals_d)
       deallocate(min_vals)
    end if

    if (allocated(min_temp)) deallocate(min_temp)
    allocate(this%volume_per_gll_lvl(this%n_gll_lvls))

    this%volume_per_gll_lvl = 0.0_rp

    do i = 1, n
       this%volume_per_gll_lvl(this%pt_lvl(i, 1, 1, 1)) = &
            this%volume_per_gll_lvl(this%pt_lvl(i, 1, 1, 1)) + &
            coef%B(i, 1, 1, 1)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, this%volume_per_gll_lvl, this%n_gll_lvls, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

    allocate(this%coord_per_gll_lvl(this%n_gll_lvls))

    this%coord_per_gll_lvl = 0.0_rp

    do i = 1, n
       this%coord_per_gll_lvl(this%pt_lvl(i, 1, 1, 1)) = &
            this%coord_per_gll_lvl(this%pt_lvl(i, 1, 1, 1)) + &
            line(i, 1, 1, 1) * coef%B(i, 1, 1, 1)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, this%coord_per_gll_lvl, this%n_gll_lvls, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

    this%coord_per_gll_lvl = this%coord_per_gll_lvl / this%volume_per_gll_lvl

  end subroutine map_1d_init

  subroutine map_1d_init_char(this, coef, dir, tol)
    class(map_1d_t) :: this
    type(coef_t), intent(inout), target :: coef
    character(len=*), intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    integer :: idir

    if (trim(dir) .eq. 'yz' .or. trim(dir) .eq. 'zy') then
       idir = 1
    else if (trim(dir) .eq. 'xz' .or. trim(dir) .eq. 'zx') then
       idir = 2
    else if (trim(dir) .eq. 'xy' .or. trim(dir) .eq. 'yx') then
       idir = 3
    else
       call neko_error('homogenous direction not supported')
    end if

    call this%init(coef, idir, tol)

  end subroutine map_1d_init_char

  !> Destructor.
  subroutine map_1d_free(this)
    class(map_1d_t) :: this

    call this%accumulate_free()
    if (allocated(this%dir_el)) deallocate(this%dir_el)
    if (allocated(this%el_lvl)) deallocate(this%el_lvl)
    if (allocated(this%pt_lvl)) deallocate(this%pt_lvl)
    if (associated(this%dof)) nullify(this%dof)
    if (associated(this%msh)) nullify(this%msh)
    if (associated(this%coef)) nullify(this%coef)
    if (allocated(this%volume_per_gll_lvl)) deallocate(this%volume_per_gll_lvl)
    if (allocated(this%coord_per_gll_lvl)) deallocate(this%coord_per_gll_lvl)
    this%dir = 0
    this%n_el_lvls = 0
    this%n_gll_lvls = 0

  end subroutine map_1d_free


  !> Computes average if field list in two directions and outputs matrix
  !! with averaged values
  !! avg_planes contains coordinates in first row, avg. of fields in the rest
  !! @param avg_planes output averages
  !! @param field_list list of fields to be averaged
  subroutine map_1d_average_field_list(this, avg_planes, field_list)
    class(map_1d_t), intent(in) :: this
    type(field_list_t), intent(in) :: field_list
    type(matrix_t), intent(inout) :: avg_planes
    integer :: n, j, i

    call avg_planes%free()
    call avg_planes%init(this%n_gll_lvls, field_list%size() + 1)
    avg_planes = 0.0_rp

    n = this%dof%size()
    do j = 2, field_list%size() + 1
       do i = 1, n
          avg_planes%x(this%pt_lvl(i,1,1,1), j) = &
               avg_planes%x(this%pt_lvl(i,1,1,1), j) + &
               field_list%items(j-1)%ptr%x(i,1,1,1) * this%coef%B(i,1,1,1)
       end do
    end do

    call map_1d_finalize_average(this, avg_planes, field_list%size())

  end subroutine map_1d_average_field_list

  !> Computes average of vector_pt in two directions and outputs matrix
  !! with averaged values
  !! avg_planes contains coordinates in first row, avg. of fields in the rest
  !! @param avg_planes output averages
  !! @param vector_pts to vectors to be averaged
  subroutine map_1d_average_vector_ptr(this, avg_planes, vector_ptr)
    class(map_1d_t), intent(in) :: this
    !Observe is an array...
    type(vector_ptr_t), intent(in) :: vector_ptr(:)
    type(matrix_t), intent(inout) :: avg_planes
    integer :: n, j, i

    call avg_planes%free()
    call avg_planes%init(this%n_gll_lvls, size(vector_ptr) + 1)
    avg_planes = 0.0_rp

    n = this%dof%size()
    do j = 2, size(vector_ptr) + 1
       do i = 1, n
          avg_planes%x(this%pt_lvl(i,1,1,1), j) = &
               avg_planes%x(this%pt_lvl(i,1,1,1), j) + &
               vector_ptr(j-1)%ptr%x(i) * this%coef%B(i,1,1,1)
       end do
    end do

    call map_1d_finalize_average(this, avg_planes, size(vector_ptr))

  end subroutine map_1d_average_vector_ptr

  !> Reduces the volume-weighted level sums of `n_fields` fields, stored in
  !! columns 2 to `n_fields + 1` of `avg_planes`, over all ranks, divides
  !! them by the level volumes and stores the level coordinates in the first
  !! column.
  !! @param avg_planes Level sums on input, averages on output.
  !! @param n_fields Number of fields.
  subroutine map_1d_finalize_average(this, avg_planes, n_fields)
    class(map_1d_t), intent(in) :: this
    type(matrix_t), intent(inout) :: avg_planes
    integer, intent(in) :: n_fields
    real(kind=rp), allocatable :: sums(:,:)
    integer :: ierr, j

    allocate(sums(this%n_gll_lvls, n_fields))
    sums = avg_planes%x(:, 2:n_fields + 1)
    if (pe_size .gt. 1 .and. n_fields .gt. 0) then
       call MPI_Allreduce(MPI_IN_PLACE, sums, n_fields * this%n_gll_lvls, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if

    avg_planes%x(:, 1) = this%coord_per_gll_lvl
    do j = 1, n_fields
       avg_planes%x(:, j + 1) = sums(:, j) / this%volume_per_gll_lvl
    end do
    deallocate(sums)

  end subroutine map_1d_finalize_average

  !> Prepares the accumulation of `n_fields` fields averaged over the
  !! planes, see `accumulate`.
  !! @param n_fields Number of fields accumulated per sample.
  subroutine map_1d_accumulate_init(this, n_fields)
    class(map_1d_t), intent(inout) :: this
    integer, intent(in) :: n_fields
    integer, allocatable :: cnt(:)
    integer :: str(3), e, h, d, dl, d1, d2, lvl, lx, lxy, lxyz, nelv, p, m

    call this%accumulate_free()
    this%n_acc = n_fields
    lx = this%dof%Xh%lx
    lxy = this%dof%Xh%lxy
    lxyz = this%dof%Xh%lxyz
    nelv = this%msh%nelv
    str = [1, lx, lxy]

    ! The first reduction sums over the lower of the two homogeneous local
    ! directions and keeps the other two in increasing order, the second
    ! sums over the remaining homogeneous direction and keeps the level.
    allocate(this%el_sa(nelv), this%el_sb(nelv), this%el_sh(nelv))
    allocate(this%el_sa2(nelv), this%el_sb2(nelv), this%el_sh2(nelv))
    allocate(this%el_code0(nelv))
    do e = 1, nelv
       dl = this%dir_el(e)
       d1 = 0
       d2 = 0
       do d = 1, 3
          if (d .ne. dl) then
             if (d1 .eq. 0) then
                d1 = d
             else
                d2 = d
             end if
          end if
       end do
       this%el_sh(e) = str(d1)
       if (dl .lt. d2) then
          this%el_sa(e) = str(dl)
          this%el_sb(e) = str(d2)
          this%el_sa2(e) = 1
          this%el_sh2(e) = lx
       else
          this%el_sa(e) = str(d2)
          this%el_sb(e) = str(dl)
          this%el_sa2(e) = lx
          this%el_sh2(e) = 1
       end if
       this%el_sb2(e) = 0
       this%el_code0(e) = 0
    end do
    allocate(this%tbl0_lxy(lxy), this%tbl0_lx(lx))
    do m = 1, lxy
       this%tbl0_lxy(m) = m - 1
    end do
    do m = 1, lx
       this%tbl0_lx(m) = m - 1
    end do

    ! Row lvl of the level sums gathers the level sums of the elements at
    ! that level, in element order.
    allocate(this%csr_ptr(this%n_gll_lvls + 1))
    allocate(this%csr_list(max(lx * nelv, 1)))
    allocate(cnt(this%n_gll_lvls))
    cnt = 0
    do e = 1, nelv
       do h = 0, lx - 1
          p = (e - 1) * lxyz + h * str(this%dir_el(e)) + 1
          lvl = this%pt_lvl(p, 1, 1, 1)
          cnt(lvl) = cnt(lvl) + 1
       end do
    end do
    this%csr_ptr(1) = 0
    do lvl = 1, this%n_gll_lvls
       this%csr_ptr(lvl + 1) = this%csr_ptr(lvl) + cnt(lvl)
    end do
    cnt = 0
    do e = 1, nelv
       do h = 0, lx - 1
          p = (e - 1) * lxyz + h * str(this%dir_el(e)) + 1
          lvl = this%pt_lvl(p, 1, 1, 1)
          this%csr_list(this%csr_ptr(lvl) + cnt(lvl) + 1) = h + lx * (e - 1)
          cnt(lvl) = cnt(lvl) + 1
       end do
    end do
    deallocate(cnt)

    allocate(this%tmp(max(lxy * nelv, 1)), this%tmp1(max(lx * nelv, 1)))
    allocate(this%acc(this%n_gll_lvls, 0:n_fields))

    ! The per-element tables are only mapped on ranks with elements; the
    ! other arrays always, since the level sums are gathered on every rank.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       if (nelv .gt. 0) then
          call device_map(this%el_sa, this%el_sa_d, nelv)
          call device_map(this%el_sb, this%el_sb_d, nelv)
          call device_map(this%el_sh, this%el_sh_d, nelv)
          call device_map(this%el_sa2, this%el_sa2_d, nelv)
          call device_map(this%el_sb2, this%el_sb2_d, nelv)
          call device_map(this%el_sh2, this%el_sh2_d, nelv)
          call device_map(this%el_code0, this%el_code0_d, nelv)
          call device_memcpy(this%el_sa, this%el_sa_d, nelv, HOST_TO_DEVICE, &
               sync = .false.)
          call device_memcpy(this%el_sb, this%el_sb_d, nelv, HOST_TO_DEVICE, &
               sync = .false.)
          call device_memcpy(this%el_sh, this%el_sh_d, nelv, HOST_TO_DEVICE, &
               sync = .false.)
          call device_memcpy(this%el_sa2, this%el_sa2_d, nelv, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%el_sb2, this%el_sb2_d, nelv, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%el_sh2, this%el_sh2_d, nelv, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%el_code0, this%el_code0_d, nelv, &
               HOST_TO_DEVICE, sync = .false.)
       end if
       call device_map(this%tbl0_lxy, this%tbl0_lxy_d, lxy)
       call device_map(this%tbl0_lx, this%tbl0_lx_d, lx)
       call device_map(this%csr_ptr, this%csr_ptr_d, size(this%csr_ptr))
       call device_map(this%csr_list, this%csr_list_d, size(this%csr_list))
       call device_map(this%tmp, this%tmp_d, size(this%tmp))
       call device_map(this%tmp1, this%tmp1_d, size(this%tmp1))
       call device_map(this%acc, this%acc_d, size(this%acc))
       call device_memcpy(this%tbl0_lxy, this%tbl0_lxy_d, lxy, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%tbl0_lx, this%tbl0_lx_d, lx, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%csr_ptr, this%csr_ptr_d, size(this%csr_ptr), &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%csr_list, this%csr_list_d, &
            size(this%csr_list), HOST_TO_DEVICE, sync = .false.)
       ! Zero the work and accumulator on the device first: on unified
       ! memory the device then faults the pages (device first touch), and
       ! the host must not write to them while device work is in flight.
       call device_rzero(this%tmp_d, size(this%tmp))
       call device_rzero(this%tmp1_d, size(this%tmp1))
       call device_rzero(this%acc_d, size(this%acc))
       call device_sync()
    end if
    this%tmp = 0.0_rp
    this%tmp1 = 0.0_rp
    this%acc = 0.0_rp

  end subroutine map_1d_accumulate_init

  !> Frees the accumulation data.
  subroutine map_1d_accumulate_free(this)
    class(map_1d_t), intent(inout) :: this

    if (c_associated(this%el_sa_d)) call device_unmap(this%el_sa, this%el_sa_d)
    if (c_associated(this%el_sb_d)) call device_unmap(this%el_sb, this%el_sb_d)
    if (c_associated(this%el_sh_d)) call device_unmap(this%el_sh, this%el_sh_d)
    if (c_associated(this%el_sa2_d)) then
       call device_unmap(this%el_sa2, this%el_sa2_d)
    end if
    if (c_associated(this%el_sb2_d)) then
       call device_unmap(this%el_sb2, this%el_sb2_d)
    end if
    if (c_associated(this%el_sh2_d)) then
       call device_unmap(this%el_sh2, this%el_sh2_d)
    end if
    if (c_associated(this%el_code0_d)) then
       call device_unmap(this%el_code0, this%el_code0_d)
    end if
    if (c_associated(this%tbl0_lxy_d)) then
       call device_unmap(this%tbl0_lxy, this%tbl0_lxy_d)
    end if
    if (c_associated(this%tbl0_lx_d)) then
       call device_unmap(this%tbl0_lx, this%tbl0_lx_d)
    end if
    if (c_associated(this%csr_ptr_d)) then
       call device_unmap(this%csr_ptr, this%csr_ptr_d)
    end if
    if (c_associated(this%csr_list_d)) then
       call device_unmap(this%csr_list, this%csr_list_d)
    end if
    if (c_associated(this%tmp_d)) call device_unmap(this%tmp, this%tmp_d)
    if (c_associated(this%tmp1_d)) call device_unmap(this%tmp1, this%tmp1_d)
    if (c_associated(this%acc_d)) call device_unmap(this%acc, this%acc_d)
    this%el_sa_d = C_NULL_PTR
    this%el_sb_d = C_NULL_PTR
    this%el_sh_d = C_NULL_PTR
    this%el_sa2_d = C_NULL_PTR
    this%el_sb2_d = C_NULL_PTR
    this%el_sh2_d = C_NULL_PTR
    this%el_code0_d = C_NULL_PTR
    this%tbl0_lxy_d = C_NULL_PTR
    this%tbl0_lx_d = C_NULL_PTR
    this%csr_ptr_d = C_NULL_PTR
    this%csr_list_d = C_NULL_PTR
    this%tmp_d = C_NULL_PTR
    this%tmp1_d = C_NULL_PTR
    this%acc_d = C_NULL_PTR

    if (allocated(this%el_sa)) deallocate(this%el_sa)
    if (allocated(this%el_sb)) deallocate(this%el_sb)
    if (allocated(this%el_sh)) deallocate(this%el_sh)
    if (allocated(this%el_sa2)) deallocate(this%el_sa2)
    if (allocated(this%el_sb2)) deallocate(this%el_sb2)
    if (allocated(this%el_sh2)) deallocate(this%el_sh2)
    if (allocated(this%el_code0)) deallocate(this%el_code0)
    if (allocated(this%tbl0_lxy)) deallocate(this%tbl0_lxy)
    if (allocated(this%tbl0_lx)) deallocate(this%tbl0_lx)
    if (allocated(this%csr_ptr)) deallocate(this%csr_ptr)
    if (allocated(this%csr_list)) deallocate(this%csr_list)
    if (allocated(this%tmp)) deallocate(this%tmp)
    if (allocated(this%tmp1)) deallocate(this%tmp1)
    if (allocated(this%acc)) deallocate(this%acc)
    this%n_acc = 0

  end subroutine map_1d_accumulate_free

  !> Adds the plane sums of `k * f * B` to slot `slot` of the accumulated
  !! level sums, where `B` is the mass matrix. The field is expected to be
  !! up to date on the device, or on the host for the CPU backend.
  !! @param f Field to accumulate.
  !! @param slot Slot of the field, 1 to `n_acc`.
  !! @param k Scaling of the sample, typically the time since the last one.
  subroutine map_1d_accumulate(this, f, slot, k)
    class(map_1d_t), intent(inout) :: this
    type(field_t), intent(in) :: f
    integer, intent(in) :: slot
    real(kind=rp), intent(in) :: k

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call map_1d_accumulate_device(this, f%x_d, slot, k)
    else
       call map_1d_accumulate_host(this, slot, k, f%x)
    end if

  end subroutine map_1d_accumulate

  !> Adds the plane volumes, `k * B` summed over the planes, to slot 0 of
  !! the accumulated level sums.
  !! @param k Scaling of the sample, typically the time since the last one.
  subroutine map_1d_accumulate_volume(this, k)
    class(map_1d_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call map_1d_accumulate_device(this, C_NULL_PTR, 0, k)
    else
       call map_1d_accumulate_host(this, 0, k)
    end if

  end subroutine map_1d_accumulate_volume

  !> Device implementation of the accumulation, a null `f_d` counts as one.
  !! @param f_d Device pointer of the field, or null.
  !! @param slot Slot of the field, 0 to `n_acc`.
  !! @param k Scaling of the sample.
  subroutine map_1d_accumulate_device(this, f_d, slot, k)
    class(map_1d_t), intent(inout) :: this
    type(c_ptr), intent(in) :: f_d
    integer, intent(in) :: slot
    real(kind=rp), intent(in) :: k
    integer :: lx, lxy, nelv

    lx = this%dof%Xh%lx
    lxy = this%dof%Xh%lxy
    nelv = this%msh%nelv

    call device_slab_sum(this%tmp_d, f_d, this%coef%B_d, this%el_sa_d, &
         this%el_sb_d, this%el_sh_d, this%el_code0_d, this%tbl0_lxy_d, lxy, &
         lx, nelv, this%dof%Xh%lxyz, lxy)
    call device_slab_sum(this%tmp1_d, this%tmp_d, C_NULL_PTR, this%el_sa2_d, &
         this%el_sb2_d, this%el_sh2_d, this%el_code0_d, this%tbl0_lx_d, lx, &
         lx, nelv, lxy, lx)
    call device_gather_add(this%acc_d, slot * this%n_gll_lvls, this%tmp1_d, &
         this%csr_ptr_d, this%csr_list_d, this%n_gll_lvls, k)

  end subroutine map_1d_accumulate_device

  !> Host implementation of the accumulation, a missing `f` counts as one.
  !! @param slot Slot of the field, 0 to `n_acc`.
  !! @param k Scaling of the sample.
  !! @param f Field to accumulate, optional.
  subroutine map_1d_accumulate_host(this, slot, k, f)
    class(map_1d_t), intent(inout) :: this
    integer, intent(in) :: slot
    real(kind=rp), intent(in) :: k
    real(kind=rp), intent(in), optional :: f(*)
    integer :: lx, lxy, nelv

    lx = this%dof%Xh%lx
    lxy = this%dof%Xh%lxy
    nelv = this%msh%nelv

    call slab_sum(this%tmp, this%el_sa, this%el_sb, this%el_sh, &
         this%el_code0, this%tbl0_lxy, lxy, lx, nelv, this%dof%Xh%lxyz, lxy, &
         f, this%coef%B)
    call slab_sum(this%tmp1, this%el_sa2, this%el_sb2, this%el_sh2, &
         this%el_code0, this%tbl0_lx, lx, lx, nelv, lxy, lx, f = this%tmp)
    call gather_add(this%acc, slot * this%n_gll_lvls, this%tmp1, &
         this%csr_ptr, this%csr_list, this%n_gll_lvls, k)

  end subroutine map_1d_accumulate_host

  !> Resets the accumulated level sums to zero.
  subroutine map_1d_accumulate_reset(this)
    class(map_1d_t), intent(inout) :: this

    ! Device first, see accumulate_init.
    if (c_associated(this%acc_d)) then
       call device_rzero(this%acc_d, size(this%acc))
       call device_sync()
    end if
    if (allocated(this%acc)) this%acc = 0.0_rp

  end subroutine map_1d_accumulate_reset

  !> Sums the accumulated level sums over all ranks, divides them by the
  !! accumulated level volumes and outputs the averages as a matrix with
  !! the level coordinates in the first column and the slots in the
  !! following ones.
  !! @param avg_planes Output averages.
  subroutine map_1d_accumulated_average(this, avg_planes)
    class(map_1d_t), intent(inout) :: this
    type(matrix_t), intent(inout) :: avg_planes
    real(kind=rp), allocatable :: sums(:,:)
    integer :: j, nf, ierr

    nf = this%n_acc
    if (c_associated(this%acc_d)) then
       call device_memcpy(this%acc, this%acc_d, size(this%acc), &
            DEVICE_TO_HOST, sync = .true.)
    end if

    allocate(sums(this%n_gll_lvls, 0:nf))
    sums = this%acc
    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, sums, this%n_gll_lvls * (nf + 1), &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
    ! A level without volume has received no sample and averages to zero.
    call avg_planes%free()
    call avg_planes%init(this%n_gll_lvls, nf + 1)
    avg_planes%x(:, 1) = this%coord_per_gll_lvl
    do j = 1, nf
       where (sums(:, 0) .gt. 0.0_rp)
          avg_planes%x(:, j + 1) = sums(:, j) / sums(:, 0)
       elsewhere
          avg_planes%x(:, j + 1) = 0.0_rp
       end where
    end do
    deallocate(sums)

  end subroutine map_1d_accumulated_average

end module map_1d
