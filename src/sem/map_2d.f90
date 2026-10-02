! Copyright (c) 2024-2026, The Neko Authors
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
!> Maps a 3D dofmap to a 2D spectral element grid.
!! @details
!! The 2D grid consists of the elements at the first level of a `map_1d_t`
!! in the requested direction, i.e. the elements touching the lower boundary
!! of the domain. The constructor assigns every element of the 3D mesh to the
!! column of the 2D element whose in-plane position it shares, together with
!! the permutation of its in-plane nodes that maps them onto the nodes of that
!! 2D element. The assignment is done by comparing the in-plane geometry of
!! the elements, so the mesh is assumed to consist of straight-sided elements
!! stacked in the homogeneous direction.
!!
!! With the column map in place, averaging a field in the homogeneous
!! direction is a single pass over the field: every point is weighted with the
!! mass matrix and added to its column, after which the column sums are sent
!! to the rank owning the 2D element and divided by the column volume. No
!! gather-scatter operations or device transfers are involved; the fields
!! are expected to be up to date on the host.
module map_2d
  use num_types, only : rp
  use dofmap, only : dofmap_t
  use map_1d, only : map_1d_t
  use mesh, only : mesh_t
  use field, only : field_t
  use field_list, only : field_list_t
  use coefs, only : coef_t
  use vector, only : vector_ptr_t
  use utils, only : neko_error
  use math, only : sort, glmax, glmin
  use comm, only : NEKO_COMM, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_Allgather, MPI_Allgatherv, &
       MPI_Alltoall, MPI_Alltoallv, MPI_Exscan, MPI_INTEGER, MPI_SUM
  use fld_file_data, only : fld_file_data_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, device_unmap, device_memcpy, &
       device_sync, HOST_TO_DEVICE, DEVICE_TO_HOST
  use device_math, only : device_slab_sum, device_gather_add, device_rzero
  use, intrinsic :: iso_c_binding, only : c_ptr, c_associated, C_NULL_PTR
  implicit none
  private

  !> Number of in-plane node permutations (the symmetries of a square).
  integer, parameter :: N_PERM = 8
  !> Number of reals describing the in-plane geometry of an element.
  integer, parameter :: N_GEO = 10

  !> A field viewed as a contiguous 1D array.
  type :: field_ptr_1d_t
     real(kind=rp), pointer, contiguous :: x(:) => null()
  end type field_ptr_1d_t

  type, public :: map_2d_t
     integer :: nelv_2d = 0 !< Number of elements in 2D mesh on this rank
     integer :: glb_nelv_2d = 0 !< global number of elements in 2d
     integer :: offset_el_2d = 0 !< element offset for this rank
     integer :: lxy = 0 !< number of gll points per 2D element
     integer :: n_2d = 0 !< total number of gll points (nelv_2d*lxy)
     integer, allocatable :: idx_2d(:) !< Mapping of GLL point from 3D to 2D
     integer, allocatable :: el_idx_2d(:) !< Mapping of element in 3D to 2D
     type(map_1d_t) :: map_1d !< 1D map in normal direction to 2D plane
     type(mesh_t), pointer :: msh => null() !< 3D mesh
     type(dofmap_t), pointer :: dof => null() !< 3D dofmap
     type(coef_t), pointer :: coef => null() !< 3D SEM coefs
     integer :: dir = 0 !< direction normal to 2D plane
     !> Number of distinct columns among the elements on this rank.
     integer :: n_cols = 0
     !> Column of each element. Columns are numbered in the order of the
     !! global index of their 2D element, which also orders them by owner.
     integer, allocatable :: el_col(:)
     !> Permutation code of the in-plane nodes of each element.
     integer, allocatable :: el_perm(:)
     !> In-plane node permutation tables, `(lxy, N_PERM)`.
     integer, allocatable :: perm_tbl(:,:)
     !> Number of local columns owned by each rank.
     integer, allocatable :: send_cols(:)
     !> Offset of the first local column owned by each rank.
     integer, allocatable :: send_displ(:)
     !> Number of column contributions received from each rank.
     integer, allocatable :: recv_cols(:)
     !> Offset of the first column contribution received from each rank.
     integer, allocatable :: recv_displ(:)
     !> Total number of column contributions received by this rank.
     integer :: n_recv = 0
     !> Local 2D element of each received column contribution.
     integer, allocatable :: recv_el(:)
     !> Number of fields accumulated per sample, see `accumulate`.
     integer :: n_acc = 0
     !> 0-based strides of the in-plane indices and of the homogeneous
     !! index of each element, see `element_strides`.
     integer, allocatable :: el_sa(:), el_sb(:), el_sh(:)
     !> 0-based permutation code of each element.
     integer, allocatable :: el_code0(:)
     !> 0-based in-plane node permutation tables, `(lxy * N_PERM)`.
     integer, allocatable :: tbl0(:)
     !> Rows of the column sums: `csr_ptr` holds `(lxy * n_cols + 1)` 0-based
     !! offsets into `csr_list`, the 0-based entries of `tmp` summed by each
     !! row.
     integer, allocatable :: csr_ptr(:), csr_list(:)
     !> Slab sums of the elements, `(lxy * nelv)`.
     real(kind=rp), allocatable :: tmp(:)
     !> Accumulated column sums, `(lxy, n_cols, 0:n_acc)`. Slot 0 holds the
     !! accumulated column volumes.
     real(kind=rp), allocatable :: acc(:,:,:)
     type(c_ptr) :: el_sa_d = C_NULL_PTR
     type(c_ptr) :: el_sb_d = C_NULL_PTR
     type(c_ptr) :: el_sh_d = C_NULL_PTR
     type(c_ptr) :: el_code0_d = C_NULL_PTR
     type(c_ptr) :: tbl0_d = C_NULL_PTR
     type(c_ptr) :: csr_ptr_d = C_NULL_PTR
     type(c_ptr) :: csr_list_d = C_NULL_PTR
     type(c_ptr) :: tmp_d = C_NULL_PTR
     type(c_ptr) :: acc_d = C_NULL_PTR
   contains
     procedure, pass(this) :: init_int => map_2d_init
     procedure, pass(this) :: init_char => map_2d_init_char
     procedure, pass(this) :: free => map_2d_free
     generic :: init => init_int, init_char
     procedure, pass(this) :: average_file => map_2d_average
     procedure, pass(this) :: average_list => map_2d_average_field_list
     generic :: average => average_list, average_file
     !> Prepare the accumulation of fields averaged in the homogeneous
     !! direction.
     procedure, pass(this) :: accumulate_init => map_2d_accumulate_init
     !> Add the column sums of a field, scaled, to the accumulated sums.
     procedure, pass(this) :: accumulate => map_2d_accumulate
     !> Add the column volumes, scaled, to the accumulated sums.
     procedure, pass(this) :: accumulate_volume => map_2d_accumulate_volume
     !> Reset the accumulated sums to zero.
     procedure, pass(this) :: accumulate_reset => map_2d_accumulate_reset
     !> Output the accumulated averages as a 2D field.
     procedure, pass(this) :: accumulated_average => &
          map_2d_accumulated_average
     !> Free the accumulation data.
     procedure, pass(this) :: accumulate_free => map_2d_accumulate_free
  end type map_2d_t

contains

  !> Constructor.
  !! @param coef SEM coefficients of the 3D mesh.
  !! @param dir Direction normal to the 2D plane, 1, 2 or 3 for x, y or z.
  !! @param tol Tolerance, relative to the in-plane extent of the domain,
  !! for comparing coordinates.
  subroutine map_2d_init(this, coef, dir, tol)
    class(map_2d_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    real(kind=rp), contiguous, pointer :: x_ptr(:), y_ptr(:)
    real(kind=rp), allocatable :: loc_geo(:,:), glb_geo(:,:), key(:)
    real(kind=rp) :: geo(N_GEO), atol
    integer, allocatable :: geo_cnt(:), geo_displ(:), glb_offset(:)
    integer, allocatable :: sorted(:), el_glb(:), el_order(:), col_glb(:)
    integer, allocatable :: recv_glb(:)
    integer :: i, j, k, e, a, b, ap, bp, c, r, g, lx, lxy, n, nelv, ierr
    integer :: base, sa, sb, sh

    call this%free()
    call this%map_1d%init(coef, dir, tol)
    this%msh => coef%msh
    this%coef => coef
    this%dof => coef%dof
    this%dir = dir

    n = this%dof%size()
    nelv = this%msh%nelv
    lx = this%dof%Xh%lx
    lxy = this%dof%Xh%lxy
    this%lxy = lxy

    ! The 2D mesh consists of the elements at the first level, distributed
    ! in rank order.
    this%nelv_2d = 0
    do e = 1, nelv
       if (this%map_1d%el_lvl(e) .eq. 1) this%nelv_2d = this%nelv_2d + 1
    end do
    this%glb_nelv_2d = 0
    call MPI_Allreduce(this%nelv_2d, this%glb_nelv_2d, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    this%offset_el_2d = 0
    call MPI_Exscan(this%nelv_2d, this%offset_el_2d, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    allocate(this%el_idx_2d(this%nelv_2d))
    do i = 1, this%nelv_2d
       this%el_idx_2d(i) = this%offset_el_2d + i
    end do
    this%n_2d = this%nelv_2d * lxy
    allocate(this%idx_2d(this%n_2d))
    j = 0
    do e = 1, nelv
       if (this%map_1d%el_lvl(e) .eq. 1) then
          call element_strides(this, e, base, sa, sb, sh)
          do b = 1, lx
             do a = 1, lx
                this%idx_2d(j * lxy + a + lx * (b - 1)) = &
                     base + (a - 1) * sa + (b - 1) * sb
             end do
          end do
          j = j + 1
       end if
    end do

    ! Tables of the in-plane node permutations.
    allocate(this%perm_tbl(lxy, N_PERM))
    do k = 1, N_PERM
       do b = 1, lx
          do a = 1, lx
             call dihedral_map(k, a, b, lx, ap, bp)
             this%perm_tbl(a + lx * (b - 1), k) = ap + lx * (bp - 1)
          end do
       end do
    end do

    ! Gather the in-plane geometry of all 2D elements on every rank. Entry
    ! g describes the 2D element with global index g, since the ranks
    ! contribute their elements in rank order.
    call inplane_coords(this, x_ptr, y_ptr)
    atol = tol * max(glmax(x_ptr, n) - glmin(x_ptr, n), &
         glmax(y_ptr, n) - glmin(y_ptr, n))

    allocate(loc_geo(N_GEO, max(this%nelv_2d, 1)))
    j = 0
    do e = 1, nelv
       if (this%map_1d%el_lvl(e) .eq. 1) then
          j = j + 1
          call element_geometry(this, e, x_ptr, y_ptr, loc_geo(:, j))
       end if
    end do

    allocate(geo_cnt(0:pe_size - 1), geo_displ(0:pe_size - 1))
    allocate(glb_offset(0:pe_size - 1))
    call MPI_Allgather(this%nelv_2d, 1, MPI_INTEGER, glb_offset, 1, &
         MPI_INTEGER, NEKO_COMM, ierr)
    geo_displ(0) = 0
    do r = 1, pe_size - 1
       geo_displ(r) = geo_displ(r - 1) + N_GEO * glb_offset(r - 1)
    end do
    geo_cnt = N_GEO * glb_offset
    ! Number of 2D elements on the ranks before each rank.
    glb_offset = geo_displ / N_GEO

    allocate(glb_geo(N_GEO, this%glb_nelv_2d))
    call MPI_Allgatherv(loc_geo, N_GEO * this%nelv_2d, MPI_REAL_PRECISION, &
         glb_geo, geo_cnt, geo_displ, MPI_REAL_PRECISION, NEKO_COMM, ierr)

    ! Sort the 2D elements by their first centroid coordinate for lookup.
    allocate(key(this%glb_nelv_2d), sorted(this%glb_nelv_2d))
    key = glb_geo(1, :)
    call sort(key, sorted, this%glb_nelv_2d)

    ! Assign every element to a 2D element and find the permutation of its
    ! in-plane nodes.
    allocate(el_glb(nelv), this%el_perm(nelv))
    do e = 1, nelv
       call element_geometry(this, e, x_ptr, y_ptr, geo)
       g = find_column(key, sorted, glb_geo, this%glb_nelv_2d, &
            geo(1), geo(2), atol)
       if (g .eq. 0) then
          call neko_error('map_2d: element not aligned with any 2D element' &
               // ', the mesh must be stacked in the homogeneous direction')
       end if
       el_glb(e) = g
       this%el_perm(e) = find_permutation(geo, glb_geo(:, g), lx, atol)
       if (this%el_perm(e) .eq. 0) then
          call neko_error('map_2d: element corners do not match its 2D' &
               // ' element, the mesh must be stacked in the homogeneous' &
               // ' direction')
       end if
    end do

    ! Group the elements into columns, ordered by the global index of the
    ! 2D element. Since 2D elements are distributed in rank order, this also
    ! orders the columns by owning rank.
    allocate(el_order(nelv), col_glb(max(nelv, 1)), this%el_col(nelv))
    call sort(el_glb, el_order, nelv)
    this%n_cols = 0
    do i = 1, nelv
       if (i .eq. 1) then
          this%n_cols = this%n_cols + 1
          col_glb(this%n_cols) = el_glb(i)
       else if (el_glb(i) .ne. el_glb(i - 1)) then
          this%n_cols = this%n_cols + 1
          col_glb(this%n_cols) = el_glb(i)
       end if
       this%el_col(el_order(i)) = this%n_cols
    end do

    ! Set up the exchange of column sums with the owning ranks.
    allocate(this%send_cols(0:pe_size - 1), this%send_displ(0:pe_size - 1))
    allocate(this%recv_cols(0:pe_size - 1), this%recv_displ(0:pe_size - 1))
    this%send_cols = 0
    do c = 1, this%n_cols
       r = owner_rank(col_glb(c), glb_offset, pe_size)
       this%send_cols(r) = this%send_cols(r) + 1
    end do
    call MPI_Alltoall(this%send_cols, 1, MPI_INTEGER, &
         this%recv_cols, 1, MPI_INTEGER, NEKO_COMM, ierr)
    this%send_displ(0) = 0
    this%recv_displ(0) = 0
    do r = 1, pe_size - 1
       this%send_displ(r) = this%send_displ(r - 1) + this%send_cols(r - 1)
       this%recv_displ(r) = this%recv_displ(r - 1) + this%recv_cols(r - 1)
    end do
    this%n_recv = sum(this%recv_cols)

    ! Tell the owners which 2D element each contribution belongs to.
    allocate(recv_glb(max(this%n_recv, 1)), this%recv_el(this%n_recv))
    call MPI_Alltoallv(col_glb, this%send_cols, this%send_displ, MPI_INTEGER, &
         recv_glb, this%recv_cols, this%recv_displ, MPI_INTEGER, &
         NEKO_COMM, ierr)
    do k = 1, this%n_recv
       this%recv_el(k) = recv_glb(k) - this%offset_el_2d
       if (this%recv_el(k) .lt. 1 .or. this%recv_el(k) .gt. this%nelv_2d) then
          call neko_error('map_2d: received a column for a 2D element' &
               // ' owned by another rank')
       end if
    end do

    deallocate(loc_geo, glb_geo, key, sorted, geo_cnt, geo_displ, glb_offset)
    deallocate(el_glb, el_order, col_glb, recv_glb)

  end subroutine map_2d_init

  !> Constructor from a character direction.
  !! @param coef SEM coefficients of the 3D mesh.
  !! @param dir Direction normal to the 2D plane, 'x', 'y' or 'z'.
  !! @param tol Tolerance, relative to the in-plane extent of the domain,
  !! for comparing coordinates.
  subroutine map_2d_init_char(this, coef, dir, tol)
    class(map_2d_t) :: this
    type(coef_t), intent(inout), target :: coef
    character(len=*), intent(in) :: dir
    real(kind=rp), intent(in) :: tol
    integer :: idir

    if (trim(dir) .eq. 'x') then
       idir = 1
    else if (trim(dir) .eq. 'y') then
       idir = 2
    else if (trim(dir) .eq. 'z') then
       idir = 3
    else
       call neko_error('Direction not supported, map_2d')
    end if

    call this%init(coef, idir, tol)

  end subroutine map_2d_init_char

  !> Destructor.
  subroutine map_2d_free(this)
    class(map_2d_t), intent(inout) :: this

    if (allocated(this%idx_2d)) deallocate(this%idx_2d)
    if (allocated(this%el_idx_2d)) deallocate(this%el_idx_2d)
    if (allocated(this%el_col)) deallocate(this%el_col)
    if (allocated(this%el_perm)) deallocate(this%el_perm)
    if (allocated(this%perm_tbl)) deallocate(this%perm_tbl)
    if (allocated(this%send_cols)) deallocate(this%send_cols)
    if (allocated(this%send_displ)) deallocate(this%send_displ)
    if (allocated(this%recv_cols)) deallocate(this%recv_cols)
    if (allocated(this%recv_displ)) deallocate(this%recv_displ)
    if (allocated(this%recv_el)) deallocate(this%recv_el)

    call this%accumulate_free()
    call this%map_1d%free()

    nullify(this%msh)
    nullify(this%dof)
    nullify(this%coef)

    this%nelv_2d = 0
    this%glb_nelv_2d = 0
    this%offset_el_2d = 0
    this%lxy = 0
    this%n_2d = 0
    this%dir = 0
    this%n_cols = 0
    this%n_recv = 0

  end subroutine map_2d_free

  !> Computes the average of a list of fields in the homogeneous direction
  !! and outputs a 2D field with the averaged values.
  !! The fields are expected to be up to date on the host and are not
  !! modified.
  !! @param fld_data2D output 2D averages
  !! @param fld_data3D list of fields to be averaged
  subroutine map_2d_average_field_list(this, fld_data2D, fld_data3D)
    class(map_2d_t), intent(in) :: this
    type(fld_file_data_t), intent(inout) :: fld_data2D
    type(field_list_t), intent(in) :: fld_data3D
    type(field_ptr_1d_t), allocatable :: flds(:)
    type(field_t), pointer :: f
    integer :: i, n, nf

    nf = fld_data3D%size()
    n = this%dof%size()
    allocate(flds(nf))
    do i = 1, nf
       f => fld_data3D%items(i)%ptr
       flds(i)%x(1:n) => f%x
    end do

    call map_2d_average_fields(this, fld_data2D, flds, nf)

    deallocate(flds)

  end subroutine map_2d_average_field_list

  !> Computes the average of the fields in a `fld_file_data_t` in the
  !! homogeneous direction and outputs a 2D field with the averaged values.
  !! The fields are not modified.
  !! @param fld_data2D output 2D averages
  !! @param fld_data3D fld_file_data of fields to be averaged
  subroutine map_2d_average(this, fld_data2D, fld_data3D)
    class(map_2d_t), intent(in) :: this
    type(fld_file_data_t), intent(inout) :: fld_data2D
    type(fld_file_data_t), intent(in) :: fld_data3D
    type(field_ptr_1d_t), allocatable :: flds(:)
    type(vector_ptr_t), allocatable :: fields3d(:)
    integer :: i, nf

    nf = fld_data3D%size()
    allocate(fields3d(nf), flds(nf))
    call fld_data3D%get_list(fields3d, nf)
    do i = 1, nf
       flds(i)%x => fields3d(i)%ptr%x
    end do

    call map_2d_average_fields(this, fld_data2D, flds, nf)

    deallocate(fields3d, flds)

  end subroutine map_2d_average

  !> Averages `nf` fields in the homogeneous direction and stores the result
  !! in a 2D field.
  subroutine map_2d_average_fields(this, fld_data2D, flds, nf)
    class(map_2d_t), intent(in) :: this
    type(fld_file_data_t), intent(inout) :: fld_data2D
    integer, intent(in) :: nf
    type(field_ptr_1d_t), intent(in) :: flds(nf)
    real(kind=rp), allocatable :: avg(:,:,:)

    allocate(avg(this%lxy, nf, this%nelv_2d))
    call map_2d_column_average(this, flds, nf, avg)
    call map_2d_fill_output(this, fld_data2D, avg, nf)
    deallocate(avg)

  end subroutine map_2d_average_fields

  !> Sets up a 2D field and fills it with the averages on the 2D elements
  !! of this rank.
  !! @param fld_data2D Output 2D field.
  !! @param avg Averages on the 2D elements of this rank.
  !! @param nf Number of fields.
  subroutine map_2d_fill_output(this, fld_data2D, avg, nf)
    class(map_2d_t), intent(in) :: this
    type(fld_file_data_t), intent(inout) :: fld_data2D
    integer, intent(in) :: nf
    real(kind=rp), intent(in) :: avg(this%lxy, nf, this%nelv_2d)
    type(vector_ptr_t), allocatable :: fields2d(:)
    integer :: i, e, m, lxy

    lxy = this%lxy

    call fld_data2D%init(this%nelv_2d, this%offset_el_2d)
    fld_data2D%gdim = 2
    fld_data2D%lx = this%dof%Xh%lx
    fld_data2D%ly = this%dof%Xh%ly
    fld_data2D%lz = 1
    fld_data2D%glb_nelv = this%glb_nelv_2d
    call map_2d_output_coords(this, fld_data2D)
    call fld_data2D%init_n_fields(nf, this%n_2d)

    allocate(fields2d(nf))
    call fld_data2D%get_list(fields2d, nf)
    do i = 1, nf
       do e = 1, this%nelv_2d
          do m = 1, lxy
             fields2d(i)%ptr%x(m + lxy * (e - 1)) = avg(m, i, e)
          end do
       end do
    end do

    deallocate(fields2d)

  end subroutine map_2d_fill_output

  !> Sets the element indices and the in-plane coordinates of a 2D field.
  subroutine map_2d_output_coords(this, fld_data2D)
    class(map_2d_t), intent(in) :: this
    type(fld_file_data_t), intent(inout) :: fld_data2D
    real(kind=rp), contiguous, pointer :: x_ptr(:), y_ptr(:)
    integer :: j

    call fld_data2D%x%init(this%n_2d)
    call fld_data2D%y%init(this%n_2d)
    allocate(fld_data2D%idx(this%nelv_2d))
    do j = 1, this%nelv_2d
       fld_data2D%idx(j) = this%el_idx_2d(j)
    end do

    call inplane_coords(this, x_ptr, y_ptr)
    do j = 1, this%n_2d
       fld_data2D%x%x(j) = x_ptr(this%idx_2d(j))
       fld_data2D%y%x(j) = y_ptr(this%idx_2d(j))
    end do

  end subroutine map_2d_output_coords

  !> Sums the fields over the columns, weighted with the mass matrix, and
  !! reduces the sums to the ranks owning the 2D elements, where they are
  !! divided by the column volumes.
  !! @param flds Fields to average, up to date on the host.
  !! @param nf Number of fields.
  !! @param avg Averages on the 2D elements owned by this rank.
  subroutine map_2d_column_average(this, flds, nf, avg)
    class(map_2d_t), intent(in) :: this
    integer, intent(in) :: nf
    type(field_ptr_1d_t), intent(in) :: flds(nf)
    real(kind=rp), intent(out) :: avg(this%lxy, nf, this%nelv_2d)
    real(kind=rp), allocatable :: acc(:,:,:)
    real(kind=rp), contiguous, pointer :: wt(:)
    real(kind=rp) :: s
    integer :: e, c, code, j, a, b, h, m, p, p0, lx, lxy, n, nelv
    integer :: base, sa, sb, sh

    lx = this%dof%Xh%lx
    lxy = this%lxy
    n = this%dof%size()
    nelv = this%msh%nelv
    wt(1:n) => this%coef%B

    ! Column sums of the mass matrix (index 0) and of each weighted field.
    allocate(acc(lxy, 0:nf, this%n_cols))
    acc = 0.0_rp
    do e = 1, nelv
       c = this%el_col(e)
       code = this%el_perm(e)
       call element_strides(this, e, base, sa, sb, sh)
       do b = 1, lx
          do a = 1, lx
             m = this%perm_tbl(a + lx * (b - 1), code)
             p0 = base + (a - 1) * sa + (b - 1) * sb
             s = 0.0_rp
             do h = 1, lx
                s = s + wt(p0 + (h - 1) * sh)
             end do
             acc(m, 0, c) = acc(m, 0, c) + s
          end do
       end do
       do j = 1, nf
          do b = 1, lx
             do a = 1, lx
                m = this%perm_tbl(a + lx * (b - 1), code)
                p0 = base + (a - 1) * sa + (b - 1) * sb
                s = 0.0_rp
                do h = 1, lx
                   p = p0 + (h - 1) * sh
                   s = s + flds(j)%x(p) * wt(p)
                end do
                acc(m, j, c) = acc(m, j, c) + s
             end do
          end do
       end do
    end do

    call map_2d_reduce_columns(this, acc, nf, avg)

    deallocate(acc)

  end subroutine map_2d_column_average

  !> Sends the column sums to the ranks owning the 2D elements, where they
  !! are summed and divided by the column volumes.
  !! @param sums Column sums of the local columns, slot 0 holding the
  !! column volumes.
  !! @param nf Number of fields.
  !! @param avg Averages on the 2D elements owned by this rank.
  subroutine map_2d_reduce_columns(this, sums, nf, avg)
    class(map_2d_t), intent(in) :: this
    integer, intent(in) :: nf
    real(kind=rp), intent(in) :: sums(this%lxy, 0:nf, this%n_cols)
    real(kind=rp), intent(out) :: avg(this%lxy, nf, this%nelv_2d)
    real(kind=rp), allocatable :: recv(:,:,:), vol(:,:)
    integer, allocatable :: send_cnt(:), send_dsp(:), recv_cnt(:), recv_dsp(:)
    integer :: e, j, k, lxy, slab, ierr

    lxy = this%lxy
    slab = lxy * (nf + 1)
    allocate(send_cnt(0:pe_size - 1), send_dsp(0:pe_size - 1))
    allocate(recv_cnt(0:pe_size - 1), recv_dsp(0:pe_size - 1))
    send_cnt = this%send_cols * slab
    send_dsp = this%send_displ * slab
    recv_cnt = this%recv_cols * slab
    recv_dsp = this%recv_displ * slab
    allocate(recv(lxy, 0:nf, max(this%n_recv, 1)))
    call MPI_Alltoallv(sums, send_cnt, send_dsp, MPI_REAL_PRECISION, &
         recv, recv_cnt, recv_dsp, MPI_REAL_PRECISION, NEKO_COMM, ierr)

    allocate(vol(lxy, this%nelv_2d))
    vol = 0.0_rp
    avg = 0.0_rp
    do k = 1, this%n_recv
       e = this%recv_el(k)
       vol(:, e) = vol(:, e) + recv(:, 0, k)
       do j = 1, nf
          avg(:, j, e) = avg(:, j, e) + recv(:, j, k)
       end do
    end do
    do e = 1, this%nelv_2d
       if (any(vol(:, e) .le. 0.0_rp)) then
          call neko_error('map_2d: 2D element with a non-positive volume')
       end if
       do j = 1, nf
          avg(:, j, e) = avg(:, j, e) / vol(:, e)
       end do
    end do

    deallocate(recv, vol, send_cnt, send_dsp, recv_cnt, recv_dsp)

  end subroutine map_2d_reduce_columns

  !> Prepares the accumulation of `n_fields` fields averaged in the
  !! homogeneous direction, see `accumulate`.
  !! @param n_fields Number of fields accumulated per sample.
  subroutine map_2d_accumulate_init(this, n_fields)
    class(map_2d_t), intent(inout) :: this
    integer, intent(in) :: n_fields
    integer, allocatable :: cnt(:)
    integer :: e, c, k, m, lxy, nelv, nrows, base, sa, sb, sh

    call this%accumulate_free()
    this%n_acc = n_fields
    lxy = this%lxy
    nelv = this%msh%nelv

    allocate(this%el_sa(nelv), this%el_sb(nelv), this%el_sh(nelv))
    allocate(this%el_code0(nelv))
    do e = 1, nelv
       call element_strides(this, e, base, sa, sb, sh)
       this%el_sa(e) = sa
       this%el_sb(e) = sb
       this%el_sh(e) = sh
       this%el_code0(e) = this%el_perm(e) - 1
    end do
    allocate(this%tbl0(lxy * N_PERM))
    do k = 1, N_PERM
       do m = 1, lxy
          this%tbl0(m + lxy * (k - 1)) = this%perm_tbl(m, k) - 1
       end do
    end do

    ! Row (m, c) of the column sums gathers node m of the slab sums of the
    ! elements in column c, in element order.
    nrows = lxy * this%n_cols
    allocate(this%csr_ptr(nrows + 1), this%csr_list(max(lxy * nelv, 1)))
    allocate(cnt(max(this%n_cols, 1)))
    cnt = 0
    do e = 1, nelv
       cnt(this%el_col(e)) = cnt(this%el_col(e)) + 1
    end do
    this%csr_ptr(1) = 0
    do c = 1, this%n_cols
       do m = 1, lxy
          k = m + lxy * (c - 1)
          this%csr_ptr(k + 1) = this%csr_ptr(k) + cnt(c)
       end do
    end do
    cnt = 0
    do e = 1, nelv
       c = this%el_col(e)
       do m = 1, lxy
          k = m + lxy * (c - 1)
          this%csr_list(this%csr_ptr(k) + cnt(c) + 1) = &
               (m - 1) + lxy * (e - 1)
       end do
       cnt(c) = cnt(c) + 1
    end do
    deallocate(cnt)

    allocate(this%tmp(max(lxy * nelv, 1)))
    allocate(this%acc(lxy, max(this%n_cols, 1), 0:n_fields))

    if (NEKO_BCKND_DEVICE .eq. 1 .and. nelv .gt. 0) then
       call device_map(this%el_sa, this%el_sa_d, nelv)
       call device_map(this%el_sb, this%el_sb_d, nelv)
       call device_map(this%el_sh, this%el_sh_d, nelv)
       call device_map(this%el_code0, this%el_code0_d, nelv)
       call device_map(this%tbl0, this%tbl0_d, size(this%tbl0))
       call device_map(this%csr_ptr, this%csr_ptr_d, size(this%csr_ptr))
       call device_map(this%csr_list, this%csr_list_d, size(this%csr_list))
       call device_map(this%tmp, this%tmp_d, size(this%tmp))
       call device_map(this%acc, this%acc_d, size(this%acc))
       call device_memcpy(this%el_sa, this%el_sa_d, nelv, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%el_sb, this%el_sb_d, nelv, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%el_sh, this%el_sh_d, nelv, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(this%el_code0, this%el_code0_d, nelv, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%tbl0, this%tbl0_d, size(this%tbl0), &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%csr_ptr, this%csr_ptr_d, size(this%csr_ptr), &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%csr_list, this%csr_list_d, &
            size(this%csr_list), HOST_TO_DEVICE, sync = .false.)
       ! Zero the work and accumulator on the device first: on unified
       ! memory the device then faults the pages (device first touch), and
       ! the host must not write to them while device work is in flight.
       call device_rzero(this%tmp_d, size(this%tmp))
       call device_rzero(this%acc_d, size(this%acc))
       call device_sync()
    end if
    this%tmp = 0.0_rp
    this%acc = 0.0_rp

  end subroutine map_2d_accumulate_init

  !> Frees the accumulation data.
  subroutine map_2d_accumulate_free(this)
    class(map_2d_t), intent(inout) :: this

    if (c_associated(this%el_sa_d)) call device_unmap(this%el_sa, this%el_sa_d)
    if (c_associated(this%el_sb_d)) call device_unmap(this%el_sb, this%el_sb_d)
    if (c_associated(this%el_sh_d)) call device_unmap(this%el_sh, this%el_sh_d)
    if (c_associated(this%el_code0_d)) then
       call device_unmap(this%el_code0, this%el_code0_d)
    end if
    if (c_associated(this%tbl0_d)) call device_unmap(this%tbl0, this%tbl0_d)
    if (c_associated(this%csr_ptr_d)) then
       call device_unmap(this%csr_ptr, this%csr_ptr_d)
    end if
    if (c_associated(this%csr_list_d)) then
       call device_unmap(this%csr_list, this%csr_list_d)
    end if
    if (c_associated(this%tmp_d)) call device_unmap(this%tmp, this%tmp_d)
    if (c_associated(this%acc_d)) call device_unmap(this%acc, this%acc_d)
    this%el_sa_d = C_NULL_PTR
    this%el_sb_d = C_NULL_PTR
    this%el_sh_d = C_NULL_PTR
    this%el_code0_d = C_NULL_PTR
    this%tbl0_d = C_NULL_PTR
    this%csr_ptr_d = C_NULL_PTR
    this%csr_list_d = C_NULL_PTR
    this%tmp_d = C_NULL_PTR
    this%acc_d = C_NULL_PTR

    if (allocated(this%el_sa)) deallocate(this%el_sa)
    if (allocated(this%el_sb)) deallocate(this%el_sb)
    if (allocated(this%el_sh)) deallocate(this%el_sh)
    if (allocated(this%el_code0)) deallocate(this%el_code0)
    if (allocated(this%tbl0)) deallocate(this%tbl0)
    if (allocated(this%csr_ptr)) deallocate(this%csr_ptr)
    if (allocated(this%csr_list)) deallocate(this%csr_list)
    if (allocated(this%tmp)) deallocate(this%tmp)
    if (allocated(this%acc)) deallocate(this%acc)
    this%n_acc = 0

  end subroutine map_2d_accumulate_free

  !> Adds the column sums of `k * f * B` to slot `slot` of the accumulated
  !! column sums, where `B` is the mass matrix. The field is expected to be
  !! up to date on the device, or on the host for the CPU backend.
  !! @param f Field to accumulate.
  !! @param slot Slot of the field, 1 to `n_acc`.
  !! @param k Scaling of the sample, typically the time since the last one.
  subroutine map_2d_accumulate(this, f, slot, k)
    class(map_2d_t), intent(inout) :: this
    type(field_t), intent(in) :: f
    integer, intent(in) :: slot
    real(kind=rp), intent(in) :: k

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call map_2d_accumulate_device(this, f%x_d, slot, k)
    else
       call map_2d_accumulate_host(this, slot, k, f%x)
    end if

  end subroutine map_2d_accumulate

  !> Adds the column volumes, `k * B` summed over the columns, to slot 0 of
  !! the accumulated column sums.
  !! @param k Scaling of the sample, typically the time since the last one.
  subroutine map_2d_accumulate_volume(this, k)
    class(map_2d_t), intent(inout) :: this
    real(kind=rp), intent(in) :: k

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call map_2d_accumulate_device(this, C_NULL_PTR, 0, k)
    else
       call map_2d_accumulate_host(this, 0, k)
    end if

  end subroutine map_2d_accumulate_volume

  !> Device implementation of the accumulation, a null `f_d` counts as one.
  subroutine map_2d_accumulate_device(this, f_d, slot, k)
    class(map_2d_t), intent(inout) :: this
    type(c_ptr), intent(in) :: f_d
    integer, intent(in) :: slot
    real(kind=rp), intent(in) :: k
    integer :: lx, lxy, nelv, nrows

    lx = this%dof%Xh%lx
    lxy = this%lxy
    nelv = this%msh%nelv
    nrows = lxy * this%n_cols

    call device_slab_sum(this%tmp_d, f_d, this%coef%B_d, this%el_sa_d, &
         this%el_sb_d, this%el_sh_d, this%el_code0_d, this%tbl0_d, lxy, lx, &
         nelv, this%dof%Xh%lxyz, lxy)
    call device_gather_add(this%acc_d, slot * nrows, this%tmp_d, &
         this%csr_ptr_d, this%csr_list_d, nrows, k)

  end subroutine map_2d_accumulate_device

  !> Host implementation of the accumulation, a missing `f` counts as one.
  subroutine map_2d_accumulate_host(this, slot, k, f)
    class(map_2d_t), intent(inout) :: this
    integer, intent(in) :: slot
    real(kind=rp), intent(in) :: k
    real(kind=rp), intent(in), optional :: f(*)
    real(kind=rp), contiguous, pointer :: wt(:)
    real(kind=rp) :: s
    integer :: e, c, m, a, b, h, p, p0, r, j, lx, lxy, lxyz, n, nelv, code

    lx = this%dof%Xh%lx
    lxy = this%lxy
    lxyz = this%dof%Xh%lxyz
    n = this%dof%size()
    nelv = this%msh%nelv
    wt(1:n) => this%coef%B

    do e = 1, nelv
       code = this%el_code0(e)
       do m = 0, lxy - 1
          a = mod(m, lx)
          b = m / lx
          p0 = (e - 1) * lxyz + a * this%el_sa(e) + b * this%el_sb(e)
          s = 0.0_rp
          if (present(f)) then
             do h = 0, lx - 1
                p = p0 + h * this%el_sh(e) + 1
                s = s + f(p) * wt(p)
             end do
          else
             do h = 0, lx - 1
                p = p0 + h * this%el_sh(e) + 1
                s = s + wt(p)
             end do
          end if
          this%tmp((e - 1) * lxy + this%tbl0(m + 1 + lxy * code) + 1) = s
       end do
    end do

    do c = 1, this%n_cols
       do m = 1, lxy
          r = m + lxy * (c - 1)
          s = 0.0_rp
          do j = this%csr_ptr(r) + 1, this%csr_ptr(r + 1)
             s = s + this%tmp(this%csr_list(j) + 1)
          end do
          this%acc(m, c, slot) = this%acc(m, c, slot) + k * s
       end do
    end do

  end subroutine map_2d_accumulate_host

  !> Resets the accumulated column sums to zero.
  subroutine map_2d_accumulate_reset(this)
    class(map_2d_t), intent(inout) :: this

    ! Device first, see accumulate_init.
    if (c_associated(this%acc_d)) then
       call device_rzero(this%acc_d, size(this%acc))
       call device_sync()
    end if
    if (allocated(this%acc)) this%acc = 0.0_rp

  end subroutine map_2d_accumulate_reset

  !> Divides the accumulated column sums by the accumulated column volumes
  !! and outputs the averages as a 2D field, in the order of the slots.
  !! @param fld_data2D Output 2D averages.
  subroutine map_2d_accumulated_average(this, fld_data2D)
    class(map_2d_t), intent(inout) :: this
    type(fld_file_data_t), intent(inout) :: fld_data2D
    real(kind=rp), allocatable :: sums(:,:,:), avg(:,:,:)
    integer :: c, j, nf

    nf = this%n_acc
    if (c_associated(this%acc_d)) then
       call device_memcpy(this%acc, this%acc_d, size(this%acc), &
            DEVICE_TO_HOST, sync = .true.)
    end if

    allocate(sums(this%lxy, 0:nf, this%n_cols))
    do c = 1, this%n_cols
       do j = 0, nf
          sums(:, j, c) = this%acc(:, c, j)
       end do
    end do

    allocate(avg(this%lxy, nf, this%nelv_2d))
    call map_2d_reduce_columns(this, sums, nf, avg)
    call map_2d_fill_output(this, fld_data2D, avg, nf)

    deallocate(sums, avg)

  end subroutine map_2d_accumulated_average

  !> Pointers to the two in-plane coordinates of the dofmap, in the order
  !! they are written to the 2D output.
  subroutine inplane_coords(this, x_ptr, y_ptr)
    class(map_2d_t), intent(in) :: this
    real(kind=rp), contiguous, pointer, intent(out) :: x_ptr(:), y_ptr(:)
    integer :: n

    n = this%dof%size()
    if (this%dir .eq. 1) then
       x_ptr(1:n) => this%dof%z%x
       y_ptr(1:n) => this%dof%y%x
    else if (this%dir .eq. 2) then
       x_ptr(1:n) => this%dof%x%x
       y_ptr(1:n) => this%dof%z%x
    else
       x_ptr(1:n) => this%dof%x%x
       y_ptr(1:n) => this%dof%y%x
    end if

  end subroutine inplane_coords

  !> Index offsets of the nodes of element `e`, such that the node with
  !! in-plane indices `(a, b)` and index `h` in the homogeneous direction is
  !! `base + (a - 1) * sa + (b - 1) * sb + (h - 1) * sh`. The in-plane
  !! indices are the two local directions that are not homogeneous, in
  !! increasing order.
  subroutine element_strides(this, e, base, sa, sb, sh)
    class(map_2d_t), intent(in) :: this
    integer, intent(in) :: e
    integer, intent(out) :: base, sa, sb, sh
    integer :: lx, lxy

    lx = this%dof%Xh%lx
    lxy = this%dof%Xh%lxy
    base = (e - 1) * this%dof%Xh%lxyz + 1

    select case (this%map_1d%dir_el(e))
    case (1) ! r is homogeneous, (a, b) = (s, t)
       sa = lx
       sb = lxy
       sh = 1
    case (2) ! s is homogeneous, (a, b) = (r, t)
       sa = 1
       sb = lxy
       sh = lx
    case default ! t is homogeneous, (a, b) = (r, s)
       sa = 1
       sb = lx
       sh = lxy
    end select

  end subroutine element_strides

  !> In-plane geometry of element `e`: `geo(1:2)` is the mean of the four
  !! in-plane corner nodes and `geo(2c+1:2c+2)` the coordinates of corner
  !! `c`, where the corners are numbered (1,1), (lx,1), (1,lx) and (lx,lx)
  !! in the in-plane node indices.
  subroutine element_geometry(this, e, x_ptr, y_ptr, geo)
    class(map_2d_t), intent(in) :: this
    integer, intent(in) :: e
    real(kind=rp), intent(in) :: x_ptr(:), y_ptr(:)
    real(kind=rp), intent(out) :: geo(N_GEO)
    integer :: base, sa, sb, sh, c, a, b, p, lx

    lx = this%dof%Xh%lx
    call element_strides(this, e, base, sa, sb, sh)

    geo(1:2) = 0.0_rp
    do c = 1, 4
       a = 1 + (lx - 1) * mod(c - 1, 2)
       b = 1 + (lx - 1) * ((c - 1) / 2)
       p = base + (a - 1) * sa + (b - 1) * sb
       geo(2 * c + 1) = x_ptr(p)
       geo(2 * c + 2) = y_ptr(p)
       geo(1) = geo(1) + 0.25_rp * x_ptr(p)
       geo(2) = geo(2) + 0.25_rp * y_ptr(p)
    end do

  end subroutine element_geometry

  !> Finds the 2D element whose centroid is within `atol` of `(cx, cy)`,
  !! using the elements sorted by their first centroid coordinate.
  !! Returns 0 if there is none.
  function find_column(key, sorted, glb_geo, n, cx, cy, atol) result(g)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: key(n), glb_geo(N_GEO, n), cx, cy, atol
    integer, intent(in) :: sorted(n)
    integer :: g, i, lo, hi, mid

    ! First element with key >= cx - atol.
    lo = 1
    hi = n + 1
    do while (lo .lt. hi)
       mid = (lo + hi) / 2
       if (key(mid) .lt. cx - atol) then
          lo = mid + 1
       else
          hi = mid
       end if
    end do

    g = 0
    do i = lo, n
       if (key(i) .gt. cx + atol) exit
       if (abs(glb_geo(2, sorted(i)) - cy) .le. atol) then
          g = sorted(i)
          return
       end if
    end do

  end function find_column

  !> Finds the permutation of the in-plane node indices that maps the
  !! corners of an element onto the corners of its 2D element.
  !! Returns 0 if none matches.
  function find_permutation(geo, ref, lx, atol) result(code)
    real(kind=rp), intent(in) :: geo(N_GEO), ref(N_GEO), atol
    integer, intent(in) :: lx
    integer :: code, c, cp, a, b, ap, bp
    logical :: match

    do code = 1, N_PERM
       match = .true.
       do c = 1, 4
          a = 1 + (lx - 1) * mod(c - 1, 2)
          b = 1 + (lx - 1) * ((c - 1) / 2)
          call dihedral_map(code, a, b, lx, ap, bp)
          cp = 1
          if (ap .eq. lx) cp = cp + 1
          if (bp .eq. lx) cp = cp + 2
          if (abs(geo(2 * c + 1) - ref(2 * cp + 1)) .gt. atol .or. &
               abs(geo(2 * c + 2) - ref(2 * cp + 2)) .gt. atol) then
             match = .false.
             exit
          end if
       end do
       if (match) return
    end do
    code = 0

  end function find_permutation

  !> The eight symmetries of a square applied to in-plane node indices.
  pure subroutine dihedral_map(code, a, b, lx, ap, bp)
    integer, intent(in) :: code, a, b, lx
    integer, intent(out) :: ap, bp

    select case (code)
    case (1)
       ap = a
       bp = b
    case (2)
       ap = lx + 1 - a
       bp = b
    case (3)
       ap = a
       bp = lx + 1 - b
    case (4)
       ap = lx + 1 - a
       bp = lx + 1 - b
    case (5)
       ap = b
       bp = a
    case (6)
       ap = lx + 1 - b
       bp = a
    case (7)
       ap = b
       bp = lx + 1 - a
    case default
       ap = lx + 1 - b
       bp = lx + 1 - a
    end select

  end subroutine dihedral_map

  !> Rank owning the 2D element with global index `g`, given the number of
  !! 2D elements on the ranks before each rank.
  pure function owner_rank(g, glb_offset, nranks) result(r)
    integer, intent(in) :: g, nranks
    integer, intent(in) :: glb_offset(0:nranks - 1)
    integer :: r, lo, hi, mid

    ! Last rank with glb_offset < g.
    lo = 0
    hi = nranks - 1
    do while (lo .lt. hi)
       mid = (lo + hi + 1) / 2
       if (glb_offset(mid) .lt. g) then
          lo = mid
       else
          hi = mid - 1
       end if
    end do
    r = lo

  end function owner_rank

end module map_2d
