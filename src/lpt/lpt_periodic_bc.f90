! Copyright (c) 2026, The Neko Authors
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
!> Periodic and cyclic boundary-condition support for LPT.
module lpt_periodic_bc
  use num_types, only : rp
  use neko_config, only : NEKO_BCKND_DEVICE
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use mesh, only : mesh_t
  use comm, only : pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_Allgather, MPI_INTEGER, &
       MPI_MIN, MPI_MAX
  use vector, only : vector_t
  use lpt_periodic_bc_cpu, only : lpt_periodic_bc_wrap_rotational_cpu, &
       lpt_periodic_bc_wrap_translational_cpu
  use lpt_periodic_bc_device, only : lpt_periodic_bc_wrap_rotational_device, &
       lpt_periodic_bc_wrap_translational_device
  use device, only : DEVICE_TO_HOST, device_sync
  implicit none
  private

  real(kind=rp), parameter :: LPT_PERIODIC_TOL = 1.0e-8_rp

  type, public :: lpt_periodic_bc_t
     logical :: periodic_enabled = .false.
     logical :: rotational_periodic_enabled = .false.
     integer :: n_periodic_dirs = 0
     real(kind=rp) :: periodic_shift_x1 = 0.0_rp
     real(kind=rp) :: periodic_shift_y1 = 0.0_rp
     real(kind=rp) :: periodic_shift_z1 = 0.0_rp
     real(kind=rp) :: periodic_shift_x2 = 0.0_rp
     real(kind=rp) :: periodic_shift_y2 = 0.0_rp
     real(kind=rp) :: periodic_shift_z2 = 0.0_rp
     real(kind=rp) :: periodic_shift_x3 = 0.0_rp
     real(kind=rp) :: periodic_shift_y3 = 0.0_rp
     real(kind=rp) :: periodic_shift_z3 = 0.0_rp
     real(kind=rp) :: periodic_dir_x1 = 0.0_rp
     real(kind=rp) :: periodic_dir_y1 = 0.0_rp
     real(kind=rp) :: periodic_dir_z1 = 0.0_rp
     real(kind=rp) :: periodic_dir_x2 = 0.0_rp
     real(kind=rp) :: periodic_dir_y2 = 0.0_rp
     real(kind=rp) :: periodic_dir_z2 = 0.0_rp
     real(kind=rp) :: periodic_dir_x3 = 0.0_rp
     real(kind=rp) :: periodic_dir_y3 = 0.0_rp
     real(kind=rp) :: periodic_dir_z3 = 0.0_rp
     real(kind=rp) :: periodic_min1 = 0.0_rp
     real(kind=rp) :: periodic_min2 = 0.0_rp
     real(kind=rp) :: periodic_min3 = 0.0_rp
     real(kind=rp) :: periodic_max1 = 0.0_rp
     real(kind=rp) :: periodic_max2 = 0.0_rp
     real(kind=rp) :: periodic_max3 = 0.0_rp
     real(kind=rp) :: periodic_len1 = 0.0_rp
     real(kind=rp) :: periodic_len2 = 0.0_rp
     real(kind=rp) :: periodic_len3 = 0.0_rp
     real(kind=rp) :: rotational_theta_min = 0.0_rp
     real(kind=rp) :: rotational_theta_max = 0.0_rp
     real(kind=rp) :: rotational_theta_len = 0.0_rp
   contains
     procedure, pass(this) :: free => lpt_periodic_bc_free
     procedure, pass(this) :: init => lpt_periodic_bc_init
     procedure, pass(this) :: wrap => lpt_periodic_bc_wrap
     procedure, private, pass(this) :: reset => lpt_periodic_bc_reset
     procedure, private, pass(this) :: init_rotational => &
          lpt_periodic_bc_init_rotational
     procedure, private, pass(this) :: init_translational => &
          lpt_periodic_bc_init_translational
  end type lpt_periodic_bc_t

contains

  !> Free periodic-boundary wrapping metadata.
  subroutine lpt_periodic_bc_free(this)
    class(lpt_periodic_bc_t), intent(inout) :: this

    call this%reset()
  end subroutine lpt_periodic_bc_free

  !> Initialise rotational or translational periodic wrapping from the mesh.
  !! @param msh Mesh containing periodic facet pairing.
  !! @param dm_Xh Coordinate dofmap used for rotational extents.
  !! @param coef Coefficients used to identify cyclic geometry.
  subroutine lpt_periodic_bc_init(this, msh, dm_Xh, coef)
    class(lpt_periodic_bc_t), intent(inout) :: this
    type(mesh_t), intent(in), target :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef

    call this%reset()

    if (msh%periodic%size .eq. 0) return

    call this%init_rotational(msh, dm_Xh, coef)
    if (this%rotational_periodic_enabled) return

    call this%init_translational(msh)
  end subroutine lpt_periodic_bc_init

  !> Wrap particle positions and optional vector histories across boundaries.
  !! @param x Particle x coordinates.
  !! @param y Particle y coordinates.
  !! @param z Particle z coordinates.
  !! @param n Number of local particles.
  subroutine lpt_periodic_bc_wrap(this, x, y, z, n, u, v, w, u_lag, v_lag, &
       w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, &
       acc_xlaglag, acc_ylaglag, acc_zlaglag)
    class(lpt_periodic_bc_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x, y, z
    integer, intent(in) :: n
    type(vector_t), intent(inout), optional :: u, v, w
    type(vector_t), intent(inout), optional :: u_lag, v_lag, w_lag
    type(vector_t), intent(inout), optional :: u_laglag, v_laglag, w_laglag
    type(vector_t), intent(inout), optional :: acc_xlag, acc_ylag, acc_zlag
    type(vector_t), intent(inout), optional :: acc_xlaglag
    type(vector_t), intent(inout), optional :: acc_ylaglag
    type(vector_t), intent(inout), optional :: acc_zlaglag

    if (n .eq. 0) return
    if (.not. this%periodic_enabled .and. &
         .not. this%rotational_periodic_enabled) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       if (this%rotational_periodic_enabled) then
          call lpt_periodic_bc_wrap_rotational_device(x, y, z, n, &
               this%rotational_theta_min, this%rotational_theta_max, &
               this%rotational_theta_len, u, v, w, u_lag, v_lag, w_lag, &
               u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, &
               acc_xlaglag, acc_ylaglag, acc_zlaglag)
       else
          call lpt_periodic_bc_wrap_translational_device(x, y, z, n, &
               this%n_periodic_dirs, this%periodic_dir_x1, &
               this%periodic_dir_y1, this%periodic_dir_z1, &
               this%periodic_dir_x2, this%periodic_dir_y2, &
               this%periodic_dir_z2, this%periodic_dir_x3, &
               this%periodic_dir_y3, this%periodic_dir_z3, &
               this%periodic_min1, this%periodic_min2, this%periodic_min3, &
               this%periodic_max1, this%periodic_max2, this%periodic_max3, &
               this%periodic_shift_x1, this%periodic_shift_y1, &
               this%periodic_shift_z1, this%periodic_shift_x2, &
               this%periodic_shift_y2, this%periodic_shift_z2, &
               this%periodic_shift_x3, this%periodic_shift_y3, &
               this%periodic_shift_z3, this%periodic_len1, &
               this%periodic_len2, this%periodic_len3)
       end if
       call lpt_periodic_bc_sync_from_device(x, y, z, u, v, w, u_lag, v_lag, &
            w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, &
            acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag)
       return
    end if

    if (this%rotational_periodic_enabled) then
       call lpt_periodic_bc_wrap_rotational_cpu(x, y, z, n, &
            this%rotational_theta_min, this%rotational_theta_max, &
            this%rotational_theta_len, u, v, w, u_lag, v_lag, w_lag, &
            u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, &
            acc_xlaglag, acc_ylaglag, acc_zlaglag)
       return
    end if

    call lpt_periodic_bc_wrap_translational_cpu(x, y, z, n, &
         this%periodic_enabled, this%n_periodic_dirs, this%periodic_dir_x1, &
         this%periodic_dir_y1, this%periodic_dir_z1, this%periodic_dir_x2, &
         this%periodic_dir_y2, this%periodic_dir_z2, this%periodic_dir_x3, &
         this%periodic_dir_y3, this%periodic_dir_z3, this%periodic_min1, &
         this%periodic_min2, this%periodic_min3, this%periodic_max1, &
         this%periodic_max2, this%periodic_max3, this%periodic_shift_x1, &
         this%periodic_shift_y1, this%periodic_shift_z1, &
         this%periodic_shift_x2, this%periodic_shift_y2, &
         this%periodic_shift_z2, this%periodic_shift_x3, &
         this%periodic_shift_y3, this%periodic_shift_z3, &
         this%periodic_len1, this%periodic_len2, this%periodic_len3)
  end subroutine lpt_periodic_bc_wrap

  !> Copy wrapped particle coordinates and optional histories back to host.
  subroutine lpt_periodic_bc_sync_from_device(x, y, z, u, v, w, u_lag, &
       v_lag, w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, &
       acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag)
    type(vector_t), intent(inout) :: x, y, z
    type(vector_t), intent(inout), optional :: u, v, w
    type(vector_t), intent(inout), optional :: u_lag, v_lag, w_lag
    type(vector_t), intent(inout), optional :: u_laglag, v_laglag, w_laglag
    type(vector_t), intent(inout), optional :: acc_xlag, acc_ylag, acc_zlag
    type(vector_t), intent(inout), optional :: acc_xlaglag
    type(vector_t), intent(inout), optional :: acc_ylaglag
    type(vector_t), intent(inout), optional :: acc_zlaglag

    call x%copy_from(DEVICE_TO_HOST, .false.)
    call y%copy_from(DEVICE_TO_HOST, .false.)
    call z%copy_from(DEVICE_TO_HOST, .false.)
    if (present(u)) call u%copy_from(DEVICE_TO_HOST, .false.)
    if (present(v)) call v%copy_from(DEVICE_TO_HOST, .false.)
    if (present(w)) call w%copy_from(DEVICE_TO_HOST, .false.)
    if (present(u_lag)) call u_lag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(v_lag)) call v_lag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(w_lag)) call w_lag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(u_laglag)) call u_laglag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(v_laglag)) call v_laglag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(w_laglag)) call w_laglag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(acc_xlag)) call acc_xlag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(acc_ylag)) call acc_ylag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(acc_zlag)) call acc_zlag%copy_from(DEVICE_TO_HOST, .false.)
    if (present(acc_xlaglag)) then
       call acc_xlaglag%copy_from(DEVICE_TO_HOST, .false.)
    end if
    if (present(acc_ylaglag)) then
       call acc_ylaglag%copy_from(DEVICE_TO_HOST, .false.)
    end if
    if (present(acc_zlaglag)) then
       call acc_zlaglag%copy_from(DEVICE_TO_HOST, .false.)
    end if
    call device_sync()
  end subroutine lpt_periodic_bc_sync_from_device

  !> Reset all periodic wrapping flags, directions, and extents.
  subroutine lpt_periodic_bc_reset(this)
    class(lpt_periodic_bc_t), intent(inout) :: this

    this%periodic_enabled = .false.
    this%rotational_periodic_enabled = .false.
    this%n_periodic_dirs = 0
    call lpt_periodic_bc_clear_translational(this)
    this%rotational_theta_min = 0.0_rp
    this%rotational_theta_max = 0.0_rp
    this%rotational_theta_len = 0.0_rp
  end subroutine lpt_periodic_bc_reset

  !> Detect and initialise rotational periodic wrapping.
  !! @param msh Mesh used to determine dimensionality.
  !! @param dm_Xh Coordinate dofmap used to compute angular extents.
  !! @param coef Coefficients carrying cyclic-geometry metadata.
  subroutine lpt_periodic_bc_init_rotational(this, msh, dm_Xh, coef)
    class(lpt_periodic_bc_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp) :: theta_min
    real(kind=rp) :: theta_max
    real(kind=rp) :: pi
    integer :: ierr

    if (.not. coef%cyclic .or. msh%gdim .lt. 2) return

    pi = acos(-1.0_rp)
    if (msh%nelv .gt. 0) then
       theta_min = minval(modulo(atan2(dm_Xh%y, dm_Xh%x) + &
            2.0_rp * pi, 2.0_rp * pi))
       theta_max = maxval(modulo(atan2(dm_Xh%y, dm_Xh%x) + &
            2.0_rp * pi, 2.0_rp * pi))
    else
       theta_min = huge(0.0_rp)
       theta_max = -huge(0.0_rp)
    end if

    call MPI_Allreduce(theta_min, this%rotational_theta_min, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
    call MPI_Allreduce(theta_max, this%rotational_theta_max, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)

    this%rotational_theta_len = this%rotational_theta_max - &
         this%rotational_theta_min
    this%rotational_periodic_enabled = &
         this%rotational_theta_len .gt. LPT_PERIODIC_TOL
  end subroutine lpt_periodic_bc_init_rotational

  !> Detect unique translational periodic directions from facet pairs.
  !! @param msh Mesh containing periodic facet-pair metadata.
  subroutine lpt_periodic_bc_init_translational(this, msh)
    class(lpt_periodic_bc_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    integer :: i
    integer :: j
    integer :: idx
    integer :: ierr
    integer :: n_local
    integer :: n_global
    integer :: max_local
    integer :: match_idx
    integer, allocatable :: counts(:)
    integer, allocatable :: padded_meta(:)
    integer, allocatable :: gathered_meta(:)
    integer, allocatable :: global_meta(:)
    real(kind=rp), allocatable :: padded_center(:)
    real(kind=rp), allocatable :: gathered_center(:)
    real(kind=rp), allocatable :: global_center(:)
    real(kind=rp) :: src_center(3)
    real(kind=rp) :: tgt_center(3)
    real(kind=rp) :: shift(3)
    real(kind=rp) :: proj_src
    real(kind=rp) :: proj_tgt
    real(kind=rp) :: dir(3)

    if (msh%periodic%size .eq. 0) return

    n_local = msh%periodic%size
    allocate(counts(pe_size))
    call MPI_Allgather(n_local, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, &
         NEKO_COMM, ierr)
    n_global = sum(counts)
    max_local = max(1, maxval(counts))

    allocate(padded_meta(2 * max_local))
    allocate(padded_center(3 * max_local))
    padded_meta = 0
    padded_center = 0.0_rp

    ! collect meta data from periodic facets
    do i = 1, n_local
       call lpt_get_periodic_center(msh, i, src_center)

       padded_meta(2 * (i - 1) + 1) = msh%periodic%facet_el(i)%x(1)
       padded_meta(2 * (i - 1) + 2) = msh%elements( &
            msh%periodic%facet_el(i)%x(2))%e%id()
       padded_center(3 * (i - 1) + 1:3 * i) = src_center
    end do

    ! broadcast meta data to all ranks
    allocate(gathered_meta(2 * max_local * pe_size))
    allocate(gathered_center(3 * max_local * pe_size))
    call MPI_Allgather(padded_meta, 2 * max_local, MPI_INTEGER, gathered_meta, &
         2 * max_local, MPI_INTEGER, NEKO_COMM, ierr)
    call MPI_Allgather(padded_center, 3 * max_local, MPI_REAL_PRECISION, &
         gathered_center, 3 * max_local, MPI_REAL_PRECISION, NEKO_COMM, ierr)

    ! reorganise the data to global arrays of size n_global
    allocate(global_meta(2 * n_global))
    allocate(global_center(3 * n_global))
    idx = 0
    do i = 1, pe_size
       if (counts(i) .gt. 0) then
          global_meta(2 * idx + 1:2 * (idx + counts(i))) = &
               gathered_meta(2 * max_local * (i - 1) + 1: &
               2 * max_local * (i - 1) + 2 * counts(i))
          global_center(3 * idx + 1:3 * (idx + counts(i))) = &
               gathered_center(3 * max_local * (i - 1) + 1: &
               3 * max_local * (i - 1) + 3 * counts(i))
          idx = idx + counts(i)
       end if
    end do

    ! cycle through the periodic facets and identify unique periodic directions
    ! and their extents by matching the facet pairs and comparing their centers
    do i = 1, msh%periodic%size
       ! get the center of the facet (averaged over the pair)
       call lpt_get_periodic_center(msh, i, src_center)

       ! find whether there is a match facet from the global list
       match_idx = 0
       do j = 1, n_global
          if (global_meta(2 * (j - 1) + 1) .eq. msh%periodic%p_facet_el(i)%x(1)&
               .and. global_meta(2 * (j - 1) + 2) .eq. &
               msh%periodic%p_facet_el(i)%x(2)) then
             match_idx = j
             exit
          end if
       end do
       if (match_idx .eq. 0) cycle

       ! fetch the target facet center and compute the shift
       tgt_center = global_center(3 * (match_idx - 1) + 1:3 * match_idx)
       shift = tgt_center - src_center
       if (norm2(shift) .le. LPT_PERIODIC_TOL) cycle

       ! compute the direction and projection
       dir = shift
       call lpt_normalize(dir)
       proj_src = dot_product(src_center, dir)
       proj_tgt = dot_product(tgt_center, dir)
       idx = 0

       ! identfy whether the direction is a unique one
       do j = 1, this%n_periodic_dirs
          if (abs(dot_product(dir, &
               lpt_periodic_bc_get_dir(this, j))) .gt. &
               1.0_rp - 1.0e-6_rp) then
             idx = j
             exit
          end if
       end do

       if (idx .ne. 0) cycle

       ! if unique, add it to the list of periodic directions
       ! with its extent and shift
       if (this%n_periodic_dirs .ge. 3) cycle
       idx = this%n_periodic_dirs + 1
       this%n_periodic_dirs = idx
       call lpt_periodic_bc_set_translational(this, idx, dir, &
            min(proj_src, proj_tgt), max(proj_src, proj_tgt), shift, &
            norm2(shift))
    end do

    deallocate(global_center)
    deallocate(global_meta)
    deallocate(gathered_center)
    deallocate(gathered_meta)
    deallocate(padded_center)
    deallocate(padded_meta)
    deallocate(counts)

    this%periodic_enabled = this%n_periodic_dirs .gt. 0
  end subroutine lpt_periodic_bc_init_translational

  !> Gather the corner points of a mesh facet.
  !! @param msh Mesh containing the element.
  !! @param el Local element index.
  !! @param facet Facet index on the element.
  !! @param pts Output corner coordinates.
  subroutine lpt_get_facet_points(msh, el, facet, pts)
    type(mesh_t), intent(in) :: msh
    integer, intent(in) :: el
    integer, intent(in) :: facet
    real(kind=rp), intent(out) :: pts(3, 4)
    integer, dimension(4, 6) :: face_nodes = reshape([ &
         1,5,7,3, &
         2,6,8,4, &
         1,2,6,5, &
         3,4,8,7, &
         1,2,4,3, &
         5,6,8,7], [4, 6])
    integer, dimension(2, 4) :: edge_nodes = reshape([ &
         1,3, &
         2,4, &
         1,2, &
         3,4], [2, 4])
    integer :: i

    pts = 0.0_rp
    if (msh%gdim .eq. 3) then
       do i = 1, 4
          pts(:, i) = real(msh%elements(el)%e%pts(face_nodes(i, facet))%p%x, rp)
       end do
    else
       do i = 1, 2
          pts(:, i) = real(msh%elements(el)%e%pts(edge_nodes(i, facet))%p%x, rp)
       end do
    end if
  end subroutine lpt_get_facet_points

  !> Compute the geometric center of a periodic source facet.
  !! @param msh Mesh containing periodic metadata.
  !! @param i_periodic Periodic-facet entry index.
  !! @param pt Output facet-center coordinate.
  subroutine lpt_get_periodic_center(msh, i_periodic, pt)
    type(mesh_t), intent(in) :: msh
    integer, intent(in) :: i_periodic
    real(kind=rp), intent(out) :: pt(3)
    integer :: npts
    real(kind=rp) :: pts(3, 4)

    npts = merge(4, 2, msh%gdim .eq. 3)
    call lpt_get_facet_points(msh, msh%periodic%facet_el(i_periodic)%x(2), &
         msh%periodic%facet_el(i_periodic)%x(1), pts)
    pt = sum(pts(:, 1:npts), dim = 2) / real(npts, rp)
  end subroutine lpt_get_periodic_center

  !> Normalise a 3-vector when its norm is above LPT periodic tolerance.
  !! @param v Vector to normalise in place.
  subroutine lpt_normalize(v)
    real(kind=rp), intent(inout) :: v(3)
    real(kind=rp) :: vnorm

    vnorm = norm2(v)
    if (vnorm .gt. LPT_PERIODIC_TOL) v = v / vnorm
  end subroutine lpt_normalize

  !> Clear translational periodic direction, extent, and shift storage.
  subroutine lpt_periodic_bc_clear_translational(this)
    class(lpt_periodic_bc_t), intent(inout) :: this

    this%periodic_shift_x1 = 0.0_rp
    this%periodic_shift_y1 = 0.0_rp
    this%periodic_shift_z1 = 0.0_rp
    this%periodic_shift_x2 = 0.0_rp
    this%periodic_shift_y2 = 0.0_rp
    this%periodic_shift_z2 = 0.0_rp
    this%periodic_shift_x3 = 0.0_rp
    this%periodic_shift_y3 = 0.0_rp
    this%periodic_shift_z3 = 0.0_rp
    this%periodic_dir_x1 = 0.0_rp
    this%periodic_dir_y1 = 0.0_rp
    this%periodic_dir_z1 = 0.0_rp
    this%periodic_dir_x2 = 0.0_rp
    this%periodic_dir_y2 = 0.0_rp
    this%periodic_dir_z2 = 0.0_rp
    this%periodic_dir_x3 = 0.0_rp
    this%periodic_dir_y3 = 0.0_rp
    this%periodic_dir_z3 = 0.0_rp
    this%periodic_min1 = 0.0_rp
    this%periodic_min2 = 0.0_rp
    this%periodic_min3 = 0.0_rp
    this%periodic_max1 = 0.0_rp
    this%periodic_max2 = 0.0_rp
    this%periodic_max3 = 0.0_rp
    this%periodic_len1 = 0.0_rp
    this%periodic_len2 = 0.0_rp
    this%periodic_len3 = 0.0_rp
  end subroutine lpt_periodic_bc_clear_translational

  pure function lpt_periodic_bc_get_dir(this, idx) result(dir)
    class(lpt_periodic_bc_t), intent(in) :: this
    integer, intent(in) :: idx
    real(kind=rp) :: dir(3)

    select case (idx)
    case (1)
       dir = [this%periodic_dir_x1, this%periodic_dir_y1, &
            this%periodic_dir_z1]
    case (2)
       dir = [this%periodic_dir_x2, this%periodic_dir_y2, &
            this%periodic_dir_z2]
    case (3)
       dir = [this%periodic_dir_x3, this%periodic_dir_y3, &
            this%periodic_dir_z3]
    case default
       dir = 0.0_rp
    end select
  end function lpt_periodic_bc_get_dir

  !> Store one translational periodic wrapping direction.
  !! @param idx Direction slot in the fixed three-direction storage.
  !! @param dir Unit direction of periodic wrapping.
  !! @param periodic_min Minimum coordinate projected on `dir`.
  !! @param periodic_max Maximum coordinate projected on `dir`.
  !! @param shift Translation vector between paired facets.
  !! @param periodic_len Periodic length along `dir`.
  subroutine lpt_periodic_bc_set_translational(this, idx, dir, periodic_min, &
       periodic_max, shift, periodic_len)
    class(lpt_periodic_bc_t), intent(inout) :: this
    integer, intent(in) :: idx
    real(kind=rp), intent(in) :: dir(3)
    real(kind=rp), intent(in) :: periodic_min
    real(kind=rp), intent(in) :: periodic_max
    real(kind=rp), intent(in) :: shift(3)
    real(kind=rp), intent(in) :: periodic_len

    select case (idx)
    case (1)
       this%periodic_dir_x1 = dir(1)
       this%periodic_dir_y1 = dir(2)
       this%periodic_dir_z1 = dir(3)
       this%periodic_min1 = periodic_min
       this%periodic_max1 = periodic_max
       this%periodic_shift_x1 = shift(1)
       this%periodic_shift_y1 = shift(2)
       this%periodic_shift_z1 = shift(3)
       this%periodic_len1 = periodic_len
    case (2)
       this%periodic_dir_x2 = dir(1)
       this%periodic_dir_y2 = dir(2)
       this%periodic_dir_z2 = dir(3)
       this%periodic_min2 = periodic_min
       this%periodic_max2 = periodic_max
       this%periodic_shift_x2 = shift(1)
       this%periodic_shift_y2 = shift(2)
       this%periodic_shift_z2 = shift(3)
       this%periodic_len2 = periodic_len
    case (3)
       this%periodic_dir_x3 = dir(1)
       this%periodic_dir_y3 = dir(2)
       this%periodic_dir_z3 = dir(3)
       this%periodic_min3 = periodic_min
       this%periodic_max3 = periodic_max
       this%periodic_shift_x3 = shift(1)
       this%periodic_shift_y3 = shift(2)
       this%periodic_shift_z3 = shift(3)
       this%periodic_len3 = periodic_len
    end select
  end subroutine lpt_periodic_bc_set_translational

end module lpt_periodic_bc
