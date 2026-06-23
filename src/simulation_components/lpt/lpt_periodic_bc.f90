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
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use mesh, only : mesh_t
  use comm, only : pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_Allgather, MPI_INTEGER, &
       MPI_MIN, MPI_MAX
  implicit none
  private

  real(kind=rp), parameter :: LPT_PERIODIC_TOL = 1.0e-8_rp

  type, public :: lpt_periodic_bc_t
     logical :: periodic_enabled = .false.
     logical :: rotational_periodic_enabled = .false.
     integer :: n_periodic_dirs = 0
     real(kind=rp) :: periodic_shift(3, 3) = 0.0_rp
     real(kind=rp) :: periodic_dir(3, 3) = 0.0_rp
     real(kind=rp) :: periodic_min(3) = 0.0_rp
     real(kind=rp) :: periodic_max(3) = 0.0_rp
     real(kind=rp) :: periodic_len(3) = 0.0_rp
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
     procedure, private, pass(this) :: wrap_rotational => &
          lpt_periodic_bc_wrap_rotational
     procedure, private, pass(this) :: wrap_translational => &
          lpt_periodic_bc_wrap_translational
  end type lpt_periodic_bc_t

contains

  subroutine lpt_periodic_bc_free(this)
    class(lpt_periodic_bc_t), intent(inout) :: this

    call this%reset()
  end subroutine lpt_periodic_bc_free

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

  subroutine lpt_periodic_bc_wrap(this, x, y, z, n, u, v, w, u_lag, v_lag, &
       w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, &
       acc_xlaglag, acc_ylaglag, acc_zlaglag)
    class(lpt_periodic_bc_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    integer, intent(in) :: n
    real(kind=rp), intent(inout), optional :: u(:), v(:), w(:)
    real(kind=rp), intent(inout), optional :: u_lag(:), v_lag(:), w_lag(:)
    real(kind=rp), intent(inout), optional :: u_laglag(:), v_laglag(:)
    real(kind=rp), intent(inout), optional :: w_laglag(:)
    real(kind=rp), intent(inout), optional :: acc_xlag(:), acc_ylag(:)
    real(kind=rp), intent(inout), optional :: acc_zlag(:)
    real(kind=rp), intent(inout), optional :: acc_xlaglag(:)
    real(kind=rp), intent(inout), optional :: acc_ylaglag(:)
    real(kind=rp), intent(inout), optional :: acc_zlaglag(:)

    if (n .eq. 0) return

    if (this%rotational_periodic_enabled) then
       call this%wrap_rotational(x, y, z, n, u, v, w, u_lag, v_lag, &
            w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, acc_ylag, &
            acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag)
       return
    end if

    call this%wrap_translational(x, y, z, n)
  end subroutine lpt_periodic_bc_wrap

  subroutine lpt_periodic_bc_reset(this)
    class(lpt_periodic_bc_t), intent(inout) :: this

    this%periodic_enabled = .false.
    this%rotational_periodic_enabled = .false.
    this%n_periodic_dirs = 0
    this%periodic_shift = 0.0_rp
    this%periodic_dir = 0.0_rp
    this%periodic_min = 0.0_rp
    this%periodic_max = 0.0_rp
    this%periodic_len = 0.0_rp
    this%rotational_theta_min = 0.0_rp
    this%rotational_theta_max = 0.0_rp
    this%rotational_theta_len = 0.0_rp
  end subroutine lpt_periodic_bc_reset

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
          if (global_meta(2 * (j - 1) + 1) .eq. msh%periodic%p_facet_el(i)%x(1) &
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
          if (abs(dot_product(dir, this%periodic_dir(:, j))) .gt. &
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
       this%periodic_dir(:, idx) = dir
       this%periodic_min(idx) = min(proj_src, proj_tgt)
       this%periodic_max(idx) = max(proj_src, proj_tgt)
       this%periodic_shift(:, idx) = shift
       this%periodic_len(idx) = norm2(shift)
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

  subroutine lpt_periodic_bc_wrap_rotational(this, x, y, z, n, u, v, w, &
       u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, acc_xlag, &
       acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, acc_zlaglag)
    class(lpt_periodic_bc_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    integer, intent(in) :: n
    real(kind=rp), intent(inout), optional :: u(:), v(:), w(:)
    real(kind=rp), intent(inout), optional :: u_lag(:), v_lag(:), w_lag(:)
    real(kind=rp), intent(inout), optional :: u_laglag(:), v_laglag(:)
    real(kind=rp), intent(inout), optional :: w_laglag(:)
    real(kind=rp), intent(inout), optional :: acc_xlag(:), acc_ylag(:)
    real(kind=rp), intent(inout), optional :: acc_zlag(:)
    real(kind=rp), intent(inout), optional :: acc_xlaglag(:)
    real(kind=rp), intent(inout), optional :: acc_ylaglag(:)
    real(kind=rp), intent(inout), optional :: acc_zlaglag(:)
    integer :: i
    real(kind=rp) :: radius
    real(kind=rp) :: theta_old
    real(kind=rp) :: theta
    real(kind=rp) :: dtheta
    real(kind=rp) :: pi

    pi = acos(-1.0_rp)
    do i = 1, n
       radius = sqrt(x(i) * x(i) + y(i) * y(i))
       theta_old = modulo(atan2(y(i), x(i)) + 2.0_rp * pi, &
            2.0_rp * pi)
       theta = theta_old

       do while (theta .lt. this%rotational_theta_min - LPT_PERIODIC_TOL)
          theta = theta + this%rotational_theta_len
       end do

       do while (theta .gt. this%rotational_theta_max + LPT_PERIODIC_TOL)
          theta = theta - this%rotational_theta_len
       end do

       dtheta = theta - theta_old
       x(i) = radius * cos(theta)
       y(i) = radius * sin(theta)
       if (abs(dtheta) .le. LPT_PERIODIC_TOL) cycle

       if (present(u) .and. present(v)) call lpt_rotate_xy(u(i), v(i), dtheta)
       if (present(u_lag) .and. present(v_lag)) &
            call lpt_rotate_xy(u_lag(i), v_lag(i), dtheta)
       if (present(u_laglag) .and. present(v_laglag)) &
            call lpt_rotate_xy(u_laglag(i), v_laglag(i), dtheta)
       if (present(acc_xlag) .and. present(acc_ylag)) &
            call lpt_rotate_xy(acc_xlag(i), acc_ylag(i), dtheta)
       if (present(acc_xlaglag) .and. present(acc_ylaglag)) &
            call lpt_rotate_xy(acc_xlaglag(i), acc_ylaglag(i), dtheta)
    end do
  end subroutine lpt_periodic_bc_wrap_rotational

  subroutine lpt_periodic_bc_wrap_translational(this, x, y, z, n)
    class(lpt_periodic_bc_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: x(:), y(:), z(:)
    integer, intent(in) :: n
    integer :: i
    integer :: j
    real(kind=rp) :: point(3)
    real(kind=rp) :: coord

    if (.not. this%periodic_enabled) return

    do i = 1, n
       point = [x(i), y(i), z(i)]
       do j = 1, this%n_periodic_dirs
          coord = dot_product(point, this%periodic_dir(:, j))

          do while (coord .lt. this%periodic_min(j) - LPT_PERIODIC_TOL)
             point = point + this%periodic_shift(:, j)
             coord = coord + this%periodic_len(j)
          end do

          do while (coord .gt. this%periodic_max(j) + LPT_PERIODIC_TOL)
             point = point - this%periodic_shift(:, j)
             coord = coord - this%periodic_len(j)
          end do
       end do
       x(i) = point(1)
       y(i) = point(2)
       z(i) = point(3)
    end do
  end subroutine lpt_periodic_bc_wrap_translational

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

  subroutine lpt_normalize(v)
    real(kind=rp), intent(inout) :: v(3)
    real(kind=rp) :: vnorm

    vnorm = norm2(v)
    if (vnorm .gt. LPT_PERIODIC_TOL) v = v / vnorm
  end subroutine lpt_normalize

  subroutine lpt_rotate_vector_xy(vec, theta)
    real(kind=rp), intent(inout) :: vec(3)
    real(kind=rp), intent(in) :: theta
    real(kind=rp) :: x_old
    real(kind=rp) :: y_old
    real(kind=rp) :: cos_theta
    real(kind=rp) :: sin_theta

    x_old = vec(1)
    y_old = vec(2)
    cos_theta = cos(theta)
    sin_theta = sin(theta)

    vec(1) = cos_theta * x_old - sin_theta * y_old
    vec(2) = sin_theta * x_old + cos_theta * y_old
  end subroutine lpt_rotate_vector_xy

  subroutine lpt_rotate_xy(x, y, theta)
    real(kind=rp), intent(inout) :: x
    real(kind=rp), intent(inout) :: y
    real(kind=rp), intent(in) :: theta
    real(kind=rp) :: x_old
    real(kind=rp) :: y_old
    real(kind=rp) :: cos_theta
    real(kind=rp) :: sin_theta

    x_old = x
    y_old = y
    cos_theta = cos(theta)
    sin_theta = sin(theta)

    x = cos_theta * x_old - sin_theta * y_old
    y = sin_theta * x_old + cos_theta * y_old
  end subroutine lpt_rotate_xy

end module lpt_periodic_bc
