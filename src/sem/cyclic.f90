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
!> Coefficients
!> Cyclic (rotationally periodic) boundary conditions
module cyclic
  use num_types, only : rp
  use neko_config, only : NEKO_BCKND_DEVICE
  use math, only : NEKO_EPS
  use logger, only : neko_log, LOG_SIZE
  use mesh, only : mesh_t
  use space, only : space_t
  use dofmap, only : dofmap_t
  use device, only : device_map, device_memcpy, device_unmap, HOST_TO_DEVICE
  use utils, only : index_is_on_facet, linear_index, neko_error
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_size
  use mpi_f08, only : MPI_Allreduce, MPI_Allgather, MPI_INTEGER, MPI_SUM, &
       MPI_IN_PLACE
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  !> Cyclic boundaries of a sector domain.
  !! @details The two periodic boundaries must be planes through a common
  !! rotation axis, and periodic points are related by a rotation about it.
  !! Vectors at the cyclic points are rotated into a local frame shared by the
  !! matching points before gather-scatter, see rotate_cyc in operators.
  !!
  !! The rotation axis is the intersection line of the two cyclic planes,
  !! \f$ a = n_1 \times n_2 / |n_1 \times n_2| \f$. At each cyclic point the
  !! frame is \f$ (e_1, a \times e_1, a) \f$, with \f$ e_1 = n_1 \f$ on the
  !! first boundary and \f$ e_1 = -n_2 \f$ on the second. Since the outward
  !! normals of a sector point in opposite angular directions, the rotation
  !! mapping one boundary onto the other maps \f$ e_1 \f$ onto itself, so
  !! matching periodic points share the same components in this frame.
  type, public :: cyclic_t
     !> Whether the periodic boundaries are cyclic. A setting, kept by free()
     logical :: enabled = .false.
     !> Linear indices of the cyclic points, msk(0) holds their number + 1
     integer, allocatable :: msk(:)
     !> Rotation matrix from Cartesian into the cyclic frame at each cyclic
     !! point, stored row-wise: R(1:3,i) is the boundary normal, R(4:6,i)
     !! the tangent axis x normal and R(7:9,i) the axis
     real(kind=rp), allocatable :: R(:,:)
     !> Unit outward normals of the two planar cyclic boundaries, one per
     !! column, identical on all ranks
     real(kind=rp) :: normal(3,2) = 0.0_rp
     !> Unit rotation axis, parallel to normal(:,1) x normal(:,2) and
     !! oriented so that its largest component is positive
     real(kind=rp) :: axis(3) = 0.0_rp
     !> Point on the rotation axis, the one closest to the origin
     real(kind=rp) :: origin(3) = 0.0_rp
     !> Sector angle swept by the domain about axis, in radians
     real(kind=rp) :: angle = 0.0_rp
     type(c_ptr) :: msk_d = C_NULL_PTR
     type(c_ptr) :: R_d = C_NULL_PTR
   contains
     procedure, pass(this) :: init => cyclic_init
     procedure, pass(this) :: free => cyclic_free
     procedure, pass(this) :: generate => cyclic_generate
  end type cyclic_t

contains

  !> Allocate the cyclic mask and frames, sized for all periodic facets
  !! and initialised to no rotation.
  !! @param msh Mesh holding the periodic facets.
  !! @param Xh Function space.
  subroutine cyclic_init(this, msh, Xh)
    class(cyclic_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    integer :: ncyc

    call this%free()

    ncyc = msh%periodic%size * Xh%lx * Xh%lx
    allocate(this%msk(0:ncyc))
    this%msk(0) = ncyc + 1
    if (ncyc .gt. 0) then
       allocate(this%R(9, ncyc))

       ! Default values correspond to no rotation
       this%R = 0.0_rp
       this%R(1,:) = 1.0_rp
       this%R(5,:) = 1.0_rp
       this%R(9,:) = 1.0_rp

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(this%msk, this%msk_d, ncyc+1)
          call device_map(this%R, this%R_d, 9 * ncyc)

          call device_memcpy(this%msk, this%msk_d, ncyc+1, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(this%R, this%R_d, 9 * ncyc, &
               HOST_TO_DEVICE, sync = .false.)
       end if
    end if

  end subroutine cyclic_init

  !> Deallocate the cyclic mask and frames. Keeps @a enabled.
  subroutine cyclic_free(this)
    class(cyclic_t), intent(inout) :: this

    if (allocated(this%msk)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(this%msk, this%msk_d)
       deallocate(this%msk)
    end if

    if (allocated(this%R)) then
       if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(this%R, this%R_d)
       deallocate(this%R)
    end if

    this%normal = 0.0_rp
    this%axis = 0.0_rp
    this%origin = 0.0_rp
    this%angle = 0.0_rp

  end subroutine cyclic_free

  !> Find the cyclic geometry and build the frames used by rotate_cyc.
  !! Does nothing unless @a enabled.
  !! @param msh Mesh holding the periodic facets.
  !! @param Xh Function space.
  !! @param dof Dofmap with the GLL coordinates.
  !! @param nx, ny, nz Facet normals, as in coef_t.
  !! @note Collective over NEKO_COMM.
  subroutine cyclic_generate(this, msh, Xh, dof, nx, ny, nz)
    class(cyclic_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(dofmap_t), intent(in) :: dof
    real(kind=rp), intent(in) :: nx(:,:,:,:), ny(:,:,:,:), nz(:,:,:,:)
    real(kind=rp) :: un(3), len, e1(3), e2(3), a(3)
    integer :: lx, ly, lz, np, np_glb, ierr
    integer :: i, j, k, pf, pe, n, nc, ncyc
    character(len=LOG_SIZE) :: log_buf

    if (.not. this%enabled) return

    np = msh%periodic%size
    call MPI_Allreduce(np, np_glb, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

    if (np_glb .eq. 0) then
       call neko_error("There are no periodic boundaries. " // &
            "Switch cyclic off in the case file.")
    end if

    ! Collective, every rank must call it before an early return
    call cyclic_find_normals(this, msh, Xh, nx, ny, nz)

    a = cross(this%normal(:,1), this%normal(:,2))
    len = norm2(a)
    if (len .lt. 1.0e-3_rp) then
       call neko_error("The two cyclic boundaries are parallel, the " // &
            "rotation axis is undefined. Use plain periodic boundaries " // &
            "for translations. Sectors of 180 degrees are not supported.")
    end if
    a = a / len
    ! Fix the orientation independently of the order the normals were found
    ! in, which depends on the partitioning
    if (a(maxloc(abs(a), 1)) .lt. 0.0_rp) a = -a
    this%axis = a

    call cyclic_find_origin_and_angle(this, msh, Xh, dof, nx, ny, nz)

    write(log_buf, '(A,3ES13.5)') 'Cyclic axis     : ', this%axis
    call neko_log%message(log_buf)
    write(log_buf, '(A,3ES13.5)') 'Cyclic origin   : ', this%origin
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.5)') 'Cyclic angle    : ', this%angle
    call neko_log%message(log_buf)

    if (np .eq. 0) return

    lx = Xh%lx
    ly = Xh%ly
    lz = Xh%lz
    ncyc = this%msk(0) - 1
    nc = 1
    do n = 1, np
       pf = msh%periodic%facet_el(n)%x(1)
       pe = msh%periodic%facet_el(n)%x(2)
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                if (index_is_on_facet(i, j, k, lx, ly, lz, pf)) then
                   un = facet_normal(nx, ny, nz, i, j, k, pe, pf)

                   if (cyclic_side(this, un) .eq. 1) then
                      e1 = this%normal(:,1)
                   else
                      e1 = -this%normal(:,2)
                   end if
                   e2 = cross(a, e1)

                   this%msk(nc) = linear_index(i, j, k, pe, lx, ly, lz)
                   this%R(1:3, nc) = e1
                   this%R(4:6, nc) = e2
                   this%R(7:9, nc) = a
                   nc = nc + 1
                end if
             end do
          end do
       end do
    end do

    if (nc - 1 /= ncyc) then
       call neko_error("The number of cyclic GLL points were " // &
            "not estimated correctly.")
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%msk, this%msk_d, ncyc+1, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%R, this%R_d, 9 * ncyc, &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine cyclic_generate

  !> Find the normals of the two planar cyclic boundaries.
  !! Every GLL point on a periodic facet is checked against the normals found
  !! so far, which also verifies that each cyclic boundary is a plane. The
  !! local normals are then gathered and merged in rank order, so that
  !! @a normal is identical on all ranks.
  !! @note Collective over NEKO_COMM.
  subroutine cyclic_find_normals(this, msh, Xh, nx, ny, nz)
    class(cyclic_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    real(kind=rp), intent(in) :: nx(:,:,:,:), ny(:,:,:,:), nz(:,:,:,:)
    !> Max distance between two unit normals regarded as the same direction
    real(kind=rp), parameter :: tol = 1.0e-4_rp
    real(kind=rp) :: un(3), len, nrm(3,2)
    real(kind=rp) :: buf(7)
    real(kind=rp), allocatable :: buf_glb(:,:)
    integer :: i, j, k, n, r, m, pf, pe, lx, ly, lz, cnt, ierr
    character(len=LOG_SIZE) :: log_buf

    lx = Xh%lx
    ly = Xh%ly
    lz = Xh%lz

    ! Distinct normals on this rank
    cnt = 0
    nrm = 0.0_rp
    do n = 1, msh%periodic%size
       pf = msh%periodic%facet_el(n)%x(1)
       pe = msh%periodic%facet_el(n)%x(2)
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                if (.not. index_is_on_facet(i, j, k, lx, ly, lz, pf)) cycle
                un = facet_normal(nx, ny, nz, i, j, k, pe, pf)
                len = norm2(un)
                if (len .le. NEKO_EPS) then
                   call neko_error("Zero surface normal on a cyclic " // &
                        "boundary.")
                end if
                call add_normal(nrm, cnt, un / len)
             end do
          end do
       end do
    end do

    ! Merge the normals of all ranks, in rank order
    buf(1:3) = nrm(:,1)
    buf(4:6) = nrm(:,2)
    buf(7) = real(cnt, rp)
    allocate(buf_glb(7, pe_size))
    call MPI_Allgather(buf, 7, MPI_REAL_PRECISION, buf_glb, 7, &
         MPI_REAL_PRECISION, NEKO_COMM, ierr)

    cnt = 0
    nrm = 0.0_rp
    do r = 1, pe_size
       do m = 1, nint(buf_glb(7, r))
          call add_normal(nrm, cnt, buf_glb(3*m-2:3*m, r))
       end do
    end do
    deallocate(buf_glb)

    if (cnt .ne. 2) then
       call neko_error("Cyclic boundaries must consist of exactly two " // &
            "planes.")
    end if

    this%normal = nrm

    do m = 1, 2
       write(log_buf, '(A,I1,A,3ES13.5)') 'Cyclic normal ', m, ' : ', &
            this%normal(:,m)
       call neko_log%message(log_buf)
    end do

  contains

    !> Add @a v to the list of distinct normals @a nrm, unless already there
    subroutine add_normal(nrm, cnt, v)
      real(kind=rp), intent(inout) :: nrm(3,2)
      integer, intent(inout) :: cnt
      real(kind=rp), intent(in) :: v(3)
      integer :: l

      do l = 1, cnt
         if (norm2(v - nrm(:,l)) .lt. tol) return
      end do

      if (cnt .eq. 2) then
         call neko_error("More than two distinct surface normals on the " // &
              "cyclic boundaries. Each cyclic boundary must be planar.")
      end if

      cnt = cnt + 1
      nrm(:,cnt) = v
    end subroutine add_normal

  end subroutine cyclic_find_normals

  !> Find a point on the rotation axis and the sector angle.
  !! The axis is the intersection of the planes \f$ n_k \cdot x = c_k \f$,
  !! where \f$ c_k \f$ is averaged over the points of boundary k. The point
  !! closest to the origin is picked by also requiring \f$ a \cdot x = 0 \f$.
  !! The angle is measured from the radial direction of boundary 1 to that of
  !! boundary 2, turning into the domain, i.e. along \f$ -n_1 \f$.
  !! @note Collective over NEKO_COMM. Needs @a normal and @a axis.
  subroutine cyclic_find_origin_and_angle(this, msh, Xh, dof, nx, ny, nz)
    class(cyclic_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(dofmap_t), intent(in) :: dof
    real(kind=rp), intent(in) :: nx(:,:,:,:), ny(:,:,:,:), nz(:,:,:,:)
    real(kind=rp) :: n1(3), n2(3), a(3), x(3), d(3), er(3,2), b(3)
    real(kind=rp) :: buf(8), c(2), det, ang, pi
    integer :: i, j, k, n, pf, pe, lx, ly, lz, m, ierr

    lx = Xh%lx
    ly = Xh%ly
    lz = Xh%lz
    n1 = this%normal(:,1)
    n2 = this%normal(:,2)
    a = this%axis
    pi = acos(-1.0_rp)

    ! Plane offsets c_k, as sums and counts in buf(1:4)
    buf = 0.0_rp
    do n = 1, msh%periodic%size
       pf = msh%periodic%facet_el(n)%x(1)
       pe = msh%periodic%facet_el(n)%x(2)
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                if (.not. index_is_on_facet(i, j, k, lx, ly, lz, pf)) cycle
                m = cyclic_side(this, facet_normal(nx, ny, nz, i, j, k, pe, pf))
                x = [dof%x%x(i, j, k, pe), dof%y%x(i, j, k, pe), &
                     dof%z%x(i, j, k, pe)]
                buf(m) = buf(m) + dot_product(this%normal(:,m), x)
                buf(m + 2) = buf(m + 2) + 1.0_rp
             end do
          end do
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, buf, 4, MPI_REAL_PRECISION, MPI_SUM, &
         NEKO_COMM, ierr)
    if (buf(3) .eq. 0.0_rp .or. buf(4) .eq. 0.0_rp) then
       call neko_error("A cyclic boundary has no points.")
    end if
    c = buf(1:2) / buf(3:4)

    ! Intersection of n1.x = c1, n2.x = c2 and a.x = 0
    det = dot_product(n1, cross(n2, a))
    this%origin = (c(1) * cross(n2, a) + c(2) * cross(a, n1)) / det

    ! Mean radial directions of the two boundaries, in buf(1:6)
    buf = 0.0_rp
    do n = 1, msh%periodic%size
       pf = msh%periodic%facet_el(n)%x(1)
       pe = msh%periodic%facet_el(n)%x(2)
       do k = 1, lz
          do j = 1, ly
             do i = 1, lx
                if (.not. index_is_on_facet(i, j, k, lx, ly, lz, pf)) cycle
                m = cyclic_side(this, facet_normal(nx, ny, nz, i, j, k, pe, pf))
                d = [dof%x%x(i, j, k, pe), dof%y%x(i, j, k, pe), &
                     dof%z%x(i, j, k, pe)] - this%origin
                d = d - dot_product(d, a) * a
                buf(3*m-2:3*m) = buf(3*m-2:3*m) + d
             end do
          end do
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, buf, 6, MPI_REAL_PRECISION, MPI_SUM, &
         NEKO_COMM, ierr)
    do m = 1, 2
       er(:,m) = buf(3*m-2:3*m)
       if (norm2(er(:,m)) .le. NEKO_EPS) then
          call neko_error("A cyclic boundary is centred on the rotation " // &
               "axis, cannot determine the sector angle.")
       end if
       er(:,m) = er(:,m) / norm2(er(:,m))
    end do

    ! Angle from er1 to er2 counter-clockwise about a, in [0, 2 pi)
    ang = atan2(dot_product(cross(er(:,1), er(:,2)), a), &
         dot_product(er(:,1), er(:,2)))
    if (ang .lt. 0.0_rp) ang = ang + 2.0_rp * pi

    ! The domain lies along -n1 from boundary 1. If that is clockwise about a,
    ! the sector is the complement.
    b = cross(a, er(:,1))
    if (dot_product(b, -n1) .lt. 0.0_rp) ang = 2.0_rp * pi - ang
    this%angle = ang

  end subroutine cyclic_find_origin_and_angle

  !> Which cyclic boundary, 1 or 2, a facet normal belongs to
  pure function cyclic_side(this, un) result(side)
    class(cyclic_t), intent(in) :: this
    real(kind=rp), intent(in) :: un(3)
    integer :: side

    if (dot_product(un, this%normal(:,1)) .ge. &
         dot_product(un, this%normal(:,2))) then
       side = 1
    else
       side = 2
    end if
  end function cyclic_side

  !> Facet normal at GLL point (i, j, k) of facet @a facet in element @a e,
  !! from normal arrays laid out as in coef_t, see coef_get_normal()
  pure function facet_normal(nx, ny, nz, i, j, k, e, facet) result(normal)
    real(kind=rp), intent(in) :: nx(:,:,:,:), ny(:,:,:,:), nz(:,:,:,:)
    integer, intent(in) :: i, j, k, e, facet
    real(kind=rp) :: normal(3)

    select case (facet)
    case (1, 2)
       normal = [nx(j, k, facet, e), ny(j, k, facet, e), nz(j, k, facet, e)]
    case (3, 4)
       normal = [nx(i, k, facet, e), ny(i, k, facet, e), nz(i, k, facet, e)]
    case default
       normal = [nx(i, j, facet, e), ny(i, j, facet, e), nz(i, j, facet, e)]
    end select
  end function facet_normal

  !> Cross product of two 3-vectors
  pure function cross(u, v) result(w)
    real(kind=rp), intent(in) :: u(3), v(3)
    real(kind=rp) :: w(3)

    w(1) = u(2) * v(3) - u(3) * v(2)
    w(2) = u(3) * v(1) - u(1) * v(3)
    w(3) = u(1) * v(2) - u(2) * v(1)
  end function cross

end module cyclic
