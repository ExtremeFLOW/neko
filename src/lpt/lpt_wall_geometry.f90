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
!> High-order wall geometry and trajectory-based collision detection for
!! Lagrangian particle tracking.
!!
!! Each wall facet is represented by its tensor-product GLL nodes,
!! \f$ X_w(s,t) = \sum_{a,b} X_{ab} \ell_a(s) \ell_b(t) \f$ on
!! \f$ [-1,1]^2 \f$, i.e. exactly the spectral-element geometry of the
!! fluid solver, with the outward normal
!! \f$ n_w = \pm X_s \times X_t / |X_s \times X_t| \f$.
!!
!! A collision is the earliest outgoing intersection of the straight
!! particle-centre path \f$ x(\alpha) = x^n + \alpha \Delta x \f$,
!! \f$ \alpha \in [0,1] \f$, with the offset wall
!! \f$ X_a = X_w - a n_w \f$, found by Newton's method on
!! \f$ F(\alpha,s,t) = x^n + \alpha \Delta x - X_w + a n_w = 0 \f$ with
!! \f$ J = [\Delta x, -X_s + a n_s, -X_t + a n_t] \f$. Uniqueness of the
!! root in a (sub)patch is certified with a normal cone: if all normals of
!! the patch make an angle below 90 degrees with \f$ \Delta x \f$, the
!! offset patch is a graph along \f$ \Delta x \f$ and the path crosses it
!! at most once. Patches failing the test, or where Newton fails, are
!! subdivided. Finite-radius particles may also touch facet edges
!! (path against a tube of radius a around the edge curve) and corners.
!!
!! Assumptions:
!! - The particle centre moves on the straight segment between its
!!   positions before and after the step; the curvature of the path within
!!   the step is not resolved.
!! - The particle radius is smaller than the wall's radius of curvature, so
!!   that the offset surface is regular.
!! - The displacement in one step plus the particle radius does not exceed
!!   the bounding radius of the element that owned the particle at the
!!   start of the step (a CFL-type condition). Violations are counted.
!! - Bounding spheres and normal cones of patches are computed from samples
!!   on a grid finer than the GLL grid, enlarged by the largest change
!!   possible between samples estimated from the sampled derivatives.
module lpt_wall_geometry
  use num_types, only : rp, i8
  use math, only : NEKO_EPS
  use fast3d, only : fd_weights_full
  use mesh, only : mesh_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use utils, only : neko_error
  use comm, only : NEKO_COMM, pe_rank, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allgather, MPI_Alltoall, MPI_Alltoallv, &
       MPI_INTEGER
  use device, only : device_map, device_memcpy, device_unmap, &
       device_alloc, device_free, HOST_TO_DEVICE, DEVICE_TO_HOST
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated, &
       c_size_t
  implicit none
  private

  !> Maximum number of wall collisions resolved for one particle per step.
  integer, public, parameter :: LPT_WALL_MAX_BOUNCES = 8
  !> Normal coefficient of restitution of the elastic wall model.
  real(kind=rp), public, parameter :: LPT_WALL_RESTITUTION = 1.0_rp
  !> Maximum number of Newton iterations of an intersection solve.
  integer, parameter :: NEWTON_MAX_ITER = 25
  !> Maximum number of step halvings in a Newton iteration.
  integer, parameter :: NEWTON_MAX_HALVINGS = 10
  !> Maximum depth of the patch subdivision.
  integer, parameter :: PATCH_MAX_DEPTH = 4
  !> Size of the patch stack (three siblings pending per level).
  integer, parameter :: PATCH_STACK = 3 * PATCH_MAX_DEPTH + 2
  !> Number of reals in the packed per-facet device record, see to_device
  !! and lpt_wall_intersect.h.
  integer, public, parameter :: LPT_WALL_NFR = 64
  !> Number of diagnostic counters kept on the device, in the order of
  !! lpt_wall_stats_t (and LPT_ST_* in lpt_wall_intersect.h).
  integer, public, parameter :: LPT_WALL_NSTATS = 2
  !> Largest lx supported by the device kernels.
  integer, public, parameter :: LPT_WALL_MAX_LX_DEVICE = 16
  !> Newton solve outcomes.
  integer, parameter :: NEWTON_CONVERGED = 0, NEWTON_SINGULAR = 1, &
       NEWTON_NOT_CONVERGED = 2, NEWTON_OUTSIDE = 3

  !> Counters of degraded collision handling, accumulated over the run.
  type, public :: lpt_wall_stats_t
     !> Patches left unresolved at the finest level.
     integer(kind=i8) :: n_unresolved = 0_i8
     !> Particles exceeding the supported per-step displacement.
     integer(kind=i8) :: n_reach_exceeded = 0_i8
  end type lpt_wall_stats_t

  !> High-order wall facets near the local mesh partition, including copies
  !! of wall facets owned by neighbouring ranks and periodic images, and
  !! per-element candidate lists.
  type, public :: lpt_wall_geometry_t
     !> Number of GLL points per direction.
     integer :: lx = 0
     !> Number of facets stored (local, halo and images).
     integer :: n_facets = 0
     !> Number of facets owned by this rank.
     integer :: n_local_facets = 0
     !> GLL points on [-1, 1].
     real(kind=rp), allocatable :: zg(:)
     !> Facet nodes relative to the facet origin, (3, lx, lx, n_facets).
     !! All intersection computations use these facet-local coordinates,
     !! so that round-off scales with the facet size rather than with the
     !! distance from the coordinate origin.
     real(kind=rp), allocatable :: xf(:,:,:,:)
     !> Facet origin (node centroid), (3, n_facets).
     real(kind=rp), allocatable :: origin(:,:)
     !> Orientation of \f$ X_s \times X_t \f$ relative to the outward normal.
     real(kind=rp), allocatable :: orient(:)
     !> Local element of a facet owned by this rank, 0 otherwise.
     integer, allocatable :: fel(:)
     !> Element facet number (1..6) of a facet owned by this rank, or 0.
     integer, allocatable :: ffacet(:)
     !> Centre of the bounding sphere of each facet.
     real(kind=rp), allocatable :: center(:,:)
     !> Radius of the bounding sphere of each facet.
     real(kind=rp), allocatable :: radius(:)
     !> Axis of the normal cone of each facet.
     real(kind=rp), allocatable :: cone_axis(:,:)
     !> Half-angle of the normal cone of each facet.
     real(kind=rp), allocatable :: cone_angle(:)
     !> Range of (X - center) . cone_axis over each facet (a slab enclosing
     !! the surface), (2, n_facets).
     real(kind=rp), allocatable :: slab(:,:)
     !> Length of the chord of each edge, (4, n_facets).
     real(kind=rp), allocatable :: edge_len(:,:)
     !> Corner points and outward normals, (3, 4, n_facets).
     real(kind=rp), allocatable :: corner_x(:,:,:), corner_n(:,:,:)
     !> Bounding sphere of each edge and its largest deviation from the
     !! polyline through its GLL nodes.
     real(kind=rp), allocatable :: edge_c(:,:,:), edge_r(:,:)
     real(kind=rp), allocatable :: edge_bulge(:,:)
     !> Coordinate scale used for round-off aware tolerances.
     real(kind=rp), allocatable :: scale(:)
     !> Whether all facet nodes lie in one plane (exactly planar surface).
     logical, allocatable :: planar(:)
     !> Bounding-sphere centre of each local element.
     real(kind=rp), allocatable :: el_center(:,:)
     !> Bounding-sphere radius of each local element, also the maximum
     !! supported per-step displacement plus particle radius.
     real(kind=rp), allocatable :: el_radius(:)
     !> CSR offsets of the candidate wall facets of each local element.
     integer, allocatable :: el_ptr(:)
     !> CSR list of candidate wall facets.
     integer, allocatable :: el_facets(:)
     !> Counters.
     type(lpt_wall_stats_t) :: stats
     !> Packed per-facet data for the device kernels, (LPT_WALL_NFR, nf)
     !! reals and (2, nf) integers (owning element, element facet).
     real(kind=rp), allocatable :: dev_fr(:,:)
     integer, allocatable :: dev_fi(:,:)
     !> Per-launch device counters, same order as lpt_wall_stats_t.
     integer, allocatable :: dev_stats(:)
     !> Per-particle scratch of the OpenCL and Metal kernels: hit counts
     !! and contact normals (3 x LPT_WALL_MAX_BOUNCES), for up to
     !! scratch_cap particles.
     type(c_ptr) :: hits_d = C_NULL_PTR
     type(c_ptr) :: normals_d = C_NULL_PTR
     integer :: scratch_cap = 0
     !> Device copies.
     type(c_ptr) :: xf_d = C_NULL_PTR
     type(c_ptr) :: fr_d = C_NULL_PTR
     type(c_ptr) :: fi_d = C_NULL_PTR
     type(c_ptr) :: zg_d = C_NULL_PTR
     type(c_ptr) :: el_ptr_d = C_NULL_PTR
     type(c_ptr) :: el_facets_d = C_NULL_PTR
     type(c_ptr) :: el_radius_d = C_NULL_PTR
     type(c_ptr) :: stats_d = C_NULL_PTR
   contains
     !> Initialise from the mesh and a mask of wall facets.
     procedure, pass(this) :: init => lpt_wall_geometry_init
     !> Release all data.
     procedure, pass(this) :: free => lpt_wall_geometry_free
     !> Earliest outgoing intersection of a path with a facet.
     procedure, pass(this) :: intersect => lpt_wall_facet_intersect
     !> Resolve all wall collisions of one particle in one step.
     procedure, pass(this) :: collide => lpt_wall_collide
     !> Pack the geometry for the device kernels and copy it to the device.
     procedure, pass(this) :: to_device => lpt_wall_to_device
     !> Pack the geometry for the device kernels on the host only.
     procedure, pass(this) :: pack => lpt_wall_pack
     !> Size the per-particle device scratch of the two-kernel backends.
     procedure, pass(this) :: reserve_scratch => lpt_wall_reserve_scratch
     !> Add the per-launch device counters to the statistics.
     procedure, pass(this) :: gather_stats => lpt_wall_gather_stats
  end type lpt_wall_geometry_t

contains

  !> Initialise the wall geometry from the mesh.
  !! @param msh Mesh of the fluid.
  !! @param dm_Xh Dofmap holding the high-order node coordinates.
  !! @param coef Coefficients providing the outward facet normals.
  !! @param wall_facet_mask Mask of facets that are walls, (2*gdim, nelv).
  !! @param images Periodic image transforms, (4, m): a rotation about the
  !! z axis by images(4, k) followed by the translation images(1:3, k).
  !! Wall facets are also stored at their images, so that collisions of
  !! particles crossing a periodic boundary in the same step are found
  !! before the particle is wrapped.
  subroutine lpt_wall_geometry_init(this, msh, dm_Xh, coef, wall_facet_mask, &
       images)
    class(lpt_wall_geometry_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: wall_facet_mask(:,:)
    real(kind=rp), intent(in), optional :: images(:,:)
    real(kind=rp), allocatable :: xloc(:,:,:,:), oloc(:), transforms(:,:)
    integer, allocatable :: eloc(:), floc(:)
    integer :: e, facet, nloc, lx, m

    call this%free()
    if (msh%gdim .ne. 3) then
       call neko_error("lpt curved wall collisions require a 3D mesh")
    end if
    call coef%require_facets('lpt wall collisions')

    lx = dm_Xh%Xh%lx
    nloc = count(wall_facet_mask)
    allocate(xloc(3, lx, lx, max(nloc, 1)), oloc(max(nloc, 1)))
    allocate(eloc(max(nloc, 1)), floc(max(nloc, 1)))

    call init_basis(this, dm_Xh%Xh%zg(:,1), lx)

    nloc = 0
    do e = 1, msh%nelv
       do facet = 1, 6
          if (.not. wall_facet_mask(facet, e)) cycle
          nloc = nloc + 1
          call extract_facet(dm_Xh, e, facet, lx, xloc(:,:,:,nloc))
          oloc(nloc) = facet_orientation(this, xloc(:,:,:,nloc), coef, &
               e, facet)
          eloc(nloc) = e
          floc(nloc) = facet
       end do
    end do

    ! The identity first, then the periodic images.
    m = 0
    if (present(images)) m = size(images, 2)
    allocate(transforms(4, m + 1))
    transforms(:, 1) = 0.0_rp
    if (m .gt. 0) transforms(:, 2:m + 1) = images

    call element_spheres(this, msh, dm_Xh)
    call exchange_halo(this, xloc, oloc, nloc, transforms)
    allocate(this%fel(max(this%n_facets, 1)), this%ffacet(max(this%n_facets, 1)))
    this%fel = 0
    this%ffacet = 0
    if (nloc .gt. 0) then
       this%fel(1:nloc) = eloc(1:nloc)
       this%ffacet(1:nloc) = floc(1:nloc)
    end if
    call finalize_facets(this)
    call build_candidates(this, msh%nelv)

    deallocate(xloc, oloc)

  end subroutine lpt_wall_geometry_init

  !> Release all data.
  subroutine lpt_wall_geometry_free(this)
    class(lpt_wall_geometry_t), intent(inout) :: this

    if (c_associated(this%xf_d)) call device_unmap(this%xf, this%xf_d)
    if (c_associated(this%fr_d)) call device_unmap(this%dev_fr, this%fr_d)
    if (c_associated(this%fi_d)) call device_unmap(this%dev_fi, this%fi_d)
    if (c_associated(this%zg_d)) call device_unmap(this%zg, this%zg_d)
    if (c_associated(this%el_ptr_d)) then
       call device_unmap(this%el_ptr, this%el_ptr_d)
    end if
    if (c_associated(this%el_facets_d)) then
       call device_unmap(this%el_facets, this%el_facets_d)
    end if
    if (c_associated(this%el_radius_d)) then
       call device_unmap(this%el_radius, this%el_radius_d)
    end if
    if (c_associated(this%stats_d)) then
       call device_unmap(this%dev_stats, this%stats_d)
    end if
    if (c_associated(this%hits_d)) call device_free(this%hits_d)
    if (c_associated(this%normals_d)) call device_free(this%normals_d)
    this%scratch_cap = 0
    if (allocated(this%dev_fr)) deallocate(this%dev_fr)
    if (allocated(this%dev_fi)) deallocate(this%dev_fi)
    if (allocated(this%dev_stats)) deallocate(this%dev_stats)

    if (allocated(this%zg)) deallocate(this%zg)
    if (allocated(this%xf)) deallocate(this%xf)
    if (allocated(this%origin)) deallocate(this%origin)
    if (allocated(this%orient)) deallocate(this%orient)
    if (allocated(this%center)) deallocate(this%center)
    if (allocated(this%radius)) deallocate(this%radius)
    if (allocated(this%scale)) deallocate(this%scale)
    if (allocated(this%planar)) deallocate(this%planar)
    if (allocated(this%fel)) deallocate(this%fel)
    if (allocated(this%ffacet)) deallocate(this%ffacet)
    if (allocated(this%cone_axis)) deallocate(this%cone_axis)
    if (allocated(this%cone_angle)) deallocate(this%cone_angle)
    if (allocated(this%slab)) deallocate(this%slab)
    if (allocated(this%corner_x)) deallocate(this%corner_x)
    if (allocated(this%corner_n)) deallocate(this%corner_n)
    if (allocated(this%edge_len)) deallocate(this%edge_len)
    if (allocated(this%edge_c)) deallocate(this%edge_c)
    if (allocated(this%edge_r)) deallocate(this%edge_r)
    if (allocated(this%edge_bulge)) deallocate(this%edge_bulge)
    if (allocated(this%el_center)) deallocate(this%el_center)
    if (allocated(this%el_radius)) deallocate(this%el_radius)
    if (allocated(this%el_ptr)) deallocate(this%el_ptr)
    if (allocated(this%el_facets)) deallocate(this%el_facets)
    this%lx = 0
    this%n_facets = 0
    this%n_local_facets = 0
    this%stats = lpt_wall_stats_t()

  end subroutine lpt_wall_geometry_free

  !> Store the GLL points.
  subroutine init_basis(this, zg, lx)
    type(lpt_wall_geometry_t), intent(inout) :: this
    real(kind=rp), intent(in) :: zg(:)
    integer, intent(in) :: lx

    if (lx .lt. 2) call neko_error("lpt wall geometry needs lx >= 2")
    this%lx = lx
    allocate(this%zg(lx))
    this%zg = zg(1:lx)

  end subroutine init_basis

  !> Copy the GLL nodes of one element facet. The facet parameters (s, t)
  !! follow the element's (r, s, t) directions: (s, t) for facets 1 and 2,
  !! (r, t) for facets 3 and 4, (r, s) for facets 5 and 6.
  subroutine extract_facet(dm_Xh, e, facet, lx, xf)
    type(dofmap_t), intent(in) :: dm_Xh
    integer, intent(in) :: e, facet, lx
    real(kind=rp), intent(out) :: xf(3, lx, lx)
    integer :: a, b, i, j, k

    do b = 1, lx
       do a = 1, lx
          call facet_node_index(facet, a, b, lx, i, j, k)
          xf(1, a, b) = dm_Xh%x%x(i, j, k, e)
          xf(2, a, b) = dm_Xh%y%x(i, j, k, e)
          xf(3, a, b) = dm_Xh%z%x(i, j, k, e)
       end do
    end do

  end subroutine extract_facet

  !> Volume index (i, j, k) of facet node (a, b).
  pure subroutine facet_node_index(facet, a, b, lx, i, j, k)
    integer, intent(in) :: facet, a, b, lx
    integer, intent(out) :: i, j, k

    select case (facet)
    case (1)
       i = 1
       j = a
       k = b
    case (2)
       i = lx
       j = a
       k = b
    case (3)
       i = a
       j = 1
       k = b
    case (4)
       i = a
       j = lx
       k = b
    case (5)
       i = a
       j = b
       k = 1
    case default
       i = a
       j = b
       k = lx
    end select

  end subroutine facet_node_index

  !> Sign of \f$ X_s \times X_t \f$ relative to Neko's outward facet normal,
  !! accumulated over all facet nodes for robustness.
  function facet_orientation(this, xf, coef, e, facet) result(orient)
    type(lpt_wall_geometry_t), intent(in) :: this
    real(kind=rp), intent(in) :: xf(:,:,:)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e, facet
    real(kind=rp) :: orient
    real(kind=rp) :: x(3), xs(3), xt(3), xss(3), xst(3), xtt(3), ncoef(3)
    real(kind=rp) :: acc
    integer :: a, b, i, j, k

    acc = 0.0_rp
    do b = 1, this%lx
       do a = 1, this%lx
          call eval_nodes(this, xf, [this%zg(a), this%zg(b)], x, xs, xt, &
               xss, xst, xtt)
          call facet_node_index(facet, a, b, this%lx, i, j, k)
          ncoef = coef%get_normal(i, j, k, e, facet)
          acc = acc + dot_product(cross(xs, xt), ncoef)
       end do
    end do

    if (acc .ge. 0.0_rp) then
       orient = 1.0_rp
    else
       orient = -1.0_rp
    end if

  end function facet_orientation

  !> Bounding spheres of the local elements. The radius is the largest node
  !! distance from the node centroid, inflated by 10% to cover the
  !! polynomial geometry between the nodes.
  subroutine element_spheres(this, msh, dm_Xh)
    type(lpt_wall_geometry_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(dofmap_t), intent(in) :: dm_Xh
    integer :: e, i, j, k, lx, n
    real(kind=rp) :: c(3), r

    lx = this%lx
    n = lx**3
    allocate(this%el_center(3, msh%nelv), this%el_radius(msh%nelv))
    do e = 1, msh%nelv
       c(1) = sum(dm_Xh%x%x(:,:,:,e)) / n
       c(2) = sum(dm_Xh%y%x(:,:,:,e)) / n
       c(3) = sum(dm_Xh%z%x(:,:,:,e)) / n
       r = 0.0_rp
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                r = max(r, norm2([dm_Xh%x%x(i,j,k,e), dm_Xh%y%x(i,j,k,e), &
                     dm_Xh%z%x(i,j,k,e)] - c))
             end do
          end do
       end do
       this%el_center(:, e) = c
       this%el_radius(e) = 1.1_rp * r
    end do

  end subroutine element_spheres

  !> Send copies of local wall facets, and of their periodic images, to
  !! every rank whose search region they intersect (images also to this
  !! rank), and store local facets followed by the received ones.
  !! @param transforms Identity followed by the periodic images, (4, m).
  subroutine exchange_halo(this, xloc, oloc, nloc, transforms)
    type(lpt_wall_geometry_t), intent(inout) :: this
    real(kind=rp), intent(in) :: xloc(:,:,:,:), oloc(:)
    integer, intent(in) :: nloc
    real(kind=rp), intent(in) :: transforms(:,:)
    real(kind=rp) :: box(6), c(3), r
    real(kind=rp), allocatable :: boxes(:,:), sendbuf(:), recvbuf(:)
    real(kind=rp), allocatable :: ximg(:,:,:)
    integer, allocatable :: send_n(:), recv_n(:), sdisp(:), rdisp(:)
    integer, allocatable :: send_pos(:)
    integer :: q, f, e, nrecv, ierr, lx, rec, pos, k

    lx = this%lx
    rec = 3 * lx * lx + 1

    ! Search region of this rank: elements grown by their reach.
    box(1:3) = huge(1.0_rp)
    box(4:6) = -huge(1.0_rp)
    do e = 1, size(this%el_radius)
       box(1:3) = min(box(1:3), this%el_center(:, e) - 3.0_rp * &
            this%el_radius(e))
       box(4:6) = max(box(4:6), this%el_center(:, e) + 3.0_rp * &
            this%el_radius(e))
    end do

    allocate(boxes(6, 0:pe_size - 1))
    call MPI_Allgather(box, 6, MPI_REAL_PRECISION, boxes, 6, &
         MPI_REAL_PRECISION, NEKO_COMM, ierr)

    allocate(send_n(0:pe_size - 1), recv_n(0:pe_size - 1))
    allocate(sdisp(0:pe_size - 1), rdisp(0:pe_size - 1))
    allocate(ximg(3, lx, lx))
    send_n = 0
    do f = 1, nloc
       do k = 1, size(transforms, 2)
          call transform_nodes(xloc(:,:,:,f), transforms(:, k), ximg)
          call node_sphere(ximg, lx, c, r)
          r = 2.0_rp * r
          do q = 0, pe_size - 1
             if (q .eq. pe_rank .and. k .eq. 1) cycle
             if (sphere_box_overlap(c, r, boxes(:, q))) then
                send_n(q) = send_n(q) + 1
             end if
          end do
       end do
    end do

    call MPI_Alltoall(send_n, 1, MPI_INTEGER, recv_n, 1, MPI_INTEGER, &
         NEKO_COMM, ierr)

    send_n = send_n * rec
    recv_n = recv_n * rec
    sdisp(0) = 0
    rdisp(0) = 0
    do q = 1, pe_size - 1
       sdisp(q) = sdisp(q - 1) + send_n(q - 1)
       rdisp(q) = rdisp(q - 1) + recv_n(q - 1)
    end do

    allocate(sendbuf(max(sum(send_n), 1)), recvbuf(max(sum(recv_n), 1)))
    allocate(send_pos(0:pe_size - 1))
    send_pos = sdisp
    do f = 1, nloc
       do k = 1, size(transforms, 2)
          call transform_nodes(xloc(:,:,:,f), transforms(:, k), ximg)
          call node_sphere(ximg, lx, c, r)
          r = 2.0_rp * r
          do q = 0, pe_size - 1
             if (q .eq. pe_rank .and. k .eq. 1) cycle
             if (.not. sphere_box_overlap(c, r, boxes(:, q))) cycle
             pos = send_pos(q)
             sendbuf(pos + 1:pos + rec - 1) = reshape(ximg, [rec - 1])
             ! A rigid motion keeps the orientation of the parametrisation.
             sendbuf(pos + rec) = oloc(f)
             send_pos(q) = pos + rec
          end do
       end do
    end do

    call MPI_Alltoallv(sendbuf, send_n, sdisp, MPI_REAL_PRECISION, &
         recvbuf, recv_n, rdisp, MPI_REAL_PRECISION, NEKO_COMM, ierr)

    nrecv = sum(recv_n) / rec
    this%n_local_facets = nloc
    this%n_facets = nloc + nrecv
    allocate(this%xf(3, lx, lx, max(this%n_facets, 1)))
    allocate(this%orient(max(this%n_facets, 1)))
    if (nloc .gt. 0) then
       this%xf(:,:,:,1:nloc) = xloc(:,:,:,1:nloc)
       this%orient(1:nloc) = oloc(1:nloc)
    end if
    do f = 1, nrecv
       pos = (f - 1) * rec
       this%xf(:,:,:,nloc + f) = reshape(recvbuf(pos + 1:pos + rec - 1), &
            [3, lx, lx])
       this%orient(nloc + f) = recvbuf(pos + rec)
    end do

  end subroutine exchange_halo

  !> Apply a rotation about the z axis by t(4) followed by the translation
  !! t(1:3) to facet nodes.
  pure subroutine transform_nodes(xf, t, xout)
    real(kind=rp), intent(in) :: xf(:,:,:), t(4)
    real(kind=rp), intent(out) :: xout(:,:,:)
    real(kind=rp) :: ct, st

    ct = cos(t(4))
    st = sin(t(4))
    xout(1,:,:) = ct * xf(1,:,:) - st * xf(2,:,:) + t(1)
    xout(2,:,:) = st * xf(1,:,:) + ct * xf(2,:,:) + t(2)
    xout(3,:,:) = xf(3,:,:) + t(3)

  end subroutine transform_nodes

  !> Node centroid and largest node distance of a facet.
  pure subroutine node_sphere(xf, lx, c, r)
    real(kind=rp), intent(in) :: xf(:,:,:)
    integer, intent(in) :: lx
    real(kind=rp), intent(out) :: c(3), r
    integer :: a, b

    c(1) = sum(xf(1,:,:)) / (lx * lx)
    c(2) = sum(xf(2,:,:)) / (lx * lx)
    c(3) = sum(xf(3,:,:)) / (lx * lx)
    r = 0.0_rp
    do b = 1, lx
       do a = 1, lx
          r = max(r, norm2(xf(:, a, b) - c))
       end do
    end do

  end subroutine node_sphere

  !> Whether a sphere overlaps an axis-aligned box (lo(3), hi(3)).
  pure logical function sphere_box_overlap(c, r, box) result(overlap)
    real(kind=rp), intent(in) :: c(3), r, box(6)
    real(kind=rp) :: q(3)

    q = max(box(1:3), min(c, box(4:6)))
    overlap = norm2(q - c) .le. r

  end function sphere_box_overlap

  !> Bounding spheres, normal cones, scales and planarity of all facets.
  subroutine finalize_facets(this)
    type(lpt_wall_geometry_t), intent(inout) :: this
    integer :: f, a, b, nf, k
    real(kind=rp) :: x(3), xs(3), xt(3), xss(3), xst(3), xtt(3), n0(3)
    real(kind=rp) :: dev, nrm, cloc(3)
    logical :: ok

    nf = max(this%n_facets, 1)
    allocate(this%center(3, nf), this%radius(nf), this%scale(nf))
    allocate(this%origin(3, nf))
    allocate(this%cone_axis(3, nf), this%cone_angle(nf), this%planar(nf))
    allocate(this%slab(2, nf))
    allocate(this%corner_x(3, 4, nf), this%corner_n(3, 4, nf))
    allocate(this%edge_len(4, nf), this%edge_c(3, 4, nf))
    allocate(this%edge_r(4, nf), this%edge_bulge(4, nf))
    this%planar = .false.
    this%edge_len = 0.0_rp

    do f = 1, this%n_facets
       ! Move to facet-local coordinates.
       this%origin(1, f) = sum(this%xf(1,:,:,f)) / (this%lx * this%lx)
       this%origin(2, f) = sum(this%xf(2,:,:,f)) / (this%lx * this%lx)
       this%origin(3, f) = sum(this%xf(3,:,:,f)) / (this%lx * this%lx)
       do b = 1, this%lx
          do a = 1, this%lx
             this%xf(:, a, b, f) = this%xf(:, a, b, f) - this%origin(:, f)
          end do
       end do

       do k = 1, 4
          call surface_frame(this, f, corner_st(k), this%corner_x(:, k, f), &
               xs, xt, this%corner_n(:, k, f), ok)
          this%edge_len(k, f) = norm2(edge_node(this, f, k, this%lx) - &
               edge_node(this, f, k, 1))
          if (this%edge_len(k, f) .le. 0.0_rp) cycle
          call edge_sphere(this, f, k, this%edge_c(:, k, f), &
               this%edge_r(k, f), this%edge_bulge(k, f))
       end do

       call patch_bounds(this, f, [-1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp], &
            this%center(:, f), this%radius(f), this%cone_axis(:, f), &
            this%cone_angle(f), this%slab(:, f))
       this%scale(f) = max(maxval(abs(this%center(:, f))), this%radius(f))
       cloc = this%center(:, f)
       ! Bounding sphere centre in global coordinates, for the candidate
       ! search; everything else stays facet-local.
       this%center(:, f) = cloc + this%origin(:, f)

       ! A surface spanned by nodes lying in one plane is exactly planar,
       ! because the Lagrange basis sums to one.
       call eval_local(this, f, [0.0_rp, 0.0_rp], x, xs, xt, xss, xst, xtt)
       n0 = cross(xs, xt)
       nrm = norm2(n0)
       if (nrm .le. 0.0_rp) cycle
       n0 = this%orient(f) * n0 / nrm
       dev = 0.0_rp
       do b = 1, this%lx
          do a = 1, this%lx
             dev = max(dev, abs(dot_product(this%xf(:, a, b, f) - x, n0)))
          end do
       end do
       if (dev .le. 64.0_rp * NEKO_EPS * this%scale(f)) then
          this%planar(f) = .true.
          this%cone_axis(:, f) = n0
          this%cone_angle(f) = 0.0_rp
          this%slab(:, f) = dot_product(x - cloc, n0) + [-dev, dev]
       end if
    end do

  end subroutine finalize_facets

  !> Bounding sphere and normal cone of the patch box = [s0,s1] x [t0,t1]
  !! of facet f. The surface is sampled on (2 lx - 1)^2 points of the box;
  !! the sphere radius is the largest sampled distance from the sample mean
  !! plus half a sample spacing times the largest sampled tangent in each
  !! direction, and the cone half-angle is the largest sampled angle from
  !! the mean normal plus the largest angle between neighbouring samples,
  !! both enlarged by 1%. The slab is the range of (X - c) . axis over the
  !! samples, grown by the same derivative margin as the radius.
  subroutine patch_bounds(this, f, box, c, r, axis, theta, slab)
    type(lpt_wall_geometry_t), intent(in) :: this
    integer, intent(in) :: f
    real(kind=rp), intent(in) :: box(4)
    real(kind=rp), intent(out) :: c(3), r, axis(3), theta, slab(2)
    real(kind=rp) :: xs_(3), xt_(3), nv(3), margin, h
    real(kind=rp), allocatable :: xp(:,:,:), np(:,:,:)
    real(kind=rp) :: hs, ht, gs, gt, nrm, dmax
    integer :: m, a, b

    m = max(3, 2 * this%lx - 1)
    allocate(xp(3, m, m), np(3, m, m))
    hs = (box(2) - box(1)) / (m - 1)
    ht = (box(4) - box(3)) / (m - 1)
    gs = 0.0_rp
    gt = 0.0_rp
    do b = 1, m
       do a = 1, m
          call eval_first(this, f, [box(1) + (a - 1) * hs, &
               box(3) + (b - 1) * ht], xp(:, a, b), xs_, xt_)
          gs = max(gs, norm2(xs_))
          gt = max(gt, norm2(xt_))
          nv = cross(xs_, xt_)
          nrm = norm2(nv)
          if (nrm .gt. 0.0_rp) then
             np(:, a, b) = this%orient(f) * nv / nrm
          else
             np(:, a, b) = 0.0_rp
          end if
       end do
    end do

    c(1) = sum(xp(1,:,:)) / (m * m)
    c(2) = sum(xp(2,:,:)) / (m * m)
    c(3) = sum(xp(3,:,:)) / (m * m)
    r = 0.0_rp
    do b = 1, m
       do a = 1, m
          r = max(r, norm2(xp(:, a, b) - c))
       end do
    end do
    margin = 0.5_rp * (hs * gs + ht * gt)
    r = 1.01_rp * (r + margin)

    axis(1) = sum(np(1,:,:))
    axis(2) = sum(np(2,:,:))
    axis(3) = sum(np(3,:,:))
    nrm = norm2(axis)
    if (nrm .le. 0.0_rp .or. any(norm2(np, dim = 1) .le. 0.0_rp)) then
       ! Degenerate normals: no cone, forces subdivision.
       axis = [0.0_rp, 0.0_rp, 1.0_rp]
       theta = huge(1.0_rp)
       slab = [-r, r]
       return
    end if
    axis = axis / nrm
    slab = [huge(1.0_rp), -huge(1.0_rp)]
    do b = 1, m
       do a = 1, m
          h = dot_product(xp(:, a, b) - c, axis)
          slab = [min(slab(1), h), max(slab(2), h)]
       end do
    end do
    slab = slab + 1.01_rp * margin * [-1.0_rp, 1.0_rp]
    theta = 0.0_rp
    dmax = 0.0_rp
    do b = 1, m
       do a = 1, m
          theta = max(theta, angle(np(:, a, b), axis))
          if (a .gt. 1) dmax = max(dmax, angle(np(:, a, b), np(:, a - 1, b)))
          if (b .gt. 1) dmax = max(dmax, angle(np(:, a, b), np(:, a, b - 1)))
       end do
    end do
    theta = 1.01_rp * (theta + dmax)

  end subroutine patch_bounds

  !> Angle between two unit vectors.
  pure real(kind=rp) function angle(u, v)
    real(kind=rp), intent(in) :: u(3), v(3)

    angle = atan2(norm2(cross(u, v)), dot_product(u, v))

  end function angle

  !> Parameter interval [a_in, a_out] (clipped to [0, 1]) in which the path
  !! x0 + alpha dx lies inside the sphere (c, r).
  pure subroutine segment_sphere_interval(x0, dx, c, r, hit, a_in, a_out)
    real(kind=rp), intent(in) :: x0(3), dx(3), c(3), r
    logical, intent(out) :: hit
    real(kind=rp), intent(out) :: a_in, a_out
    real(kind=rp) :: qa, qb, qc, disc, sq

    hit = .false.
    a_in = 0.0_rp
    a_out = 0.0_rp
    qa = dot_product(dx, dx)
    qb = dot_product(dx, x0 - c)
    qc = dot_product(x0 - c, x0 - c) - r * r
    if (qa .le. 0.0_rp) then
       hit = qc .le. 0.0_rp
       a_out = 1.0_rp
       return
    end if
    disc = qb * qb - qa * qc
    if (disc .lt. 0.0_rp) return
    sq = sqrt(disc)
    a_in = max(0.0_rp, (-qb - sq) / qa)
    a_out = min(1.0_rp, (-qb + sq) / qa)
    hit = a_in .le. a_out
  end subroutine segment_sphere_interval

  !> Surface point, tangents, unit normal and, optionally, the normal
  !! derivatives of facet f.
  subroutine surface_frame(this, f, st, x, xs, xt, nrm, ok, ns, nt)
    type(lpt_wall_geometry_t), intent(in) :: this
    integer, intent(in) :: f
    real(kind=rp), intent(in) :: st(2)
    real(kind=rp), intent(out) :: x(3), xs(3), xt(3), nrm(3)
    logical, intent(out) :: ok
    real(kind=rp), intent(out), optional :: ns(3), nt(3)
    real(kind=rp) :: xss(3), xst(3), xtt(3), nv(3), nu(3), ds(3), dt(3)
    real(kind=rp) :: nlen

    if (present(ns)) then
       call eval_local(this, f, st, x, xs, xt, xss, xst, xtt)
    else
       call eval_first(this, f, st, x, xs, xt)
    end if
    nv = cross(xs, xt)
    nlen = norm2(nv)
    ok = nlen .gt. 0.0_rp
    if (.not. ok) then
       nrm = 0.0_rp
       if (present(ns)) ns = 0.0_rp
       if (present(nt)) nt = 0.0_rp
       return
    end if
    nu = nv / nlen
    nrm = this%orient(f) * nu
    if (present(ns) .and. present(nt)) then
       ! d(N/|N|) = (I - n n^T) dN / |N|
       ds = cross(xss, xt) + cross(xs, xst)
       dt = cross(xst, xt) + cross(xs, xtt)
       ns = this%orient(f) * (ds - nu * dot_product(nu, ds)) / nlen
       nt = this%orient(f) * (dt - nu * dot_product(nu, dt)) / nlen
    end if

  end subroutine surface_frame

  !> Newton's method for \f$ F(\alpha,s,t) = x_0 + \alpha \Delta x - X_w +
  !! a n_w = 0 \f$. With exact = .false. the normal derivatives are left out
  !! of the Jacobian, which is exact for a = 0 and for planar facets. Steps
  !! are halved until |F| decreases; iterates are kept in the patch grown by
  !! half its width (and in [-1.25, 1.25]^2) and alpha in [-0.5, 1.5].
  !! Convergence: |F| at round-off level, or no further decrease of |F|
  !! once it is below a few hundred round-off units. An iteration that
  !! stalls with the iterate held at the boundary of the search region
  !! returns NEWTON_OUTSIDE: the root lies outside the region, beyond the
  !! patch by at least half its width, or at alpha outside [-0.5, 1.5].
  subroutine newton_intersect(this, f, x0, dx, a, box, st0, alpha0, exact, &
       st, alpha, nrm, iters, status)
    type(lpt_wall_geometry_t), intent(in) :: this
    integer, intent(in) :: f
    real(kind=rp), intent(in) :: x0(3), dx(3), a, box(4), st0(2), alpha0
    logical, intent(in) :: exact
    real(kind=rp), intent(out) :: st(2), alpha, nrm(3)
    integer, intent(out) :: iters, status
    real(kind=rp) :: lo(2), hi(2), ext, tolf, tolh, fv(3), fn, c1(3), c2(3)
    real(kind=rp) :: c3(3)
    real(kind=rp) :: det, dscale, d(3), x(3), xs(3), xt(3), ns(3), nt(3)
    real(kind=rp) :: sttry(2), atry, ftry(3), fntry, xtry(3), ntry(3), lam
    real(kind=rp) :: xst_(3), xtt_(3), nstry(3), nttry(3), c23(3), stepn
    logical :: ok, oktry
    integer :: it, ih

    ext = 0.5_rp * max(box(2) - box(1), box(4) - box(3))
    lo = max(-1.25_rp, [box(1), box(3)] - ext)
    hi = min(1.25_rp, [box(2), box(4)] + ext)
    tolf = 64.0_rp * NEKO_EPS * (this%scale(f) + maxval(abs(x0)))
    tolh = 4.0_rp * NEKO_EPS * (this%scale(f) + maxval(abs(x0)))
    status = NEWTON_NOT_CONVERGED
    st = max(lo, min(hi, st0))
    alpha = max(-0.5_rp, min(1.5_rp, alpha0))
    iters = 0

    ns = 0.0_rp
    nt = 0.0_rp
    nstry = 0.0_rp
    nttry = 0.0_rp
    if (exact) then
       call surface_frame(this, f, st, x, xs, xt, nrm, ok, ns, nt)
    else
       call surface_frame(this, f, st, x, xs, xt, nrm, ok)
    end if
    if (.not. ok) then
       status = NEWTON_SINGULAR
       return
    end if
    fv = x0 + alpha * dx - x + a * nrm
    fn = norm2(fv)

    do it = 1, NEWTON_MAX_ITER
       iters = it
       ! Iterate down to round-off: near grazing incidence a residual of
       ! size tolf moves the contact along the path by tolf / sin(angle).
       if (fn .le. tolh) then
          status = NEWTON_CONVERGED
          return
       end if

       c1 = dx
       c2 = -xs
       c3 = -xt
       if (exact .and. a .gt. 0.0_rp) then
          c2 = c2 + a * ns
          c3 = c3 + a * nt
       end if
       c23 = cross(c2, c3)
       det = dot_product(c1, c23)
       dscale = norm2(c1) * norm2(c2) * norm2(c3)
       if (abs(det) .le. 1.0e3_rp * NEKO_EPS * dscale) then
          status = NEWTON_SINGULAR
          if (fn .le. 16.0_rp * tolf) status = NEWTON_CONVERGED
          return
       end if
       ! Cramer's rule for J d = -F
       d(1) = -dot_product(fv, c23) / det
       d(2) = -dot_product(c1, cross(fv, c3)) / det
       d(3) = -dot_product(c1, cross(c2, fv)) / det

       lam = 1.0_rp
       do ih = 1, NEWTON_MAX_HALVINGS
          atry = max(-0.5_rp, min(1.5_rp, alpha + lam * d(1)))
          sttry = max(lo, min(hi, st + lam * d(2:3)))
          if (exact) then
             call surface_frame(this, f, sttry, xtry, xst_, xtt_, ntry, &
                  oktry, nstry, nttry)
          else
             call surface_frame(this, f, sttry, xtry, xst_, xtt_, ntry, oktry)
          end if
          if (oktry) then
             ftry = x0 + atry * dx - xtry + a * ntry
             fntry = norm2(ftry)
             if (fntry .lt. fn) exit
          end if
          lam = 0.5_rp * lam
       end do
       if (.not. oktry .or. fntry .ge. fn) then
          ! No decrease: converged at round-off level, or stuck.
          if (fn .le. 16.0_rp * tolf) then
             status = NEWTON_CONVERGED
          else if (at_bound(st, alpha, lo, hi)) then
             status = NEWTON_OUTSIDE
          end if
          return
       end if

       stepn = max(abs(atry - alpha), maxval(abs(sttry - st)))
       alpha = atry
       st = sttry
       x = xtry
       xs = xst_
       xt = xtt_
       nrm = ntry
       ns = nstry
       nt = nttry
       fv = ftry
       fn = fntry
       if (stepn .le. 100.0_rp * NEKO_EPS .and. fn .le. 16.0_rp * tolf) then
          status = NEWTON_CONVERGED
          return
       end if
    end do
    if (fn .le. 16.0_rp * tolf) then
       status = NEWTON_CONVERGED
    else if (at_bound(st, alpha, lo, hi)) then
       status = NEWTON_OUTSIDE
    end if

  end subroutine newton_intersect

  !> Whether a Newton iterate is held at the boundary of its search region.
  pure logical function at_bound(st, alpha, lo, hi)
    real(kind=rp), intent(in) :: st(2), alpha, lo(2), hi(2)

    at_bound = any(st .le. lo) .or. any(st .ge. hi) .or. &
         alpha .le. -0.5_rp .or. alpha .ge. 1.5_rp

  end function at_bound

  !> Earliest outgoing intersection, with alpha < alpha_max, of the path
  !! x0 + alpha dx (alpha in [0, 1]) with the offset facet X_w - a n_w.
  !!
  !! Patches, starting from the whole facet, are processed from a stack:
  !! - a patch whose bounding sphere, grown by a, misses the path before
  !!   the best alpha so far is dropped;
  !! - a patch whose normals all point against dx (normal cone) can only be
  !!   entered, never left, and is dropped;
  !! - a patch whose normals all make an angle below 90 degrees with dx is a
  !!   graph along dx, so the path crosses it at most once: Newton gives
  !!   that crossing, a converged root outside the patch proves there is
  !!   none, and a failed solve leads to subdivision;
  !! - other patches are subdivided into four, down to PATCH_MAX_DEPTH,
  !!   where Newton is applied without the uniqueness guarantee.
  !! A failed Newton solve thus never by itself implies "no collision".
  !! @param f Facet index.
  !! @param x0g Path start (global coordinates).
  !! @param dx Path displacement.
  !! @param a Particle radius.
  !! @param st_guess Initial reference coordinates for the whole facet.
  !! @param have_guess Whether st_guess is valid.
  !! @param alpha On input the largest alpha of interest, on output the
  !! earliest intersection if found.
  !! @param found Whether an intersection before the input alpha was found.
  !! @param nc Outward wall normal at the intersection.
  !! @param is_rim Whether the intersection is an edge or corner contact.
  !! @param alpha_tie An edge or corner contact is preferred to a surface
  !! intersection only if it is earlier by more than this (default 0). At
  !! smooth facet junctions both coincide, and the surface normal is the
  !! better conditioned one.
  subroutine lpt_wall_facet_intersect(this, f, x0g, dx, a, st_guess, &
       have_guess, alpha, found, nc, is_rim, alpha_tie)
    class(lpt_wall_geometry_t), intent(inout) :: this
    integer, intent(in) :: f
    real(kind=rp), intent(in) :: x0g(3), dx(3), a, st_guess(2)
    logical, intent(in) :: have_guess
    real(kind=rp), intent(inout) :: alpha
    logical, intent(out) :: found
    real(kind=rp), intent(out) :: nc(3)
    logical, intent(out), optional :: is_rim
    real(kind=rp), intent(in), optional :: alpha_tie
    real(kind=rp) :: al_r, nc_r(3), tie, x0(3), cl(3)
    logical :: fd_r
    real(kind=rp) :: boxes(4, PATCH_STACK), box(4), c(3), r, axis(3), theta
    real(kind=rp) :: slab(2), hmax
    real(kind=rp) :: a_in, a_out, dxn(3), cosang, sth, st0(2), al0, xg(3)
    real(kind=rp) :: xs(3), xt(3), ng(3), den, st(2), al, nrm(3), tol_st
    real(kind=rp) :: tol_al, sm, tm
    integer :: depth(PATCH_STACK), top, dep, iters, status, k
    logical :: hit, graph, ok, subdivide, valid

    found = .false.
    nc = 0.0_rp
    if (present(is_rim)) is_rim = .false.
    if (norm2(dx) .le. 0.0_rp) return
    ! Facet-local coordinates.
    x0 = x0g - this%origin(:, f)
    cl = this%center(:, f) - this%origin(:, f)
    dxn = dx / norm2(dx)
    tol_st = 1.0e3_rp * NEKO_EPS
    tol_al = 64.0_rp * NEKO_EPS * this%scale(f) / norm2(dx)

    call segment_sphere_interval(x0, dx, cl, this%radius(f) + a, hit, &
         a_in, a_out)
    if (.not. hit .or. a_in .gt. alpha) return
    hmax = max(dot_product(x0 - cl, this%cone_axis(:, f)), &
         dot_product(x0 + dx - cl, this%cone_axis(:, f)))
    if (hmax .lt. this%slab(1, f) - a) return

    ! Edge and corner contacts of a finite sphere; they lie within a of
    ! the facet and so pass the same filters.
    tie = 0.0_rp
    if (present(alpha_tie)) tie = alpha_tie
    al_r = alpha
    fd_r = .false.
    nc_r = 0.0_rp
    if (a .gt. 0.0_rp) call rim_intersect(this, f, x0, dx, a, al_r, fd_r, nc_r)

    top = 1
    boxes(:, 1) = [-1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp]
    depth(1) = 0
    do while (top .ge. 1)
       box = boxes(:, top)
       dep = depth(top)
       top = top - 1

       if (dep .eq. 0) then
          c = cl
          r = this%radius(f)
          axis = this%cone_axis(:, f)
          theta = this%cone_angle(f)
          slab = this%slab(:, f)
       else
          call patch_bounds(this, f, box, c, r, axis, theta, slab)
       end if
       call segment_sphere_interval(x0, dx, c, r + a, hit, a_in, a_out)
       if (.not. hit .or. a_in .gt. alpha) cycle
       ! Slab test: the offset surface lies above slab(1) - a along the
       ! axis; a path staying below that cannot reach it.
       hmax = max(dot_product(x0 - c, axis), dot_product(x0 + dx - c, axis))
       if (hmax .lt. slab(1) - a) cycle

       cosang = dot_product(dxn, axis)
       graph = .false.
       if (theta .lt. 0.5_rp * acos(-1.0_rp)) then
          sth = sin(theta)
          if (cosang .lt. -sth) cycle
          graph = cosang .gt. sth
       end if

       subdivide = .false.
       if (graph .or. dep .eq. PATCH_MAX_DEPTH) then
          ! Initial guess: the given coordinates for the whole facet, the
          ! patch centre otherwise; alpha from the offset tangent plane.
          sm = 0.5_rp * (box(1) + box(2))
          tm = 0.5_rp * (box(3) + box(4))
          st0 = [sm, tm]
          if (dep .eq. 0 .and. have_guess) st0 = max(-1.0_rp, &
               min(1.0_rp, st_guess))
          call surface_frame(this, f, st0, xg, xs, xt, ng, ok)
          al0 = 0.5_rp * (a_in + a_out)
          if (ok) then
             den = dot_product(dx, ng)
             if (abs(den) .gt. NEKO_EPS * norm2(dx)) then
                al0 = dot_product(xg - a * ng - x0, ng) / den
             end if
          end if

          call newton_intersect(this, f, x0, dx, a, box, st0, al0, .false., &
               st, al, nrm, iters, status)
          if ((status .eq. NEWTON_SINGULAR .or. &
               status .eq. NEWTON_NOT_CONVERGED) .and. a .gt. 0.0_rp .and. &
               .not. this%planar(f)) then
             call newton_intersect(this, f, x0, dx, a, box, st0, al0, &
                  .true., st, al, nrm, iters, status)
          end if

          valid = .false.
          if (status .eq. NEWTON_CONVERGED) then
             call accept_root(st, al, nrm, valid)
          end if
          if (.not. valid .and. .not. graph .and. &
               dep .eq. PATCH_MAX_DEPTH) then
             ! Several roots possible (e.g. the contact point of a previous
             ! bounce and the next contact): retry from the far end.
             call newton_intersect(this, f, x0, dx, a, box, [sm, tm], a_out, &
                  a .gt. 0.0_rp, st, al, nrm, iters, status)
             if (status .eq. NEWTON_CONVERGED) then
                call accept_root(st, al, nrm, valid)
             end if
             if (.not. valid) then
                this%stats%n_unresolved = this%stats%n_unresolved + 1_i8
             end if
          else if (.not. valid) then
             if (status .eq. NEWTON_CONVERGED .or. &
                  (status .eq. NEWTON_OUTSIDE .and. graph)) then
                ! The unique crossing is outside this patch or this step.
             else if (dep .lt. PATCH_MAX_DEPTH) then
                subdivide = .true.
             else
                this%stats%n_unresolved = this%stats%n_unresolved + 1_i8
             end if
          end if
       else
          subdivide = .true.
       end if

       if (subdivide) then
          sm = 0.5_rp * (box(1) + box(2))
          tm = 0.5_rp * (box(3) + box(4))
          do k = 1, 4
             top = top + 1
             depth(top) = dep + 1
             select case (k)
             case (1)
                boxes(:, top) = [box(1), sm, box(3), tm]
             case (2)
                boxes(:, top) = [sm, box(2), box(3), tm]
             case (3)
                boxes(:, top) = [box(1), sm, tm, box(4)]
             case default
                boxes(:, top) = [sm, box(2), tm, box(4)]
             end select
          end do
       end if
    end do

    if (fd_r .and. (.not. found .or. al_r .lt. alpha - tie)) then
       alpha = al_r
       nc = nc_r
       found = .true.
       if (present(is_rim)) is_rim = .true.
    end if

  contains

    !> Record a converged root if it is a valid outgoing intersection: on
    !! the facet (any patch), alpha in [0, 1] up to round-off, moving into
    !! the wall, and earlier than the best so far. Accepting a genuine
    !! root found outside the current patch can only improve the result.
    subroutine accept_root(st_r, al_r, n_r, ok_r)
      real(kind=rp), intent(in) :: st_r(2), al_r, n_r(3)
      logical, intent(out) :: ok_r

      ok_r = all(abs(st_r) .le. 1.0_rp + tol_st) .and. &
           al_r .ge. -tol_al .and. al_r .le. 1.0_rp .and. &
           dot_product(dx, n_r) .gt. 0.0_rp
      if (ok_r .and. max(0.0_rp, al_r) .lt. alpha) then
         alpha = max(0.0_rp, al_r)
         found = .true.
         nc = n_r
      end if
    end subroutine accept_root

  end subroutine lpt_wall_facet_intersect

  !> Earliest contact, with alpha < alpha_max, of a sphere of radius a > 0
  !! moving on x0 + alpha dx with the four edges and four corners of facet
  !! f. An edge contact is an intersection of the path with the tube of
  !! radius a around the edge curve E(u), found by Newton's method on
  !! \f$ (|x - E|^2 - a^2)/(2a) = 0 \f$, \f$ (x - E) \cdot E'/|E'| = 0 \f$
  !! from the entry point of the tube around the edge chord; a corner
  !! contact is the entry into the sphere of radius a around the corner.
  !! Only contacts of a sphere approaching the edge from the fluid side of
  !! the facet count; contacts over the facet interior are always preceded
  !! by the offset-surface intersection.
  subroutine rim_intersect(this, f, x0, dx, a, alpha, found, nc)
    type(lpt_wall_geometry_t), intent(inout) :: this
    integer, intent(in) :: f
    real(kind=rp), intent(in) :: x0(3), dx(3), a
    real(kind=rp), intent(inout) :: alpha
    logical, intent(inout) :: found
    real(kind=rp), intent(inout) :: nc(3)
    real(kind=rp) :: e0(3), e1(3), ed(3), x(3), xs(3), xt(3), nf(3)
    real(kind=rp) :: dperp(3), qa, qb, qc, disc, al, u, tolf, rr(3)
    real(kind=rp) :: g(2), jm(2,2), det, du, dal, ep(3), epp(3)
    real(kind=rp) :: b1, v, dm, gn, lam
    real(kind=rp) :: altry, utry, gtry(2), rtry(3), eptry(3), epptry(3)
    real(kind=rp) :: al_close, al_entry, u_close
    integer :: k, it, ih, status, attempt
    logical :: ok, conv

    tolf = 64.0_rp * NEKO_EPS * (this%scale(f) + maxval(abs(x0)))

    ! Corners: |x0 + alpha dx - V|^2 = a^2
    do k = 1, 4
       e0 = this%corner_x(:, k, f)
       nf = this%corner_n(:, k, f)
       qa = dot_product(dx, dx)
       qb = 2.0_rp * dot_product(dx, x0 - e0)
       qc = dot_product(x0 - e0, x0 - e0) - a * a
       if (qb .ge. 0.0_rp .or. qa .le. 0.0_rp) cycle
       disc = qb * qb - 4.0_rp * qa * qc
       if (disc .lt. 0.0_rp) cycle
       al = max(0.0_rp, (-qb - sqrt(disc)) / (2.0_rp * qa))
       if (al .gt. 1.0_rp .or. al .ge. alpha) cycle
       x = x0 + al * dx
       if (dot_product(x - e0, nf) .gt. tolf) cycle
       alpha = al
       found = .true.
       nc = (e0 - x) / max(norm2(e0 - x), tiny(1.0_rp))
    end do

    ! Edges: k = 1, 2 at s = -1, 1 (parameter t); k = 3, 4 at t = -1, 1.
    do k = 1, 4
       if (this%edge_len(k, f) .le. 0.0_rp) cycle
       call segment_sphere_interval(x0, dx, this%edge_c(:, k, f), &
            this%edge_r(k, f) + a, conv, qa, qb)
       if (.not. conv .or. qa .gt. alpha) cycle

       ! Closest approach of the path to the polyline through the edge's GLL
       ! nodes; the edge curve deviates from it by at most edge_bulge.
       dm = huge(1.0_rp)
       al_close = 0.0_rp
       u_close = 0.0_rp
       ed = 0.0_rp
       do it = 1, this%lx - 1
          e0 = edge_node(this, f, k, it)
          e1 = edge_node(this, f, k, it + 1)
          call segment_segment(x0, dx, e0, e1 - e0, al, v, b1)
          if (b1 .lt. dm) then
             dm = b1
             al_close = al
             u_close = this%zg(it) + v * (this%zg(it + 1) - this%zg(it))
             ed = (e1 - e0) / max(norm2(e1 - e0), tiny(1.0_rp))
          end if
       end do
       if (dm .gt. a + this%edge_bulge(k, f)) cycle
       dperp = dx - dot_product(dx, ed) * ed
       if (norm2(dperp) .le. NEKO_EPS * norm2(dx)) cycle
       al_entry = al_close - sqrt(max(a * a - dm * dm, 0.0_rp)) / norm2(dperp)

       ! First from the tube entry estimated with the polyline, then, for
       ! paths grazing the tube, from the closest approach.
       do attempt = 1, 2
          u = u_close
          if (attempt .eq. 1) then
             al = max(-0.5_rp, min(1.5_rp, al_entry))
          else
             al = max(-0.5_rp, min(1.5_rp, al_close))
          end if
          conv = .false.
          status = NEWTON_NOT_CONVERGED
          call edge_point(this, f, k, u, e0, ep, epp)
          rr = x0 + al * dx - e0
          g = [(dot_product(rr, rr) - a * a) / (2.0_rp * a), &
               dot_product(rr, ep) / norm2(ep)]
          gn = norm2(g)
          do it = 1, NEWTON_MAX_ITER
             if (gn .le. tolf) then
                conv = .true.
                exit
             end if
             if (al .le. -0.5_rp .or. al .ge. 1.5_rp) then
                ! Tube entry far outside this step.
                status = NEWTON_OUTSIDE
                exit
             end if
             jm(1,1) = dot_product(rr, dx) / a
             jm(1,2) = -dot_product(rr, ep) / a
             jm(2,1) = dot_product(dx, ep) / norm2(ep)
             jm(2,2) = (-dot_product(ep, ep) + dot_product(rr, epp)) / norm2(ep)
             det = jm(1,1) * jm(2,2) - jm(1,2) * jm(2,1)
             if (abs(det) .le. 1.0e3_rp * NEKO_EPS * (abs(jm(1,1) * jm(2,2)) + &
                  abs(jm(1,2) * jm(2,1)))) then
                ! Path tangent to the tube while outside it: no entry.
                if (g(1) .gt. 0.0_rp) status = NEWTON_OUTSIDE
                exit
             end if
             dal = -(jm(2,2) * g(1) - jm(1,2) * g(2)) / det
             du = -(jm(1,1) * g(2) - jm(2,1) * g(1)) / det
             lam = 1.0_rp
             do ih = 1, NEWTON_MAX_HALVINGS
                altry = max(-0.5_rp, min(1.5_rp, al + lam * dal))
                utry = max(-1.25_rp, min(1.25_rp, u + lam * du))
                call edge_point(this, f, k, utry, e1, eptry, epptry)
                rtry = x0 + altry * dx - e1
                gtry = [(dot_product(rtry, rtry) - a * a) / (2.0_rp * a), &
                     dot_product(rtry, eptry) / norm2(eptry)]
                if (norm2(gtry) .lt. gn) exit
                lam = 0.5_rp * lam
             end do
             if (norm2(gtry) .ge. gn) then
                conv = gn .le. 16.0_rp * tolf
                ! Stalled outside the tube (|x - E| > a): the path does not
                ! come within a of the edge, no contact.
                if (.not. conv .and. (g(1) .gt. 0.0_rp .or. &
                     abs(u) .ge. 1.25_rp)) status = NEWTON_OUTSIDE
                exit
             end if
             al = altry
             u = utry
             e0 = e1
             ep = eptry
             epp = epptry
             rr = rtry
             g = gtry
             gn = norm2(g)
          end do
          if (conv) status = NEWTON_CONVERGED
          if (status .ne. NEWTON_NOT_CONVERGED) exit
       end do
       if (status .eq. NEWTON_NOT_CONVERGED) then
          this%stats%n_unresolved = this%stats%n_unresolved + 1_i8
       end if
       if (.not. conv) cycle
       if (abs(u) .gt. 1.0_rp .or. al .lt. 0.0_rp .or. al .gt. 1.0_rp) cycle
       if (al .ge. alpha) cycle
       ! Approaching the edge, from the fluid side of the facet.
       if (dot_product(dx, rr) .ge. 0.0_rp) cycle
       call surface_frame(this, f, edge_st(k, u), x, xs, xt, nf, ok)
       if (.not. ok .or. dot_product(rr, nf) .gt. tolf) cycle
       alpha = al
       found = .true.
       nc = -rr / norm2(rr)
    end do

  end subroutine rim_intersect

  !> Facet coordinates of corner k.
  pure function corner_st(k) result(st)
    integer, intent(in) :: k
    real(kind=rp) :: st(2)

    st = [merge(-1.0_rp, 1.0_rp, mod(k, 2) .eq. 1), &
         merge(-1.0_rp, 1.0_rp, k .le. 2)]

  end function corner_st

  !> Facet coordinates of edge k at parameter u.
  pure function edge_st(k, u) result(st)
    integer, intent(in) :: k
    real(kind=rp), intent(in) :: u
    real(kind=rp) :: st(2)

    select case (k)
    case (1)
       st = [-1.0_rp, u]
    case (2)
       st = [1.0_rp, u]
    case (3)
       st = [u, -1.0_rp]
    case default
       st = [u, 1.0_rp]
    end select

  end function edge_st

  !> Point, first and second derivative of edge k of facet f at u.
  subroutine edge_point(this, f, k, u, e, ep, epp)
    type(lpt_wall_geometry_t), intent(in) :: this
    integer, intent(in) :: f, k
    real(kind=rp), intent(in) :: u
    real(kind=rp), intent(out) :: e(3), ep(3), epp(3)
    real(kind=rp) :: xs(3), xt(3), xss(3), xst(3), xtt(3)

    call eval_local(this, f, edge_st(k, u), e, xs, xt, xss, xst, xtt)
    if (k .le. 2) then
       ep = xt
       epp = xtt
    else
       ep = xs
       epp = xss
    end if

  end subroutine edge_point

  !> Bounding sphere of edge k of facet f, and the largest distance of the
  !! edge curve from the polyline through its GLL nodes. Both come from
  !! samples, four intervals per node interval; the sphere radius is grown
  !! by half a sample spacing times the largest tangent, the deviation by
  !! the largest distance of the curve from the chord between two samples,
  !! h^2/8 times the largest second derivative.
  subroutine edge_sphere(this, f, k, c, r, bulge)
    type(lpt_wall_geometry_t), intent(in) :: this
    integer, intent(in) :: f, k
    real(kind=rp), intent(out) :: c(3), r, bulge
    real(kind=rp) :: e(3), ep(3), epp(3), g, h, u, p0(3), p1(3), q(3), d(3)
    real(kind=rp) :: t, hmax, g2
    integer :: i, j, m

    m = 0
    c = 0.0_rp
    g = 0.0_rp
    g2 = 0.0_rp
    bulge = 0.0_rp
    hmax = 0.0_rp
    do i = 1, this%lx - 1
       p0 = edge_node(this, f, k, i)
       p1 = edge_node(this, f, k, i + 1)
       d = p1 - p0
       h = 0.25_rp * (this%zg(i + 1) - this%zg(i))
       hmax = max(hmax, h)
       do j = 0, 4
          u = this%zg(i) + j * h
          call edge_point(this, f, k, u, e, ep, epp)
          g = max(g, norm2(ep))
          g2 = max(g2, norm2(epp))
          c = c + e
          m = m + 1
          q = e - p0
          t = max(0.0_rp, min(1.0_rp, dot_product(q, d) / &
               max(dot_product(d, d), tiny(1.0_rp))))
          bulge = max(bulge, norm2(q - t * d))
       end do
    end do
    c = c / m
    r = 0.0_rp
    do i = 1, this%lx - 1
       h = 0.25_rp * (this%zg(i + 1) - this%zg(i))
       do j = 0, 4
          call edge_point(this, f, k, this%zg(i) + j * h, e, ep, epp)
          r = max(r, norm2(e - c))
       end do
    end do
    r = 1.01_rp * (r + 0.5_rp * hmax * g)
    bulge = 1.01_rp * (bulge + 0.125_rp * hmax**2 * g2)

  end subroutine edge_sphere

  !> GLL node i of edge k of facet f.
  pure function edge_node(this, f, k, i) result(x)
    type(lpt_wall_geometry_t), intent(in) :: this
    integer, intent(in) :: f, k, i
    real(kind=rp) :: x(3)

    select case (k)
    case (1)
       x = this%xf(:, 1, i, f)
    case (2)
       x = this%xf(:, this%lx, i, f)
    case (3)
       x = this%xf(:, i, 1, f)
    case default
       x = this%xf(:, i, this%lx, f)
    end select

  end function edge_node

  !> Closest points of the segments p0 + s d0 and q0 + t d1, s, t in
  !! [0, 1], and their distance.
  pure subroutine segment_segment(p0, d0, q0, d1, s, t, dist)
    real(kind=rp), intent(in) :: p0(3), d0(3), q0(3), d1(3)
    real(kind=rp), intent(out) :: s, t, dist
    real(kind=rp) :: a, b, c, d, e, den, w(3)

    w = p0 - q0
    a = dot_product(d0, d0)
    b = dot_product(d0, d1)
    c = dot_product(d1, d1)
    d = dot_product(d0, w)
    e = dot_product(d1, w)
    den = a * c - b * b
    if (den .gt. NEKO_EPS * a * c) then
       s = max(0.0_rp, min(1.0_rp, (b * e - c * d) / den))
    else
       s = 0.0_rp
    end if
    if (c .gt. 0.0_rp) then
       t = (b * s + e) / c
    else
       t = 0.0_rp
    end if
    if (t .lt. 0.0_rp .or. t .gt. 1.0_rp) then
       t = max(0.0_rp, min(1.0_rp, t))
       if (a .gt. 0.0_rp) then
          s = max(0.0_rp, min(1.0_rp, (b * t - d) / a))
       else
          s = 0.0_rp
       end if
    end if
    dist = norm2(p0 + s * d0 - q0 - t * d1)

  end subroutine segment_segment

  !> Facet coordinates of a point from its reference coordinates in the
  !! element that owns the facet: the two coordinates tangential to the
  !! facet, in the facet's (s, t) order.
  pure function facet_st_from_rst(facet, rst) result(st)
    integer, intent(in) :: facet
    real(kind=rp), intent(in) :: rst(3)
    real(kind=rp) :: st(2)

    select case (facet)
    case (1, 2)
       st = [rst(2), rst(3)]
    case (3, 4)
       st = [rst(1), rst(3)]
    case default
       st = [rst(1), rst(2)]
    end select
    st = max(-1.0_rp, min(1.0_rp, st))

  end function facet_st_from_rst

  !> Resolve the wall collisions of one particle in one step.
  !! The earliest outgoing intersection over all candidate facets of element
  !! el (and, for a > 0, their edges and corners) is found, the particle is
  !! placed at \f$ x_c = x^n + \alpha_c \Delta x \f$, the remaining
  !! displacement \f$ D = (1 - \alpha_c) \Delta x \f$ is reflected,
  !! \f$ D^+ = D - (1 + e_n)(D \cdot n_c) n_c \f$, and the search is
  !! repeated on the reflected segment, up to LPT_WALL_MAX_BOUNCES times.
  !! The contact point is moved inward by a round-off sized distance so
  !! that the next search starts on the fluid side.
  !! @param el Local element (1-based) that owned the particle at p0.
  !! @param p0 Particle centre before the step.
  !! @param p1 Predicted particle centre, updated in place.
  !! @param a Particle radius.
  !! @param e_n Normal coefficient of restitution.
  !! @param n_hits Number of collisions applied.
  !! @param normals Contact normals of the collisions, (3, n_hits).
  !! @param rst_old Optional reference coordinates of p0 in element el, used
  !! as initial guesses on the facets of el.
  subroutine lpt_wall_collide(this, el, p0, p1, a, e_n, n_hits, normals, &
       rst_old)
    class(lpt_wall_geometry_t), intent(inout) :: this
    integer, intent(in) :: el
    real(kind=rp), intent(in) :: p0(3), a, e_n
    real(kind=rp), intent(inout) :: p1(3)
    integer, intent(out) :: n_hits
    real(kind=rp), intent(out) :: normals(3, LPT_WALL_MAX_BOUNCES)
    real(kind=rp), intent(in), optional :: rst_old(3)
    real(kind=rp) :: q0(3), q1(3), dx(3), r(3), nc(3), nbest(3), alpha
    real(kind=rp) :: rn, st_guess(2), tol_off, al_c, a_tie
    integer :: k, f, fbest, bounce
    logical :: found, have_guess, rim, rim_best

    n_hits = 0
    normals = 0.0_rp
    if (el .lt. 1 .or. el .gt. size(this%el_radius)) return
    if (this%el_ptr(el) .ge. this%el_ptr(el + 1)) return

    if (norm2(p1 - p0) + a .gt. this%el_radius(el)) then
       this%stats%n_reach_exceeded = this%stats%n_reach_exceeded + 1_i8
    end if

    q0 = p0
    q1 = p1
    do bounce = 1, LPT_WALL_MAX_BOUNCES + 1
       dx = q1 - q0
       if (norm2(dx) .le. 0.0_rp) exit
       alpha = huge(1.0_rp)
       fbest = 0
       rim_best = .false.
       a_tie = (128.0_rp * NEKO_EPS * (maxval(abs(q0)) + norm2(dx)) + &
            1.0e-3_rp * a) / norm2(dx)
       do k = this%el_ptr(el), this%el_ptr(el + 1) - 1
          f = this%el_facets(k)
          have_guess = .false.
          st_guess = 0.0_rp
          if (bounce .eq. 1 .and. present(rst_old) .and. &
               this%fel(f) .eq. el .and. this%ffacet(f) .gt. 0) then
             st_guess = facet_st_from_rst(this%ffacet(f), rst_old)
             have_guess = .true.
          end if
          ! A surface intersection replaces an edge/corner contact found on
          ! another facet unless the latter is earlier by more than a_tie.
          al_c = alpha
          if (rim_best) al_c = alpha + a_tie
          call this%intersect(f, q0, dx, a, st_guess, have_guess, al_c, &
               found, nc, rim, a_tie)
          if (found) then
             if (rim .and. .not. rim_best) then
                found = al_c .lt. alpha - a_tie
             else if (.not. rim .and. rim_best) then
                found = al_c .lt. alpha + a_tie
             else
                found = al_c .lt. alpha
             end if
          end if
          if (found) then
             alpha = al_c
             fbest = f
             nbest = nc
             rim_best = rim
          end if
       end do
       if (fbest .eq. 0) exit

       tol_off = 16.0_rp * NEKO_EPS * this%scale(fbest)
       if (n_hits .eq. LPT_WALL_MAX_BOUNCES) then
          ! Out of bounces: stop at this contact, inside the domain. The
          ! next step starts in contact and collides at alpha = 0.
          q1 = q0 + alpha * dx - tol_off * nbest
          exit
       end if
       r = (1.0_rp - alpha) * dx
       q0 = q0 + alpha * dx - tol_off * nbest
       rn = dot_product(r, nbest)
       if (rn .gt. 0.0_rp) r = r - (1.0_rp + e_n) * rn * nbest
       q1 = q0 + r
       n_hits = n_hits + 1
       normals(:, n_hits) = nbest
    end do

    p1 = q1

  end subroutine lpt_wall_collide

  !> Pack the per-facet data for the device kernels (layout in
  !! lpt_wall_intersect.h): origin(3), orientation, facet-local bounding
  !! sphere centre(3) and radius, cone axis(3) and half-angle, slab(2),
  !! scale, planar flag, edge chord lengths(4), edge sphere centres(3x4),
  !! radii(4) and polyline deviations(4), corner points(3x4) and corner
  !! normals(3x4).
  subroutine lpt_wall_pack(this)
    class(lpt_wall_geometry_t), intent(inout) :: this
    integer :: f, nf

    nf = max(this%n_facets, 1)
    if (allocated(this%dev_fr)) deallocate(this%dev_fr)
    if (allocated(this%dev_fi)) deallocate(this%dev_fi)
    if (allocated(this%dev_stats)) deallocate(this%dev_stats)
    allocate(this%dev_fr(LPT_WALL_NFR, nf), this%dev_fi(2, nf))
    allocate(this%dev_stats(LPT_WALL_NSTATS))
    this%dev_fr = 0.0_rp
    this%dev_fi = 0
    this%dev_stats = 0
    do f = 1, this%n_facets
       this%dev_fr(1:3, f) = this%origin(:, f)
       this%dev_fr(4, f) = this%orient(f)
       this%dev_fr(5:7, f) = this%center(:, f) - this%origin(:, f)
       this%dev_fr(8, f) = this%radius(f)
       this%dev_fr(9:11, f) = this%cone_axis(:, f)
       this%dev_fr(12, f) = min(this%cone_angle(f), huge(1.0_rp) / 4)
       this%dev_fr(13:14, f) = this%slab(:, f)
       this%dev_fr(15, f) = this%scale(f)
       this%dev_fr(16, f) = merge(1.0_rp, 0.0_rp, this%planar(f))
       this%dev_fr(17:20, f) = this%edge_len(:, f)
       this%dev_fr(21:32, f) = reshape(this%edge_c(:, :, f), [12])
       this%dev_fr(33:36, f) = this%edge_r(:, f)
       this%dev_fr(37:40, f) = this%edge_bulge(:, f)
       this%dev_fr(41:52, f) = reshape(this%corner_x(:, :, f), [12])
       this%dev_fr(53:64, f) = reshape(this%corner_n(:, :, f), [12])
       this%dev_fi(1, f) = this%fel(f)
       this%dev_fi(2, f) = this%ffacet(f)
    end do

  end subroutine lpt_wall_pack

  !> Pack the geometry and copy it to the device, once at setup.
  subroutine lpt_wall_to_device(this)
    class(lpt_wall_geometry_t), intent(inout) :: this
    integer, allocatable :: tmp(:)

    if (this%lx .gt. LPT_WALL_MAX_LX_DEVICE) then
       call neko_error("lpt wall collisions on the device support lx <= 16")
    end if
    call this%pack()
    ! Never map empty arrays.
    if (size(this%el_facets) .eq. 0) then
       allocate(tmp(1))
       tmp = 0
       call move_alloc(tmp, this%el_facets)
    end if
    call device_map(this%xf, this%xf_d, size(this%xf))
    call device_memcpy(this%xf, this%xf_d, size(this%xf), HOST_TO_DEVICE, &
         sync = .false.)
    call device_map(this%dev_fr, this%fr_d, size(this%dev_fr))
    call device_memcpy(this%dev_fr, this%fr_d, size(this%dev_fr), &
         HOST_TO_DEVICE, sync = .false.)
    call device_map(this%dev_fi, this%fi_d, size(this%dev_fi))
    call device_memcpy(this%dev_fi, this%fi_d, size(this%dev_fi), &
         HOST_TO_DEVICE, sync = .false.)
    call device_map(this%zg, this%zg_d, size(this%zg))
    call device_memcpy(this%zg, this%zg_d, size(this%zg), HOST_TO_DEVICE, &
         sync = .false.)
    call device_map(this%el_ptr, this%el_ptr_d, size(this%el_ptr))
    call device_memcpy(this%el_ptr, this%el_ptr_d, size(this%el_ptr), &
         HOST_TO_DEVICE, sync = .false.)
    call device_map(this%el_facets, this%el_facets_d, size(this%el_facets))
    call device_memcpy(this%el_facets, this%el_facets_d, &
         size(this%el_facets), HOST_TO_DEVICE, sync = .false.)
    call device_map(this%el_radius, this%el_radius_d, size(this%el_radius))
    call device_memcpy(this%el_radius, this%el_radius_d, &
         size(this%el_radius), HOST_TO_DEVICE, sync = .false.)
    call device_map(this%dev_stats, this%stats_d, LPT_WALL_NSTATS)
    call device_memcpy(this%dev_stats, this%stats_d, LPT_WALL_NSTATS, &
         HOST_TO_DEVICE, sync = .true.)

  end subroutine lpt_wall_to_device

  !> Make sure the per-particle device scratch of the OpenCL and Metal
  !! kernels holds n particles; grown by 50% when too small.
  subroutine lpt_wall_reserve_scratch(this, n)
    class(lpt_wall_geometry_t), intent(inout) :: this
    integer, intent(in) :: n
    integer :: cap

    if (n .le. this%scratch_cap) return
    if (c_associated(this%hits_d)) call device_free(this%hits_d)
    if (c_associated(this%normals_d)) call device_free(this%normals_d)
    cap = max(n, (3 * this%scratch_cap) / 2)
    call device_alloc(this%hits_d, int(cap, c_size_t) * 4_c_size_t)
    call device_alloc(this%normals_d, int(cap, c_size_t) * &
         int(3 * LPT_WALL_MAX_BOUNCES, c_size_t) * &
         int(storage_size(1.0_rp) / 8, c_size_t))
    this%scratch_cap = cap

  end subroutine lpt_wall_reserve_scratch

  !> Add the per-launch device counters to the statistics and reset them on
  !! the device. Called after every launch, so the 32-bit counters cannot
  !! overflow and the warnings are raised in the step they occur.
  subroutine lpt_wall_gather_stats(this)
    class(lpt_wall_geometry_t), intent(inout) :: this

    call device_memcpy(this%dev_stats, this%stats_d, LPT_WALL_NSTATS, &
         DEVICE_TO_HOST, sync = .true.)
    if (all(this%dev_stats .eq. 0)) return
    this%stats%n_unresolved = this%stats%n_unresolved + &
         int(this%dev_stats(1), i8)
    this%stats%n_reach_exceeded = this%stats%n_reach_exceeded + &
         int(this%dev_stats(2), i8)
    this%dev_stats = 0
    call device_memcpy(this%dev_stats, this%stats_d, LPT_WALL_NSTATS, &
         HOST_TO_DEVICE, sync = .false.)

  end subroutine lpt_wall_gather_stats

  !> Per-element candidate facet lists in CSR form. A facet is a candidate
  !! of element e if its bounding sphere comes within two element radii of
  !! the element centre: the particle starts within one radius and may move,
  !! together with its own radius, by at most one more. Facets are binned on
  !! a uniform grid to avoid testing all element-facet pairs.
  subroutine build_candidates(this, nelv)
    type(lpt_wall_geometry_t), intent(inout) :: this
    integer, intent(in) :: nelv
    real(kind=rp) :: lo(3), hi(3), cell, reach
    integer :: dims(3), ncell, f, e, i, j, k, ijk, pass, nnz, cnt
    integer :: ilo(3), ihi(3)
    integer, allocatable :: cell_ptr(:), cell_list(:), fill(:), mark(:)
    real(kind=rp) :: rmax

    allocate(this%el_ptr(nelv + 1))
    if (this%n_facets .eq. 0) then
       this%el_ptr = 1
       allocate(this%el_facets(0))
       return
    end if

    lo = huge(1.0_rp)
    hi = -huge(1.0_rp)
    rmax = 0.0_rp
    do f = 1, this%n_facets
       lo = min(lo, this%center(:, f) - this%radius(f))
       hi = max(hi, this%center(:, f) + this%radius(f))
       rmax = max(rmax, this%radius(f))
    end do
    cell = max(2.0_rp * rmax, tiny(1.0_rp))
    do
       dims = max(1, ceiling((hi - lo) / cell))
       if (real(dims(1), rp) * dims(2) * dims(3) .le. &
            8.0_rp * this%n_facets + 64.0_rp) exit
       cell = 2.0_rp * cell
    end do
    ncell = dims(1) * dims(2) * dims(3)

    ! Bin facets into all cells overlapped by their bounding box.
    allocate(cell_ptr(ncell + 1), fill(ncell))
    cell_ptr = 0
    do pass = 1, 2
       fill = 0
       do f = 1, this%n_facets
          call cell_range(this%center(:, f), this%radius(f), lo, cell, dims, &
               ilo, ihi)
          do k = ilo(3), ihi(3)
             do j = ilo(2), ihi(2)
                do i = ilo(1), ihi(1)
                   ijk = i + dims(1) * ((j - 1) + dims(2) * (k - 1))
                   fill(ijk) = fill(ijk) + 1
                   if (pass .eq. 2) then
                      cell_list(cell_ptr(ijk) + fill(ijk) - 1) = f
                   end if
                end do
             end do
          end do
       end do
       if (pass .eq. 1) then
          cell_ptr(1) = 1
          do ijk = 1, ncell
             cell_ptr(ijk + 1) = cell_ptr(ijk) + fill(ijk)
          end do
          allocate(cell_list(cell_ptr(ncell + 1) - 1))
       end if
    end do

    ! Element queries, counted first and stored second.
    allocate(mark(this%n_facets))
    mark = 0
    nnz = 0
    do pass = 1, 2
       nnz = 0
       do e = 1, nelv
          if (pass .eq. 2) this%el_ptr(e) = nnz + 1
          reach = 2.0_rp * this%el_radius(e)
          call cell_range(this%el_center(:, e), reach, lo, cell, dims, &
               ilo, ihi)
          cnt = 0
          do k = ilo(3), ihi(3)
             do j = ilo(2), ihi(2)
                do i = ilo(1), ihi(1)
                   ijk = i + dims(1) * ((j - 1) + dims(2) * (k - 1))
                   do f = cell_ptr(ijk), cell_ptr(ijk + 1) - 1
                      if (mark(cell_list(f)) .eq. e) cycle
                      mark(cell_list(f)) = e
                      if (norm2(this%center(:, cell_list(f)) - &
                           this%el_center(:, e)) .gt. &
                           this%radius(cell_list(f)) + reach) cycle
                      nnz = nnz + 1
                      if (pass .eq. 2) this%el_facets(nnz) = cell_list(f)
                   end do
                end do
             end do
          end do
       end do
       if (pass .eq. 1) then
          allocate(this%el_facets(nnz))
          mark = 0
       end if
    end do
    this%el_ptr(nelv + 1) = nnz + 1

  end subroutine build_candidates

  !> Range of grid cells overlapped by the bounding box of a sphere,
  !! empty (ilo > ihi) if the box misses the grid.
  pure subroutine cell_range(c, r, lo, cell, dims, ilo, ihi)
    real(kind=rp), intent(in) :: c(3), r, lo(3), cell
    integer, intent(in) :: dims(3)
    integer, intent(out) :: ilo(3), ihi(3)
    real(kind=rp) :: a(3), b(3)

    a = (c - r - lo) / cell
    b = (c + r - lo) / cell
    where (b .lt. 0.0_rp .or. a .gt. real(dims, rp))
       ilo = 1
       ihi = 0
    elsewhere
       ilo = max(1, min(dims, int(a) + 1))
       ihi = max(1, min(dims, int(b) + 1))
    end where
    if (any(ihi .lt. ilo)) then
       ilo = 1
       ihi = 0
    end if

  end subroutine cell_range

  !> Evaluate facet f in facet-local coordinates.
  subroutine eval_local(this, f, st, x, xs, xt, xss, xst, xtt)
    type(lpt_wall_geometry_t), intent(in) :: this
    integer, intent(in) :: f
    real(kind=rp), intent(in) :: st(2)
    real(kind=rp), intent(out) :: x(3), xs(3), xt(3), xss(3), xst(3), xtt(3)

    call eval_nodes(this, this%xf(:,:,:,f), st, x, xs, xt, xss, xst, xtt)

  end subroutine eval_local

  !> Evaluate a tensor-product GLL surface and its derivatives using the
  !! Lagrange weights of fd_weights_full.
  subroutine eval_nodes(this, xf, st, x, xs, xt, xss, xst, xtt)
    type(lpt_wall_geometry_t), intent(in) :: this
    real(kind=rp), intent(in) :: xf(:,:,:)
    real(kind=rp), intent(in) :: st(2)
    real(kind=rp), intent(out) :: x(3), xs(3), xt(3), xss(3), xst(3), xtt(3)
    real(kind=rp) :: cs(this%lx, 0:2), ct(this%lx, 0:2)
    real(kind=rp) :: a0(3), a1(3), a2(3)
    integer :: a, b, lx

    lx = this%lx
    call fd_weights_full(st(1), this%zg, lx - 1, 2, cs)
    call fd_weights_full(st(2), this%zg, lx - 1, 2, ct)

    x = 0.0_rp
    xs = 0.0_rp
    xt = 0.0_rp
    xss = 0.0_rp
    xst = 0.0_rp
    xtt = 0.0_rp
    do b = 1, lx
       a0 = 0.0_rp
       a1 = 0.0_rp
       a2 = 0.0_rp
       do a = 1, lx
          a0 = a0 + xf(:, a, b) * cs(a, 0)
          a1 = a1 + xf(:, a, b) * cs(a, 1)
          a2 = a2 + xf(:, a, b) * cs(a, 2)
       end do
       x = x + a0 * ct(b, 0)
       xs = xs + a1 * ct(b, 0)
       xt = xt + a0 * ct(b, 1)
       xss = xss + a2 * ct(b, 0)
       xst = xst + a1 * ct(b, 1)
       xtt = xtt + a0 * ct(b, 2)
    end do

  end subroutine eval_nodes

  !> Evaluate a stored facet and its first derivatives only.
  subroutine eval_first(this, f, st, x, xs, xt)
    type(lpt_wall_geometry_t), intent(in) :: this
    integer, intent(in) :: f
    real(kind=rp), intent(in) :: st(2)
    real(kind=rp), intent(out) :: x(3), xs(3), xt(3)
    real(kind=rp) :: cs(this%lx, 0:1), ct(this%lx, 0:1), a0(3), a1(3)
    integer :: a, b, lx

    lx = this%lx
    call fd_weights_full(st(1), this%zg, lx - 1, 1, cs)
    call fd_weights_full(st(2), this%zg, lx - 1, 1, ct)
    x = 0.0_rp
    xs = 0.0_rp
    xt = 0.0_rp
    do b = 1, lx
       a0 = 0.0_rp
       a1 = 0.0_rp
       do a = 1, lx
          a0 = a0 + this%xf(:, a, b, f) * cs(a, 0)
          a1 = a1 + this%xf(:, a, b, f) * cs(a, 1)
       end do
       x = x + a0 * ct(b, 0)
       xs = xs + a1 * ct(b, 0)
       xt = xt + a0 * ct(b, 1)
    end do

  end subroutine eval_first

  !> Cross product.
  pure function cross(a, b) result(c)
    real(kind=rp), intent(in) :: a(3), b(3)
    real(kind=rp) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end function cross

end module lpt_wall_geometry
