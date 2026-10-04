! Copyright (c) 2021-2026, The Neko Authors
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
!> Subroutines to add advection terms to the RHS of a transport equation.
!!
!! Dealiased (over-integrated) advection. The convecting velocity and the
!! convected field are interpolated onto a finer Gauss-Legendre (GL) space,
!! the product is formed there, and the result is projected back onto the
!! simulation (GLL) space with the transposed interpolation. Written in the
!! "rst form" of Nek5000's `convect_new`: the weighted contravariant
!! convecting velocity
!! \f$ c_r = w_3 J (\partial_r x_j) u_j \f$ (and \f$ c_s, c_t \f$) is formed
!! once per call and each convected field \f$ u \f$ then only needs
!! \f$ c_r \partial_r u + c_s \partial_s u + c_t \partial_t u \f$.
!!
!! The GL metrics (the cofactors \f$ J \partial r_i / \partial x_j \f$ at the
!! GL points) can be formed in two ways, see adv_dealias_t::metrics. With
!! the default `exact` choice the GLL coordinates are interpolated onto the GL
!! space (exact, since they are polynomials of the simulation order) and
!! differentiated there, which gives the exact cofactors at the GL points. The
!! `interpolated` choice interpolates the GLL cofactors instead, as earlier
!! versions did; the two coincide to rounding on affine elements and differ at
!! geometry-interpolation level on curved ones, where `exact` is the more
!! accurate. Either way the metrics are rebuilt on every call by default on
!! device backends, and built once and stored on the CPU, SX and XSMM backends
!! where the rebuild is costly, see adv_dealias_t::store_metrics.
!!
!! On the device backends the elements are processed in chunks of
!! adv_dealias_t::chunk elements, every operation being element local, so the
!! GL work arrays are proportional to the chunk size rather than to the mesh,
!! and they come from the scratch registry so that they are shared between the
!! fluid, every scalar and the OIFS scheme. The CPU, SX and XSMM backends work
!! element by element on the stack and need no work arrays at all.
module adv_dealias
  use advection, only : advection_t
  use num_types, only : rp, i8
  use math, only : sub2, add2, copy
  use mxm_wrapper, only : mxm
  use space, only : space_t, GL
  use field, only : field_t
  use coefs, only : coef_t
  use vector, only : vector_t
  use device_math, only : device_sub2, device_add2, device_col3
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use logger, only : neko_log, LOG_SIZE
  use profiler, only : profiler_start_region, profiler_end_region
  use interpolation, only : interpolator_t
  use device, only : device_mp_count, device_total_mem, device_map, &
       device_unmap
  use device_coef, only : device_coef_generate_dxydrst
  use tensor_device, only : tnsr3d_device
  use opr_device, only : opr_device_opgrad_ptr, &
       opr_device_set_convect_rst_ptr, opr_device_convect_scalar_gl
  use scratch_registry, only : neko_scratch_registry
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_intptr_t, &
       c_sizeof
  implicit none
  private

  !> Exact GL cofactors, from the coordinates differentiated on the GL space.
  integer, public, parameter :: DEALIAS_METRICS_EXACT = 1
  !> GL cofactors interpolated from the GLL cofactors.
  integer, public, parameter :: DEALIAS_METRICS_INTERPOLATED = 2

  !> Peak number of GL-sized work arrays live at once on the device for one
  !! chunk when the metrics are rebuilt per call: 3 coordinates, 9 forward
  !! derivatives, 9 cofactors and the two Jacobian arrays. Used to size the
  !! chunk.
  integer, parameter :: DEALIAS_PEAK_ARRAYS_REBUILD = 23
  !> The same with stored metrics (the ALE path, 3 mesh velocities, 1 field,
  !! 1 flux, 3 gradients and 1 accumulator).
  integer, parameter :: DEALIAS_PEAK_ARRAYS_STORED = 9
  !> Waves of blocks per launch below which chunking stops paying off: the
  !! partial last wave of a launch then costs at most 1/8 of it.
  integer, parameter :: DEALIAS_MIN_WAVES = 8
  !> Blocks per multiprocessor assumed resident at once when sizing a wave,
  !! typical of the shared memory bound tensor product kernels used here.
  integer, parameter :: DEALIAS_BLOCKS_PER_MP = 4
  !> Share of the device memory the automatic chunk size may use for the
  !! work arrays, as a divisor (50 = 2%).
  integer(kind=i8), parameter :: DEALIAS_MEM_DIVISOR = 50_i8
  !> Floor for the automatic work array budget in bytes (256 MiB).
  integer(kind=i8), parameter :: DEALIAS_MEM_FLOOR = 268435456_i8
  !> Budget used when the device memory is unknown (1 GiB).
  integer(kind=i8), parameter :: DEALIAS_MEM_DEFAULT = 1073741824_i8

  !> A scratch vector together with its registry index.
  type :: dealias_vec_t
     type(vector_t), pointer :: v => null()
     integer :: idx = 0
  end type dealias_vec_t

  !> Type encapsulating advection routines with dealiasing
  type, public, extends(advection_t) :: adv_dealias_t
     !> Coeffs of the original space in the simulation
     type(coef_t), pointer :: coef_GLL => null()
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t) :: GLL_to_GL
     !> The additional higher-order space used in dealiasing
     type(space_t) :: Xh_GL
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL => null()
     !> Elements processed per chunk on whole-field backends.
     integer :: chunk = 0
     !> Chunk size asked for by the case file: 0 automatic, -1 all elements
     !! at once, otherwise the number of elements.
     integer :: chunk_request = 0
     !> How the GL metrics are formed, DEALIAS_METRICS_EXACT or
     !! DEALIAS_METRICS_INTERPOLATED.
     integer :: metrics = DEALIAS_METRICS_EXACT
     !> Whether the GL cofactors are built once and kept in `gl_cof` (9 GL
     !! arrays of mesh size) rather than rebuilt on every call.
     logical :: store_metrics = .false.
     !> Stored GL cofactors, (lxyz_GL, nelv, 9) in the order drdx, drdy,
     !! drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz. Allocated iff
     !! `store_metrics`.
     real(kind=rp), allocatable :: gl_cof(:,:,:)
     !> Device pointer for `gl_cof`
     type(c_ptr) :: gl_cof_d = C_NULL_PTR
   contains
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS.
     procedure, pass(this) :: compute => compute_advection_dealias
     !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_scalar => compute_scalar_advection_dealias
     !> Add the advection term in ALE framework.
     procedure, pass(this) :: compute_ale => compute_ale_advection_dealias
     !> Update any metrics needed for the advection computation in ALE.
     procedure, pass(this) :: recompute_metrics => recompute_metrics_dealias
     !> Constructor
     procedure, pass(this) :: init => init_dealias
     !> Destructor
     procedure, pass(this) :: free => free_dealias
  end type adv_dealias_t

contains

  !> Constructor
  !! @param lxd The polynomial order of the space used in the dealiasing.
  !! @param coef The coefficients of the (space, mesh) pair.
  !! @param chunk Elements per chunk on device backends: 0 or absent for an
  !! automatic choice based on the device, -1 for all elements at once.
  !! @param store_metrics Build the GL cofactors once and keep them (9 GL
  !! arrays of mesh size) rather than rebuilding them on every call. Defaults
  !! to true on the CPU, SX and XSMM backends and false on device backends.
  !! @param metrics `exact` (default) or `interpolated`, see the module
  !! description.
  subroutine init_dealias(this, lxd, coef, chunk, store_metrics, metrics)
    class(adv_dealias_t), target, intent(inout) :: this
    integer, intent(in) :: lxd
    type(coef_t), intent(inout), target :: coef
    integer, intent(in), optional :: chunk
    logical, intent(in), optional :: store_metrics
    character(len=*), intent(in), optional :: metrics
    character(len=LOG_SIZE) :: log_buf
    integer :: nchunks, nelv

    call this%free()

    call this%Xh_GL%init(GL, lxd, lxd, lxd)
    this%Xh_GLL => coef%Xh
    this%coef_GLL => coef
    call this%GLL_to_GL%init(this%Xh_GL, this%Xh_GLL)
    nelv = coef%msh%nelv

    this%metrics = DEALIAS_METRICS_EXACT
    if (present(metrics)) then
       select case (trim(metrics))
       case ('exact')
          this%metrics = DEALIAS_METRICS_EXACT
       case ('interpolated')
          this%metrics = DEALIAS_METRICS_INTERPOLATED
       case default
          call neko_error("Unknown dealias_metrics '" // trim(metrics) // &
               "', expected 'exact' or 'interpolated'")
       end select
    end if

    this%store_metrics = (NEKO_BCKND_DEVICE .ne. 1)
    if (present(store_metrics)) this%store_metrics = store_metrics

    this%chunk_request = 0
    if (present(chunk)) this%chunk_request = max(chunk, -1)
    this%chunk = dealias_chunk_size(this, nelv)

    if (this%store_metrics) then
       allocate(this%gl_cof(this%Xh_GL%lxyz, max(nelv, 1), 9))
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(this%gl_cof, this%gl_cof_d, &
               this%Xh_GL%lxyz * max(nelv, 1) * 9)
       end if
       call dealias_build_stored_metrics(this, coef)
    end if

    if (this%metrics .eq. DEALIAS_METRICS_EXACT) then
       log_buf = 'Dealiasing : exact metrics'
    else
       log_buf = 'Dealiasing : interpolated metrics'
    end if
    if (this%store_metrics) then
       log_buf = trim(log_buf) // ', stored'
    else
       log_buf = trim(log_buf) // ', rebuilt per call'
    end if
    call neko_log%message(log_buf)
    if (whole_field_backend() .and. nelv .gt. 0) then
       nchunks = (nelv + this%chunk - 1) / this%chunk
       if (nchunks .gt. 1) then
          write(log_buf, '(A,I0,A,I0,A)') 'Dealiasing : ', nchunks, &
               ' chunks of ', this%chunk, ' elements'
       else
          log_buf = 'Dealiasing : all elements at once'
       end if
       call neko_log%message(log_buf)
    end if

  end subroutine init_dealias

  !> Destructor
  subroutine free_dealias(this)
    class(adv_dealias_t), intent(inout) :: this

    if (allocated(this%gl_cof)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%gl_cof, this%gl_cof_d)
       end if
       deallocate(this%gl_cof)
    end if
    this%gl_cof_d = C_NULL_PTR

    call this%GLL_to_GL%free()
    call this%Xh_GL%free()

    nullify(this%Xh_GLL)
    nullify(this%coef_GLL)
    this%chunk = 0

  end subroutine free_dealias

  !> Whether the active backend works on whole-field buffers (and hence
  !! chunks the elements), as opposed to the element-local path used on the
  !! CPU, SX and XSMM backends.
  pure function whole_field_backend() result(whole)
    logical :: whole
    whole = (NEKO_BCKND_DEVICE .eq. 1)
  end function whole_field_backend

  !> Choose the number of elements per chunk.
  !!
  !! An explicit request is honoured as is (-1 meaning all elements at
  !! once). Otherwise the chunk is sized so that the peak GL work set stays
  !! within a small share of the device memory, while every launch still
  !! covers at least DEALIAS_MIN_WAVES complete waves of blocks, where a wave
  !! is DEALIAS_BLOCKS_PER_MP blocks on each multiprocessor. Chunks are made
  !! equal in size and a multiple of a wave, so only the last wave of the last
  !! chunk can be partial. Small meshes are left unchunked.
  !! @param nelv Number of local elements.
  function dealias_chunk_size(this, nelv) result(chunk)
    class(adv_dealias_t), intent(in) :: this
    integer, intent(in) :: nelv
    integer :: chunk
    integer(kind=i8) :: budget, per_elem, total_mem, chunk_max
    integer :: nmp, wave, nchunks, min_chunk, peak_arrays

    if (.not. whole_field_backend() .or. nelv .le. 0) then
       chunk = max(nelv, 1)
       return
    end if

    if (this%chunk_request .lt. 0) then
       chunk = nelv
       return
    else if (this%chunk_request .gt. 0) then
       chunk = min(this%chunk_request, nelv)
       return
    end if

    ! Bytes of GL work arrays per element at the peak of a chunk
    if (this%store_metrics) then
       peak_arrays = DEALIAS_PEAK_ARRAYS_STORED
    else
       peak_arrays = DEALIAS_PEAK_ARRAYS_REBUILD
    end if
    per_elem = int(peak_arrays, i8) * int(this%Xh_GL%lxyz, i8) * &
         int(c_sizeof(1.0_rp), i8)

    total_mem = device_total_mem()
    if (total_mem .gt. 0_i8) then
       budget = max(total_mem / DEALIAS_MEM_DIVISOR, DEALIAS_MEM_FLOOR)
    else
       budget = DEALIAS_MEM_DEFAULT
    end if
    chunk_max = max(budget / per_elem, 1_i8)

    ! A wave fills every multiprocessor with resident blocks (one element
    ! per block in the kernels used here)
    nmp = device_mp_count()
    if (nmp .gt. 0) then
       wave = nmp * DEALIAS_BLOCKS_PER_MP
    else
       wave = 1
    end if
    min_chunk = DEALIAS_MIN_WAVES * wave

    if (chunk_max .ge. int(nelv, i8) .or. nelv .le. min_chunk) then
       chunk = nelv
    else
       nchunks = int((int(nelv, i8) + chunk_max - 1_i8) / chunk_max)
       chunk = (nelv + nchunks - 1) / nchunks
       chunk = ((chunk + wave - 1) / wave) * wave
       chunk = max(chunk, min_chunk)
       chunk = min(chunk, nelv)
    end if

  end function dealias_chunk_size

  !> Device pointer `n` reals past `ptr`.
  !! @param ptr Base device pointer.
  !! @param n Offset in reals (not bytes).
  function dev_ptr_offset(ptr, n) result(p)
    type(c_ptr), intent(in) :: ptr
    integer(kind=i8), intent(in) :: n
    type(c_ptr) :: p
    integer(c_intptr_t) :: addr

    addr = transfer(ptr, addr)
    addr = addr + int(n, c_intptr_t) * int(c_sizeof(1.0_rp), c_intptr_t)
    p = transfer(addr, p)

  end function dev_ptr_offset

  !> Device pointer to stored cofactor `k` of elements `e0` onwards.
  function stored_cof_ptr(this, k, e0) result(p)
    class(adv_dealias_t), intent(in) :: this
    integer, intent(in) :: k, e0
    type(c_ptr) :: p
    integer(kind=i8) :: off

    off = (int(k - 1, i8) * int(size(this%gl_cof, 2), i8) + int(e0 - 1, i8)) &
         * int(this%Xh_GL%lxyz, i8)
    p = dev_ptr_offset(this%gl_cof_d, off)

  end function stored_cof_ptr

  !> Request `size(vecs)` scratch vectors of `n` entries each.
  subroutine dealias_request(vecs, n)
    type(dealias_vec_t), intent(inout) :: vecs(:)
    integer, intent(in) :: n
    integer :: i

    do i = 1, size(vecs)
       call neko_scratch_registry%request_vector(vecs(i)%v, vecs(i)%idx, &
            n, .false.)
    end do

  end subroutine dealias_request

  !> Hand the scratch vectors in `vecs` back to the registry.
  subroutine dealias_relinquish(vecs)
    type(dealias_vec_t), intent(inout) :: vecs(:)
    integer :: i

    do i = 1, size(vecs)
       if (associated(vecs(i)%v)) then
          call neko_scratch_registry%relinquish(vecs(i)%idx)
          nullify(vecs(i)%v)
          vecs(i)%idx = 0
       end if
    end do

  end subroutine dealias_relinquish

  !> Build the stored GL cofactors for the whole mesh.
  subroutine dealias_build_stored_metrics(this, coef)
    class(adv_dealias_t), target, intent(inout) :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: work
    type(c_ptr) :: cof_d(9)
    integer :: e, e0, ne, k, nel

    nel = coef%msh%nelv

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do e0 = 1, nel, this%chunk
          ne = min(this%chunk, nel - e0 + 1)
          do k = 1, 9
             cof_d(k) = stored_cof_ptr(this, k, e0)
          end do
          call dealias_device_build_metrics(this, coef, e0, ne, cof_d)
       end do
    else
       !$omp parallel do private(e, work)
       do e = 1, nel
          call dealias_local_build_metrics(this, coef, e, &
               this%gl_cof(:, e, 1), this%gl_cof(:, e, 2), &
               this%gl_cof(:, e, 3), this%gl_cof(:, e, 4), &
               this%gl_cof(:, e, 5), this%gl_cof(:, e, 6), &
               this%gl_cof(:, e, 7), this%gl_cof(:, e, 8), &
               this%gl_cof(:, e, 9), work)
       end do
       !$omp end parallel do
    end if

  end subroutine dealias_build_stored_metrics

  !
  ! ---------------------------------------------------------------------
  ! Device backend: chunked, whole-field kernels
  ! ---------------------------------------------------------------------
  !

  !> Interpolate GLL fields of elements `e0 .. e0+ne-1` onto the GL space
  !! of a chunk.
  !! @param this The object.
  !! @param src_d Device pointers to the full GLL fields.
  !! @param dst_d Device pointers to chunk-sized GL arrays receiving the
  !! interpolants.
  !! @param e0 First element of the chunk.
  !! @param ne Elements in the chunk.
  subroutine dealias_device_to_gl(this, src_d, dst_d, e0, ne)
    class(adv_dealias_t), intent(in) :: this
    type(c_ptr), intent(in) :: src_d(:)
    type(c_ptr), intent(in) :: dst_d(:)
    integer, intent(in) :: e0, ne
    integer(kind=i8) :: off
    integer :: i

    off = int(e0 - 1, i8) * int(this%Xh_GLL%lxyz, i8)
    do i = 1, size(src_d)
       call tnsr3d_device(dst_d(i), this%Xh_GL%lx, &
            dev_ptr_offset(src_d(i), off), this%Xh_GLL%lx, &
            this%GLL_to_GL%Yh_Xh_d, this%GLL_to_GL%Yh_XhT_d, &
            this%GLL_to_GL%Yh_XhT_d, ne)
    end do

  end subroutine dealias_device_to_gl

  !> Build the GL cofactors \f$ J \partial r_i / \partial x_j \f$ of
  !! elements `e0 .. e0+ne-1` into the arrays `cof_d`.
  !!
  !! With exact metrics the GLL coordinates are interpolated onto the GL
  !! space, which is exact for the polynomial geometry, and differentiated
  !! there with the GL derivative matrices; the cofactors follow pointwise
  !! and the 14 arrays used along the way are released on return. With
  !! interpolated metrics the GLL cofactors are interpolated directly.
  !! @param this The object.
  !! @param coef The GLL coefficients.
  !! @param e0 First element of the chunk.
  !! @param ne Elements in the chunk.
  !! @param cof_d Device pointers to the cofactors, ordered drdx, drdy, drdz,
  !! dsdx, dsdy, dsdz, dtdx, dtdy, dtdz.
  subroutine dealias_device_build_metrics(this, coef, e0, ne, cof_d)
    class(adv_dealias_t), intent(in) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e0, ne
    type(c_ptr), intent(in) :: cof_d(9)
    type(dealias_vec_t) :: xyz(3), fwd(9), jac(2)
    type(c_ptr) :: src_d(9), dst_d(3)
    integer :: n_gl, i

    if (this%metrics .eq. DEALIAS_METRICS_INTERPOLATED) then
       src_d(1) = coef%drdx_d
       src_d(2) = coef%drdy_d
       src_d(3) = coef%drdz_d
       src_d(4) = coef%dsdx_d
       src_d(5) = coef%dsdy_d
       src_d(6) = coef%dsdz_d
       src_d(7) = coef%dtdx_d
       src_d(8) = coef%dtdy_d
       src_d(9) = coef%dtdz_d
       call dealias_device_to_gl(this, src_d, cof_d, e0, ne)
       return
    end if

    ! Work vectors are always sized for a full chunk, so that a shorter last
    ! chunk reuses them rather than adding a second set to the registry
    n_gl = this%chunk * this%Xh_GL%lxyz

    call dealias_request(xyz, n_gl)
    call dealias_request(fwd, n_gl)
    call dealias_request(jac, n_gl)

    src_d(1) = coef%dof%x%x_d
    src_d(2) = coef%dof%y%x_d
    src_d(3) = coef%dof%z%x_d
    do i = 1, 3
       dst_d(i) = xyz(i)%v%x_d
    end do
    call dealias_device_to_gl(this, src_d(1:3), dst_d, e0, ne)

    ! fwd holds dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt, dzdt,
    ! jac(1) the inverse Jacobian and jac(2) the Jacobian
    call device_coef_generate_dxydrst(cof_d(1), cof_d(2), cof_d(3), &
         cof_d(4), cof_d(5), cof_d(6), cof_d(7), cof_d(8), cof_d(9), &
         fwd(1)%v%x_d, fwd(2)%v%x_d, fwd(3)%v%x_d, &
         fwd(4)%v%x_d, fwd(5)%v%x_d, fwd(6)%v%x_d, &
         fwd(7)%v%x_d, fwd(8)%v%x_d, fwd(9)%v%x_d, &
         this%Xh_GL%dx_d, this%Xh_GL%dy_d, this%Xh_GL%dz_d, &
         xyz(1)%v%x_d, xyz(2)%v%x_d, xyz(3)%v%x_d, &
         jac(1)%v%x_d, jac(2)%v%x_d, this%Xh_GL%lx, ne)

    call dealias_relinquish(jac)
    call dealias_relinquish(fwd)
    call dealias_relinquish(xyz)

  end subroutine dealias_device_build_metrics

  !> Provide the GL cofactors of a chunk: pointers into the stored metrics,
  !! or freshly built ones in scratch vectors.
  !! @param this The object.
  !! @param coef The GLL coefficients.
  !! @param e0 First element of the chunk.
  !! @param ne Elements in the chunk.
  !! @param cof Scratch vectors holding the cofactors when they are rebuilt;
  !! untouched otherwise. Released by the caller (a no-op if untouched).
  !! @param cof_d Device pointers to the cofactors, ordered drdx, drdy, drdz,
  !! dsdx, dsdy, dsdz, dtdx, dtdy, dtdz.
  subroutine dealias_device_metrics(this, coef, e0, ne, cof, cof_d)
    class(adv_dealias_t), intent(in) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e0, ne
    type(dealias_vec_t), intent(inout) :: cof(9)
    type(c_ptr), intent(inout) :: cof_d(9)
    integer :: k

    if (this%store_metrics) then
       do k = 1, 9
          cof_d(k) = stored_cof_ptr(this, k, e0)
       end do
    else
       call dealias_request(cof, this%chunk * this%Xh_GL%lxyz)
       do k = 1, 9
          cof_d(k) = cof(k)%v%x_d
       end do
       call dealias_device_build_metrics(this, coef, e0, ne, cof_d)
    end if

  end subroutine dealias_device_metrics

  !> Form the weighted contravariant convecting velocity of a chunk.
  !! @param this The object.
  !! @param coef The GLL coefficients.
  !! @param vx_d, vy_d, vz_d Device pointers to the convecting velocity.
  !! @param e0 First element of the chunk.
  !! @param ne Elements in the chunk.
  !! @param c Receives \f$ (c_r, c_s, c_t) \f$; requested here, released by
  !! the caller.
  subroutine dealias_device_convecting(this, coef, vx_d, vy_d, vz_d, &
       e0, ne, c)
    class(adv_dealias_t), intent(in) :: this
    type(coef_t), intent(in) :: coef
    type(c_ptr), intent(in) :: vx_d, vy_d, vz_d
    integer, intent(in) :: e0, ne
    type(dealias_vec_t), intent(inout) :: c(3)
    type(dealias_vec_t) :: cof(9), u(3)
    type(c_ptr) :: cof_d(9), src_d(3), dst_d(3)
    integer :: n_gl, i

    n_gl = this%chunk * this%Xh_GL%lxyz

    call dealias_device_metrics(this, coef, e0, ne, cof, cof_d)

    call dealias_request(u, n_gl)
    src_d(1) = vx_d
    src_d(2) = vy_d
    src_d(3) = vz_d
    do i = 1, 3
       dst_d(i) = u(i)%v%x_d
    end do
    call dealias_device_to_gl(this, src_d, dst_d, e0, ne)

    call dealias_request(c, n_gl)
    call opr_device_set_convect_rst_ptr(c(1)%v%x_d, c(2)%v%x_d, c(3)%v%x_d, &
         u(1)%v%x_d, u(2)%v%x_d, u(3)%v%x_d, &
         cof_d(1), cof_d(4), cof_d(7), &
         cof_d(2), cof_d(5), cof_d(8), &
         cof_d(3), cof_d(6), cof_d(9), &
         this%Xh_GL%w3_d, ne, this%Xh_GL%lx)

    call dealias_relinquish(u)
    call dealias_relinquish(cof)

  end subroutine dealias_device_convecting

  !> Subtract the dealiased \f$ c \cdot \nabla u \f$ of one field from its
  !! right-hand side on a chunk.
  !! @param this The object.
  !! @param u_d Device pointer to the convected GLL field.
  !! @param f_d Device pointer to the GLL right-hand side to update.
  !! @param c The weighted contravariant velocity of the chunk.
  !! @param e0 First element of the chunk.
  !! @param ne Elements in the chunk.
  subroutine dealias_device_convect(this, u_d, f_d, c, e0, ne)
    class(adv_dealias_t), intent(in) :: this
    type(c_ptr), intent(in) :: u_d, f_d
    type(dealias_vec_t), intent(in) :: c(3)
    integer, intent(in) :: e0, ne
    type(dealias_vec_t) :: ug(1), du(1), tg(1)
    type(c_ptr) :: src_d(1), dst_d(1)
    integer(kind=i8) :: off
    integer :: n_gl, n_gll

    n_gl = this%chunk * this%Xh_GL%lxyz
    n_gll = ne * this%Xh_GLL%lxyz
    off = int(e0 - 1, i8) * int(this%Xh_GLL%lxyz, i8)

    call dealias_request(ug, n_gl)
    call dealias_request(du, n_gl)
    call dealias_request(tg, this%chunk * this%Xh_GLL%lxyz)

    src_d(1) = u_d
    dst_d(1) = ug(1)%v%x_d
    call dealias_device_to_gl(this, src_d, dst_d, e0, ne)

    call opr_device_convect_scalar_gl(du(1)%v%x_d, ug(1)%v%x_d, &
         c(1)%v%x_d, c(2)%v%x_d, c(3)%v%x_d, &
         this%Xh_GL%dx_d, this%Xh_GL%dy_d, this%Xh_GL%dz_d, ne, this%Xh_GL%lx)

    ! Transposed interpolation back to the GLL space
    call tnsr3d_device(tg(1)%v%x_d, this%Xh_GLL%lx, du(1)%v%x_d, &
         this%Xh_GL%lx, this%GLL_to_GL%Yh_XhT_d, this%GLL_to_GL%Yh_Xh_d, &
         this%GLL_to_GL%Yh_Xh_d, ne)

    call device_sub2(dev_ptr_offset(f_d, off), tg(1)%v%x_d, n_gll)

    call dealias_relinquish(tg)
    call dealias_relinquish(du)
    call dealias_relinquish(ug)

  end subroutine dealias_device_convect

  !> Weighted gradient on a chunk of the GL space with the chunk's cofactors,
  !! see opgrad.
  subroutine dealias_device_opgrad(this, ux_d, uy_d, uz_d, u_d, cof_d, ne)
    class(adv_dealias_t), intent(in) :: this
    type(c_ptr), intent(inout) :: ux_d, uy_d, uz_d
    type(c_ptr), intent(in) :: u_d
    type(c_ptr), intent(in) :: cof_d(9)
    integer, intent(in) :: ne

    call opr_device_opgrad_ptr(ux_d, uy_d, uz_d, u_d, &
         this%Xh_GL%dx_d, this%Xh_GL%dy_d, this%Xh_GL%dz_d, &
         cof_d(1), cof_d(4), cof_d(7), &
         cof_d(2), cof_d(5), cof_d(8), &
         cof_d(3), cof_d(6), cof_d(9), &
         this%Xh_GL%w3_d, ne, this%Xh_GL%lx)

  end subroutine dealias_device_opgrad

  !
  ! ---------------------------------------------------------------------
  ! CPU backend: element-local helpers
  ! ---------------------------------------------------------------------
  !

  !> Reference-space gradient of one element,
  !! \f$ u_r = D_r u, u_s = D_s u, u_t = D_t u \f$.
  !! @param ur, us, ut The derivatives.
  !! @param u The field on the element.
  !! @param Xh The space (for the derivative matrices).
  subroutine dealias_local_grad(ur, us, ut, u, Xh)
    type(space_t), intent(in) :: Xh
    real(kind=rp), intent(inout) :: ur(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp), intent(inout) :: us(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp), intent(inout) :: ut(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%lx, Xh%lx)
    integer :: k

    associate(lx => Xh%lx)
      call mxm(Xh%dx, lx, u, lx, ur, lx*lx)
      do k = 1, lx
         call mxm(u(1,1,k), lx, Xh%dyt, lx, us(1,1,k), lx)
      end do
      call mxm(u, lx*lx, Xh%dzt, lx, ut, lx)
    end associate

  end subroutine dealias_local_grad

  !> Build the GL cofactors of one element, see adv_dealias_t::metrics.
  !! @param this The object.
  !! @param coef The GLL coefficients.
  !! @param e The element.
  !! @param drdx, ..., dtdz The cofactors \f$ J \partial r_i / \partial x_j
  !! \f$ at the GL points.
  !! @param work Work array of GL size.
  subroutine dealias_local_build_metrics(this, coef, e, drdx, drdy, drdz, &
       dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, work)
    class(adv_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e
    real(kind=rp), dimension(this%Xh_GL%lxyz), intent(inout) :: &
         drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz
    real(kind=rp), dimension(this%Xh_GL%lxyz), intent(inout) :: work
    real(kind=rp) :: xr, xs, xt, yr, ys, yt, zr, zs, zt
    integer :: i

    if (this%metrics .eq. DEALIAS_METRICS_INTERPOLATED) then
       call this%GLL_to_GL%map(drdx, coef%drdx(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(drdy, coef%drdy(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(drdz, coef%drdz(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(dsdx, coef%dsdx(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(dsdy, coef%dsdy(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(dsdz, coef%dsdz(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(dtdx, coef%dtdx(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(dtdy, coef%dtdy(1,1,1,e), 1, this%Xh_GL)
       call this%GLL_to_GL%map(dtdz, coef%dtdz(1,1,1,e), 1, this%Xh_GL)
       return
    end if

    ! Interpolate the coordinates onto the GL space (exact for the
    ! polynomial geometry) and differentiate there: the x derivatives land
    ! in (drdx, drdy, drdz), the y derivatives in (dsdx, dsdy, dsdz) and the
    ! z derivatives in (dtdx, dtdy, dtdz), to be turned into cofactors in
    ! place below
    call this%GLL_to_GL%map(work, coef%dof%x%x(1,1,1,e), 1, this%Xh_GL)
    call dealias_local_grad(drdx, drdy, drdz, work, this%Xh_GL)
    call this%GLL_to_GL%map(work, coef%dof%y%x(1,1,1,e), 1, this%Xh_GL)
    call dealias_local_grad(dsdx, dsdy, dsdz, work, this%Xh_GL)
    call this%GLL_to_GL%map(work, coef%dof%z%x(1,1,1,e), 1, this%Xh_GL)
    call dealias_local_grad(dtdx, dtdy, dtdz, work, this%Xh_GL)

    do i = 1, this%Xh_GL%lxyz
       xr = drdx(i)
       xs = drdy(i)
       xt = drdz(i)
       yr = dsdx(i)
       ys = dsdy(i)
       yt = dsdz(i)
       zr = dtdx(i)
       zs = dtdy(i)
       zt = dtdz(i)

       drdx(i) = ys*zt - yt*zs
       drdy(i) = xt*zs - xs*zt
       drdz(i) = xs*yt - xt*ys
       dsdx(i) = yt*zr - yr*zt
       dsdy(i) = xr*zt - xt*zr
       dsdz(i) = xt*yr - xr*yt
       dtdx(i) = yr*zs - ys*zr
       dtdy(i) = xs*zr - xr*zs
       dtdz(i) = xr*ys - xs*yr
    end do

  end subroutine dealias_local_build_metrics

  !> Provide the GL cofactors of one element, from the stored metrics or
  !! freshly built, see dealias_local_build_metrics.
  subroutine dealias_local_metrics(this, coef, e, drdx, drdy, drdz, &
       dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, work)
    class(adv_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e
    real(kind=rp), dimension(this%Xh_GL%lxyz), intent(inout) :: &
         drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz
    real(kind=rp), dimension(this%Xh_GL%lxyz), intent(inout) :: work

    if (this%store_metrics) then
       associate(n => this%Xh_GL%lxyz)
         call copy(drdx, this%gl_cof(:, e, 1), n)
         call copy(drdy, this%gl_cof(:, e, 2), n)
         call copy(drdz, this%gl_cof(:, e, 3), n)
         call copy(dsdx, this%gl_cof(:, e, 4), n)
         call copy(dsdy, this%gl_cof(:, e, 5), n)
         call copy(dsdz, this%gl_cof(:, e, 6), n)
         call copy(dtdx, this%gl_cof(:, e, 7), n)
         call copy(dtdy, this%gl_cof(:, e, 8), n)
         call copy(dtdz, this%gl_cof(:, e, 9), n)
       end associate
    else
       call dealias_local_build_metrics(this, coef, e, drdx, drdy, drdz, &
            dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, work)
    end if

  end subroutine dealias_local_metrics

  !> Weighted gradient \f$ w_3 J \nabla u \f$ of one element on the GL
  !! space, given its cofactors, see opgrad.
  subroutine dealias_local_opgrad(this, ux, uy, uz, u, drdx, drdy, drdz, &
       dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, w3)
    class(adv_dealias_t), intent(in) :: this
    real(kind=rp), dimension(this%Xh_GL%lxyz), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(this%Xh_GL%lxyz), intent(in) :: u
    real(kind=rp), dimension(this%Xh_GL%lxyz), intent(in) :: &
         drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz
    real(kind=rp), dimension(this%Xh_GL%lxyz), intent(in) :: w3
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: ur, us, ut
    integer :: i

    call dealias_local_grad(ur, us, ut, u, this%Xh_GL)
    do i = 1, this%Xh_GL%lxyz
       ux(i) = w3(i) * (drdx(i)*ur(i) + dsdx(i)*us(i) + dtdx(i)*ut(i))
       uy(i) = w3(i) * (drdy(i)*ur(i) + dsdy(i)*us(i) + dtdy(i)*ut(i))
       uz(i) = w3(i) * (drdz(i)*ur(i) + dsdz(i)*us(i) + dtdz(i)*ut(i))
    end do

  end subroutine dealias_local_opgrad

  !> Weighted contravariant convecting velocity of one element on the GL
  !! space, \f$ c_r = w_3 (r_x u + r_y v + r_z w) \f$ and likewise for
  !! \f$ s, t \f$, given the cofactors.
  subroutine dealias_local_convecting(cr, cs, ct, u, v, w, drdx, drdy, drdz, &
       dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, w3, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(n), intent(in) :: u, v, w
    real(kind=rp), dimension(n), intent(in) :: &
         drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz
    real(kind=rp), dimension(n), intent(in) :: w3
    integer :: i

    do i = 1, n
       cr(i) = w3(i) * (drdx(i)*u(i) + drdy(i)*v(i) + drdz(i)*w(i))
       cs(i) = w3(i) * (dsdx(i)*u(i) + dsdy(i)*v(i) + dsdz(i)*w(i))
       ct(i) = w3(i) * (dtdx(i)*u(i) + dtdy(i)*v(i) + dtdz(i)*w(i))
    end do

  end subroutine dealias_local_convecting

  !
  ! ---------------------------------------------------------------------
  ! Advection operators
  ! ---------------------------------------------------------------------
  !

  !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
  !! the RHS.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_advection_dealias(this, vx, vy, vz, fx, fy, fz, Xh, &
       coef, n, dt)
    class(adv_dealias_t), intent(inout) :: this
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: fx, fy, fz
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tx, ty, tz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: ur, us, ut, fg
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: cr, cs, ct
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: drdx, drdy, drdz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dsdx, dsdy, dsdz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dtdx, dtdy, dtdz
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: tempx, tempy, tempz
    type(dealias_vec_t) :: c(3)
    integer :: e, i, nel, e0, ne

    nel = coef%msh%nelv
    call profiler_start_region('Dealias advection')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do e0 = 1, nel, this%chunk
          ne = min(this%chunk, nel - e0 + 1)
          call dealias_device_convecting(this, coef, vx%x_d, vy%x_d, vz%x_d, &
               e0, ne, c)
          call dealias_device_convect(this, vx%x_d, fx%x_d, c, e0, ne)
          call dealias_device_convect(this, vy%x_d, fy%x_d, c, e0, ne)
          call dealias_device_convect(this, vz%x_d, fz%x_d, c, e0, ne)
          call dealias_relinquish(c)
       end do
    else
       associate(w3 => this%Xh_GL%w3, lxyz_GL => this%Xh_GL%lxyz)
         !$omp parallel do private(e, i, tempx, tempy, tempz), &
         !$omp& private(tx, ty, tz, ur, us, ut, fg, cr, cs, ct), &
         !$omp& private(drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz)
         do e = 1, nel
            call dealias_local_metrics(this, coef, e, drdx, drdy, drdz, &
                 dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, fg)

            call this%GLL_to_GL%map(tx, vx%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(ty, vy%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(tz, vz%x(1,1,1,e), 1, this%Xh_GL)

            ! Weighted contravariant convecting velocity
            call dealias_local_convecting(cr, cs, ct, tx, ty, tz, &
                 drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
                 w3, lxyz_GL)

            call dealias_local_grad(ur, us, ut, tx, this%Xh_GL)
            do i = 1, lxyz_GL
               fg(i) = cr(i)*ur(i) + cs(i)*us(i) + ct(i)*ut(i)
            end do
            call this%GLL_to_GL%map(tempx, fg, 1, this%Xh_GLL)

            call dealias_local_grad(ur, us, ut, ty, this%Xh_GL)
            do i = 1, lxyz_GL
               fg(i) = cr(i)*ur(i) + cs(i)*us(i) + ct(i)*ut(i)
            end do
            call this%GLL_to_GL%map(tempy, fg, 1, this%Xh_GLL)

            call dealias_local_grad(ur, us, ut, tz, this%Xh_GL)
            do i = 1, lxyz_GL
               fg(i) = cr(i)*ur(i) + cs(i)*us(i) + ct(i)*ut(i)
            end do
            call this%GLL_to_GL%map(tempz, fg, 1, this%Xh_GLL)

            call sub2(fx%x(1,1,1,e), tempx, this%Xh_GLL%lxyz)
            call sub2(fy%x(1,1,1,e), tempy, this%Xh_GLL%lxyz)
            call sub2(fz%x(1,1,1,e), tempz, this%Xh_GLL%lxyz)
         end do
         !$omp end parallel do
       end associate
    end if

    call profiler_end_region('Dealias advection')

  end subroutine compute_advection_dealias

  !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to the
  !! RHS.
  !! @param this The object.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param s The scalar.
  !! @param fs The source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_scalar_advection_dealias(this, vx, vy, vz, s, fs, Xh, &
       coef, n, dt)
    class(adv_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: fs
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vx_GL, vy_GL, vz_GL, s_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dsdr, dsds, dsdt, f_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: cr, cs, ct
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: drdx, drdy, drdz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dsdx, dsdy, dsdz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dtdx, dtdy, dtdz
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: temp
    type(dealias_vec_t) :: c(3)
    integer :: e, i, nel, e0, ne

    nel = coef%msh%nelv
    call profiler_start_region('Dealias scalar advection')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do e0 = 1, nel, this%chunk
          ne = min(this%chunk, nel - e0 + 1)
          call dealias_device_convecting(this, coef, vx%x_d, vy%x_d, vz%x_d, &
               e0, ne, c)
          call dealias_device_convect(this, s%x_d, fs%x_d, c, e0, ne)
          call dealias_relinquish(c)
       end do
    else
       associate(w3 => this%Xh_GL%w3, lxyz_GL => this%Xh_GL%lxyz)
         !$omp parallel do private(e, i, vx_GL, vy_GL, vz_GL, s_GL), &
         !$omp& private(f_GL, temp, dsdr, dsds, dsdt, cr, cs, ct), &
         !$omp& private(drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz)
         do e = 1, nel
            call dealias_local_metrics(this, coef, e, drdx, drdy, drdz, &
                 dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, f_GL)

            ! Map advecting velocity and the scalar onto the higher-order
            ! space
            call this%GLL_to_GL%map(vx_GL, vx%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(vy_GL, vy%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(vz_GL, vz%x(1,1,1,e), 1, this%Xh_GL)
            call this%GLL_to_GL%map(s_GL, s%x(1,1,1,e), 1, this%Xh_GL)

            call dealias_local_convecting(cr, cs, ct, vx_GL, vy_GL, vz_GL, &
                 drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
                 w3, lxyz_GL)

            ! Reference-space gradient of s and the convective term
            call dealias_local_grad(dsdr, dsds, dsdt, s_GL, this%Xh_GL)
            do i = 1, lxyz_GL
               f_GL(i) = cr(i)*dsdr(i) + cs(i)*dsds(i) + ct(i)*dsdt(i)
            end do

            ! Map back the contructed operator to the original space
            call this%GLL_to_GL%map(temp, f_GL, 1, this%Xh_GLL)

            call sub2(fs%x(1,1,1,e), temp, this%Xh_GLL%lxyz)
         end do
         !$omp end parallel do
       end associate
    end if

    call profiler_end_region('Dealias scalar advection')

  end subroutine compute_scalar_advection_dealias


  !!> Add the advection term in ALE framework using dealiasing.
  !! @param this The object.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param wm_x The x component of mesh velocity.
  !! @param wm_y The y component of mesh velocity.
  !! @param wm_z The z component of mesh velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  !! Here, we compute: - div ( u_i * wm ).
  !! Based on Ho, L.W. A Legendre spectral element method for
  !! simulation of incompressible unsteady viscous free-surface flows.
  !! Ph.D. thesis, Massachusetts Institute of Technology, 1989.
  !! Note: In Nek5000, dealiasing is not done for this term.
  subroutine compute_ale_advection_dealias(this, vx, vy, vz, &
       wm_x, wm_y, wm_z, fx, fy, fz, Xh, coef, n, dt)
    class(adv_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: wm_x, wm_y, wm_z
    type(field_t), intent(inout) :: fx, fy, fz
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vx_GL, vy_GL, vz_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: wm_x_GL, wm_y_GL, wm_z_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: flux_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: grad_x, grad_y, grad_z
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: total_div_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: drdx, drdy, drdz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dsdx, dsdy, dsdz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dtdx, dtdy, dtdz
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: temp_x, temp_y, temp_z
    type(dealias_vec_t) :: cof(9), wm(3), ug(1), tmp(1), grad(3), acc(1), &
         tg(1)
    type(c_ptr) :: cof_d(9), src_d(3), dst_d(3), u_d(3), f_d(3)
    integer(kind=i8) :: off
    integer :: e, i, nel, e0, ne, n_gl, n_gll, comp

    nel = coef%msh%nelv
    call profiler_start_region('Dealias ALE advection')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       u_d(1) = vx%x_d
       u_d(2) = vy%x_d
       u_d(3) = vz%x_d
       f_d(1) = fx%x_d
       f_d(2) = fy%x_d
       f_d(3) = fz%x_d
       do e0 = 1, nel, this%chunk
          ne = min(this%chunk, nel - e0 + 1)
          n_gl = ne * this%Xh_GL%lxyz
          n_gll = ne * this%Xh_GLL%lxyz
          off = int(e0 - 1, i8) * int(this%Xh_GLL%lxyz, i8)

          call dealias_device_metrics(this, coef, e0, ne, cof, cof_d)

          ! Mesh velocity on the GL space (work vectors sized for a full
          ! chunk, see dealias_device_build_metrics)
          call dealias_request(wm, this%chunk * this%Xh_GL%lxyz)
          src_d(1) = wm_x%x_d
          src_d(2) = wm_y%x_d
          src_d(3) = wm_z%x_d
          do i = 1, 3
             dst_d(i) = wm(i)%v%x_d
          end do
          call dealias_device_to_gl(this, src_d, dst_d, e0, ne)

          call dealias_request(ug, this%chunk * this%Xh_GL%lxyz)
          call dealias_request(tmp, this%chunk * this%Xh_GL%lxyz)
          call dealias_request(grad, this%chunk * this%Xh_GL%lxyz)
          call dealias_request(acc, this%chunk * this%Xh_GL%lxyz)
          call dealias_request(tg, this%chunk * this%Xh_GLL%lxyz)

          do comp = 1, 3
             src_d(1) = u_d(comp)
             dst_d(1) = ug(1)%v%x_d
             call dealias_device_to_gl(this, src_d(1:1), dst_d(1:1), e0, ne)

             ! div(u_i wm) = d/dx (u_i wm_x) + d/dy (u_i wm_y)
             ! + d/dz (u_i wm_z), each term as the matching component of the
             ! weighted gradient of the flux
             call device_col3(tmp(1)%v%x_d, ug(1)%v%x_d, wm(1)%v%x_d, n_gl)
             call dealias_device_opgrad(this, acc(1)%v%x_d, grad(2)%v%x_d, &
                  grad(3)%v%x_d, tmp(1)%v%x_d, cof_d, ne)

             call device_col3(tmp(1)%v%x_d, ug(1)%v%x_d, wm(2)%v%x_d, n_gl)
             call dealias_device_opgrad(this, grad(1)%v%x_d, grad(2)%v%x_d, &
                  grad(3)%v%x_d, tmp(1)%v%x_d, cof_d, ne)
             call device_add2(acc(1)%v%x_d, grad(2)%v%x_d, n_gl)

             call device_col3(tmp(1)%v%x_d, ug(1)%v%x_d, wm(3)%v%x_d, n_gl)
             call dealias_device_opgrad(this, grad(1)%v%x_d, grad(2)%v%x_d, &
                  grad(3)%v%x_d, tmp(1)%v%x_d, cof_d, ne)
             call device_add2(acc(1)%v%x_d, grad(3)%v%x_d, n_gl)

             ! Map the divergence back to the GLL space and add to the RHS
             call tnsr3d_device(tg(1)%v%x_d, this%Xh_GLL%lx, acc(1)%v%x_d, &
                  this%Xh_GL%lx, this%GLL_to_GL%Yh_XhT_d, &
                  this%GLL_to_GL%Yh_Xh_d, this%GLL_to_GL%Yh_Xh_d, ne)
             call device_add2(dev_ptr_offset(f_d(comp), off), tg(1)%v%x_d, &
                  n_gll)
          end do

          call dealias_relinquish(tg)
          call dealias_relinquish(acc)
          call dealias_relinquish(grad)
          call dealias_relinquish(tmp)
          call dealias_relinquish(ug)
          call dealias_relinquish(wm)
          call dealias_relinquish(cof)
       end do
    else
       !$omp parallel do private (e, i, vx_GL, vy_GL, vz_GL), &
       !$omp& private (wm_x_GL, wm_y_GL, wm_z_GL, flux_GL, total_div_GL), &
       !$omp& private (grad_x, grad_y, grad_z, temp_x, temp_y, temp_z), &
       !$omp& private (drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz)
       do e = 1, nel
          call dealias_local_metrics(this, coef, e, drdx, drdy, drdz, &
               dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, flux_GL)

          ! Map advecting velocity and mesh velocity onto the
          ! higher-order space
          call this%GLL_to_GL%map(vx_GL, vx%x(1,1,1,e), 1, this%Xh_GL)
          call this%GLL_to_GL%map(vy_GL, vy%x(1,1,1,e), 1, this%Xh_GL)
          call this%GLL_to_GL%map(vz_GL, vz%x(1,1,1,e), 1, this%Xh_GL)
          call this%GLL_to_GL%map(wm_x_GL, wm_x%x(1,1,1,e), 1, this%Xh_GL)
          call this%GLL_to_GL%map(wm_y_GL, wm_y%x(1,1,1,e), 1, this%Xh_GL)
          call this%GLL_to_GL%map(wm_z_GL, wm_z%x(1,1,1,e), 1, this%Xh_GL)

          ! --------------------- X-Momentum
          ! div(u * wm_*) = d/dx (u * wm_x) + d/dy (u * wm_y) +
          ! d/dz (u * wm_z)
          flux_GL = vx_GL * wm_x_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = grad_x
          flux_GL = vx_GL * wm_y_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = total_div_GL + grad_y
          flux_GL = vx_GL * wm_z_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = total_div_GL + grad_z

          ! Map back the constructed operator to the original space
          call this%GLL_to_GL%map(temp_x, total_div_GL, 1, this%Xh_GLL)

          ! --------------------- Y-Momentum
          flux_GL = vy_GL * wm_x_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = grad_x
          flux_GL = vy_GL * wm_y_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = total_div_GL + grad_y
          flux_GL = vy_GL * wm_z_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = total_div_GL + grad_z

          call this%GLL_to_GL%map(temp_y, total_div_GL, 1, this%Xh_GLL)

          ! --------------------- Z-Momentum
          flux_GL = vz_GL * wm_x_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = grad_x
          flux_GL = vz_GL * wm_y_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = total_div_GL + grad_y
          flux_GL = vz_GL * wm_z_GL
          call dealias_local_opgrad(this, grad_x, grad_y, grad_z, flux_GL, &
               drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
               this%Xh_GL%w3)
          total_div_GL = total_div_GL + grad_z

          call this%GLL_to_GL%map(temp_z, total_div_GL, 1, this%Xh_GLL)

          ! Note we add (+) here since the ALE advection term is
          ! - div(u * wm) on the LHS. So on the RHS it will be + div(u * wm)
          call add2(fx%x(1,1,1,e), temp_x, this%Xh_GLL%lxyz)
          call add2(fy%x(1,1,1,e), temp_y, this%Xh_GLL%lxyz)
          call add2(fz%x(1,1,1,e), temp_z, this%Xh_GLL%lxyz)
       end do
       !$omp end parallel do
    end if

    call profiler_end_region('Dealias ALE advection')

  end subroutine compute_ale_advection_dealias

  !> Refresh the stored GL metrics after the mesh has moved. With metrics
  !! rebuilt per call there is nothing to do, as they derive from the current
  !! coordinates (or GLL cofactors) every time.
  !! @param coef The updated GLL coefficients.
  !! @param moving_boundary Whether the mesh actually moved.
  subroutine recompute_metrics_dealias(this, coef, moving_boundary)
    class(adv_dealias_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    logical, intent(in) :: moving_boundary

    if (.not. moving_boundary) return
    if (.not. this%store_metrics) return

    call dealias_build_stored_metrics(this, coef)

  end subroutine recompute_metrics_dealias

end module adv_dealias
