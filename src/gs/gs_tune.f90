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
!> Runtime autotuning of the gather-scatter communication
submodule (gather_scatter) gs_tune
  use gs_shmem, only : GS_SHMEM_AVAIL
  use gs_caf, only : GS_CAF_AVAIL, gs_caf_usable, gs_caf_signal_auto, &
       gs_caf_signal_modes, gs_caf_mode_name, gs_caf_mode_get, &
       gs_caf_set_mode
  use gs_utofu, only : GS_UTOFU_AVAIL
  use gs_mpi_rma, only : GS_MPI_RMA_AVAIL
  use gs_device_nccl, only : GS_DEVICE_NCCL_AVAIL
  use gs_device_shmem, only : GS_DEVICE_SHMEM_AVAIL
  use neko_config, only : NEKO_DEVICE_MPI
  implicit none

  !> The comm. backends the runtime autotuning can benchmark, in the order
  !! they are benchmarked and reported: first the host backends, which on a
  !! device build stage the halo through the host mirror of the shared
  !! buffer, then the device-resident ones, which keep it on the device and
  !! are candidates on CUDA and HIP builds only (see gs_comm_tunable).
  integer, parameter :: GS_TUNE_BCKND(11) = [GS_COMM_MPI, GS_COMM_NEIGHBOUR, &
       GS_COMM_MPIRMA, GS_COMM_OPENSHMEM, GS_COMM_CAF, GS_COMM_UTOFU, &
       GS_COMM_CRYSTAL, GS_COMM_MPIGPU, GS_COMM_NCCL, GS_COMM_CRYSTALGPU, &
       GS_COMM_NVSHMEM]

  !> Which of GS_TUNE_BCKND are benchmarked when NEKO_GS_TUNE names no
  !! candidates of its own: all of them this build and run supports but CAF
  !! and NVSHMEM, which have to be asked for (see gs_tune_select). A host
  !! build drops every device backend in gs_comm_tunable, so the one set
  !! serves both.
  !!
  !! @note This is the set a comparison starts from, not a statement that a
  !! comparison happens. On a build with device-aware MPI none does unless
  !! NEKO_GS_TUNE is set at all (see gs_tune_requested and gs_init); the
  !! device backends are in the set so that device MPI, the default there,
  !! defends its place once a comparison is asked for.
  logical, parameter :: GS_TUNE_DEFAULT(11) = [.true., .true., .true., &
       .true., .false., .true., .true., .true., .true., .true., .false.]

  !> Whether this build has what the device-resident comm. backends need to
  !! pack and unpack the halo in device memory. Their kernels are written
  !! for CUDA and HIP only and abort on the OpenCL and Metal backends, so
  !! those builds tune over the host backends alone, staging the halo
  !! through the host mirror of the shared buffer as they always have.
  logical, parameter :: GS_TUNE_DEV_AVAIL = (NEKO_BCKND_HIP .eq. 1) .or. &
       (NEKO_BCKND_CUDA .eq. 1)

  !> Whether the coarray signaling mode has already been selected by the
  !! autotuner (NEKO_GS_CAF_SIGNALING=auto). The mode is a program-wide
  !! binding shared by every gs_caf_t instance, so it is benchmarked once,
  !! on the first gs_t that tunes the CAF backend, and kept: re-binding it
  !! later would change the protocol under instances already running.
  logical, save :: caf_signal_tuned = .false.

contains

  !> Whether the autotuning has been asked for explicitly, that is whether
  !! NEKO_GS_TUNE is set to anything at all. What it is set to selects the
  !! candidates (see gs_tune_select); merely being set is what turns the
  !! benchmark on where benchmarking is not the default, namely on a build
  !! with device-aware MPI (see gs_init).
  !! @return whether NEKO_GS_TUNE is set
  module function gs_tune_requested() result(requested)
    logical :: requested
    character(len=255) :: env_val
    integer :: env_len

    call get_environment_variable("NEKO_GS_TUNE", env_val, env_len)
    requested = (env_len .gt. 0)

  end function gs_tune_requested

  !> Comm. backends to benchmark in the runtime autotuning: the ones
  !! selected by NEKO_GS_TUNE (see gs_tune_select) that this build and run
  !! can actually use (see gs_comm_tunable), in GS_TUNE_BCKND order. A
  !! selected backend the run cannot use is dropped, and reported in the log
  !! if it was named explicitly rather than merely part of the default set.
  !! @return the candidate backends, as GS_COMM_* constants. May be empty,
  !! or hold a single entry, in which case gs_tune_comm benchmarks nothing
  !! @note The list says nothing about which backend the gs kernel holds
  !! when the tuning starts; gs_tune_comm is told that separately.
  function gs_comm_cand() result(cand)
    integer, allocatable :: cand(:)
    character(len=LOG_SIZE) :: log_buf
    character(len=13) :: label
    logical :: sel(size(GS_TUNE_BCKND)), named(size(GS_TUNE_BCKND))
    integer :: c(size(GS_TUNE_BCKND)), i, n

    call gs_tune_select(sel, named)

    n = 0
    do i = 1, size(GS_TUNE_BCKND)
       if (.not. sel(i)) cycle
       if (gs_comm_tunable(GS_TUNE_BCKND(i))) then
          n = n + 1
          c(n) = GS_TUNE_BCKND(i)
       else if (named(i)) then
          ! Asked for by name, but this build or run cannot drive it: say so
          ! rather than leaving the candidate quietly missing from the
          ! comparison below
          label = adjustl(gs_comm_name(GS_TUNE_BCKND(i)))
          write(log_buf, '(A,A,A12)') label, ': ', 'unavailable'
          call neko_log%message(log_buf)
       end if
    end do

    allocate(cand(n))
    cand = c(1:n)

  end function gs_comm_cand

  !> Whether the comm. backend @a comm_bcknd can be driven by this build and
  !! this run, and may therefore be benchmarked. MPI and the neighbourhood
  !! collective need nothing but MPI-3; the rest abort in their init, or
  !! exchange nothing, unless the corresponding support was built in. The
  !! device-resident backends additionally need a build whose gather-scatter
  !! kernels can pack the halo in device memory (GS_TUNE_DEV_AVAIL).
  !! @param comm_bcknd comm. backend to check, one of the GS_TUNE_BCKND
  !! entries
  !! @return whether the backend is usable as a candidate
  function gs_comm_tunable(comm_bcknd) result(tunable)
    integer, intent(in) :: comm_bcknd
    logical :: tunable

    select case (comm_bcknd)
    case (GS_COMM_MPI, GS_COMM_NEIGHBOUR)
       tunable = .true.
    case (GS_COMM_CRYSTAL)
       ! Nothing beyond MPI-1, and it addresses its partners by their rank
       ! in NEKO_COMM, so a communicator split does not rule it out either.
       ! It aggregates the halo into one message per routing stage instead
       ! of one per peer, which is a win only where per-message overhead
       ! dominates; the benchmark is what decides that
       tunable = .true.
    case (GS_COMM_MPIRMA)
       ! Nothing beyond MPI-3, and it addresses its peers by their rank in
       ! NEKO_COMM rather than by a global PE or image number, so unlike the
       ! other one-sided backends a communicator split (NEKO_COMM_ID) does
       ! not rule it out. What it does assume is that the implementation
       ! makes RMA progress without the target entering MPI; on one that
       ! does not (Open MPI osc/pt2pt) it is slow rather than wrong, so it
       ! loses the benchmark and is dropped, which is the outcome the
       ! autotuning exists to produce. Rule it out with NEKO_GS_TUNE=-MPIRMA
       ! if paying for that measurement is not worth it.
       tunable = GS_MPI_RMA_AVAIL
    case (GS_COMM_OPENSHMEM)
       ! The one-sided backends that address their peers by global PE or
       ! image number (OpenSHMEM PEs, coarray images) are only correct when
       ! NEKO_COMM spans every process, so skip them when the run has been
       ! split into several communicators (NEKO_COMM_ID). uTofu exchanges
       ! its addresses over NEKO_COMM itself and is unaffected.
       tunable = GS_SHMEM_AVAIL .and. (pe_size .eq. global_pe_size)
    case (GS_COMM_CAF)
       ! Coarray support at configure time says nothing about the job
       ! actually running more than one image, see gs_caf_usable
       tunable = GS_CAF_AVAIL .and. gs_caf_usable() .and. &
            (pe_size .eq. global_pe_size)
    case (GS_COMM_UTOFU)
       tunable = GS_UTOFU_AVAIL
    case (GS_COMM_MPIGPU, GS_COMM_CRYSTALGPU)
       ! Both hand device pointers straight to MPI, which only holds on a
       ! build configured for device-aware MPI (--enable-device-mpi). They
       ! address their peers by their rank in NEKO_COMM, so a communicator
       ! split does not rule them out.
       tunable = GS_TUNE_DEV_AVAIL .and. NEKO_DEVICE_MPI
    case (GS_COMM_NCCL)
       ! The NCCL communicator is built once at startup
       ! (neko_comm_nccl_init) over the ranks of NEKO_COMM, but from a
       ! unique id broadcast over MPI_COMM_WORLD, so it is only well defined
       ! when NEKO_COMM spans every process, as for the PE addressed
       ! backends above
       tunable = GS_TUNE_DEV_AVAIL .and. GS_DEVICE_NCCL_AVAIL .and. &
            (pe_size .eq. global_pe_size)
    case (GS_COMM_NVSHMEM)
       ! Addresses its peers by NVSHMEM PE, as OpenSHMEM does, so the same
       ! restriction on communicator splits applies
       tunable = GS_TUNE_DEV_AVAIL .and. GS_DEVICE_SHMEM_AVAIL .and. &
            (pe_size .eq. global_pe_size)
    case default
       tunable = .false.
    end select

  end function gs_comm_tunable

  !> Which of GS_TUNE_BCKND the autotuning should consider, as requested by
  !! the NEKO_GS_TUNE environment variable. Its value is a list of backend
  !! names (comma or space separated, case insensitive) spelled as for
  !! NEKO_GS_COMM, in either of two forms:
  !!
  !! - plain names, which replace the default set outright
  !!   (`NEKO_GS_TUNE=MPI,NEIGHBOUR`), or
  !! - names prefixed with `+` or `-`, which add to and remove from it
  !!   (`NEKO_GS_TUNE=+CAF`).
  !!
  !! The two forms cannot be mixed. When the variable is unset the default
  !! set GS_TUNE_DEFAULT is used, which is everything but CAF and NVSHMEM.
  !! A compiler may accept coarrays and still give the job a single image,
  !! and a backend that cannot communicate is at best wasted setup, so CAF
  !! has to be asked for; NVSHMEM aborts outright unless the peer lists come
  !! out symmetric and aligned, which is too much to risk on a run that
  !! never asked for it. Selecting a backend outright (NEKO_GS_COMM)
  !! bypasses all of this.
  !!
  !! @param sel whether each entry of GS_TUNE_BCKND was selected
  !! @param named whether each entry was named explicitly, i.e. is a request
  !! rather than a default, which is what makes it worth reporting when the
  !! run cannot use it
  subroutine gs_tune_select(sel, named)
    logical, intent(out) :: sel(:), named(:)
    character(len=255) :: env_val
    character(len=32) :: tok, name
    character(len=1) :: op
    integer :: env_len, i, j, k
    logical :: delta, first

    sel = GS_TUNE_DEFAULT
    named = .false.

    call get_environment_variable("NEKO_GS_TUNE", env_val, env_len)
    if (env_len .eq. 0) return
    ! env_len is the length of the value, which may well exceed the buffer
    env_len = min(env_len, len(env_val))

    delta = .false.
    first = .true.
    i = 1
    do while (i .le. env_len)
       if (scan(env_val(i:i), ', ') .ne. 0) then
          i = i + 1
          cycle
       end if

       j = i
       do while (j .le. env_len)
          if (scan(env_val(j:j), ', ') .ne. 0) exit
          j = j + 1
       end do
       tok = gs_tune_upcase(env_val(i:j-1))
       i = j

       op = ' '
       name = tok
       if (tok(1:1) .eq. '+' .or. tok(1:1) .eq. '-') then
          op = tok(1:1)
          name = tok(2:)
       end if

       if (first) then
          ! The first name decides how the list is read; a plain one starts
          ! the set from nothing rather than from the default
          delta = (op .ne. ' ')
          if (.not. delta) sel = .false.
          first = .false.
       else if ((op .ne. ' ') .neqv. delta) then
          call neko_error('NEKO_GS_TUNE: plain backend names and +/- ' // &
               'prefixed ones cannot be mixed')
       end if

       k = gs_tune_index(name)
       if (k .eq. 0) then
          call neko_error('NEKO_GS_TUNE: not a tunable Gather-scatter ' // &
               'comm. backend: ' // trim(name))
       end if
       sel(k) = (op .ne. '-')
       named(k) = .true.
    end do

  end subroutine gs_tune_select

  !> Index into GS_TUNE_BCKND of the comm. backend named @a name, using the
  !! spellings NEKO_GS_COMM accepts, or 0 if it names no backend that takes
  !! part in the autotuning. Host and device OpenSHMEM are one name to
  !! NEKO_GS_COMM, which picks between them by build; here they can both be
  !! candidates of the same run, so they are spelled apart: SHMEM is the
  !! host backend and NVSHMEM the device one.
  !! @param name upper case backend name to look up
  !! @return the index into GS_TUNE_BCKND, or 0
  function gs_tune_index(name) result(idx)
    character(len=*), intent(in) :: name
    integer :: idx, bcknd, i

    select case (trim(name))
    case ('MPI')
       bcknd = GS_COMM_MPI
    case ('NEIGHBOUR', 'NEIGHBOR')
       bcknd = GS_COMM_NEIGHBOUR
    case ('SHMEM', 'OPENSHMEM')
       bcknd = GS_COMM_OPENSHMEM
    case ('CAF')
       bcknd = GS_COMM_CAF
    case ('UTOFU')
       bcknd = GS_COMM_UTOFU
    case ('MPIRMA', 'RMA')
       bcknd = GS_COMM_MPIRMA
    case ('CRYSTAL')
       bcknd = GS_COMM_CRYSTAL
    case ('MPIGPU')
       bcknd = GS_COMM_MPIGPU
    case ('NCCL', 'RCCL')
       bcknd = GS_COMM_NCCL
    case ('NVSHMEM')
       bcknd = GS_COMM_NVSHMEM
    case ('CRYSTALGPU')
       bcknd = GS_COMM_CRYSTALGPU
    case default
       bcknd = 0
    end select

    idx = 0
    do i = 1, size(GS_TUNE_BCKND)
       if (GS_TUNE_BCKND(i) .eq. bcknd) idx = i
    end do

  end function gs_tune_index

  !> Upper case copy of @a str, so that the backend names in NEKO_GS_TUNE
  !! may be given in any case
  !! @param str string to convert
  !! @return the converted string
  function gs_tune_upcase(str) result(upper)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: upper
    integer :: i, c

    do i = 1, len(str)
       c = iachar(str(i:i))
       if (c .ge. iachar('a') .and. c .le. iachar('z')) then
          upper(i:i) = achar(c - (iachar('a') - iachar('A')))
       else
          upper(i:i) = str(i:i)
       end if
    end do

  end function gs_tune_upcase

  !> Switch the comm. backend of @a gs to @a comm_bcknd, retaining the
  !! schedule computed in gs_schedule (the dof lists are moved over to the
  !! new backend, not recomputed). The old backend is freed before the new
  !! one is set up, so backends holding scarce resources (uTofu VCQs,
  !! symmetric memory, coarrays) never overlap.
  !!
  !! The gather-scatter backend is told where to leave the gathered shared
  !! dofs to match: in the host mirror of the shared buffer for a host
  !! backend, in device memory for a device-resident one.
  !! @param gs gather-scatter kernel whose comm. backend is replaced
  !! @param comm_bcknd comm. backend to switch to, one of the GS_COMM_*
  !! constants
  !! @note Collective, every rank must switch to the same backend.
  subroutine gs_comm_switch(gs, comm_bcknd)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: comm_bcknd
    class(gs_comm_t), allocatable :: comm_new

    call gs_comm_alloc(comm_new, comm_bcknd)
    call comm_new%take_schedule(gs%comm)
    call gs%comm%free()
    deallocate(gs%comm)
    call move_alloc(comm_new, gs%comm)
    call gs%comm%init_schedule()

    gs%bcknd%shared_on_host = .not. gs_comm_on_device(comm_bcknd)

  end subroutine gs_comm_switch

  !> Time @a ntrials gather-scatter operations on @a u, returning the
  !! average wall time (s) per operation. The backend is warmed up and all
  !! ranks are synchronised before the timing window is opened.
  !! @param gs gather-scatter kernel to time
  !! @param u working vector to operate on
  !! @param n length of @a u
  !! @param op gather-scatter operation, one of the GS_OP_* constants
  !! @param ntrials number of timed operations to average over
  !! @return the average wall time (s) per operation on this rank
  module function gs_time_ops(gs, u, n, op, ntrials) result(t)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer, intent(in) :: op, ntrials
    real(kind=dp) :: t
    integer :: i

    do i = 1, GS_TUNE_NWARMUP
       call gs_op_vector(gs, u, n, op)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) call device_sync
    call MPI_Barrier(NEKO_COMM)

    t = MPI_Wtime()
    do i = 1, ntrials
       call gs_op_vector(gs, u, n, op)
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) call device_sync
    t = (MPI_Wtime() - t) / real(ntrials, dp)

  end function gs_time_ops

  !> Select the fastest of the candidate comm. backends at runtime.
  !!
  !! Each candidate returned by gs_comm_cand is benchmarked in turn on
  !! a dummy vector of length @a n, reusing the schedule computed by
  !! gs_schedule (handed over from one candidate to the next, see
  !! gs_comm_t%take_schedule), and @a gs is left with the fastest backend
  !! allocated. Fewer than two candidates leaves nothing to compare: @a gs
  !! is switched to the single candidate, if it is not already running it,
  !! and no benchmark is run.
  !!
  !! @param gs gather-scatter kernel to tune, left holding the fastest
  !! backend
  !! @param n length of the dummy vector to benchmark on (the dofmap size)
  !! @param comm_bcknd the comm. backend @a gs holds on entry, i.e. the one
  !! the schedule was built with, which need not be a candidate itself
  !!
  !! A candidate with a sub-choice of its own has it benchmarked here too,
  !! and contributes the winning variant's time: the device MPI
  !! synchronisation strategy, and the coarray signaling mode.
  !!
  !! @note Setting up and running a comm. backend is collective over
  !! NEKO_COMM (the neighbourhood backend builds a graph communicator, the
  !! one-sided backends allocate symmetric memory), so every rank has to
  !! arrive at the same decision: the per-rank timings are averaged over all
  !! ranks before the winner is picked.
  module subroutine gs_tune_comm(gs, n, comm_bcknd)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    integer, intent(in) :: comm_bcknd
    character(len=LOG_SIZE) :: log_buf
    character(len=13) :: label
    integer, allocatable :: cand(:)
    real(kind=dp), allocatable :: cand_time(:)
    real(kind=rp), allocatable :: tmp(:)
    type(c_ptr) :: tmp_d
    integer :: i, best, nmin, cur, dev_strtgy

    allocate(cand, source = gs_comm_cand())

    ! Track what gs is running as the sweep proceeds. The candidate list
    ! need neither start with the backend the schedule was built with nor
    ! contain it at all, and switching is what costs, so every candidate is
    ! switched to exactly once
    cur = comm_bcknd

    ! A rank without any dofs skips the halo exchange in gs_op_vector, which
    ! would leave every other rank hanging in a collective based backend.
    ! Keep the current backend in that case. NEKO_GS_TUNE may also have left
    ! us with nothing to compare, in which case a lone candidate is simply
    ! switched to, as if it had been requested with NEKO_GS_COMM.
    nmin = n
    call MPI_Allreduce(MPI_IN_PLACE, nmin, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM)
    if (nmin .eq. 0 .or. size(cand) .lt. 2) then
       call neko_log%message('Comm tuning  :      skipped')
       if (nmin .gt. 0 .and. size(cand) .eq. 1) then
          if (cand(1) .ne. cur) call gs_comm_switch(gs, cand(1))
          cur = cand(1)
          ! Nothing to compare the backend against, but device MPI still has
          ! a strategy to pick
          if (cur .eq. GS_COMM_MPIGPU) call gs_tune_dev_strtgy(gs, n)
       end if
       call neko_log%message('Tuned comm   : ' // gs_comm_name(cur))
       return
    end if

    allocate(cand_time(size(cand)))

    tmp_d = C_NULL_PTR
    allocate(tmp(n))
    tmp = 1.0_rp
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(tmp, tmp_d, n)
       call device_memcpy(tmp, tmp_d, n, HOST_TO_DEVICE, sync = .false.)
    end if

    ! GS_OP_MIN leaves the working vector untouched, so the benchmark can
    ! run for any number of trials without the values growing out of range
    dev_strtgy = -1
    do i = 1, size(cand)
       if (cand(i) .eq. GS_COMM_CAF .and. gs_caf_signal_auto() .and. &
            .not. caf_signal_tuned) then
          ! Benchmark the signaling modes as well; this leaves the CAF
          ! backend allocated in the winning mode
          cand_time(i) = gs_tune_caf_signal(gs, tmp, n)
          caf_signal_tuned = .true.
       else
          if (cand(i) .ne. cur) call gs_comm_switch(gs, cand(i))
          if (cand(i) .eq. GS_COMM_MPIGPU) then
             ! Benchmark the synchronisation strategies as well; this leaves
             ! the fastest one bound, and it has to be remembered in case the
             ! sweep switches away from device MPI and back again
             call gs_tune_strtgy(gs, tmp, n, cand_time(i))
             dev_strtgy = gs_get_strtgy(gs)
          else
             cand_time(i) = gs_time_ops(gs, tmp, n, GS_OP_MIN, &
                  GS_TUNE_NTRIALS)
          end if
       end if
       cur = cand(i)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) call device_unmap(tmp, tmp_d)
    deallocate(tmp)

    call MPI_Allreduce(MPI_IN_PLACE, cand_time, size(cand), &
         MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM)
    cand_time = cand_time / pe_size

    best = minloc(cand_time, 1)

    do i = 1, size(cand)
       label = adjustl(gs_comm_name(cand(i)))
       ! What the coarray backend costs depends on the signaling mode in
       ! force, so report the timing under the mode it was measured in
       if (cand(i) .eq. GS_COMM_CAF .and. gs_caf_mode_get() .ne. 0) then
          label = 'CAF (' // &
               trim(adjustl(gs_caf_mode_name(gs_caf_mode_get()))) // ')'
       end if
       ! ES10.3 + ' s' fills the same 12 columns as the right-adjusted
       ! backend names, so the unit lines up with their last letter
       write(log_buf, '(A,A,ES10.3,A)') label, ': ', cand_time(i), ' s'
       call neko_log%message(log_buf)
    end do

    if (cand(best) .ne. cur) call gs_comm_switch(gs, cand(best))

    ! A switch back to device MPI allocates a fresh backend, which comes up
    ! with the default strategy rather than the one benchmarked above
    if (dev_strtgy .ge. 0) call gs_set_strtgy(gs, dev_strtgy)

    call neko_log%message('Tuned comm   : ' // gs_comm_name(cand(best)))

    deallocate(cand_time)

  end subroutine gs_tune_comm

  !> Bind the non-blocking synchronisation strategy of the device MPI comm.
  !! backend held by @a gs, benchmarking the strategies on a dummy vector.
  !!
  !! Called by gs_init when device MPI was requested outright; when the
  !! comm. backend is left to the autotuning, gs_tune_comm benchmarks the
  !! strategies on its own working vector instead (see gs_tune_strtgy).
  !!
  !! @param gs gather-scatter kernel to tune, holding a device MPI comm.
  !! backend
  !! @param n length of the dummy vector to benchmark on (the dofmap size)
  module subroutine gs_tune_dev_strtgy(gs, n)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), allocatable :: tmp(:)
    type(c_ptr) :: tmp_d

    tmp_d = C_NULL_PTR
    allocate(tmp(n))
    tmp = 1.0_rp
    call device_map(tmp, tmp_d, n)
    call device_memcpy(tmp, tmp_d, n, HOST_TO_DEVICE, sync = .false.)

    call gs_tune_strtgy(gs, tmp, n)

    call device_unmap(tmp, tmp_d)
    deallocate(tmp)

  end subroutine gs_tune_dev_strtgy

  !> Bind the non-blocking synchronisation strategy of the device MPI comm.
  !! backend held by @a gs, and optionally return the average time per
  !! gather-scatter operation under it.
  !!
  !! The strategy is the pair of independent choices of packing the halo per
  !! peer or in one go, and of waiting for the incoming messages one by one
  !! or all at once. All four are benchmarked on @a u and the fastest one is
  !! bound, unless NEKO_GS_STRTGY names one (1-4), which is then bound and
  !! merely timed.
  !!
  !! Unlike the comm. backend itself this is a rank-local choice -- it
  !! changes how a rank drives its own streams, not what goes on the wire --
  !! so every rank keeps its own winner and only the reported strategy is
  !! averaged over the ranks.
  !!
  !! @param gs gather-scatter kernel to tune, holding a device MPI comm.
  !! backend. Doing this to any other backend binds nothing and simply times
  !! it four times over
  !! @param u working vector to operate on
  !! @param n length of @a u
  !! @param t the average wall time (s) per operation under the bound
  !! strategy
  subroutine gs_tune_strtgy(gs, u, n, t)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    real(kind=dp), intent(out), optional :: t
    character(len=LOG_SIZE) :: log_buf
    integer, parameter :: strtgy(4) = [int(B'00'), int(B'01'), int(B'10'), &
         int(B'11')]
    real(kind=dp) :: strtgy_time(size(strtgy))
    character(len=255) :: env_strtgy
    integer :: i, env_len, best, avg

    call get_environment_variable("NEKO_GS_STRTGY", env_strtgy, env_len)

    if (env_len .eq. 0) then
       do i = 1, size(strtgy)
          call gs_set_strtgy(gs, strtgy(i))
          strtgy_time(i) = gs_time_ops(gs, u, n, GS_OP_MIN, GS_TUNE_NTRIALS)
       end do

       best = minloc(strtgy_time, 1)
       call gs_set_strtgy(gs, strtgy(best))
       if (present(t)) t = strtgy_time(best)

       ! Each rank picked for itself, so report the average of what they
       ! picked rather than rank 0's own choice
       avg = best
       call MPI_Allreduce(MPI_IN_PLACE, avg, 1, MPI_INTEGER, MPI_SUM, &
            NEKO_COMM)
       avg = avg / pe_size

       write(log_buf, '(A,B0.2,A)') 'Avg. strtgy  :         [', &
            strtgy(avg), ']'
    else
       read(env_strtgy(1:env_len), *) best

       if (best .lt. 1 .or. best .gt. size(strtgy)) then
          call neko_error('Invalid gs sync strtgy')
       end if

       call gs_set_strtgy(gs, strtgy(best))
       if (present(t)) then
          t = gs_time_ops(gs, u, n, GS_OP_MIN, GS_TUNE_NTRIALS)
       end if

       write(log_buf, '(A,B0.2,A)') 'Env. strtgy  :         [', &
            strtgy(best), ']'
    end if

    call neko_log%message(log_buf)

  end subroutine gs_tune_strtgy

  !> Bind the non-blocking synchronisation strategy @a strtgy of the device
  !! MPI comm. backend held by @a gs. A no-op for any other backend, so that
  !! the caller need not know which one is in place.
  !! @param gs gather-scatter kernel whose comm. backend is bound
  !! @param strtgy the strategy, a two bit mask (see gs_device_mpi_t)
  subroutine gs_set_strtgy(gs, strtgy)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: strtgy

    select type (c => gs%comm)
    type is (gs_device_mpi_t)
       c%nb_strtgy = strtgy
    end select

  end subroutine gs_set_strtgy

  !> The non-blocking synchronisation strategy bound in the device MPI comm.
  !! backend held by @a gs, or -1 if @a gs holds another backend.
  !! @param gs gather-scatter kernel to read the strategy of
  !! @return the strategy, a two bit mask (see gs_device_mpi_t), or -1
  function gs_get_strtgy(gs) result(strtgy)
    type(gs_t), intent(in) :: gs
    integer :: strtgy

    strtgy = -1

    select type (c => gs%comm)
    type is (gs_device_mpi_t)
       strtgy = c%nb_strtgy
    end select

  end function gs_get_strtgy

  !> Benchmark the coarray signaling modes and bind the fastest one, leaving
  !! @a gs holding a CAF backend that runs in that mode. Returns the winning
  !! mode's average time per gather-scatter operation, so the caller can use
  !! it as the CAF entry of the comm. backend benchmark.
  !!
  !! The mode is a program-wide binding shared by every gs_caf_t instance
  !! and the per-instance state depends on it, so every mode is measured on
  !! a freshly built backend. As in gs_tune_comm the timings are averaged
  !! over all ranks before the winner is picked, so every rank binds the
  !! same mode; the returned time is that rank-invariant average, which the
  !! caller's own averaging leaves unchanged.
  !!
  !! @param gs gather-scatter kernel to tune, left holding a CAF backend
  !! running in the winning mode
  !! @param u working vector to operate on
  !! @param n length of @a u
  !! @return the winning mode's average wall time (s) per operation
  function gs_tune_caf_signal(gs, u, n) result(t)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    real(kind=dp) :: t
    character(len=LOG_SIZE) :: log_buf
    character(len=13) :: label
    integer, allocatable :: mode(:)
    real(kind=dp), allocatable :: mode_time(:)
    integer :: i, best

    allocate(mode, source = gs_caf_signal_modes())
    allocate(mode_time(size(mode)))

    do i = 1, size(mode)
       ! Bind the mode before building the backend: gs_caf_init sets up its
       ! per-instance signaling state according to the mode in force
       call gs_caf_set_mode(mode(i))
       call gs_comm_switch(gs, GS_COMM_CAF)
       mode_time(i) = gs_time_ops(gs, u, n, GS_OP_MIN, GS_TUNE_NTRIALS)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, mode_time, size(mode), &
         MPI_DOUBLE_PRECISION, MPI_SUM, NEKO_COMM)
    mode_time = mode_time / pe_size

    best = minloc(mode_time, 1)

    do i = 1, size(mode)
       label = 'CAF ' // adjustl(gs_caf_mode_name(mode(i)))
       write(log_buf, '(A,A,ES10.3,A)') label, ': ', mode_time(i), ' s'
       call neko_log%message(log_buf)
    end do

    ! The last mode is the one currently allocated
    if (best .ne. size(mode)) then
       call gs_caf_set_mode(mode(best))
       call gs_comm_switch(gs, GS_COMM_CAF)
    end if

    call neko_log%message('Tuned CAF sig: ' // gs_caf_mode_name(mode(best)))

    t = mode_time(best)

    deallocate(mode, mode_time)

  end function gs_tune_caf_signal

end submodule gs_tune
