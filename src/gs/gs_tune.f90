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
  implicit none

  !> The host comm. backends the runtime autotuning can benchmark, in the
  !! order they are benchmarked and reported. The device backends are not
  !! among them: the autotuning only runs on host builds (see gs_init).
  integer, parameter :: GS_TUNE_BCKND(6) = [GS_COMM_MPI, GS_COMM_NEIGHBOUR, &
       GS_COMM_MPIRMA, GS_COMM_OPENSHMEM, GS_COMM_CAF, GS_COMM_UTOFU]

  !> Which of GS_TUNE_BCKND are benchmarked when NEKO_GS_TUNE is unset: all
  !! of them the build supports but CAF, which has to be asked for (see
  !! gs_tune_select).
  logical, parameter :: GS_TUNE_DEFAULT(6) = [.true., .true., .true., &
       .true., .false., .true.]

  !> Whether the coarray signaling mode has already been selected by the
  !! autotuner (NEKO_GS_CAF_SIGNALING=auto). The mode is a program-wide
  !! binding shared by every gs_caf_t instance, so it is benchmarked once,
  !! on the first gs_t that tunes the CAF backend, and kept: re-binding it
  !! later would change the protocol under instances already running.
  logical, save :: caf_signal_tuned = .false.

contains

  !> Host comm. backends to benchmark in the runtime autotuning: the ones
  !! selected by NEKO_GS_TUNE (see gs_tune_select) that this run can
  !! actually use (see gs_comm_tunable), in GS_TUNE_BCKND order. A selected
  !! backend the run cannot use is dropped, and reported in the log if it
  !! was named explicitly rather than merely part of the default set.
  !! @return the candidate backends, as GS_COMM_* constants. May be empty,
  !! or hold a single entry, in which case gs_tune_comm benchmarks nothing
  !! @note The list says nothing about which backend the gs kernel holds
  !! when the tuning starts; gs_tune_comm is told that separately.
  function gs_comm_host_cand() result(cand)
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

  end function gs_comm_host_cand

  !> Whether the comm. backend @a comm_bcknd can be driven by this build and
  !! this run, and may therefore be benchmarked. MPI and the neighbourhood
  !! collective need nothing but MPI-3; the rest abort in their init unless
  !! the corresponding support was built in.
  !! @param comm_bcknd comm. backend to check, one of the GS_TUNE_BCKND
  !! entries
  !! @return whether the backend is usable as a candidate
  function gs_comm_tunable(comm_bcknd) result(tunable)
    integer, intent(in) :: comm_bcknd
    logical :: tunable

    select case (comm_bcknd)
    case (GS_COMM_MPI, GS_COMM_NEIGHBOUR)
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
  !! set GS_TUNE_DEFAULT is used, which is everything but CAF -- a compiler
  !! may accept coarrays and still give the job a single image, and a
  !! backend that cannot communicate is at best wasted setup, so CAF has to
  !! be asked for. Selecting a backend outright (NEKO_GS_COMM) bypasses all
  !! of this.
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
  !! part in the autotuning (the device backends among them: the autotuning
  !! only runs on host builds).
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
  !! Each candidate returned by gs_comm_host_cand is benchmarked in turn on
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
    integer :: i, best, nmin, cur

    allocate(cand, source = gs_comm_host_cand())

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
    do i = 1, size(cand)
       if (cand(i) .eq. GS_COMM_CAF .and. gs_caf_signal_auto() .and. &
            .not. caf_signal_tuned) then
          ! Benchmark the signaling modes as well; this leaves the CAF
          ! backend allocated in the winning mode
          cand_time(i) = gs_tune_caf_signal(gs, tmp, n)
          caf_signal_tuned = .true.
       else
          if (cand(i) .ne. cur) call gs_comm_switch(gs, cand(i))
          cand_time(i) = gs_time_ops(gs, tmp, n, GS_OP_MIN, GS_TUNE_NTRIALS)
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

    call neko_log%message('Tuned comm   : ' // gs_comm_name(cand(best)))

    deallocate(cand_time)

  end subroutine gs_tune_comm

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
