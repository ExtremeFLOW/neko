! Copyright (c) 2020-2026, The Neko Authors
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
!> Gather-scatter
module gather_scatter
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_HIP, &
       NEKO_BCKND_CUDA, NEKO_BCKND_OPENCL, NEKO_BCKND_METAL, NEKO_DEVICE_MPI
  use gs_bcknd, only : gs_bcknd_t, GS_BCKND_CPU, GS_BCKND_SX, GS_BCKND_DEV
  use gs_device, only : gs_device_t
  use gs_sx, only : gs_sx_t
  use gs_cpu, only : gs_cpu_t
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use gs_comm, only : gs_comm_t, GS_COMM_MPI, GS_COMM_MPIGPU, GS_COMM_NCCL, &
       GS_COMM_NVSHMEM, GS_COMM_OPENSHMEM, GS_COMM_CAF, GS_COMM_NEIGHBOUR, &
       GS_COMM_UTOFU, GS_VEC_NC
  use gs_mpi, only : gs_mpi_t
  use gs_neighbour, only : gs_neighbour_t
  use gs_shmem, only : gs_shmem_t, GS_SHMEM_AVAIL
  use gs_caf, only : gs_caf_t, GS_CAF_AVAIL, gs_caf_signal_auto, &
       gs_caf_signal_modes, gs_caf_mode_name, gs_caf_mode_get, gs_caf_set_mode
  use gs_utofu, only : gs_utofu_t, GS_UTOFU_AVAIL
  use gs_device_mpi, only : gs_device_mpi_t
  use gs_device_nccl, only : gs_device_nccl_t
  use gs_device_shmem, only : gs_device_shmem_t
  use mesh, only : mesh_t
  use comm, only : pe_rank, pe_size, global_pe_size, NEKO_COMM
  use mpi_f08, only : MPI_Reduce, MPI_Allreduce, MPI_Barrier, MPI_IN_PLACE, &
       MPI_Wtime, MPI_SUM, MPI_MIN, MPI_INTEGER, MPI_INTEGER8, &
       MPI_DOUBLE_PRECISION
  use dofmap, only : dofmap_t
  use field, only : field_t
  use num_types, only : rp, dp, i2, i8, c_rp
  use htable, only : htable_i8_t, htable_iter_i8_t
  use stack, only : stack_i4_t, stack_i8_t
  use crystal_router, only : crystal_router_transfer, crystal_router_pack
  use math, only : sort
  use utils, only : neko_error, linear_index
  use logger, only : neko_log, LOG_SIZE
  use profiler, only : profiler_start_region, profiler_end_region
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_DEVICE, &
       device_sync, device_map, device_unmap, device_get_ptr, &
       device_event_record, device_deassociate
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_intptr_t, &
       c_sizeof, c_associated, c_size_t
  !$ use omp_lib, only : omp_get_thread_num
  implicit none
  private

  !> Gather-scatter kernel
  !! @note A gs_t instance is not re-entrant: each object owns mutable
  !! scratch buffers, backend device state and communication state that are
  !! reused by every gather-scatter operation. Concurrent calls to @a op on
  !! the same object must be serialised by the caller (or use one object
  !! per thread).
  type, public :: gs_t
     real(kind=rp), allocatable :: local_gs(:) !< Buffer for local gs-ops
     integer, allocatable :: local_dof_gs(:) !< Local dof to gs mapping
     integer, allocatable :: local_gs_dof(:) !< Local gs to dof mapping
     integer, allocatable :: local_blk_len(:) !< Local non-facet blocks
     integer, allocatable :: local_blk_off(:) !< Local non-facet blocks offset
     real(kind=rp), allocatable :: shared_gs(:) !< Buffer for shared gs-op
     !> Compact multi-component shared buffer for the fused vector gs (gs_op_r3),
     !! sized GS_VEC_NC*nshared, component-outer. On device it is device-mapped
     !! and the fused exchange operates on its device pointer.
     real(kind=rp), allocatable :: shared_gs_v(:)
     type(c_ptr) :: shared_gs_v_d = C_NULL_PTR !< Dev. ptr for shared_gs_v
     integer, allocatable :: shared_dof_gs(:) !< Shared dof to gs map.
     integer, allocatable :: shared_gs_dof(:) !< Shared gs to dof map.
     integer, allocatable :: shared_blk_len(:) !< Shared non-facet blocks
     integer, allocatable :: shared_blk_off(:) !< Shared non-facet blocks offset
     type(dofmap_t), pointer ::dofmap !< Dofmap for gs-ops
     type(htable_i8_t) :: shared_dofs !< Htable of shared dofs
     integer :: nlocal !< Local gs-ops
     integer :: nshared !< Shared gs-ops
     integer :: nlocal_blks !< Number of local blks
     integer :: nshared_blks !< Number of shared blks
     integer :: local_facet_offset !< offset for loc. facets
     integer :: shared_facet_offset !< offset for shr. facets
     class(gs_bcknd_t), allocatable :: bcknd !< Gather-scatter backend
     class(gs_comm_t), allocatable :: comm !< Comm. method
   contains
     procedure, private, pass(gs) :: gs_op_fld
     procedure, private, pass(gs) :: gs_op_r4
     procedure, pass(gs) :: gs_op_vector
     procedure, pass(gs) :: gs_op_r3
     procedure, pass(gs) :: gs_op_vector3
     procedure, pass(gs) :: init => gs_init
     procedure, pass(gs) :: free => gs_free
     generic :: op => gs_op_fld, gs_op_r4, gs_op_vector, gs_op_r3, &
          gs_op_vector3
  end type gs_t

  ! Expose available gather-scatter operation
  public :: GS_OP_ADD, GS_OP_MUL, GS_OP_MIN, GS_OP_MAX

  ! Expose available gather-scatter backends
  public :: GS_BCKND_CPU, GS_BCKND_SX, GS_BCKND_DEV

  ! Expose available gather-scatter comm. backends
  public :: GS_COMM_MPI, GS_COMM_MPIGPU, GS_COMM_NCCL, GS_COMM_NVSHMEM, &
       GS_COMM_OPENSHMEM, GS_COMM_CAF, GS_COMM_NEIGHBOUR, GS_COMM_UTOFU

  !> Number of timed (and untimed, warm-up) gather-scatter operations per
  !! candidate in the runtime autotuning of comm. backends and device
  !! synchronisation strategies
  integer, parameter :: GS_TUNE_NTRIALS = 100
  integer, parameter :: GS_TUNE_NWARMUP = 2

  !> Whether the coarray signaling mode has already been selected by the
  !! autotuner (NEKO_GS_CAF_SIGNALING=auto). The mode is a program-wide
  !! binding shared by every gs_caf_t instance, so it is benchmarked once,
  !! on the first gs_t that tunes the CAF backend, and kept: re-binding it
  !! later would change the protocol under instances already running.
  logical, save :: caf_signal_tuned = .false.

contains

  !> Initialize a gather-scatter kernel
  !> @param dofmap, global numbering of points and connectivity to base gs on
  !> @param bcknd, backend for executing the gs_ops
  !> @param comm_bcknd, backend for excuting the communication with
  subroutine gs_init(gs, dofmap, bcknd, comm_bcknd)
    class(gs_t), intent(inout) :: gs
    type(dofmap_t), target, intent(inout) :: dofmap
    character(len=LOG_SIZE) :: log_buf
    character(len=20) :: bcknd_str
    integer, optional :: bcknd, comm_bcknd
    integer :: i, ierr, bcknd_, comm_bcknd_
    integer(i8) :: glb_nshared, glb_nlocal
    logical :: use_device_mpi, use_device_nccl, use_device_shmem, use_host_mpi
    logical :: use_host_shmem
    logical :: use_caf
    logical :: use_neighbour
    logical :: use_utofu
    logical :: tune_comm
    real(kind=rp), allocatable :: tmp(:)
    type(c_ptr) :: tmp_d = C_NULL_PTR
    integer :: strtgy(4) = [int(B'00'), int(B'01'), int(B'10'), int(B'11')]
    integer :: avg_strtgy, env_len
    character(len=255) :: env_strtgy, env_gscomm
    real(kind=dp) :: strtgy_time(4)

    call gs%free()

    call neko_log%section('Gather-Scatter')
    ! Currently this uses the dofmap which also contains geometric information
    ! Only connectivity/numbering of points is technically necessary for gs
    gs%dofmap => dofmap

    use_device_mpi = .false.
    use_device_nccl = .false.
    use_device_shmem = .false.
    use_host_mpi = .false.
    use_host_shmem = .false.
    use_caf = .false.
    use_neighbour = .false.
    use_utofu = .false.
    tune_comm = .false.

    ! Check if a comm-backend is requested via env. variables
    call get_environment_variable("NEKO_GS_COMM", env_gscomm, env_len)
    if (env_len .gt. 0) then
       if (env_gscomm(1:env_len) .eq. "MPI") then
          use_host_mpi = .true.
       else if (env_gscomm(1:env_len) .eq. "MPIGPU") then
          use_device_mpi = .true.
       else if (env_gscomm(1:env_len) .eq. "NCCL") then
          use_device_nccl = .true.
       else if (env_gscomm(1:env_len) .eq. "SHMEM") then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             use_device_shmem = .true.
          else
             use_host_shmem = .true.
          end if
       else if (env_gscomm(1:env_len) .eq. "CAF") then
          use_caf = .true.
       else if (env_gscomm(1:env_len) .eq. "NEIGHBOUR" .or. &
            env_gscomm(1:env_len) .eq. "NEIGHBOR") then
          use_neighbour = .true.
       else if (env_gscomm(1:env_len) .eq. "UTOFU") then
          use_utofu = .true.
       else
          call neko_error('Unknown Gather-scatter comm. backend')
       end if
    end if


    if (present(comm_bcknd)) then
       comm_bcknd_ = comm_bcknd
    else if (use_host_mpi) then
       comm_bcknd_ = GS_COMM_MPI
    else if (use_device_mpi) then
       comm_bcknd_ = GS_COMM_MPIGPU
    else if (use_device_nccl) then
       comm_bcknd_ = GS_COMM_NCCL
    else if (use_device_shmem) then
       comm_bcknd_ = GS_COMM_NVSHMEM
    else if (use_host_shmem) then
       comm_bcknd_ = GS_COMM_OPENSHMEM
    else if (use_caf) then
       comm_bcknd_ = GS_COMM_CAF
    else if (use_neighbour) then
       comm_bcknd_ = GS_COMM_NEIGHBOUR
    else if (use_utofu) then
       comm_bcknd_ = GS_COMM_UTOFU
    else
       if (NEKO_DEVICE_MPI) then
          comm_bcknd_ = GS_COMM_MPIGPU
          use_device_mpi = .true.
       else
          ! No backend requested, benchmark the host backends once the
          ! schedule is known and keep the fastest one (see gs_tune_comm)
          comm_bcknd_ = GS_COMM_MPI
          tune_comm = (pe_size .gt. 1)
       end if
    end if

    call gs_comm_alloc(gs%comm, comm_bcknd_)

    if (tune_comm) then
       call neko_log%message('Comm         :         auto')
    else
       call neko_log%message('Comm         : ' // gs_comm_name(comm_bcknd_))
    end if
    ! Initialize a stack for each rank containing which dofs to send/recv at
    ! that rank
    call gs%comm%init_dofs()
    ! Initialize mapping between local ids and gather-scatter ids
    ! based on the global numbering in dofmap
    call gs_init_mapping(gs)
    ! Setup buffers and which ranks to send/recv data from based on mapping
    ! and initializes gs%comm (sets up gs%comm%send_dof and gs%comm%recv_dof and
    ! recv_pe/send_pe)
    call gs_schedule(gs)
    ! Global number of points not needing to be sent over mpi for gs operations
    ! "Internal points"
    glb_nlocal = int(gs%nlocal, i8)
    ! Global number of points needing to be communicated with other pes/ranks
    ! "external points"
    glb_nshared = int(gs%nshared, i8)
    ! Can be thought of a measure of the volume of this rank (glb_nlocal) and
    ! the surface area (glb_nshared) that is shared with other ranks
    ! Lots of internal volume compared to surface that needs communication is
    ! good

    if (pe_rank .eq. 0) then
       call MPI_Reduce(MPI_IN_PLACE, glb_nlocal, 1, &
            MPI_INTEGER8, MPI_SUM, 0, NEKO_COMM, ierr)

       call MPI_Reduce(MPI_IN_PLACE, glb_nshared, 1, &
            MPI_INTEGER8, MPI_SUM, 0, NEKO_COMM, ierr)
    else
       call MPI_Reduce(glb_nlocal, glb_nlocal, 1, &
            MPI_INTEGER8, MPI_SUM, 0, NEKO_COMM, ierr)

       call MPI_Reduce(glb_nshared, glb_nshared, 1, &
            MPI_INTEGER8, MPI_SUM, 0, NEKO_COMM, ierr)
    end if

    write(log_buf, '(A,I12)') 'Avg. internal: ', glb_nlocal/pe_size
    call neko_log%message(log_buf)
    write(log_buf, '(A,I12)') 'Avg. external: ', glb_nshared/pe_size
    call neko_log%message(log_buf)

    if (present(bcknd)) then
       bcknd_ = bcknd
    else
       if (NEKO_BCKND_SX .eq. 1) then
          bcknd_ = GS_BCKND_SX
       else if (NEKO_BCKND_DEVICE .eq. 1) then
          bcknd_ = GS_BCKND_DEV
       else
          bcknd_ = GS_BCKND_CPU
       end if
    end if

    ! Setup Gather-scatter backend
    select case (bcknd_)
    case (GS_BCKND_CPU)
       allocate(gs_cpu_t::gs%bcknd)
       bcknd_str = '         std'
    case (GS_BCKND_DEV)
       allocate(gs_device_t::gs%bcknd)
       if (NEKO_BCKND_HIP .eq. 1) then
          bcknd_str = '         hip'
       else if (NEKO_BCKND_CUDA .eq. 1) then
          bcknd_str = '        cuda'
       else if (NEKO_BCKND_OPENCL .eq. 1) then
          bcknd_str = '      opencl'
       else if (NEKO_BCKND_METAL .eq. 1) then
          bcknd_str = '       metal'
       end if
    case (GS_BCKND_SX)
       allocate(gs_sx_t::gs%bcknd)
       bcknd_str = '          sx'
    case default
       call neko_error('Unknown Gather-scatter backend')
    end select

    write(log_buf, '(A)') 'Backend      : ' // trim(bcknd_str)
    call neko_log%message(log_buf)


    call gs%bcknd%init(gs%nlocal, gs%nshared, gs%nlocal_blks, gs%nshared_blks)

    ! Plain base-type assignment; setting this through a
    ! select type (gs_device_t) miscompiles with CCE 21 at -O2/-O3,
    ! silently leaving shared points on the host so that the scatter
    ! overwrites the unpacked halo data with the stale host buffer
    if (use_device_mpi .or. use_device_nccl .or. use_device_shmem) then
       gs%bcknd%shared_on_host = .false.
    end if

    if (use_device_mpi) then
       if (pe_size .gt. 1) then
          ! Select fastest device MPI strategy at runtime
          select type (c => gs%comm)
          type is (gs_device_mpi_t)
             call get_environment_variable("NEKO_GS_STRTGY", env_strtgy, &
                  env_len)
             if (env_len .eq. 0) then
                allocate(tmp(dofmap%size()))
                call device_map(tmp, tmp_d, dofmap%size())
                tmp = 1.0_rp
                call device_memcpy(tmp, tmp_d, dofmap%size(), &
                     HOST_TO_DEVICE, sync = .false.)

                do i = 1, size(strtgy)
                   c%nb_strtgy = strtgy(i)
                   strtgy_time(i) = gs_time_ops(gs, tmp, dofmap%size(), &
                        GS_OP_ADD, GS_TUNE_NTRIALS)
                end do

                call device_unmap(tmp, tmp_d)
                deallocate(tmp)

                c%nb_strtgy = strtgy(minloc(strtgy_time, 1))

                avg_strtgy = minloc(strtgy_time, 1)
                call MPI_Allreduce(MPI_IN_PLACE, avg_strtgy, 1, &
                     MPI_INTEGER, MPI_SUM, NEKO_COMM)
                avg_strtgy = avg_strtgy / pe_size

                write(log_buf, '(A,B0.2,A)') 'Avg. strtgy  :         [', &
                     strtgy(avg_strtgy), ']'

             else
                read(env_strtgy(1:env_len), *) i

                if (i .lt. 1 .or. i .gt. 4) then
                   call neko_error('Invalid gs sync strtgy')
                end if

                c%nb_strtgy = strtgy(i)
                avg_strtgy = i

                write(log_buf, '(A,B0.2,A)') 'Env. strtgy  :         [', &
                     strtgy(avg_strtgy), ']'
             end if

             call neko_log%message(log_buf)

          end select
       end if
    end if

    ! Select the fastest host comm. backend at runtime
    if (tune_comm) then
       call gs_tune_comm(gs, dofmap%size(), gs_comm_host_cand())
    end if

    call neko_log%end_section()

  end subroutine gs_init

  !> Allocate a gather-scatter comm. backend of type @a comm_bcknd
  !! @param comm the allocated backend. Must not hold a live backend on
  !! entry: intent(out) releases it without running its free, which would
  !! leak whatever resources it holds (symmetric memory, coarrays, uTofu
  !! VCQs). Both callers free the previous backend themselves, see gs_free
  !! and gs_comm_switch.
  !! @param comm_bcknd comm. backend to allocate, one of the GS_COMM_*
  !! constants
  subroutine gs_comm_alloc(comm, comm_bcknd)
    class(gs_comm_t), allocatable, intent(out) :: comm
    integer, intent(in) :: comm_bcknd

    select case (comm_bcknd)
    case (GS_COMM_MPI)
       allocate(gs_mpi_t::comm)
    case (GS_COMM_MPIGPU)
       allocate(gs_device_mpi_t::comm)
    case (GS_COMM_NCCL)
       allocate(gs_device_nccl_t::comm)
    case (GS_COMM_NVSHMEM)
       allocate(gs_device_shmem_t::comm)
    case (GS_COMM_OPENSHMEM)
       allocate(gs_shmem_t::comm)
    case (GS_COMM_CAF)
       allocate(gs_caf_t::comm)
    case (GS_COMM_NEIGHBOUR)
       allocate(gs_neighbour_t::comm)
    case (GS_COMM_UTOFU)
       allocate(gs_utofu_t::comm)
    case default
       call neko_error('Unknown Gather-scatter comm. backend')
    end select

  end subroutine gs_comm_alloc

  !> Host comm. backends to consider in the runtime autotuning, in the order
  !! they are benchmarked. MPI and the neighbourhood collective need nothing
  !! but MPI-3; the remaining ones are only included if the corresponding
  !! support was built in, as their init aborts otherwise.
  !! @note The first entry is the backend the gs schedule is built with in
  !! gs_init, and must stay that way (see gs_tune_comm).
  !! @return the candidate backends, as GS_COMM_* constants
  function gs_comm_host_cand() result(cand)
    integer, allocatable :: cand(:)
    integer :: c(5), n

    n = 2
    c(1) = GS_COMM_MPI
    c(2) = GS_COMM_NEIGHBOUR

    ! The one-sided backends that address their peers by global PE or image
    ! number (OpenSHMEM PEs, coarray images) are only correct when NEKO_COMM
    ! spans every process, so skip them when the run has been split into
    ! several communicators (NEKO_COMM_ID). uTofu exchanges its addresses
    ! over NEKO_COMM itself and is unaffected.
    if (GS_SHMEM_AVAIL .and. (pe_size .eq. global_pe_size)) then
       n = n + 1
       c(n) = GS_COMM_OPENSHMEM
    end if

    if (GS_CAF_AVAIL .and. (pe_size .eq. global_pe_size)) then
       n = n + 1
       c(n) = GS_COMM_CAF
    end if

    if (GS_UTOFU_AVAIL) then
       n = n + 1
       c(n) = GS_COMM_UTOFU
    end if

    allocate(cand(n))
    cand = c(1:n)

  end function gs_comm_host_cand

  !> Name of the gather-scatter comm. backend @a comm_bcknd, right-adjusted
  !! for the log
  !! @param comm_bcknd comm. backend to name, one of the GS_COMM_* constants
  function gs_comm_name(comm_bcknd) result(name)
    integer, intent(in) :: comm_bcknd
    character(len=12) :: name

    select case (comm_bcknd)
    case (GS_COMM_MPI)
       name = '         MPI'
    case (GS_COMM_MPIGPU)
       name = '  Device MPI'
    case (GS_COMM_NCCL)
       name = '        NCCL'
    case (GS_COMM_NVSHMEM)
       name = '     NVSHMEM'
    case (GS_COMM_OPENSHMEM)
       name = '   OpenSHMEM'
    case (GS_COMM_CAF)
       name = '         CAF'
    case (GS_COMM_NEIGHBOUR)
       name = '  MPI neigh.'
    case (GS_COMM_UTOFU)
       name = '       uTofu'
    case default
       name = '     unknown'
       call neko_error('Unknown Gather-scatter comm. backend')
    end select

  end function gs_comm_name

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
  function gs_time_ops(gs, u, n, op, ntrials) result(t)
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

  !> Select the fastest of the comm. backends in @a cand at runtime.
  !!
  !! Each candidate is benchmarked in turn on a dummy vector of length @a n,
  !! reusing the schedule computed by gs_schedule (handed over from one
  !! candidate to the next, see gs_comm_t%take_schedule), and @a gs is left
  !! with the fastest backend allocated. On entry @a gs%comm must be the
  !! first candidate, i.e. the backend the schedule was built with.
  !!
  !! @param gs gather-scatter kernel to tune, left holding the fastest
  !! backend
  !! @param n length of the dummy vector to benchmark on (the dofmap size)
  !! @param cand comm. backends to consider, as GS_COMM_* constants, in the
  !! order they are benchmarked (see gs_comm_host_cand)
  !!
  !! @note Setting up and running a comm. backend is collective over
  !! NEKO_COMM (the neighbourhood backend builds a graph communicator, the
  !! one-sided backends allocate symmetric memory), so every rank has to
  !! arrive at the same decision: the per-rank timings are averaged over all
  !! ranks before the winner is picked.
  subroutine gs_tune_comm(gs, n, cand)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    integer, intent(in) :: cand(:)
    character(len=LOG_SIZE) :: log_buf
    character(len=13) :: label
    real(kind=dp), allocatable :: cand_time(:)
    real(kind=rp), allocatable :: tmp(:)
    type(c_ptr) :: tmp_d
    integer :: i, best, nmin

    ! A rank without any dofs skips the halo exchange in gs_op_vector, which
    ! would leave every other rank hanging in a collective based backend.
    ! Keep the default backend in that case.
    nmin = n
    call MPI_Allreduce(MPI_IN_PLACE, nmin, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM)
    if (nmin .eq. 0) then
       call neko_log%message('Comm tuning  :      skipped')
       call neko_log%message('Tuned comm   : ' // gs_comm_name(cand(1)))
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
          if (i .gt. 1) call gs_comm_switch(gs, cand(i))
          cand_time(i) = gs_time_ops(gs, tmp, n, GS_OP_MIN, GS_TUNE_NTRIALS)
       end if
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

    ! The last candidate is the one currently allocated
    if (best .ne. size(cand)) call gs_comm_switch(gs, cand(best))

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

  !> Deallocate a gather-scatter kernel
  subroutine gs_free(gs)
    class(gs_t), intent(inout) :: gs

    nullify(gs%dofmap)

    ! The device backend lazily maps these arrays in gather/scatter,
    ! and its free only releases the device side; remove any stale
    ! address table entries before deallocating the host side
    if (allocated(gs%local_gs)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%local_gs)
       end if
       deallocate(gs%local_gs)
    end if

    if (allocated(gs%local_dof_gs)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%local_dof_gs)
       end if
       deallocate(gs%local_dof_gs)
    end if

    if (allocated(gs%local_gs_dof)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%local_gs_dof)
       end if
       deallocate(gs%local_gs_dof)
    end if

    if (allocated(gs%local_blk_len)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%local_blk_len)
       end if
       deallocate(gs%local_blk_len)
    end if

    if (allocated(gs%local_blk_off)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%local_blk_off)
       end if
       deallocate(gs%local_blk_off)
    end if

    if (allocated(gs%shared_gs)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%shared_gs)
       end if
       deallocate(gs%shared_gs)
    end if

    if (allocated(gs%shared_gs_v)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(gs%shared_gs_v_d)) then
          call device_unmap(gs%shared_gs_v, gs%shared_gs_v_d)
       end if
       deallocate(gs%shared_gs_v)
    end if

    if (allocated(gs%shared_dof_gs)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%shared_dof_gs)
       end if
       deallocate(gs%shared_dof_gs)
    end if

    if (allocated(gs%shared_gs_dof)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%shared_gs_dof)
       end if
       deallocate(gs%shared_gs_dof)
    end if

    if (allocated(gs%shared_blk_len)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%shared_blk_len)
       end if
       deallocate(gs%shared_blk_len)
    end if

    if (allocated(gs%shared_blk_off)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_deassociate(gs%shared_blk_off)
       end if
       deallocate(gs%shared_blk_off)
    end if

    gs%nlocal = 0
    gs%nshared = 0
    gs%nlocal_blks = 0
    gs%nshared_blks = 0

    call gs%shared_dofs%free()

    if (allocated(gs%bcknd)) then
       call gs%bcknd%free()
       deallocate(gs%bcknd)
    end if

    if (allocated(gs%comm)) then
       call gs%comm%free()
       deallocate(gs%comm)
    end if

  end subroutine gs_free

  !> Setup mapping of dofs to gather-scatter operations
  subroutine gs_init_mapping(gs)
    type(gs_t), target, intent(inout) :: gs
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dofmap
    type(stack_i4_t), target :: local_dof, dof_local, shared_dof, dof_shared
    type(stack_i4_t), target :: local_face_dof, face_dof_local
    type(stack_i4_t), target :: shared_face_dof, face_dof_shared
    integer :: i, j, k, l, lx, ly, lz, max_id, max_sid, id, lid, dm_size
    type(htable_i8_t) :: dm !>
    type(htable_i8_t), pointer :: sdm

    dofmap => gs%dofmap
    msh => dofmap%msh
    sdm => gs%shared_dofs

    lx = dofmap%Xh%lx
    ly = dofmap%Xh%ly
    lz = dofmap%Xh%lz
    dm_size = dofmap%size()/lx

    call dm%init(dm_size, i)
    !>@note this might be a bit overkill,
    !!but having many collisions makes the init take too long.
    !!This is really critical to performance of the init
    call sdm%init(dofmap%size(), i)


    call local_dof%init()
    call dof_local%init()

    call local_face_dof%init()
    call face_dof_local%init()

    call shared_dof%init()
    call dof_shared%init()

    call shared_face_dof%init()
    call face_dof_shared%init()

    !
    ! Setup mapping for dofs points
    !

    max_id = 0
    max_sid = 0
    do i = 1, msh%nelv
       ! Local id of vertices
       lid = linear_index(1, 1, 1, i, lx, ly, lz)
       ! Check if this dof is shared among ranks or not
       if (dofmap%shared_dof(1, 1, 1, i)) then
          id = gs_mapping_add_dof(sdm, dofmap%dof(1, 1, 1, i), max_sid)
          !If add unique gather-scatter id to shared_dof stack
          call shared_dof%push(id)
          !If add local id to dof_shared stack
          call dof_shared%push(lid)
          !Now we have the mapping of local id <-> gather scatter id!
       else
          ! Same here, only here we know the point is local
          ! It will as such not need to be sent to other ranks later
          id = gs_mapping_add_dof(dm, dofmap%dof(1, 1, 1, i), max_id)
          call local_dof%push(id)
          call dof_local%push(lid)
       end if
       ! This procedure is then repeated for all vertices and edges
       ! Facets can be treated a little bit differently since they only have one
       ! neighbor

       lid = linear_index(lx, 1, 1, i, lx, ly, lz)
       if (dofmap%shared_dof(lx, 1, 1, i)) then
          id = gs_mapping_add_dof(sdm, dofmap%dof(lx, 1, 1, i), max_sid)
          call shared_dof%push(id)
          call dof_shared%push(lid)
       else
          id = gs_mapping_add_dof(dm, dofmap%dof(lx, 1, 1, i), max_id)
          call local_dof%push(id)
          call dof_local%push(lid)
       end if

       lid = linear_index(1, ly, 1, i, lx, ly, lz)
       if (dofmap%shared_dof(1, ly, 1, i)) then
          id = gs_mapping_add_dof(sdm, dofmap%dof(1, ly, 1, i), max_sid)
          call shared_dof%push(id)
          call dof_shared%push(lid)
       else
          id = gs_mapping_add_dof(dm, dofmap%dof(1, ly, 1, i), max_id)
          call local_dof%push(id)
          call dof_local%push(lid)
       end if

       lid = linear_index(lx, ly, 1, i, lx, ly, lz)
       if (dofmap%shared_dof(lx, ly, 1, i)) then
          id = gs_mapping_add_dof(sdm, dofmap%dof(lx, ly, 1, i), max_sid)
          call shared_dof%push(id)
          call dof_shared%push(lid)
       else
          id = gs_mapping_add_dof(dm, dofmap%dof(lx, ly, 1, i), max_id)
          call local_dof%push(id)
          call dof_local%push(lid)
       end if
       if (lz .gt. 1) then
          lid = linear_index(1, 1, lz, i, lx, ly, lz)
          if (dofmap%shared_dof(1, 1, lz, i)) then
             id = gs_mapping_add_dof(sdm, dofmap%dof(1, 1, lz, i), max_sid)
             call shared_dof%push(id)
             call dof_shared%push(lid)
          else
             id = gs_mapping_add_dof(dm, dofmap%dof(1, 1, lz, i), max_id)
             call local_dof%push(id)
             call dof_local%push(lid)
          end if

          lid = linear_index(lx, 1, lz, i, lx, ly, lz)
          if (dofmap%shared_dof(lx, 1, lz, i)) then
             id = gs_mapping_add_dof(sdm, dofmap%dof(lx, 1, lz, i), max_sid)
             call shared_dof%push(id)
             call dof_shared%push(lid)
          else
             id = gs_mapping_add_dof(dm, dofmap%dof(lx, 1, lz, i), max_id)
             call local_dof%push(id)
             call dof_local%push(lid)
          end if

          lid = linear_index(1, ly, lz, i, lx, ly, lz)
          if (dofmap%shared_dof(1, ly, lz, i)) then
             id = gs_mapping_add_dof(sdm, dofmap%dof(1, ly, lz, i), max_sid)
             call shared_dof%push(id)
             call dof_shared%push(lid)
          else
             id = gs_mapping_add_dof(dm, dofmap%dof(1, ly, lz, i), max_id)
             call local_dof%push(id)
             call dof_local%push(lid)
          end if

          lid = linear_index(lx, ly, lz, i, lx, ly, lz)
          if (dofmap%shared_dof(lx, ly, lz, i)) then
             id = gs_mapping_add_dof(sdm, dofmap%dof(lx, ly, lz, i), max_sid)
             call shared_dof%push(id)
             call dof_shared%push(lid)
          else
             id = gs_mapping_add_dof(dm, dofmap%dof(lx, ly, lz, i), max_id)
             call local_dof%push(id)
             call dof_local%push(lid)
          end if
       end if
    end do

    ! Clear local dofmap table
    call dm%clear()
    ! Get gather scatter ids and local ids of edges
    if (lz .gt. 1) then
       !
       ! Setup mapping for dofs on edges
       !
       do i = 1, msh%nelv

          !
          ! dofs on edges in x-direction
          !
          if (dofmap%shared_dof(2, 1, 1, i)) then
             do j = 2, lx - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(j, 1, 1, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(j, 1, 1, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do j = 2, lx - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(j, 1, 1, i), max_id)
                call local_dof%push(id)
                id = linear_index(j, 1, 1, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if
          if (dofmap%shared_dof(2, 1, lz, i)) then
             do j = 2, lx - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(j, 1, lz, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(j, 1, lz, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do j = 2, lx - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(j, 1, lz, i), max_id)
                call local_dof%push(id)
                id = linear_index(j, 1, lz, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if

          if (dofmap%shared_dof(2, ly, 1, i)) then
             do j = 2, lx - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(j, ly, 1, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(j, ly, 1, i, lx, ly, lz)
                call dof_shared%push(id)
             end do

          else
             do j = 2, lx - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, 1, i), max_id)
                call local_dof%push(id)
                id = linear_index(j, ly, 1, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if
          if (dofmap%shared_dof(2, ly, lz, i)) then
             do j = 2, lx - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(j, ly, lz, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(j, ly, lz, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do j = 2, lx - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, lz, i), max_id)
                call local_dof%push(id)
                id = linear_index(j, ly, lz, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if

          !
          ! dofs on edges in y-direction
          !
          if (dofmap%shared_dof(1, 2, 1, i)) then
             do k = 2, ly - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(1, k, 1, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(1, k, 1, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do k = 2, ly - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(1, k, 1, i), max_id)
                call local_dof%push(id)
                id = linear_index(1, k, 1, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if
          if (dofmap%shared_dof(1, 2, lz, i)) then
             do k = 2, ly - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(1, k, lz, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(1, k, lz, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do k = 2, ly - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(1, k, lz, i), max_id)
                call local_dof%push(id)
                id = linear_index(1, k, lz, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if

          if (dofmap%shared_dof(lx, 2, 1, i)) then
             do k = 2, ly - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(lx, k, 1, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(lx, k, 1, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do k = 2, ly - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, 1, i), max_id)
                call local_dof%push(id)
                id = linear_index(lx, k, 1, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if
          if (dofmap%shared_dof(lx, 2, lz, i)) then
             do k = 2, ly - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(lx, k, lz, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(lx, k, lz, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do k = 2, ly - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, lz, i), max_id)
                call local_dof%push(id)
                id = linear_index(lx, k, lz, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if
          !
          ! dofs on edges in z-direction
          !
          if (dofmap%shared_dof(1, 1, 2, i)) then
             do l = 2, lz - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(1, 1, l, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(1, 1, l, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do l = 2, lz - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(1, 1, l, i), max_id)
                call local_dof%push(id)
                id = linear_index(1, 1, l, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if

          if (dofmap%shared_dof(lx, 1, 2, i)) then
             do l = 2, lz - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(lx, 1, l, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(lx, 1, l, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do l = 2, lz - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(lx, 1, l, i), max_id)
                call local_dof%push(id)
                id = linear_index(lx, 1, l, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if

          if (dofmap%shared_dof(1, ly, 2, i)) then
             do l = 2, lz - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(1, ly, l, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(1, ly, l, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do l = 2, lz - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(1, ly, l, i), max_id)
                call local_dof%push(id)
                id = linear_index(1, ly, l, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if

          if (dofmap%shared_dof(lx, ly, 2, i)) then
             do l = 2, lz - 1
                id = gs_mapping_add_dof(sdm, dofmap%dof(lx, ly, l, i), max_sid)
                call shared_dof%push(id)
                id = linear_index(lx, ly, l, i, lx, ly, lz)
                call dof_shared%push(id)
             end do
          else
             do l = 2, lz - 1
                id = gs_mapping_add_dof(dm, dofmap%dof(lx, ly, l, i), max_id)
                call local_dof%push(id)
                id = linear_index(lx, ly, l, i, lx, ly, lz)
                call dof_local%push(id)
             end do
          end if
       end do
    end if

    ! Clear local dofmap table
    call dm%clear()

    !
    ! Setup mapping for dofs on facets
    !
    ! This is for 2d
    if (lz .eq. 1) then
       do i = 1, msh%nelv

          !
          ! dofs on edges in x-direction
          !
          if (msh%facet_neigh(3, i) .ne. 0) then
             if (dofmap%shared_dof(2, 1, 1, i)) then
                do j = 2, lx - 1
                   id = gs_mapping_add_dof(sdm, dofmap%dof(j, 1, 1, i), max_sid)
                   call shared_face_dof%push(id)
                   id = linear_index(j, 1, 1, i, lx, ly, lz)
                   call face_dof_shared%push(id)
                end do
             else
                do j = 2, lx - 1
                   id = gs_mapping_add_dof(dm, dofmap%dof(j, 1, 1, i), max_id)
                   call local_face_dof%push(id)
                   id = linear_index(j, 1, 1, i, lx, ly, lz)
                   call face_dof_local%push(id)
                end do
             end if
          end if

          if (msh%facet_neigh(4, i) .ne. 0) then
             if (dofmap%shared_dof(2, ly, 1, i)) then
                do j = 2, lx - 1
                   id = gs_mapping_add_dof(sdm, dofmap%dof(j, ly, 1, i), &
                        max_sid)
                   call shared_face_dof%push(id)
                   id = linear_index(j, ly, 1, i, lx, ly, lz)
                   call face_dof_shared%push(id)
                end do

             else
                do j = 2, lx - 1
                   id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, 1, i), &
                        max_id)
                   call local_face_dof%push(id)
                   id = linear_index(j, ly, 1, i, lx, ly, lz)
                   call face_dof_local%push(id)
                end do
             end if
          end if

          !
          ! dofs on edges in y-direction
          !
          if (msh%facet_neigh(1, i) .ne. 0) then
             if (dofmap%shared_dof(1, 2, 1, i)) then
                do k = 2, ly - 1
                   id = gs_mapping_add_dof(sdm, dofmap%dof(1, k, 1, i), max_sid)
                   call shared_face_dof%push(id)
                   id = linear_index(1, k, 1, i, lx, ly, lz)
                   call face_dof_shared%push(id)
                end do
             else
                do k = 2, ly - 1
                   id = gs_mapping_add_dof(dm, dofmap%dof(1, k, 1, i), max_id)
                   call local_face_dof%push(id)
                   id = linear_index(1, k, 1, i, lx, ly, lz)
                   call face_dof_local%push(id)
                end do
             end if
          end if

          if (msh%facet_neigh(2, i) .ne. 0) then
             if (dofmap%shared_dof(lx, 2, 1, i)) then
                do k = 2, ly - 1
                   id = gs_mapping_add_dof(sdm, dofmap%dof(lx, k, 1, i), &
                        max_sid)
                   call shared_face_dof%push(id)
                   id = linear_index(lx, k, 1, i, lx, ly, lz)
                   call face_dof_shared%push(id)
                end do
             else
                do k = 2, ly - 1
                   id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, 1, i), &
                        max_id)
                   call local_face_dof%push(id)
                   id = linear_index(lx, k, 1, i, lx, ly, lz)
                   call face_dof_local%push(id)
                end do
             end if
          end if
       end do
    else
       do i = 1, msh%nelv

          ! Facets in x-direction (s, t)-plane
          if (msh%facet_neigh(1, i) .ne. 0) then
             if (dofmap%shared_dof(1, 2, 2, i)) then
                do l = 2, lz - 1
                   do k = 2, ly - 1
                      id = gs_mapping_add_dof(sdm, dofmap%dof(1, k, l, i), &
                           max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(1, k, l, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do l = 2, lz - 1
                   do k = 2, ly - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(1, k, l, i), &
                           max_id)
                      call local_face_dof%push(id)
                      id = linear_index(1, k, l, i, lx, ly, lz)
                      call face_dof_local%push(id)
                   end do
                end do
             end if
          end if

          if (msh%facet_neigh(2, i) .ne. 0) then
             if (dofmap%shared_dof(lx, 2, 2, i)) then
                do l = 2, lz - 1
                   do k = 2, ly - 1
                      id = gs_mapping_add_dof(sdm, dofmap%dof(lx, k, l, i), &
                           max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(lx, k, l, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do l = 2, lz - 1
                   do k = 2, ly - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, l, i), &
                           max_id)
                      call local_face_dof%push(id)
                      id = linear_index(lx, k, l, i, lx, ly, lz)
                      call face_dof_local%push(id)
                   end do
                end do
             end if
          end if

          ! Facets in y-direction (r, t)-plane
          if (msh%facet_neigh(3, i) .ne. 0) then
             if (dofmap%shared_dof(2, 1, 2, i)) then
                do l = 2, lz - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(sdm, dofmap%dof(j, 1, l, i), &
                           max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(j, 1, l, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do l = 2, lz - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(j, 1, l, i), &
                           max_id)
                      call local_face_dof%push(id)
                      id = linear_index(j, 1, l, i, lx, ly, lz)
                      call face_dof_local%push(id)
                   end do
                end do
             end if
          end if

          if (msh%facet_neigh(4, i) .ne. 0) then
             if (dofmap%shared_dof(2, ly, 2, i)) then
                do l = 2, lz - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(sdm, dofmap%dof(j, ly, l, i), &
                           max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(j, ly, l, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do l = 2, lz - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, l, i), &
                           max_id)
                      call local_face_dof%push(id)
                      id = linear_index(j, ly, l, i, lx, ly, lz)
                      call face_dof_local%push(id)
                   end do
                end do
             end if
          end if

          ! Facets in z-direction (r, s)-plane
          if (msh%facet_neigh(5, i) .ne. 0) then
             if (dofmap%shared_dof(2, 2, 1, i)) then
                do k = 2, ly - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(sdm, dofmap%dof(j, k, 1, i), &
                           max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(j, k, 1, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do k = 2, ly - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(j, k, 1, i), &
                           max_id)
                      call local_face_dof%push(id)
                      id = linear_index(j, k, 1, i, lx, ly, lz)
                      call face_dof_local%push(id)
                   end do
                end do
             end if
          end if

          if (msh%facet_neigh(6, i) .ne. 0) then
             if (dofmap%shared_dof(2, 2, lz, i)) then
                do k = 2, ly - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(sdm, dofmap%dof(j, k, lz, i), &
                           max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(j, k, lz, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do k = 2, ly - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(j, k, lz, i), &
                           max_id)
                      call local_face_dof%push(id)
                      id = linear_index(j, k, lz, i, lx, ly, lz)
                      call face_dof_local%push(id)
                   end do
                end do
             end if
          end if
       end do
    end if


    call dm%free()

    gs%nlocal = local_dof%size() + local_face_dof%size()
    gs%local_facet_offset = local_dof%size() + 1

    ! Finalize local dof to gather-scatter index
    allocate(gs%local_dof_gs(gs%nlocal))

    ! Add dofs on points and edges

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (dof_array => local_dof%data)
    type is (integer)
       j = local_dof%size()
       do i = 1, j
          gs%local_dof_gs(i) = dof_array(i)
       end do
    end select
    call local_dof%free()

    ! Add dofs on faces

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (dof_array => local_face_dof%data)
    type is (integer)
       do i = 1, local_face_dof%size()
          gs%local_dof_gs(i + j) = dof_array(i)
       end do
    end select
    call local_face_dof%free()

    ! Finalize local gather-scatter index to dof
    allocate(gs%local_gs_dof(gs%nlocal))

    ! Add gather-scatter index on points and edges

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (dof_array => dof_local%data)
    type is (integer)
       j = dof_local%size()
       do i = 1, j
          gs%local_gs_dof(i) = dof_array(i)
       end do
    end select
    call dof_local%free()

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (dof_array => face_dof_local%data)
    type is (integer)
       do i = 1, face_dof_local%size()
          gs%local_gs_dof(i+j) = dof_array(i)
       end do
    end select
    call face_dof_local%free()

    call gs_qsort_dofmap(gs%local_dof_gs, gs%local_gs_dof, &
         gs%nlocal, 1, gs%nlocal)

    call gs_find_blks(gs%local_dof_gs, gs%local_blk_len, &
         gs%local_blk_off, gs%nlocal_blks, gs%nlocal, gs%local_facet_offset)

    ! Allocate buffer for local gs-ops
    allocate(gs%local_gs(gs%nlocal))

    gs%nshared = shared_dof%size() + shared_face_dof%size()
    gs%shared_facet_offset = shared_dof%size() + 1

    ! Finalize shared dof to gather-scatter index
    allocate(gs%shared_dof_gs(gs%nshared))

    ! Add shared dofs on points and edges

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (dof_array => shared_dof%data)
    type is (integer)
       j = shared_dof%size()
       do i = 1, j
          gs%shared_dof_gs(i) = dof_array(i)
       end do
    end select
    call shared_dof%free()

    ! Add shared dofs on faces

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (dof_array => shared_face_dof%data)
    type is (integer)
       do i = 1, shared_face_dof%size()
          gs%shared_dof_gs(i + j) = dof_array(i)
       end do
    end select
    call shared_face_dof%free()

    ! Finalize shared gather-scatter index to dof
    allocate(gs%shared_gs_dof(gs%nshared))

    ! Add dofs on points and edges

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (dof_array => dof_shared%data)
    type is (integer)
       j = dof_shared%size()
       do i = 1, j
          gs%shared_gs_dof(i) = dof_array(i)
       end do
    end select
    call dof_shared%free()

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type (dof_array => face_dof_shared%data)
    type is (integer)
       do i = 1, face_dof_shared%size()
          gs%shared_gs_dof(i + j) = dof_array(i)
       end do
    end select
    call face_dof_shared%free()

    ! Allocate buffer for shared gs-ops
    allocate(gs%shared_gs(gs%nshared))

    ! Compact multi-component shared buffer for the fused vector gs. On the
    ! device it is mapped so the fused exchange can use its device pointer.
    allocate(gs%shared_gs_v(max(1, GS_VEC_NC * gs%nshared)))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(gs%shared_gs_v, gs%shared_gs_v_d, &
            max(1, GS_VEC_NC * gs%nshared))
    end if

    if (gs%nshared .gt. 0) then
       call gs_qsort_dofmap(gs%shared_dof_gs, gs%shared_gs_dof, &
            gs%nshared, 1, gs%nshared)

       call gs_find_blks(gs%shared_dof_gs, gs%shared_blk_len, &
            gs%shared_blk_off, gs%nshared_blks, gs%nshared, &
            gs%shared_facet_offset)
    end if

  contains

    !> Register a unique dof
    !! Takes the unique id dof and checks if it is in the htable map_
    !! If it is we return the gather-scatter id this global dof has been
    !! assigned to. This is done as the global id can be very large
    !! max(integer8), but the number of local points is at most max(integer4)
    !! @param map_, htable of global unique id to local unique id
    !! @param dof, global unique id of dof
    !! @param max_id, current number of entries in map_
    function gs_mapping_add_dof(map_, dof, max_id) result(id)
      type(htable_i8_t), intent(inout) :: map_
      integer(kind=i8), intent(inout) :: dof
      integer, intent(inout) :: max_id
      integer :: id

      if (map_%get(dof, id) .gt. 0) then
         max_id = max_id + 1
         call map_%set(dof, max_id)
         id = max_id
      end if

    end function gs_mapping_add_dof

    !> Sort the dof lists based on the dof to gather-scatter list
    recursive subroutine gs_qsort_dofmap(dg, gd, n, lo, hi)
      integer, intent(inout) :: n
      integer, dimension(n), intent(inout) :: dg
      integer, dimension(n), intent(inout) :: gd
      integer :: lo, hi
      integer :: tmp, i, j, pivot

      i = lo - 1
      j = hi + 1
      pivot = dg((lo + hi) / 2)
      do
         do
            i = i + 1
            if (dg(i) .ge. pivot) exit
         end do

         do
            j = j - 1
            if (dg(j) .le. pivot) exit
         end do

         if (i .lt. j) then
            tmp = dg(i)
            dg(i) = dg(j)
            dg(j) = tmp

            tmp = gd(i)
            gd(i) = gd(j)
            gd(j) = tmp
         else if (i .eq. j) then
            i = i + 1
            exit
         else
            exit
         end if
      end do
      if (lo .lt. j) call gs_qsort_dofmap(dg, gd, n, lo, j)
      if (i .lt. hi) call gs_qsort_dofmap(dg, gd, n, i, hi)

    end subroutine gs_qsort_dofmap

    !> Find blocks sharing dofs in non-facet data
    subroutine gs_find_blks(dg, blk_len, blk_off, nblks, n, m)
      integer, intent(in) :: n
      integer, intent(in) :: m
      integer, dimension(n), intent(inout) :: dg
      integer, allocatable, intent(inout) :: blk_len(:)
      integer, allocatable, intent(inout) :: blk_off(:)
      integer, intent(inout) :: nblks
      integer :: i, j
      integer :: id, count
      type(stack_i4_t), target :: blks

      call blks%init()
      i = 1
      do while (i .lt. m)
         id = dg(i)
         count = 1
         j = i
         do while ( j+1 .le. n .and. dg(j+1) .eq. id)
            j = j + 1
            count = count + 1
         end do
         call blks%push(count)
         i = j + 1
      end do

      select type (blk_array => blks%data)
      type is (integer)
         nblks = blks%size()
         allocate(blk_len(nblks))
         do i = 1, nblks
            blk_len(i) = blk_array(i)
         end do
         allocate(blk_off(nblks))
         blk_off(1) = 0
         do i = 2, nblks
            blk_off(i) = blk_off(i - 1) + blk_len(i - 1)
         end do
      end select
      call blks%free()

    end subroutine gs_find_blks

  end subroutine gs_init_mapping

  !> Schedule shared gather-scatter operations
  !! @details Discovers, for every shared dof, the set of ranks that also hold
  !! it, using a canonical-owner rendezvous routed through the crystal router
  !! instead of the previous O(P) shifted neighbour exchange. Each dof global
  !! id is hashed to an owner rank, `mod(gid, P)`; the owner gathers all holders
  !! and reflects, to each holder, the other holders. Each holder then registers
  !! the dof for both send and receive with every peer. The per-peer dof lists
  !! are ordered by dof global id, which gives the same ordering on both ranks
  !! of a pair (the matching invariant the gather-scatter exchange relies on)
  !! without ever transmitting a single sequenced key list. Scales to >1e5 ranks
  !! without O(P) buffers or unexpected-message buffer exhaustion.
  subroutine gs_schedule(gs)
    type(gs_t), target, intent(inout) :: gs
    type(htable_iter_i8_t) :: it
    type(stack_i4_t) :: send_pe, recv_pe
    type(stack_i8_t) :: cr_buf
    integer(i8), allocatable :: buf(:)
    integer(i8), pointer :: cr_data(:)
    integer(i8), allocatable :: rgid(:), gtmp(:)
    integer, allocatable :: rpeer(:), rgsid(:), rperm(:), gperm(:)
    integer(i8) :: gid
    integer :: i, j, n, owner, nrec, peer, shared_gs_id, tmp
    integer :: a, b, cnt, t

    call send_pe%init()
    call recv_pe%init()

    !
    ! Phase 1: route every local shared dof to its canonical owner.
    !   record = [dest=owner, len=2, gid, origin]
    !
    call cr_buf%init(max(gs%shared_dofs%num_entries(), 1) * 4)
    call it%init(gs%shared_dofs)
    do while (it%next())
       gid = it%key()
       owner = int(modulo(gid, int(pe_size, i8)))
       call crystal_router_pack(cr_buf, owner, [gid, int(pe_rank, i8)])
    end do

    n = cr_buf%size()
    allocate(buf(max(n, 1)))
    if (n .gt. 0) then
       cr_data => cr_buf%array()
       buf(1:n) = cr_data(1:n)
    end if
    call cr_buf%free()

    call crystal_router_transfer(buf, n)

    !
    ! Phase 2: at the owner, group holders by gid and reflect, to each holder,
    !   every other holder of the same dof.
    !   reply = [dest=holder, len=2, gid, peer]
    !
    nrec = n / 4 ! every record here has the fixed form [me, 2, gid, origin]
    allocate(rgid(max(nrec, 1)), rgsid(max(nrec, 1)), gperm(max(nrec, 1)))
    do i = 1, nrec
       rgid(i) = buf((i - 1) * 4 + 3) ! gid
       rgsid(i) = int(buf((i - 1) * 4 + 4)) ! origin rank (reuse array)
    end do
    if (nrec .gt. 0) call gs_sort_i8(rgid, gperm, nrec)

    call cr_buf%init(max(n, 1))
    i = 1
    do while (i .le. nrec)
       j = i
       do while (j .le. nrec)
          if (rgid(j) .ne. rgid(i)) exit
          j = j + 1
       end do
       ! Reflect, to each holder, every other holder of this dof.
       if (j - i .gt. 1) then
          do a = i, j - 1 ! recipient holder
             do b = i, j - 1 ! the other holder
                if (a .eq. b) cycle
                call crystal_router_pack(cr_buf, rgsid(gperm(a)), &
                     [rgid(i), int(rgsid(gperm(b)), i8)])
             end do
          end do
       end if
       i = j
    end do
    deallocate(rgid, rgsid, gperm)

    n = cr_buf%size()
    if (allocated(buf)) deallocate(buf)
    allocate(buf(max(n, 1)))
    if (n .gt. 0) then
       cr_data => cr_buf%array()
       buf(1:n) = cr_data(1:n)
    end if
    call cr_buf%free()

    call crystal_router_transfer(buf, n)

    !
    ! Phase 3: register each (dof, peer) for both send and receive. Order each
    !   peer's dof list by gid so both ranks of a pair agree on the order.
    !
    nrec = n / 4 ! replies are [me, 2, gid, peer]
    allocate(rgid(max(nrec, 1)), rpeer(max(nrec, 1)), rgsid(max(nrec, 1)), &
         rperm(max(nrec, 1)))
    do i = 1, nrec
       gid = buf((i - 1) * 4 + 3)
       rgid(i) = gid
       rpeer(i) = int(buf((i - 1) * 4 + 4))
       tmp = gs%shared_dofs%get(gid, shared_gs_id)
       rgsid(i) = shared_gs_id
    end do

    ! Sort by peer; within each peer run, sort by gid and register in that order.
    if (nrec .gt. 0) call sort(rpeer, rperm, nrec)
    a = 1
    do while (a .le. nrec)
       b = a
       do while (b .le. nrec)
          if (rpeer(b) .ne. rpeer(a)) exit
          b = b + 1
       end do
       peer = rpeer(a)
       cnt = b - a
       allocate(gtmp(cnt), gperm(cnt))
       do t = 1, cnt
          gtmp(t) = rgid(rperm(a + t - 1))
       end do
       call gs_sort_i8(gtmp, gperm, cnt)
       do t = 1, cnt
          shared_gs_id = rgsid(rperm(a + gperm(t) - 1))
          call gs%comm%send_dof(peer)%push(shared_gs_id)
          call gs%comm%recv_dof(peer)%push(shared_gs_id)
       end do
       deallocate(gtmp, gperm)
       call send_pe%push(peer)
       call recv_pe%push(peer)
       a = b
    end do
    deallocate(rgid, rpeer, rgsid, rperm)
    if (allocated(buf)) deallocate(buf)

    call gs%comm%init(send_pe, recv_pe)

    call send_pe%free()
    call recv_pe%free()

    !This arrays seems to take massive amounts of memory...
    call gs%shared_dofs%free()

  end subroutine gs_schedule

  !> Heap sort for 64-bit integer arrays, returning the permutation @a ind.
  !! Local helper for gs_schedule (the generic math sort has no i8 variant).
  subroutine gs_sort_i8(a, ind, n)
    integer, intent(in) :: n
    integer(i8), intent(inout) :: a(n)
    integer, intent(out) :: ind(n)
    integer(i8) :: aa
    integer :: j, ir, i, ii, l

    do j = 1, n
       ind(j) = j
    end do

    if (n .le. 1) return

    l = n/2 + 1
    ir = n
    do while (.true.)
       if (l .gt. 1) then
          l = l - 1
          aa = a(l)
          ii = ind(l)
       else
          aa = a(ir)
          ii = ind(ir)
          a(ir) = a(1)
          ind(ir) = ind(1)
          ir = ir - 1
          if (ir .eq. 1) then
             a(1) = aa
             ind(1) = ii
             return
          end if
       end if
       i = l
       j = l + l
       do while (j .le. ir)
          if (j .lt. ir) then
             if (a(j) .lt. a(j + 1)) j = j + 1
          end if
          if (aa .lt. a(j)) then
             a(i) = a(j)
             ind(i) = ind(j)
             i = j
             j = j + j
          else
             j = ir + 1
          end if
       end do
       a(i) = aa
       ind(i) = ii
    end do
  end subroutine gs_sort_i8

  !> Gather-scatter operation on a field @a u with op @a op
  subroutine gs_op_fld(gs, u, op, event)
    class(gs_t), intent(inout) :: gs
    type(field_t), intent(inout) :: u
    type(c_ptr), optional, intent(inout) :: event
    integer :: n, op

    n = u%msh%nelv * u%Xh%lx * u%Xh%ly * u%Xh%lz
    if (present(event)) then
       call gs_op_vector(gs, u%x, n, op, event)
    else
       call gs_op_vector(gs, u%x, n, op)
    end if

  end subroutine gs_op_fld

  !> Gather-scatter operation on a rank 4 array
  subroutine gs_op_r4(gs, u, n, op, event)
    class(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), contiguous, dimension(:,:,:,:), intent(inout) :: u
    type(c_ptr), optional, intent(inout) :: event
    integer :: op

    if (present(event)) then
       call gs_op_vector(gs, u, n, op, event)
    else
       call gs_op_vector(gs, u, n, op)
    end if

  end subroutine gs_op_r4

  !> Gather-scatter operation on a vector @a u with op @a op
  subroutine gs_op_vector(gs, u, n, op, event)
    class(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), optional, intent(inout) :: event
    integer :: m, l, op, lo, so, tid
    type(c_ptr) :: scatter_event

    lo = gs%local_facet_offset
    so = -gs%shared_facet_offset
    m = gs%nlocal
    l = gs%nshared

    ! Capture the calling thread id before opening any parallel region; it
    ! is used as the MPI tag so concurrent gs ops driven from different
    ! threads (device path) don't collide.
    tid = 0
    !$ tid = omp_get_thread_num()

    ! Resolve the optional event into a non-optional local before opening
    ! the parallel region. An absent optional dummy must not be captured by
    ! the region's data-sharing, otherwise the outlined region prologue
    ! dereferences a null descriptor (segfaults on CCE).
    scatter_event = C_NULL_PTR
    if (present(event)) scatter_event = event

    !$omp parallel if (NEKO_BCKND_DEVICE .eq. 0)
    call profiler_start_region("gather_scatter", 5)
    ! Gather shared dofs
    if (pe_size .gt. 1 .and. n .gt. 0) then
       call profiler_start_region("gs_nbrecv", 13)
       call gs%comm%nbrecv(tid)
       call profiler_end_region("gs_nbrecv", 13)
       call profiler_start_region("gs_gather_shared", 14)
       call gs%bcknd%gather(gs%shared_gs, l, so, gs%shared_dof_gs, u, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
            gs%shared_blk_off, op, .true.)
       call profiler_end_region("gs_gather_shared", 14)
       call profiler_start_region("gs_nbsend", 6)
       call gs%comm%nbsend(gs%shared_gs, l, tid, &
            gs%bcknd%gather_event, gs%bcknd%gs_stream)
       call profiler_end_region("gs_nbsend", 6)

    end if

    ! Gather-scatter local dofs
    call profiler_start_region("gs_local", 12)
    call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, gs%local_blk_off, &
         op, .false.)
    call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, gs%local_blk_off, &
         .false., C_NULL_PTR)
    call profiler_end_region("gs_local", 12)
    ! Scatter shared dofs
    if (pe_size .gt. 1 .and. n .gt. 0) then
       call profiler_start_region("gs_nbwait", 7)
       call gs%comm%nbwait(gs%shared_gs, l, op, gs%bcknd%gs_stream)
       call profiler_end_region("gs_nbwait", 7)
       call profiler_start_region("gs_scatter_shared", 15)
       call gs%bcknd%scatter(gs%shared_gs, l,&
            gs%shared_dof_gs, u, n, &
            gs%shared_gs_dof, gs%nshared_blks, &
            gs%shared_blk_len, gs%shared_blk_off, .true., scatter_event)
       call profiler_end_region("gs_scatter_shared", 15)
    end if

    call profiler_end_region("gather_scatter", 5)
    !$omp end parallel
  end subroutine gs_op_vector

  !> Gather-scatter operation on a 3-component vector of rank-4 arrays
  !! (@a u1, @a u2, @a u3) with op @a op; see gs_op_vector3.
  subroutine gs_op_r3(gs, u1, u2, u3, n, op, event)
    class(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), contiguous, dimension(:,:,:,:), intent(inout) :: u1, u2, u3
    type(c_ptr), optional, intent(inout) :: event
    integer :: op

    if (present(event)) then
       call gs_op_vector3(gs, u1, u2, u3, n, op, event)
    else
       call gs_op_vector3(gs, u1, u2, u3, n, op)
    end if

  end subroutine gs_op_r3

  !> Gather-scatter operation on a 3-component vector (@a u1, @a u2, @a u3)
  !! with op @a op. When the comm backend supports a fused vector halo
  !! exchange, the three components are communicated in a single round;
  !! otherwise this falls back to three independent scalar gs_op_vector
  !! calls (identical result, 3 rounds). The on-node gather/scatter is
  !! always done per component (scalar), so only the communication cost is
  !! reduced. On device backends the caller's event is recorded once, after
  !! the last component's scatter (or superseded by hard synchronisation on
  !! host-mirrored comms), so a single event sync covers all components.
  subroutine gs_op_vector3(gs, u1, u2, u3, n, op, event)
    class(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u1, u2, u3
    type(c_ptr), optional, intent(inout) :: event
    integer :: m, l, op, lo, so, tid
    integer, parameter :: nc = 3
    type(c_ptr) :: scatter_event

    ! Fall back to nc independent scalar exchanges when the comm backend has
    ! no fused vector path.
    if (.not. gs%comm%vec_supported) then
       if (present(event)) then
          call gs_op_vector(gs, u1, n, op, event)
          call gs_op_vector(gs, u2, n, op, event)
          call gs_op_vector(gs, u3, n, op, event)
       else
          call gs_op_vector(gs, u1, n, op)
          call gs_op_vector(gs, u2, n, op)
          call gs_op_vector(gs, u3, n, op)
       end if
       return
    end if

    lo = gs%local_facet_offset
    so = -gs%shared_facet_offset
    m = gs%nlocal
    l = gs%nshared

    tid = 0
    !$ tid = omp_get_thread_num()

    scatter_event = C_NULL_PTR
    if (present(event)) scatter_event = event

    if (NEKO_BCKND_DEVICE .eq. 0) then

       !$omp parallel
       call profiler_start_region("gather_scatter", 5)

       ! Gather each component's shared dofs directly into its column of
       ! shared_gs_v (the host backends write the actual argument), then
       ! launch ONE fused exchange covering all nc components.
       if (pe_size .gt. 1 .and. n .gt. 0) then
          call gs%comm%nbrecv_vec(tid, nc)
          call gs%bcknd%gather(gs%shared_gs_v(1), l, so, gs%shared_dof_gs, &
               u1, n, gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
               gs%shared_blk_off, op, .true.)
          call gs%bcknd%gather(gs%shared_gs_v(l + 1), l, so, &
               gs%shared_dof_gs, u2, n, gs%shared_gs_dof, gs%nshared_blks, &
               gs%shared_blk_len, gs%shared_blk_off, op, .true.)
          call gs%bcknd%gather(gs%shared_gs_v(2*l + 1), l, so, &
               gs%shared_dof_gs, u3, n, gs%shared_gs_dof, gs%nshared_blks, &
               gs%shared_blk_len, gs%shared_blk_off, op, .true.)
          call gs%comm%nbsend_vec(gs%shared_gs_v, l, nc, tid, &
               gs%bcknd%gather_event, gs%bcknd%gs_stream)
       end if

       ! Local gather-scatter, one scalar pass per component (reuses local_gs;
       ! the internal barriers make the sequential reuse safe).
       call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u1, n, &
            gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, &
            gs%local_blk_off, op, .false.)
       call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u1, n, &
            gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, &
            gs%local_blk_off, .false., C_NULL_PTR)
       call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u2, n, &
            gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, &
            gs%local_blk_off, op, .false.)
       call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u2, n, &
            gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, &
            gs%local_blk_off, .false., C_NULL_PTR)
       call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u3, n, &
            gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, &
            gs%local_blk_off, op, .false.)
       call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u3, n, &
            gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, &
            gs%local_blk_off, .false., C_NULL_PTR)

       ! Wait for the fused exchange and scatter each component back.
       if (pe_size .gt. 1 .and. n .gt. 0) then
          call gs%comm%nbwait_vec(gs%shared_gs_v, l, nc, op, &
               gs%bcknd%gs_stream)
          call gs%bcknd%scatter(gs%shared_gs_v(1), l, gs%shared_dof_gs, u1, &
               n, gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
               gs%shared_blk_off, .true., scatter_event)
          call gs%bcknd%scatter(gs%shared_gs_v(l + 1), l, gs%shared_dof_gs, &
               u2, n, gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
               gs%shared_blk_off, .true., scatter_event)
          call gs%bcknd%scatter(gs%shared_gs_v(2*l + 1), l, &
               gs%shared_dof_gs, u3, n, gs%shared_gs_dof, gs%nshared_blks, &
               gs%shared_blk_len, gs%shared_blk_off, .true., scatter_event)
       end if

       call profiler_end_region("gather_scatter", 5)
       !$omp end parallel

    else

       call gs_op_r3_device(gs, u1, u2, u3, n, op, nc, lo, so, m, l, tid, &
            scatter_event)

    end if

  end subroutine gs_op_vector3

  !> Device path for the fused 3-component gs. The device backend's shared
  !! gather/scatter always operate on its internal shared buffer, so each
  !! component is staged between that buffer and its column of the compact
  !! vector buffer gs%shared_gs_v: through the host mirror when the comm
  !! backend is a host backend (shared_on_host, e.g. OpenCL/Metal or a device
  !! build with host MPI), or with device-to-device copies when the comm is
  !! device-resident (device MPI/NCCL/NVSHMEM). Apart from the column staging
  !! the two modes are identical; the comm backends transparently pick the
  !! host array or its device pointer, as in gs_op_vector.
  subroutine gs_op_r3_device(gs, u1, u2, u3, n, op, nc, lo, so, m, l, tid, &
       scatter_event)
    class(gs_t), intent(inout) :: gs
    integer, intent(in) :: n, op, nc, lo, so, m, l, tid
    real(kind=rp), dimension(n), intent(inout) :: u1, u2, u3
    type(c_ptr), intent(inout) :: scatter_event
    type(c_ptr) :: sgs_d, col_d, col_event
    integer(c_intptr_t) :: sv_addr, off_bytes
    integer(c_size_t) :: colbytes
    real(c_rp) :: rp_dummy
    logical :: on_host

    on_host = .true.
    sgs_d = C_NULL_PTR
    select type (b => gs%bcknd)
    type is (gs_device_t)
       on_host = b%shared_on_host
       sgs_d = b%shared_gs_d
    end select

    sv_addr = transfer(gs%shared_gs_v_d, sv_addr)
    colbytes = c_sizeof(rp_dummy) * int(l, c_size_t)
    off_bytes = int(l, c_intptr_t) * int(c_sizeof(rp_dummy), c_intptr_t)

    ! With a host-mirrored shared buffer, each scatter below issues an
    ! asynchronous host-to-device copy of shared_gs; a null event makes the
    ! scatter sync so the next column may safely overwrite the host buffer.
    ! Device-resident staging is stream-ordered and carries the caller's
    ! event.
    if (on_host) then
       col_event = C_NULL_PTR
    else
       col_event = scatter_event
    end if

    if (pe_size .gt. 1 .and. n .gt. 0) then
       call gs%comm%nbrecv_vec(tid, nc)

       ! Gather each component into the backend's shared buffer, then stage
       ! it into its column of shared_gs_v.
       call gs%bcknd%gather(gs%shared_gs, l, so, gs%shared_dof_gs, u1, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
            gs%shared_blk_off, op, .true.)
       if (on_host) then
          ! The gather mirrored the shared buffer to the host (synchronous).
          gs%shared_gs_v(1:l) = gs%shared_gs(1:l)
       else
          ! shared_gs_d is created lazily on the first gather.
          if (.not. c_associated(sgs_d)) then
             select type (b => gs%bcknd)
             type is (gs_device_t)
                sgs_d = b%shared_gs_d
             end select
          end if
          col_d = transfer(sv_addr, col_d)
          call device_memcpy(col_d, sgs_d, colbytes, DEVICE_TO_DEVICE, &
               sync = .false., strm = gs%bcknd%gs_stream)
       end if

       call gs%bcknd%gather(gs%shared_gs, l, so, gs%shared_dof_gs, u2, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
            gs%shared_blk_off, op, .true.)
       if (on_host) then
          gs%shared_gs_v(l + 1:2*l) = gs%shared_gs(1:l)
       else
          col_d = transfer(sv_addr + off_bytes, col_d)
          call device_memcpy(col_d, sgs_d, colbytes, DEVICE_TO_DEVICE, &
               sync = .false., strm = gs%bcknd%gs_stream)
       end if

       call gs%bcknd%gather(gs%shared_gs, l, so, gs%shared_dof_gs, u3, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
            gs%shared_blk_off, op, .true.)
       if (on_host) then
          gs%shared_gs_v(2*l + 1:3*l) = gs%shared_gs(1:l)
       else
          col_d = transfer(sv_addr + 2_c_intptr_t*off_bytes, col_d)
          call device_memcpy(col_d, sgs_d, colbytes, DEVICE_TO_DEVICE, &
               sync = .false., strm = gs%bcknd%gs_stream)
          ! Re-record the gather event so it covers the column copies above.
          ! Comm backends that order their per-peer packing streams on this
          ! event (NCCL, NVSHMEM) would otherwise race with the copies; the
          ! device MPI backend packs on gs_stream itself and is unaffected.
          call device_event_record(gs%bcknd%gather_event, gs%bcknd%gs_stream)
       end if

       call gs%comm%nbsend_vec(gs%shared_gs_v, l, nc, tid, &
            gs%bcknd%gather_event, gs%bcknd%gs_stream)
    end if

    ! Local gather-scatter per component.
    call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u1, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, gs%local_blk_off, &
         op, .false.)
    call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u1, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, gs%local_blk_off, &
         .false., C_NULL_PTR)
    call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u2, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, gs%local_blk_off, &
         op, .false.)
    call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u2, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, gs%local_blk_off, &
         .false., C_NULL_PTR)
    call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u3, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, gs%local_blk_off, &
         op, .false.)
    call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u3, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, gs%local_blk_off, &
         .false., C_NULL_PTR)

    ! Wait for the fused exchange (reduces into shared_gs_v), then stage each
    ! column back into the shared buffer and scatter.
    if (pe_size .gt. 1 .and. n .gt. 0) then
       call gs%comm%nbwait_vec(gs%shared_gs_v, l, nc, op, gs%bcknd%gs_stream)

       if (on_host) then
          gs%shared_gs(1:l) = gs%shared_gs_v(1:l)
       else
          col_d = transfer(sv_addr, col_d)
          call device_memcpy(sgs_d, col_d, colbytes, DEVICE_TO_DEVICE, &
               sync = .false., strm = gs%bcknd%gs_stream)
       end if
       call gs%bcknd%scatter(gs%shared_gs, l, gs%shared_dof_gs, u1, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
            gs%shared_blk_off, .true., col_event)

       if (on_host) then
          gs%shared_gs(1:l) = gs%shared_gs_v(l + 1:2*l)
       else
          col_d = transfer(sv_addr + off_bytes, col_d)
          call device_memcpy(sgs_d, col_d, colbytes, DEVICE_TO_DEVICE, &
               sync = .false., strm = gs%bcknd%gs_stream)
       end if
       call gs%bcknd%scatter(gs%shared_gs, l, gs%shared_dof_gs, u2, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
            gs%shared_blk_off, .true., col_event)

       if (on_host) then
          gs%shared_gs(1:l) = gs%shared_gs_v(2*l + 1:3*l)
       else
          col_d = transfer(sv_addr + 2_c_intptr_t*off_bytes, col_d)
          call device_memcpy(sgs_d, col_d, colbytes, DEVICE_TO_DEVICE, &
               sync = .false., strm = gs%bcknd%gs_stream)
       end if
       call gs%bcknd%scatter(gs%shared_gs, l, gs%shared_dof_gs, u3, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
            gs%shared_blk_off, .true., col_event)
    end if

  end subroutine gs_op_r3_device

end module gather_scatter
