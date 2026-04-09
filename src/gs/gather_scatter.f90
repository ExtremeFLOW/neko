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
       NEKO_BCKND_CUDA, NEKO_BCKND_OPENCL, NEKO_DEVICE_MPI
  use gs_bcknd, only : gs_bcknd_t, GS_BCKND_CPU, GS_BCKND_SX, GS_BCKND_DEV
  use gs_device, only : gs_device_t
  use gs_sx, only : gs_sx_t
  use gs_cpu, only : gs_cpu_t
  use gs_ops, only : GS_OP_ADD, GS_OP_MAX, GS_OP_MIN, GS_OP_MUL
  use gs_comm, only : gs_comm_t, GS_COMM_MPI, GS_COMM_MPIGPU, GS_COMM_NCCL, &
       GS_COMM_NVSHMEM
  use gs_mpi, only : gs_mpi_t
  use gs_interp, only : gs_interp_t
  use gs_interp_cpu, only : gs_interp_cpu_t
  use gs_device_mpi, only : gs_device_mpi_t
  use gs_device_nccl, only : gs_device_nccl_t
  use gs_device_shmem, only : gs_device_shmem_t
  use math, only : swap, sort, sort_tuple
  use mesh_conn, only : mesh_conn_obj_t
  use comm, only : pe_rank, pe_size, NEKO_COMM
  use mpi_f08, only : MPI_Reduce, MPI_Allreduce, MPI_Barrier, MPI_IN_PLACE, &
       MPI_Wait, MPI_Irecv, MPI_Isend, MPI_Wtime, MPI_SUM, MPI_MAX, &
       MPI_INTEGER, MPI_INTEGER2, MPI_INTEGER8, MPI_Request, MPI_Status, &
       MPI_STATUS_IGNORE, MPI_Get_Count
  use dofmap, only : dofmap_t
  use field, only : field_t
  use num_types, only : rp, dp, i2, i8
  use htable, only : htable_i8_t, htable_iter_i8_t
  use stack, only : stack_i4_t
  use utils, only : neko_error, linear_index
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use profiler, only : profiler_start_region, profiler_end_region
  use device, only : device_memcpy, HOST_TO_DEVICE, device_sync, device_free, &
       device_map, device_deassociate
  use amr_reconstruct, only : amr_reconstruct_t
  use amr_restart_component, only : amr_restart_component_t
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  type, public, extends(amr_restart_component_t) :: gs_t
     real(kind=rp), allocatable :: local_gs(:) !< Buffer for local gs-ops
     integer, allocatable :: local_dof_gs(:) !< Local dof to gs mapping
     integer, allocatable :: local_gs_dof(:) !< Local gs to dof mapping
     integer, allocatable :: local_blk_len(:) !< Local non-facet blocks
     integer, allocatable :: local_blk_off(:) !< Local non-facet blocks offset
     real(kind=rp), allocatable :: shared_gs(:) !< Buffer for shared gs-op
     integer, allocatable :: shared_dof_gs(:) !< Shared dof to gs map.
     integer, allocatable :: shared_gs_dof(:) !< Shared gs to dof map.
     integer, allocatable :: shared_blk_len(:) !< Shared non-facet blocks
     integer, allocatable :: shared_blk_off(:) !< Shared non-facet blocks offset
     type(dofmap_t), pointer ::dofmap !< Dofmap for gs-ops
     integer :: nlocal !< Local gs-ops
     integer :: nshared !< Shared gs-ops
     integer :: nlocal_blks !< Number of local blks
     integer :: nshared_blks !< Number of shared blks
     integer :: local_facet_offset !< offset for loc. facets
     integer :: shared_facet_offset !< offset for shr. facets
     class(gs_bcknd_t), allocatable :: bcknd !< Gather-scatter backend
     class(gs_comm_t), allocatable :: comm !< Comm. method
     class(gs_interp_t), allocatable :: interp !< Face/edge interpolation
     ! no longer needed for AMR version
     type(htable_i8_t) :: shared_dofs !< Htable of shared dofs
   contains
     procedure, private, pass(gs) :: gs_op_fld
     procedure, private, pass(gs) :: gs_op_r4
     procedure, pass(gs) :: gs_op_vector
     procedure, pass(gs) :: init => gs_init
     procedure, pass(gs) :: free => gs_free
     generic :: op => gs_op_fld, gs_op_r4, gs_op_vector
     !> AMR restart
     procedure, pass(this) :: amr_restart => gs_amr_restart
  end type gs_t

  ! Expose available gather-scatter operation
  public :: GS_OP_ADD, GS_OP_MUL, GS_OP_MIN, GS_OP_MAX

  ! Expose available gather-scatter backends
  public :: GS_BCKND_CPU, GS_BCKND_SX, GS_BCKND_DEV

  ! Expose available gather-scatter comm. backends
  public :: GS_COMM_MPI, GS_COMM_MPIGPU, GS_COMM_NCCL, GS_COMM_NVSHMEM


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
    integer :: i, j, ierr, bcknd_, comm_bcknd_
    integer(i8) :: glb_nshared, glb_nlocal
    logical :: use_device_mpi, use_device_nccl, use_device_shmem, use_host_mpi
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
          use_device_shmem = .true.
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
    else
       if (NEKO_DEVICE_MPI) then
          comm_bcknd_ = GS_COMM_MPIGPU
          use_device_mpi = .true.
       else
          comm_bcknd_ = GS_COMM_MPI
       end if
    end if

    select case (comm_bcknd_)
    case (GS_COMM_MPI)
       call neko_log%message('Comm         :          MPI')
       allocate(gs_mpi_t::gs%comm)
    case (GS_COMM_MPIGPU)
       call neko_log%message('Comm         :   Device MPI')
       allocate(gs_device_mpi_t::gs%comm)
    case (GS_COMM_NCCL)
       call neko_log%message('Comm         :         NCCL')
       allocate(gs_device_nccl_t::gs%comm)
    case (GS_COMM_NVSHMEM)
       call neko_log%message('Comm         :      NVSHMEM')
       allocate(gs_device_shmem_t::gs%comm)
    case default
       call neko_error('Unknown Gather-scatter comm. backend')
    end select
    ! Initialize a stack for each rank containing which dofs to send/recv at
    ! that rank
    call gs%comm%init_dofs()

    ! Initialise mapping/scheduling using connectivitu information
    ! assumption lx = ly = lz
    if (gs%dofmap%Xh%lx .ne. gs%dofmap%Xh%ly .or. &
         gs%dofmap%Xh%lx .ne. gs%dofmap%Xh%lz) &
         call neko_error('gs_init: inconsistent polynomial order')
    call gs_init_mapping_schedule(gs)

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

    ! Setup Gather-scatter backend and interpolation
    select case (bcknd_)
    case (GS_BCKND_CPU)
       allocate(gs_cpu_t::gs%bcknd)
       allocate(gs_interp_cpu_t::gs%interp)
       bcknd_str = '         std'
    case (GS_BCKND_DEV)
       allocate(gs_device_t::gs%bcknd)
       call neko_error('gs_init: device interpolation not done yet')
       if (NEKO_BCKND_HIP .eq. 1) then
          bcknd_str = '         hip'
       else if (NEKO_BCKND_CUDA .eq. 1) then
          bcknd_str = '        cuda'
       else if (NEKO_BCKND_OPENCL .eq. 1) then
          bcknd_str = '      opencl'
       end if
    case (GS_BCKND_SX)
       allocate(gs_sx_t::gs%bcknd)
       call neko_error('gs_init: sx interpolation not done yet')
       bcknd_str = '          sx'
    case default
       call neko_error('Unknown Gather-scatter backend/interpolation')
    end select

    write(log_buf, '(A)') 'Backend/interpolation : ' // trim(bcknd_str)
    call neko_log%message(log_buf)

    call gs%bcknd%init(gs%nlocal, gs%nshared, gs%nlocal_blks, gs%nshared_blks)
    call gs%interp%init(gs%dofmap%Xh%lx, gs%dofmap%msh%conn)

    if (use_device_mpi .or. use_device_nccl .or. use_device_shmem) then
       select type (b => gs%bcknd)
       type is (gs_device_t)
          b%shared_on_host = .false.
       end select
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
                call gs_op_vector(gs, tmp, dofmap%size(), GS_OP_ADD)

                do i = 1, size(strtgy)
                   c%nb_strtgy = strtgy(i)
                   call device_sync
                   call MPI_Barrier(NEKO_COMM)
                   strtgy_time(i) = MPI_Wtime()
                   do j = 1, 100
                      call gs_op_vector(gs, tmp, dofmap%size(), GS_OP_ADD)
                   end do
                   strtgy_time(i) = (MPI_Wtime() - strtgy_time(i)) / 100d0
                end do

                call device_deassociate(tmp)
                call device_free(tmp_d)
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

    call neko_log%end_section()

  end subroutine gs_init

  !> Deallocate a gather-scatter kernel
  subroutine gs_free(gs)
    class(gs_t), intent(inout) :: gs

    nullify(gs%dofmap)

    if (allocated(gs%local_gs)) then
       deallocate(gs%local_gs)
    end if

    if (allocated(gs%local_dof_gs)) then
       deallocate(gs%local_dof_gs)
    end if

    if (allocated(gs%local_gs_dof)) then
       deallocate(gs%local_gs_dof)
    end if

    if (allocated(gs%local_blk_len)) then
       deallocate(gs%local_blk_len)
    end if

    if (allocated(gs%local_blk_off)) then
       deallocate(gs%local_blk_off)
    end if

    if (allocated(gs%shared_gs)) then
       deallocate(gs%shared_gs)
    end if

    if (allocated(gs%shared_dof_gs)) then
       deallocate(gs%shared_dof_gs)
    end if

    if (allocated(gs%shared_gs_dof)) then
       deallocate(gs%shared_gs_dof)
    end if

    if (allocated(gs%shared_blk_len)) then
       deallocate(gs%shared_blk_len)
    end if

    if (allocated(gs%shared_blk_off)) then
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

    if (allocated(gs%interp)) then
       call gs%interp%free()
       deallocate(gs%interp)
    end if

    call gs%free_amr_base()

  end subroutine gs_free

  !> Initialise mapping/scheduling using connectivity
  subroutine gs_init_mapping_schedule(gs)
    type(gs_t), target, intent(inout) :: gs
    integer :: lx
    integer :: nvrt_loc, nfcs_loc, nedg_loc, nlvrt_dof, nlfcs_dof, nledg_dof, &
         nsvrt_dof, nsfcs_dof, nsedg_dof, nlfcs_ncon, nsfcs_ncon, &
         nlfcs_ncon_dof, nsfcs_ncon_dof
    integer, dimension(:), allocatable :: vrt_mult, fcs_mult, edg_mult, &
         fcs_mult_glb, vrt_shr, fcs_shr, edg_shr, vrt_loc, fcs_loc, edg_loc
    integer :: il ,jl, kl, ll, ml, id
    integer, allocatable, dimension(:, :) :: ngh_src, ngh_dst
    integer, parameter :: lda1 = 2 ! tuple length
    integer, dimension(lda1) :: aa1 ! tmp array for sorting
    integer, dimension(:), allocatable :: ind_src, ind_dst
    integer, parameter :: nkey = 1 ! max. number of keys
    integer, dimension(nkey) :: key
    type(stack_i4_t) :: send_pe, recv_pe

    ! Number of degrees of freedom; assumption lx=ly=lz
    if (gs%dofmap%Xh%lx .ne. gs%dofmap%Xh%ly .or. &
         gs%dofmap%Xh%lx .ne. gs%dofmap%Xh%lz) &
         call neko_error('gs_init_mapping_schedule: inconsistent polynomial &
         &order')
    lx = gs%dofmap%Xh%lx

    associate (vrt => gs%dofmap%msh%conn%vrt, fcs => gs%dofmap%msh%conn%fcs, &
         edg => gs%dofmap%msh%conn%edg, conn => gs%dofmap%msh%conn)

      ! local multiplicity, shared/local object lists and object count including
      ! multiplicity
      ! shared objects lists exists already but is not sorted with respect to
      ! their global index
      ! vertices
      call gs_size_list_get(vrt, vrt_mult, vrt_loc, nvrt_loc,  nlvrt_dof, &
           vrt_shr, nsvrt_dof)

      ! edges
      call gs_size_list_get(edg, edg_mult, edg_loc, nedg_loc, &
           nledg_dof, edg_shr, nsedg_dof)

      ! faces; added global multiplicity
      call gs_size_list_get(fcs, fcs_mult, fcs_loc, nfcs_loc, &
           nlfcs_dof, fcs_shr, nsfcs_dof, fcs_mult_glb)

      ! count of all degrees of freedom including multiplicity
      gs%nlocal = nlvrt_dof + nledg_dof * (lx - 2) + nlfcs_dof * (lx - 2) * &
           (lx - 2)
      gs%nshared = nsvrt_dof + nsedg_dof * (lx - 2) + nsfcs_dof * (lx - 2) * &
           (lx - 2)

      ! Multiplicity of vertices and edges is neglected, but for faces it has
      ! to be taken into account to set up properly facet offset, that is used
      ! to flag array section that can be used for vectorisation. In conforming
      ! case all local faces have multiplicity 2 and shared one multiplicity 1,
      ! so are automatically marked for vectorisation. In AMR case hanging faces
      ! have multiplicity higher than 2. Flag those objects using the global
      ! multiplicity information and reorder shared/local lists.
      ! This influences block lists.
      ! conforming case; vertices and edges only
      gs%local_facet_offset = nlvrt_dof + nledg_dof * (lx - 2) + 1
      gs%shared_facet_offset = nsvrt_dof + nsedg_dof * (lx - 2) + 1
      ! numbers of blocks per degree of freedom
      gs%nlocal_blks = nvrt_loc + nedg_loc * (lx - 2)
      gs%nshared_blks = vrt%nshare + edg%nshare * (lx - 2)
      ! number of nonconforming faces
      nlfcs_ncon = 0
      nsfcs_ncon = 0
      ! number of nonconforming faces including multiplicity
      nlfcs_ncon_dof = 0
      nsfcs_ncon_dof = 0
      ! AMR case; find nonconforming faces
      ! local
      ! find max multiplicity for local nodes
      id = 0
      do il = 1, nfcs_loc
         id = max(id, fcs_mult(fcs_loc(il)))
      end do
      if (id .gt. 2) then
         ! local
!         call neko_error('gs_init_mapping_schedule: not done yet')
            

         gs%local_facet_offset = gs%local_facet_offset + nlfcs_ncon_dof * &
              (lx - 2) * (lx - 2)
         gs%nlocal_blks = gs%nlocal_blks + nlfcs_ncon * (lx - 2) * (lx - 2)
      end if
      ! shared
      ! For shared dof it is necessary to use global information not to
      ! destroy unique communication order.
      if (fcs%nshare .gt. 0) then
         if (maxval(fcs_mult_glb) .gt. 2) then
!            call neko_error('gs_init_mapping_schedule: not done yet')
            

            gs%shared_facet_offset = gs%shared_facet_offset + nsfcs_ncon_dof * &
                 (lx - 2) * (lx - 2)
            gs%nshared_blks = gs%nshared_blks + nsfcs_ncon * (lx - 2) * (lx - 2)
         end if
      end if

      ! allocate and fill arrays
      ! local
      if (gs%nlocal .gt. 0) then
         allocate(gs%local_dof_gs(gs%nlocal), gs%local_gs_dof(gs%nlocal), &
              gs%local_gs(gs%nlocal), gs%local_blk_len(gs%nlocal_blks), &
              gs%local_blk_off(gs%nlocal_blks))
         call gs_fill_arrays(lx, nvrt_loc, nedg_loc, nfcs_loc, nlfcs_ncon, &
              vrt_mult, edg_mult, fcs_mult, vrt_loc, edg_loc, fcs_loc, &
              vrt%lmap, vrt%lmapoff, edg%lmap, edg%lmapoff, edg%algn, &
              fcs%lmap, fcs%lmapoff, fcs%algn, gs%local_blk_len, &
              gs%local_blk_off, gs%local_dof_gs, gs%local_gs_dof, &
              gs%nlocal_blks, gs%local_facet_offset)
      end if

      ! shared
      if (gs%nshared .gt. 0) then
         allocate(gs%shared_dof_gs(gs%nshared), gs%shared_gs_dof(gs%nshared), &
              gs%shared_gs(gs%nshared), gs%shared_blk_len(gs%nshared_blks), &
              gs%shared_blk_off(gs%nshared_blks))
         call gs_fill_arrays(lx, vrt%nshare, edg%nshare, fcs%nshare, &
              nsfcs_ncon, vrt_mult, edg_mult, fcs_mult, vrt_shr, edg_shr, &
              fcs_shr, vrt%lmap, vrt%lmapoff, edg%lmap, edg%lmapoff, edg%algn, &
              fcs%lmap, fcs%lmapoff, fcs%algn, gs%shared_blk_len, &
              gs%shared_blk_off, gs%shared_dof_gs, gs%shared_gs_dof, &
              gs%nshared_blks, gs%shared_facet_offset)

         ! sanity check; gs%shared_dof_gs should be sorted inside the sections
         ! 1 : gs%shared_facet_offset - 1, and
         ! gs%shared_facet_offset : gs%nshared
         ! There may be no faces for coarse grid solver
         if (gs%shared_facet_offset .gt. gs%nshared) then
            do il = 1, gs%nshared - 1
               if (gs%shared_dof_gs(il) .gt. gs%shared_dof_gs(il + 1)) &
                    call neko_error('gs%shared_dof_gs not ordered 1')
            end do
         else
            do il = 1, gs%shared_facet_offset - 2
               if (gs%shared_dof_gs(il) .gt. gs%shared_dof_gs(il + 1)) &
                    call neko_error('gs%shared_dof_gs not ordered 2')
            end do
            do il = gs%shared_facet_offset, gs%nshared - 1
               if (gs%shared_dof_gs(il) .gt. gs%shared_dof_gs(il + 1)) &
                    call neko_error('gs%shared_dof_gs not ordered 3')
            end do
         end if

         ! Initialise communication arrays
         ! This must be done consistently on all the ranks to avoid node
         ! mismatch. To get it right all the shared objects were sorted
         ! according to the global object id and the faces were reshuffled
         ! according to the global multiplicity number. As the operation
         ! relies on object data not dof, the point numbering must be
         ! consistent with the one in gs_fill_arrays, but object alignment
         ! shouldn't matter this time. As object lists are sorted, send_dof
         ! and recv_dof arrays are identical.
         ! loop order matters
         ! vertices
         if (vrt%nshare .gt. 0) then
            do il = 1, vrt%nshare ! all shared
               id = il
               do jl = 1, vrt%nrank ! ranks
                  do kl = vrt%rankoff(jl), vrt%rankoff(jl + 1) - 1 ! rank shared
                     ! rank shared list can have different order than sorted one
                     if (vrt%rankshare(kl) .eq. vrt_shr(il)) then
                        call gs%comm%send_dof(vrt%rank(jl))%push(id)
                        call gs%comm%recv_dof(vrt%rank(jl))%push(id)
                        exit
                     end if
                  end do
               end do
            end do
         end if

         ! edges
         if (edg%nshare .gt. 0) then
            do il = 1, edg%nshare ! all shared
               id = vrt%nshare + (il - 1) * (lx - 2)
               do jl = 1, edg%nrank ! ranks
                  do kl = edg%rankoff(jl), edg%rankoff(jl + 1) - 1 ! rank shared
                     ! rank shared list can have different order than sorted one
                     if (edg%rankshare(kl) .eq. edg_shr(il)) then
                        do ll = 1, lx - 2
                           call gs%comm%send_dof(edg%rank(jl))%push(id + ll)
                           call gs%comm%recv_dof(edg%rank(jl))%push(id + ll)
                        end do
                        exit
                     end if
                  end do
               end do
            end do
         end if

         ! faces
         if (fcs%nshare .gt. 0) then
            do il = 1, fcs%nshare ! all shared
               id = vrt%nshare + edg%nshare * (lx - 2) + &
                    (il - 1) * (lx - 2) * (lx - 2)
               do jl = 1, fcs%nrank ! ranks
                  do kl = fcs%rankoff(jl), fcs%rankoff(jl + 1) - 1 ! rank shared
                     ! rank shared list can have different order than sorted one
                     if (fcs%rankshare(kl) .eq. fcs_shr(il)) then
                        do ll = 1, lx - 2
                           do ml = 1, lx - 2
                              call gs%comm%send_dof(fcs%rank(jl))%push(id + &
                                   (ll - 1) * (lx - 2) + ml)
                              call gs%comm%recv_dof(fcs%rank(jl))%push(id + &
                                   (ll - 1) * (lx - 2) + ml)
                           end do
                        end do
                        exit
                     end if
                  end do
               end do
            end do
         end if

         ! rank list
         ! Build send/receive list based on the vertex neighbour list.
         ! What follows recreates the send/receive order from the ring,
         ! but it may be not needed.
         if (vrt%nshare .gt. 0) then
            allocate(ngh_src(lda1, vrt%nrank), ind_src(vrt%nrank), &
                 ngh_dst(lda1, vrt%nrank), ind_dst(vrt%nrank))
            do il = 1, vrt%nrank
               ngh_src(1, il) = mod(pe_rank - vrt%rank(il) + pe_size, pe_size)
               ngh_src(2, il) = il
               ngh_dst(1, il) = mod(vrt%rank(il) - pe_rank + pe_size, pe_size)
               ngh_dst(2, il) = il
            end do
            il = vrt%nrank
            key(1) = 1
            call sort_tuple(ngh_src, lda1, il, key, 1, ind_src, aa1)
            call sort_tuple(ngh_dst, lda1, il, key, 1, ind_dst, aa1)

            call send_pe%init(vrt%nrank)
            call recv_pe%init(vrt%nrank)
            do il = 1, vrt%nrank
               call recv_pe%push(vrt%rank(ngh_src(2, il)))
               call send_pe%push(vrt%rank(ngh_dst(2, il)))
            end do
            deallocate(ngh_src, ind_src, ngh_dst, ind_dst)

            call gs%comm%init(send_pe, recv_pe)

            call send_pe%free()
            call recv_pe%free()
         end if
      end if

      deallocate(vrt_mult, fcs_mult, edg_mult, vrt_loc, fcs_loc, edg_loc)
      if (allocated(vrt_shr)) deallocate(vrt_shr)
      if (allocated(edg_shr)) deallocate(edg_shr)
      if (allocated(fcs_shr)) deallocate(fcs_shr)
      if (allocated(fcs_mult_glb)) deallocate(fcs_mult_glb)

    end associate

  end subroutine gs_init_mapping_schedule

  !> Get local object multiplicity and list of shared/local objects
  !! @param[in]    obj       object connectivity info
  !! @param[out]   loc_mult  object's local multiplicity
  !! @param[out]   loc_list  list of local objects
  !! @param[out]   nlist     number of local objects
  !! @param[out]   nlocal    number of local objects including multiplicity
  !! @param[out]   shr_list  sorted list of shared objects
  !! @param[out]   nshared   number of shared objects including multiplicity
  !! @param[out]   glb_mult  object's global multiplicity
  subroutine gs_size_list_get(obj, loc_mult, loc_list, nlist, nlocal, &
       shr_list, nshared, glb_mult)
    type(mesh_conn_obj_t), intent(in) :: obj
    integer, dimension(:), allocatable, intent(out):: loc_mult, shr_list, &
         loc_list
    integer, dimension(:), allocatable, intent(out), optional:: glb_mult
    integer, intent(out) :: nlist, nlocal, nshared
    integer :: il, jl
    integer(i8), dimension(:), allocatable :: gidx
    integer, dimension(:), allocatable :: ind

    nlocal = 0
    nshared = 0
    ! local
    nlist = obj%lnum - obj%nshare
    allocate(loc_mult(obj%lnum), loc_list(nlist))
    jl = 0
    do il = 1, obj%lnum
       loc_mult(il) = obj%lmapoff(il + 1) - obj%lmapoff(il)
       ! local object
       if (.not. obj%lshare(il)) then
          ! local only object with no neighbour can be skipped
          if (loc_mult(il) .gt. 1) then
             nlocal = nlocal + loc_mult(il) ! multiplicity count
             jl = jl + 1 ! object count
             loc_list(jl) = il
          end if
       end if
    end do

    ! sanity check
    if (jl .gt. nlist) &
         call neko_error('gs_size_list_get: inconsistent local object number')
    nlist = jl

    ! shared
    ! for communication consistency shared objects must be sorted according
    ! to their global id
    if (obj%nshare .gt. 0) then
       allocate(gidx(obj%nshare), ind(obj%nshare), shr_list(obj%nshare))
       shr_list(:) = obj%sharelist(:)
       do il = 1, obj%nshare
          jl = shr_list(il)
          gidx(il) = obj%gidx(jl)
       end do
       !
       call sort(gidx, ind, obj%nshare)
       call swap(shr_list, ind, obj%nshare)
       !
       deallocate(gidx, ind)
       ! all shared objects must be included
       do il = 1, obj%nshare
          jl = shr_list(il)
          nshared = nshared + loc_mult(jl) ! multiplicity count
       end do
       ! provide global object multiplicity if requested
       if (present(glb_mult)) then
          allocate(glb_mult(obj%lnum))
          glb_mult(:) = loc_mult(:)
          do il = 1, obj%nshare
             glb_mult(obj%sharelist(il)) = glb_mult(obj%sharelist(il)) + &
                  obj%gmapoff(il + 1) - obj%gmapoff(il)
          end do
       end if
    else
       nshared = 0
    end if
  end subroutine gs_size_list_get

  !> Get local object multiplicity and list of shared/local objects
  !! @param[in]    lx        number of points in 1D
  !! @param[in]    nvrt      number of vertices
  !! @param[in]    nedg      number of edges
  !! @param[in]    nfcs      number of faces
  !! @param[in]    nfcsn     number of nonconforming faces
  !! @param[in]    vrt_mult  vertex multiplicity
  !! @param[in]    edg_mult  edge multiplicity
  !! @param[in]    fcs_mult  face multiplicity
  !! @param[in]    vrt       vertex list
  !! @param[in]    edg       edge list
  !! @param[in]    fcs       face list
  !! @param[in]    vrt_map   vertex to element mapping
  !! @param[in]    vrt_off   vertex to element mapping offset
  !! @param[in]    edg_map   edge to element mapping
  !! @param[in]    edg_off   edge to element mapping offset
  !! @param[in]    edg_algn  edge alignment
  !! @param[in]    fcs_map   face to element mapping
  !! @param[in]    fcs_off   face to element mapping offset
  !! @param[in]    fcs_algn  face alignment
  !! @param[out]   blk_len   block length
  !! @param[out]   blk_off   block offset
  !! @param[out]   mdg       dof to gs mapping
  !! @param[out]   mgd       gs to dof mapping
  !! @param[in]    nblks     number of blocks
  !! @param[in]    cfcs_off  conforming facet offset
  subroutine gs_fill_arrays(lx, nvrt, nedg, nfcs, nfcsn, vrt_mult, edg_mult, &
       fcs_mult, vrt, edg, fcs, vrt_map, vrt_off, edg_map, edg_off, edg_algn, &
       fcs_map, fcs_off, fcs_algn, blk_len, blk_off, mdg, mgd, nblks, cfcs_off)
    integer, intent(in) :: lx, nvrt, nedg, nfcs, nfcsn, nblks, cfcs_off
    integer, dimension(:), intent(in) :: vrt_mult, edg_mult, fcs_mult, &
         vrt, edg, fcs, vrt_off, edg_off, fcs_off
    integer, dimension(:, :), intent(in) :: vrt_map, edg_map, edg_algn, &
         fcs_map, fcs_algn
    integer, dimension(:), intent(out) :: blk_len, blk_off, mdg, mgd
    integer :: lxl, il, jl, kl, ll, ml, itmp, mult, iel, ipos, algn, el_off, &
         edg_id, fcs_id, id, lid
    integer :: ix, iy, iz, nx, ny, nz, nxyz
    integer, dimension(3) :: edg_stride
    integer, dimension(12) :: edg_start
    integer, dimension(6) :: fcs_stride, fcs_strider, fcs_start

    ! block length
    lxl = lx -2
    ! vertex
    do il = 1, nvrt
       blk_len(il) = vrt_mult(vrt(il))
    end do
    itmp = nvrt
    ! edge
    do il = 1, nedg
       mult = edg_mult(edg(il))
       do jl = 1, lxl
          itmp = itmp + 1
          blk_len(itmp) = mult
       end do
    end do
    ! faces; only those with multiplicity bigger than 2; nonconforming only
    do il = 1, nfcsn
       mult = fcs_mult(fcs(il))
       do jl = 1, lxl
          do kl = 1, lxl
             itmp = itmp + 1
             blk_len(itmp) = mult
          end do
       end do
    end do

    ! sanity checks
    if (itmp .ne. nblks) &
         call neko_error('gs_fill_arrays: inconsistent block number')

    ! offset
    if (itmp .gt. 0) then
       blk_off(1) = 0
       do il = 2, itmp
          blk_off(il) = blk_off(il - 1) + blk_len(il - 1)
       end do

       ! sanity checks
       if (blk_off(itmp) + blk_len(itmp) .ne. cfcs_off - 1) &
            call neko_error('gs_fill_arrays: inconsistent facet offset')
    end if

    ! mappings
    itmp = 0
    lxl = lx - 1
    ! vertex
    do il = 1, nvrt
       id = il
       ! offset array
       do jl = vrt_off(vrt(il)), vrt_off(vrt(il) + 1) - 1
          ! local element and vertex position
          iel = vrt_map(1, jl)
          ipos = vrt_map(2, jl)
          ! position in the element
          ix = mod(ipos - 1, 2) * lxl + 1
          iy = (mod(ipos - 1, 4)/2) * lxl + 1
          iz = ((ipos - 1)/4) * lxl + 1
          ! linear id
          lid = linear_index(ix, iy, iz, iel, lx, lx, lx)
          itmp = itmp + 1
          ! dof => gs
          mdg(itmp) = id
          ! gs => dof
          mgd(itmp) = lid
       end do
    end do

    ! edge
    ! data position in the element; assuming lx = ly = lz
    nx = lx
    ny = lx
    nz = lx
    nxyz = nx * ny * nz
    ! data stride
    edg_stride(1)  = 1
    edg_stride(2)  = nx
    edg_stride(3)  = nx * ny
    ! data start
    edg_start(1)  = 1
    edg_start(2)  = nx * (ny - 1) + 1
    edg_start(3)  = nx * ny * (nz - 1) + 1
    edg_start(4)  = nx * (ny * nz - 1) + 1
    edg_start(5)  = 1
    edg_start(6)  = nx
    edg_start(7)  = nx * ny * (nz - 1) + 1
    edg_start(8)  = nx * ny * (nz - 1) + nx
    edg_start(9)  = 1
    edg_start(10) = nx
    edg_start(11) = nx * (ny - 1) + 1
    edg_start(12) = nx * ny

    ! For a global communication it is important to order array entrances in
    ! such a way, that "id" is always growing withing the element despite of
    ! the alignment. Combining it with the proper object order gives global
    ! unique dof ordering.
    do il = 1, nedg
       id = nvrt + (il - 1) * (nx - 2)
       ! multiplicity
       mult = edg_mult(edg(il))
       do jl = 1, mult
          ! position in map array
          ml = jl + edg_off(edg(il)) - 1
          ! local element and vertex position
          iel = edg_map(1, ml)
          ipos = edg_map(2, ml)
          el_off = (iel - 1) * nxyz + edg_start(ipos)
          edg_id = (ipos - 1)/4 + 1
          ! alignment
          algn = edg_algn(ipos, iel)
          select case (algn)
          case (0) ! identity
             ! points at the edge
             do concurrent (kl = 1: nx - 2)
                ! linear id
                lid = el_off + kl * edg_stride(edg_id)
                ! dof => gs
                mdg(itmp + (kl - 1) * mult + jl) = id + kl
                ! gs => dof
                mgd(itmp + (kl - 1) * mult + jl) = lid
             end do
          case (1) ! permutation
             ! points at the edge
             do concurrent (kl = 1: nx - 2)
                ! linear id
                lid = el_off + kl * edg_stride(edg_id)
                ! dof => gs
                mdg(itmp + (nx - 1 - kl - 1) * mult + jl) = id + nx - 1 - kl
                ! gs => dof
                mgd(itmp + (nx - 1 - kl - 1) * mult + jl) = lid
             end do
          end select
       end do
       itmp = itmp + (nx -2) * mult
    end do

    ! face
    ! data stride; within 1D row
    fcs_stride(1) = nx
    fcs_stride(2) = nx
    fcs_stride(3) = 1
    fcs_stride(4) = 1
    fcs_stride(5) = 1
    fcs_stride(6) = 1
    ! data stride; between 1D rows
    fcs_strider(1) = nx * ny
    fcs_strider(2) = nx * ny
    fcs_strider(3) = nx * ny
    fcs_strider(4) = nx * ny
    fcs_strider(5) = nx
    fcs_strider(6) = nx
    ! data start
    fcs_start(1) = nx * ny + 1
    fcs_start(2) = nx * (ny + 1)
    fcs_start(3) = nx * ny + 1
    fcs_start(4) = nx * (2 * ny - 1) + 1
    fcs_start(5) = nx + 1
    fcs_start(6) = nx * ny * (nz - 1) + nx + 1

    ! For a global communication it is important to order array entrances in
    ! such a way, that "id" is always growing withing the element despite of
    ! the alignment. Combining it with the proper object order gives global
    ! unique dof ordering.
    do il = 1, nfcs
       id = nvrt + nedg * (nx - 2) + (il - 1) * (nx - 2) * (nx - 2)
       ! multiplicity
       mult = fcs_mult(fcs(il))
       do jl = 1, mult
          ! position in map array
          ml = jl + fcs_off(fcs(il)) - 1
          ! local element and vertex position
          iel = fcs_map(1, ml)
          ipos = fcs_map(2, ml)
          el_off = (iel - 1) * nxyz + fcs_start(ipos)
          fcs_id = (ipos - 1)/2 + 1
          ! alignment
          algn = fcs_algn(ipos, iel)
          select case (algn)
          case (0) ! identity
             ! points at the face
             do concurrent (kl = 1: nx - 2, ll = 1: nx - 2)
                ! linear id
                lid = el_off + (kl - 1) * fcs_strider(ipos) + &
                     ll * fcs_stride(ipos)
                ! dof => gs
                mdg(itmp + ((kl - 1) * (nx - 2) + &
                     (ll - 1))* mult + jl) = id + (kl - 1) * (nx - 2) + ll
                ! gs => dof
                mgd(itmp + ((kl - 1) * (nx - 2) + &
                     (ll - 1))* mult + jl) = lid
             end do
          case (1) ! transpose
             ! points at the face
             do concurrent (kl = 1: nx - 2, ll = 1: nx - 2)
                ! linear id
                lid = el_off + (kl - 1) * fcs_strider(ipos) + &
                     ll * fcs_stride(ipos)
                ! dof => gs
                mdg(itmp + ((ll - 1) * (nx - 2) + &
                     (kl - 1))* mult + jl) = id + (ll - 1) * (nx - 2) + kl
                ! gs => dof
                mgd(itmp + ((ll - 1) * (nx - 2) + &
                     (kl - 1))* mult + jl) = lid
             end do
          case (2) ! permutation in X
             ! points at the face
             do concurrent (kl = 1: nx - 2, ll = 1: nx - 2)
                ! linear id
                lid = el_off + (kl - 1) * fcs_strider(ipos) + &
                     ll * fcs_stride(ipos)
                ! dof => gs
                mdg(itmp + ((kl - 1) * (nx - 2) + &
                     (nx - 2 - ll))* mult + jl) = &
                     id + (kl - 1) * (nx - 2) + nx - 1 - ll
                ! gs => dof
                mgd(itmp + ((kl - 1) * (nx - 2) + &
                     (nx - 2 - ll))* mult + jl) = lid
             end do
          case (3) ! permutation in X; transpose; inverse of 4
             ! Be careful, as depending on perspective (element realisation
             ! or a reference one) this can turn into P_Y T. Moreover,
             ! take int account that P_X T = T P_Y.
             ! points at the face
             do concurrent (kl = 1: nx - 2, ll = 1: nx - 2)
                ! linear id
                lid = el_off + (kl - 1) * fcs_strider(ipos) + &
                     ll * fcs_stride(ipos)
                ! dof => gs
                mdg(itmp + ((ll - 1) * (nx - 2) + &
                     (nx - 2 - kl))* mult + jl) = &
                     id + (ll - 1) * (nx - 2) + nx - 1 - kl
                ! gs => dof
                mgd(itmp + ((ll - 1) * (nx - 2) + &
                     (nx - 2 - kl))* mult + jl) = lid
             end do
          case (4) ! permutation in Y; transpose; inverse of 3
             ! Be careful, as depending on perspective (element realisation
             ! or a reference one) this can turn into P_X T. Moreover,
             ! take int account that P_Y T = T P_X.
             ! points at the face
             do concurrent (kl = 1: nx - 2, ll = 1: nx - 2)
                ! linear id
                lid = el_off + (kl - 1) * fcs_strider(ipos) + &
                     ll * fcs_stride(ipos)
                ! dof => gs
                mdg(itmp + ((nx - 2 - ll) * (nx - 2) + &
                     (kl - 1))* mult + jl) = &
                     id + (nx - 2 - ll) * (nx - 2) + kl
                ! gs => dof
                mgd(itmp + ((nx - 2 - ll) * (nx - 2) + &
                     (kl - 1))* mult + jl) = lid
             end do
          case (5) ! permutation in Y
             ! points at the face
             do concurrent (kl = 1: nx - 2, ll = 1: nx - 2)
                ! linear id
                lid = el_off + (kl - 1) * fcs_strider(ipos) + &
                     ll * fcs_stride(ipos)
                ! dof => gs
                mdg(itmp + ((nx - 2 - kl) * (nx - 2) + &
                     (ll - 1))* mult + jl) = id + (nx - 2 - kl) * (nx - 2) + ll
                ! gs => dof
                mgd(itmp + ((nx - 2 - kl) * (nx - 2) + &
                     (ll - 1))* mult + jl) = lid
             end do
          case (6) ! permutation in Y; permutation in X; transpose
             ! points at the face
             do concurrent (kl = 1: nx - 2, ll = 1: nx - 2)
                ! linear id
                lid = el_off + (kl - 1) * fcs_strider(ipos) + &
                     ll * fcs_stride(ipos)
                ! dof => gs
                mdg(itmp + ((nx - 2 - ll) * (nx - 2) + &
                     (nx - 2 - kl))* mult + jl) = &
                     id + (nx - 2 - ll) * (nx - 2) + nx - 1 - kl
                ! gs => dof
                mgd(itmp + ((nx - 2 - ll) * (nx - 2) + &
                     (nx - 2 - kl))* mult + jl) = lid
             end do
          case (7) ! permutation in Y; permutation in X
             ! points at the face
             do concurrent (kl = 1: nx - 2, ll = 1: nx - 2)
                ! linear id
                lid = el_off + (kl - 1) * fcs_strider(ipos) + &
                     ll * fcs_stride(ipos)
                ! dof => gs
                mdg(itmp + ((nx - 2 - kl) * (nx - 2) + &
                     (nx - 2 - ll))* mult + jl) = &
                     id + (nx - 2 - kl) * (nx - 2) + nx - 1 - ll
                ! gs => dof
                mgd(itmp + ((nx - 2 - kl) * (nx - 2) + &
                     (nx - 2 - ll))* mult + jl) = lid
             end do
          end select
       end do
       itmp = itmp + (nx -2) * (nx -2) * mult
    end do

  end subroutine gs_fill_arrays

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
    integer :: m, l, op, lo, so

    lo = gs%local_facet_offset
    so = -gs%shared_facet_offset
    m = gs%nlocal
    l = gs%nshared

    call profiler_start_region("gather_scatter", 5)
    ! Gather shared dofs
    if (pe_size .gt. 1 .and. n .gt. 0) then
       call profiler_start_region("gs_nbrecv", 13)
       call gs%comm%nbrecv()
       call profiler_end_region("gs_nbrecv", 13)
       call profiler_start_region("gs_gather_shared", 14)
       call gs%bcknd%gather(gs%shared_gs, l, so, gs%shared_dof_gs, u, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, &
            gs%shared_blk_off, op, .true.)
       call profiler_end_region("gs_gather_shared", 14)
       call profiler_start_region("gs_nbsend", 6)
       call gs%comm%nbsend(gs%shared_gs, l, &
            gs%bcknd%gather_event, gs%bcknd%gs_stream)
       call profiler_end_region("gs_nbsend", 6)

    end if

    ! Gather-scatter local dofs
    if (m .gt. 0) then
       call profiler_start_region("gs_local", 12)
       call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u, n, &
            gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, &
            gs%local_blk_off, op, .false.)
       call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u, n, &
            gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, &
            gs%local_blk_off, .false., C_NULL_PTR)
       call profiler_end_region("gs_local", 12)
    end if

    ! Scatter shared dofs
    if (pe_size .gt. 1 .and. n .gt. 0) then
       call profiler_start_region("gs_nbwait", 7)
       call gs%comm%nbwait(gs%shared_gs, l, op, gs%bcknd%gs_stream)
       call profiler_end_region("gs_nbwait", 7)
       call profiler_start_region("gs_scatter_shared", 15)
       if (present(event)) then
          call gs%bcknd%scatter(gs%shared_gs, l,&
               gs%shared_dof_gs, u, n, &
               gs%shared_gs_dof, gs%nshared_blks, &
               gs%shared_blk_len, gs%shared_blk_off, .true., event)
       else
          call gs%bcknd%scatter(gs%shared_gs, l,&
               gs%shared_dof_gs, u, n, &
               gs%shared_gs_dof, gs%nshared_blks, &
               gs%shared_blk_len, gs%shared_blk_off, .true., C_NULL_PTR)
       end if
       call profiler_end_region("gs_scatter_shared", 15)
    end if

    call profiler_end_region("gather_scatter", 5)

  end subroutine gs_op_vector

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine gs_amr_restart(this, reconstruct, counter, tstep)
    class(gs_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep
    integer :: ierr
    integer(i8) :: glb_nlocal, glb_nshared
    character(len=LOG_SIZE) :: log_buf

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    if (this%dofmap%Xh%lx .lt. 1e1) then
       write(log_buf, '(A,I2)') 'Reconstructing Gather-Scatter; lx =  ', &
            this%dofmap%Xh%lx
    else if (this%dofmap%Xh%lx .lt. 1e2) then
       write(log_buf, '(A,I2)') 'Reconstructing Gather-Scatter; lx =  ', &
            this%dofmap%Xh%lx
    end if
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

    ! reconstruct dofmap; It is safe to call it here, as AMR restart prevents
    ! recursive reconstructions
    if (associated(this%dofmap)) &
         call this%dofmap%amr_restart(reconstruct, counter, tstep)

    ! clear space
    if (this%nlocal .gt. 0) then
       deallocate(this%local_dof_gs, this%local_gs_dof, this%local_blk_len, &
            this%local_blk_off, this%local_gs)
    end if
    if (this%nshared .gt. 0) then
       deallocate(this%shared_dof_gs, this%shared_gs_dof, this%shared_blk_len, &
            this%shared_blk_off, this%shared_gs)
    end if
    this%nlocal = 0
    this%local_facet_offset = 0
    this%nlocal_blks = 0
    this%nshared = 0
    this%shared_facet_offset = 0
    this%nshared_blks = 0

    ! reconstruct comm; for now just clearing stacks and deallocating arrays,
    ! as gs_init_mapping_schedule calls comm%init
    if (allocated(this%comm)) &
         call this%comm%amr_restart(reconstruct, counter, tstep)

    ! get new mapping/scheduling
    call gs_init_mapping_schedule(this)

    ! reconstruct interpolation
    call this%interp%amr_restart(reconstruct, counter, tstep)

    ! Global number of points not needing to be sent over mpi for gs operations
    ! "Internal points"
    glb_nlocal = int(this%nlocal, i8)
    ! Global number of points needing to be communicated with other pes/ranks
    ! "external points"
    glb_nshared = int(this%nshared, i8)
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
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)
    write(log_buf, '(A,I12)') 'Avg. external: ', glb_nshared/pe_size
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)

    ! PLACE FOR RECONSTRUCTING BACK-END
    ! for now do nothing, as cpu does not store anything

  end subroutine gs_amr_restart

end module gather_scatter
