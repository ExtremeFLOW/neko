! Copyright (c) 2020-2022, The Neko Authors
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
  use neko_config
  use gs_bcknd
  use gs_device
  use gs_sx
  use gs_cpu
  use gs_ops
  use gs_comm
  use gs_mpi
  use gs_device_mpi
  use mesh
  use dofmap
  use field
  use num_types
  use mpi_f08
  use htable
  use stack
  use utils
  use logger
  use profiler
  implicit none

  type gs_t
     real(kind=rp), allocatable :: local_gs(:)        !< Buffer for local gs-ops
     integer, allocatable :: local_dof_gs(:)          !< Local dof to gs mapping
     integer, allocatable :: local_gs_dof(:)          !< Local gs to dof mapping
     integer, allocatable :: local_blk_len(:)         !< Local non-facet blocks
     real(kind=rp), allocatable :: shared_gs(:)       !< Buffer for shared gs-op
     integer, allocatable :: shared_dof_gs(:)         !< Shared dof to gs map.
     integer, allocatable :: shared_gs_dof(:)         !< Shared gs to dof map.
     integer, allocatable :: shared_blk_len(:)        !< Shared non-facet blocks
     type(dofmap_t), pointer ::dofmap                 !< Dofmap for gs-ops
     type(htable_i8_t) :: shared_dofs                 !< Htable of shared dofs
     integer :: nlocal                                !< Local gs-ops
     integer :: nshared                               !< Shared gs-ops
     integer :: nlocal_blks                           !< Number of local blks
     integer :: nshared_blks                          !< Number of shared blks
     integer :: local_facet_offset                    !< offset for loc. facets
     integer :: shared_facet_offset                   !< offset for shr. facets
     class(gs_bcknd_t), allocatable :: bcknd          !< Gather-scatter backend
     class(gs_comm_t), allocatable :: comm            !< Comm. method
  end type gs_t

  private :: gs_init_mapping, gs_schedule
  
  interface gs_op
     module procedure gs_op_fld, gs_op_r4, gs_op_vector
  end interface gs_op

contains

  !> Initialize a gather-scatter kernel
  subroutine gs_init(gs, dofmap, bcknd)
    type(gs_t), intent(inout) :: gs
    type(dofmap_t), target, intent(inout) :: dofmap
    character(len=LOG_SIZE) :: log_buf
    character(len=20) :: bcknd_str
    integer, optional :: bcknd
    integer :: i, j, ierr, bcknd_
    integer(i8) :: glb_nshared, glb_nlocal
    logical :: use_device_mpi
    real(kind=rp), allocatable :: tmp(:)
    type(c_ptr) :: tmp_d = C_NULL_PTR
    integer :: strtgy(4) = (/ int(B'00'), int(B'01'), int(B'10'), int(B'11') /)
    integer :: avg_strtgy, env_len
    character(len=255) :: env_strtgy
    real(kind=dp) :: strtgy_time(4)

    call gs_free(gs)

    call neko_log%section('Gather-Scatter')
    
    gs%dofmap => dofmap
    
    ! Here one could use some heuristic or autotuning to select comm method,
    ! such as only using device MPI when there is enough data.
    !use_device_mpi = NEKO_DEVICE_MPI .and. gs%nshared .gt. 20000
    use_device_mpi = NEKO_DEVICE_MPI

    if (use_device_mpi) then
       call neko_log%message('Comm         :   Device MPI')
       allocate(gs_device_mpi_t::gs%comm)
    else
       call neko_log%message('Comm         :          MPI')
       allocate(gs_mpi_t::gs%comm)
    end if

    call gs%comm%init_dofs()
    call gs_init_mapping(gs)

    call gs_schedule(gs)

    glb_nlocal = int(gs%nlocal, i8)
    glb_nshared = int(gs%nshared, i8)
    
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
    select case(bcknd_)
    case(GS_BCKND_CPU)
       allocate(gs_cpu_t::gs%bcknd)
       bcknd_str = '         std'
    case(GS_BCKND_DEV)
       allocate(gs_device_t::gs%bcknd)
       if (NEKO_BCKND_HIP .eq. 1) then
          bcknd_str = '         hip'
       else if (NEKO_BCKND_CUDA .eq. 1) then
          bcknd_str = '        cuda'
       else if (NEKO_BCKND_OPENCL .eq. 1) then
          bcknd_str = '      opencl'
       end if
    case(GS_BCKND_SX)
       allocate(gs_sx_t::gs%bcknd)
       bcknd_str = '          sx'
    case default
       call neko_error('Unknown Gather-scatter backend')
    end select

    write(log_buf, '(A)') 'Backend      : ' // trim(bcknd_str)
    call neko_log%message(log_buf)
    


    call gs%bcknd%init(gs%nlocal, gs%nshared, gs%nlocal_blks, gs%nshared_blks)

    if (use_device_mpi) then
       select type(b => gs%bcknd)
       type is (gs_device_t)
          b%shared_on_host = .false.
       end select

       if(pe_size .gt. 1) then
          ! Select fastest device MPI strategy at runtime
          select type(c => gs%comm)
          type is (gs_device_mpi_t)
             call get_environment_variable("NEKO_GS_STRTGY", env_strtgy, env_len)
             if (env_len .eq. 0) then             
                allocate(tmp(dofmap%size()))
                call device_map(tmp, tmp_d, dofmap%size())
                tmp = 1.0_rp
                call device_memcpy(tmp, tmp_d, dofmap%size(), HOST_TO_DEVICE)
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
                     strtgy(avg_strtgy),']'

             else
                read(env_strtgy(1:env_len), *) i

                if (i .lt. 1 .or. i .gt. 4) then
                   call neko_error('Invalid gs sync strtgy')
                end if
                
                c%nb_strtgy = strtgy(i)
                avg_strtgy = i

                write(log_buf, '(A,B0.2,A)') 'Env. strtgy  :         [', &
                     strtgy(avg_strtgy),']'
             end if
             
             call neko_log%message(log_buf)             

          end select
       end if
    end if

    call neko_log%end_section()
    
  end subroutine gs_init

  !> Deallocate a gather-scatter kernel
  subroutine gs_free(gs)
    type(gs_t), intent(inout) :: gs

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
    type(htable_i8_t) :: dm
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
       lid = linear_index(1, 1, 1, i, lx, ly, lz)
       if (dofmap%shared_dof(1, 1, 1, i)) then
          id = gs_mapping_add_dof(sdm, dofmap%dof(1, 1, 1, i), max_sid)
          call shared_dof%push(id)
          call dof_shared%push(lid)
       else
          id = gs_mapping_add_dof(dm, dofmap%dof(1, 1, 1, i), max_id)
          call local_dof%push(id)
          call dof_local%push(lid)
       end if

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
    !
    ! Setup mapping for dofs on facets
    !
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
                   id = gs_mapping_add_dof(sdm, dofmap%dof(j, ly, 1, i), max_sid)
                   call shared_face_dof%push(id)
                   id = linear_index(j, ly, 1, i, lx, ly, lz)
                   call face_dof_shared%push(id)
                end do
                
             else
                do j = 2, lx - 1
                   id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, 1, i), max_id)
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
                   id = gs_mapping_add_dof(sdm, dofmap%dof(lx, k, 1, i), max_sid)
                   call shared_face_dof%push(id)
                   id = linear_index(lx, k, 1, i, lx, ly, lz)
                   call face_dof_shared%push(id)
                end do
             else
                do k = 2, ly - 1
                   id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, 1, i), max_id)
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
                      id = gs_mapping_add_dof(sdm, dofmap%dof(1, k, l, i), max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(1, k, l, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do l = 2, lz - 1
                   do k = 2, ly - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(1, k, l, i), max_id)
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
                      id = gs_mapping_add_dof(sdm, dofmap%dof(lx, k, l,  i), max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(lx, k, l, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do l = 2, lz - 1
                   do k = 2, ly - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, l,  i), max_id)
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
                      id = gs_mapping_add_dof(sdm, dofmap%dof(j, 1, l, i), max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(j, 1, l, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do l = 2, lz - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(j, 1, l, i), max_id)
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
                      id = gs_mapping_add_dof(sdm, dofmap%dof(j, ly, l, i), max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(j, ly, l, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do l = 2, lz - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, l, i), max_id)
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
                      id = gs_mapping_add_dof(sdm, dofmap%dof(j, k, 1, i), max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(j, k, 1, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do k = 2, ly - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(j, k, 1, i), max_id)
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
                      id = gs_mapping_add_dof(sdm, dofmap%dof(j, k, lz, i), max_sid)
                      call shared_face_dof%push(id)
                      id = linear_index(j, k, lz, i, lx, ly, lz)
                      call face_dof_shared%push(id)
                   end do
                end do
             else
                do k = 2, ly - 1
                   do j = 2, lx - 1
                      id = gs_mapping_add_dof(dm, dofmap%dof(j, k, lz, i), max_id)
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
    select type(dof_array => local_dof%data)
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
    select type(dof_array => local_face_dof%data)
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
    select type(dof_array => dof_local%data)
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
    select type(dof_array => face_dof_local%data)
    type is (integer)
       do i = 1, face_dof_local%size()
          gs%local_gs_dof(i+j) = dof_array(i)
       end do
    end select
    call face_dof_local%free()
       
    call gs_qsort_dofmap(gs%local_dof_gs, gs%local_gs_dof, &
         gs%nlocal, 1, gs%nlocal)
    
    call gs_find_blks(gs%local_dof_gs, gs%local_blk_len, &
         gs%nlocal_blks, gs%nlocal, gs%local_facet_offset)
    
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
    select type(dof_array => shared_dof%data)
    type is (integer)
       j =  shared_dof%size()
       do i = 1, j
          gs%shared_dof_gs(i) = dof_array(i)
       end do
    end select
    call shared_dof%free()

    ! Add shared dofs on faces

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type(dof_array => shared_face_dof%data)
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
    select type(dof_array => dof_shared%data)
    type is(integer)
       j = dof_shared%size()
       do i = 1, j
          gs%shared_gs_dof(i) = dof_array(i)
       end do
    end select
    call dof_shared%free()

    ! We should use the %array() procedure, which works great for
    ! GNU, Intel and NEC, but it breaks horribly on Cray when using
    ! certain data types
    select type(dof_array => face_dof_shared%data)
    type is (integer)
       do i = 1, face_dof_shared%size()
          gs%shared_gs_dof(i + j) = dof_array(i)
       end do
    end select
    call face_dof_shared%free()

    ! Allocate buffer for shared gs-ops
    allocate(gs%shared_gs(gs%nshared))

    if (gs%nshared .gt. 0) then
       call gs_qsort_dofmap(gs%shared_dof_gs, gs%shared_gs_dof, &
            gs%nshared, 1, gs%nshared)

       call gs_find_blks(gs%shared_dof_gs, gs%shared_blk_len, &
            gs%nshared_blks, gs%nshared, gs%shared_facet_offset)
    end if
    
  contains
    
    !> Register a unique dof
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
    subroutine gs_find_blks(dg, blk_len, nblks, n, m)
      integer, intent(in) :: n
      integer, intent(in) :: m
      integer, dimension(n), intent(inout) :: dg
      integer, allocatable, intent(inout) :: blk_len(:)
      integer, intent(inout) :: nblks
      integer :: i, j
      integer :: id, count
      type(stack_i4_t), target :: blks
      
      call blks%init()
      i = 1
      do while( i .lt. m)
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

      select type(blk_array => blks%data)
      type is(integer)
         nblks = blks%size()
         allocate(blk_len(nblks))
         do i = 1, nblks
            blk_len(i) = blk_array(i)
         end do
      end select
      call blks%free()
      
    end subroutine gs_find_blks

  end subroutine gs_init_mapping

  !> Schedule shared gather-scatter operations
  subroutine gs_schedule(gs)
    type(gs_t), target, intent(inout) :: gs
    integer(kind=i8), allocatable :: send_buf(:), recv_buf(:)
    integer(kind=i2), allocatable :: shared_flg(:), recv_flg(:)
    type(htable_iter_i8_t) :: it
    type(stack_i4_t) :: send_pe, recv_pe
    type(MPI_Status) :: status
    type(MPI_Request) :: send_req, recv_req
    integer :: i, j, max_recv, src, dst, ierr, n_recv
    integer :: tmp, shared_gs_id
    integer :: nshared_unique

    nshared_unique = gs%shared_dofs%num_entries()
    
    call it%init(gs%shared_dofs)
    allocate(send_buf(nshared_unique))
    i = 1
    do while(it%next())
       send_buf(i) = it%key()
       i = i + 1
    end do

    call send_pe%init()
    call recv_pe%init()
    

    !
    ! Schedule exchange of shared dofs
    !

    call MPI_Allreduce(nshared_unique, max_recv, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

    allocate(recv_buf(max_recv))
    allocate(shared_flg(max_recv))
    allocate(recv_flg(max_recv))

    !> @todo Consider switching to a crystal router...
    do i = 1, size(gs%dofmap%msh%neigh_order)
       src = modulo(pe_rank - gs%dofmap%msh%neigh_order(i) + pe_size, pe_size)
       dst = modulo(pe_rank + gs%dofmap%msh%neigh_order(i), pe_size)
       
       if (gs%dofmap%msh%neigh(src)) then
          call MPI_Irecv(recv_buf, max_recv, MPI_INTEGER8, &
               src, 0, NEKO_COMM, recv_req, ierr)
       end if

       if (gs%dofmap%msh%neigh(dst)) then
          call MPI_Isend(send_buf, nshared_unique, MPI_INTEGER8, &
               dst, 0, NEKO_COMM, send_req, ierr)
       end if

       if (gs%dofmap%msh%neigh(src)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER8, n_recv, ierr)
          
          do j = 1, n_recv
             shared_flg(j) = gs%shared_dofs%get(recv_buf(j), shared_gs_id)
             if (shared_flg(j) .eq. 0) then
                !> @todo don't touch others data...
                call gs%comm%recv_dof(src)%push(shared_gs_id)
             end if
          end do
          
          if (gs%comm%recv_dof(src)%size() .gt. 0) then
             call recv_pe%push(src)
          end if
       end if

       if (gs%dofmap%msh%neigh(dst)) then
          call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
          call MPI_Irecv(recv_flg, max_recv, MPI_INTEGER2, &
               dst, 0, NEKO_COMM, recv_req, ierr)          
       end if

       if (gs%dofmap%msh%neigh(src)) then
          call MPI_Isend(shared_flg, n_recv, MPI_INTEGER2, &
               src, 0, NEKO_COMM, send_req, ierr)
       end if

       if (gs%dofmap%msh%neigh(dst)) then
          call MPI_Wait(recv_req, status, ierr)
          call MPI_Get_count(status, MPI_INTEGER2, n_recv, ierr)
       
          do j = 1, n_recv
             if (recv_flg(j) .eq. 0) then
                tmp = gs%shared_dofs%get(send_buf(j), shared_gs_id) 
                !> @todo don't touch others data...
                call gs%comm%send_dof(dst)%push(shared_gs_id)
             end if
          end do

          if (gs%comm%send_dof(dst)%size() .gt. 0) then
             call send_pe%push(dst)
          end if
       end if

       if (gs%dofmap%msh%neigh(src)) then
          call MPI_Wait(send_req, MPI_STATUS_IGNORE, ierr)
       end if
       
    end do

    call gs%comm%init(send_pe, recv_pe)
    
    call send_pe%free()
    call recv_pe%free()

    deallocate(send_buf)
    deallocate(recv_flg)
    deallocate(shared_flg)
    !This arrays seems to take massive amounts of memory... 
    call gs%shared_dofs%free()

  end subroutine gs_schedule

  !> Gather-scatter operation on a field @a u with op @a op
  subroutine gs_op_fld(gs, u, op)
    type(gs_t), intent(inout) :: gs
    type(field_t), intent(inout) :: u
    integer :: n, op
    
    n = u%msh%nelv * u%Xh%lx * u%Xh%ly * u%Xh%lz
    call gs_op_vector(gs, u%x, n, op)
    
  end subroutine gs_op_fld
  
  !> Gather-scatter operation on a rank 4 array
  subroutine gs_op_r4(gs, u, n, op)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), dimension(:,:,:,:), intent(inout) :: u
    integer :: op

    call gs_op_vector(gs, u, n, op)
    
  end subroutine gs_op_r4
  
  !> Gather-scatter operation on a vector @a u with op @a op
  subroutine gs_op_vector(gs, u, n, op)
    type(gs_t), intent(inout) :: gs
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    integer :: m, l, op, lo, so
    
    lo = gs%local_facet_offset
    so = -gs%shared_facet_offset
    m = gs%nlocal
    l = gs%nshared

    call profiler_start_region("gather-scatter")
    ! Gather shared dofs
    if (pe_size .gt. 1) then
       call profiler_start_region("gs_nbrecv")
       call gs%comm%nbrecv()
       call profiler_end_region
       call profiler_start_region("gs_gather_shared")
       call gs%bcknd%gather(gs%shared_gs, l, so, gs%shared_dof_gs, u, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, op, .true.)
       call profiler_end_region
       call profiler_start_region("gs_nbsend")
       call gs%comm%nbsend(gs%shared_gs, l, gs%bcknd%gather_event)
       call profiler_end_region
       
    end if
    
    ! Gather-scatter local dofs
    call profiler_start_region("gs_local")
    call gs%bcknd%gather(gs%local_gs, m, lo, gs%local_dof_gs, u, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, op, .false.)
    call gs%bcknd%scatter(gs%local_gs, m, gs%local_dof_gs, u, n, &
         gs%local_gs_dof, gs%nlocal_blks, gs%local_blk_len, .false.)
    call profiler_end_region
    ! Scatter shared dofs
    if (pe_size .gt. 1) then
       call profiler_start_region("gs_nbwait")
       call gs%comm%nbwait(gs%shared_gs, l, op)
       call profiler_end_region
       call profiler_start_region("gs_scatter_shared")
       call gs%bcknd%scatter(gs%shared_gs, l, gs%shared_dof_gs, u, n, &
            gs%shared_gs_dof, gs%nshared_blks, gs%shared_blk_len, .true.)
       call profiler_end_region
    end if

    call profiler_end_region
    
  end subroutine gs_op_vector
  
end module gather_scatter
