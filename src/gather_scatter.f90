!> Gather-scatter
module gather_scatter
  use mesh
  use dofmap
  use field
  use num_types
  use htable
  use stack
  use utils
  use mpi
  implicit none

  integer, parameter :: GS_OP_ADD = 1, GS_OP_MUL = 2, &
       GS_OP_MIN = 3, GS_OP_MAX = 4
  
  type, private :: gs_comm_t
     integer :: status(MPI_STATUS_SIZE)
     integer :: request
     logical :: flag
     real(kind=dp), allocatable :: data(:)
  end type gs_comm_t
  
  type gs_t
     real(kind=dp), allocatable :: local_gs(:)        !< Buffer for local gs-ops
     integer, allocatable :: local_dof_gs(:)          !< Local dof to gs mapping
     integer, allocatable :: local_gs_dof(:)          !< Local gs to dof mapping
     real(kind=dp), allocatable :: shared_gs(:)       !< Buffer for shared gs-op
     integer, allocatable :: shared_dof_gs(:)         !< Shared dof to gs map.
     integer, allocatable :: shared_gs_dof(:)         !< Shared gs to dof map.
     integer, allocatable :: send_pe(:)               !< Send order
     integer, allocatable :: recv_pe(:)               !< Recv order
     type(stack_i4_t), allocatable :: send_dof(:)     !< Send dof to shared-gs
     type(stack_i4_t), allocatable :: recv_dof(:)     !< Recv dof to shared-gs
     type(gs_comm_t), allocatable :: send_buf(:)      !< Comm. buffers
     type(gs_comm_t), allocatable :: recv_buf(:)      !< Comm. buffers
     type(dofmap_t), pointer ::dofmap                 !< Dofmap for gs-ops
     type(htable_i8_t) :: shared_dofs                 !< Htable of shared dofs
     integer :: nlocal                                !< Local gs-ops
     integer :: nshared                               !< Shared gs-ops
     integer :: local_facet_offset                    !< offset for loc. facets
     integer :: shared_facet_offset                   !< offset for shr. facets
  end type gs_t

  private :: gs_init_mapping, gs_schedule

  interface gs_op
     module procedure gs_op_fld, gs_op_vector
  end interface gs_op

contains

  !> Initialize a gather-scatter kernel
  subroutine gs_init(gs, dofmap)
    type(gs_t), intent(inout) :: gs
    type(dofmap_t), target, intent(inout) :: dofmap
    integer :: i

    call gs_free(gs)

    gs%dofmap => dofmap
    
    allocate(gs%send_dof(0:pe_size-1))
    allocate(gs%recv_dof(0:pe_size-1))

    do i = 0, pe_size -1
       call gs%send_dof(i)%init()
       call gs%recv_dof(i)%init()
    end do

    call gs_init_mapping(gs)

    call gs_schedule(gs)
    
  end subroutine gs_init

  !> Deallocate a gather-scatter kernel
  subroutine gs_free(gs)
    type(gs_t), intent(inout) :: gs
    integer :: i

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

    if (allocated(gs%shared_gs)) then
       deallocate(gs%shared_gs)
    end if
    
    if (allocated(gs%shared_dof_gs)) then
       deallocate(gs%shared_dof_gs)
    end if

    if (allocated(gs%shared_gs_dof)) then
       deallocate(gs%shared_gs_dof)
    end if

    gs%nlocal = 0
    gs%nshared = 0

    call gs%shared_dofs%free()

    if (allocated(gs%send_dof)) then
       do i = 0, pe_size - 1
          call gs%send_dof(i)%free()
       end do
       deallocate(gs%send_dof)
    end if

    if (allocated(gs%recv_dof)) then
       do i = 0, pe_size - 1
          call gs%recv_dof(i)%free()
       end do
       deallocate(gs%recv_dof)
    end if

    if (allocated(gs%send_pe)) then
       deallocate(gs%send_pe)
    end if

    if (allocated(gs%recv_pe)) then
       deallocate(gs%recv_pe)
    end if

    if (allocated(gs%send_buf)) then
       do i = 1, size(gs%send_buf)
          if (allocated(gs%send_buf(i)%data)) then
             deallocate(gs%send_buf(i)%data)
          end if
       end do
       deallocate(gs%send_buf)
    end if


    if (allocated(gs%recv_buf)) then
       do i = 1, size(gs%recv_buf)
          if (allocated(gs%recv_buf(i)%data)) then
             deallocate(gs%recv_buf(i)%data)
          end if
       end do
       deallocate(gs%recv_buf)
    end if
  end subroutine gs_free

  !> Setup mapping of dofs to gather-scatter operations
  subroutine gs_init_mapping(gs)
    type(gs_t), target, intent(inout) :: gs
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dofmap
    type(stack_i4_t) :: local_dof, dof_local, shared_dof, dof_shared
    type(stack_i4_t) :: local_face_dof, face_dof_local
    type(stack_i4_t) :: shared_face_dof, face_dof_shared
    integer :: i, j, k, l, lx, ly, lz, max_id, max_sid, id, lid
    integer, pointer :: sp(:)
    type(htable_i8_t) :: dm
    type(htable_i8_t), pointer :: sdm

    dofmap => gs%dofmap
    msh => dofmap%msh
    sdm => gs%shared_dofs

    call dm%init(msh%nelv, i)
    call sdm%init(msh%nelv, i)
    
    lx = dofmap%Xh%lx
    ly = dofmap%Xh%ly
    lz = dofmap%Xh%lz

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

    end do

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

    !
    ! Setup mapping for dofs on facets
    !
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
    
       

    call dm%free()
    
    gs%nlocal = local_dof%size() + local_face_dof%size()
    gs%local_facet_offset = local_dof%size() + 1
    
    ! Finalize local dof to gather-scatter index
    allocate(gs%local_dof_gs(gs%nlocal))

    ! Add dofs on points and edges
    sp => local_dof%array()
    j = local_dof%size()
    do i = 1, j
       gs%local_dof_gs(i) = sp(i)
    end do
    call local_dof%free()

    ! Add dofs on faces
    sp => local_face_dof%array()
    do i = 1, local_face_dof%size()
       gs%local_dof_gs(i + j) = sp(i)
    end do
    call local_face_dof%free()

    ! Finalize local gather-scatter index to dof
    allocate(gs%local_gs_dof(gs%nlocal))

    ! Add gather-scatter index on points and edges
    sp => dof_local%array()
    j = dof_local%size()
    do i = 1, j
       gs%local_gs_dof(i) = sp(i)
    end do
    call dof_local%free()

    sp => face_dof_local%array()
    do i = 1, face_dof_local%size()
       gs%local_gs_dof(i+j) = sp(i)
    end do
    call face_dof_local%free()
       
    call gs_qsort_dofmap(gs%local_dof_gs, gs%local_gs_dof, &
         gs%nlocal, 1, gs%nlocal)
    
    ! Allocate buffer for local gs-ops
    allocate(gs%local_gs(gs%nlocal))   

    gs%nshared = shared_dof%size() + shared_face_dof%size()
    gs%shared_facet_offset = shared_dof%size() + 1

    ! Finalize shared dof to gather-scatter index
    allocate(gs%shared_dof_gs(gs%nshared))

    ! Add shared dofs on points and edges
    sp => shared_dof%array()
    j =  shared_dof%size()
    do i = 1, j
       gs%shared_dof_gs(i) = sp(i)
    end do
    call shared_dof%free()

    ! Add shared dofs on faces
    sp => shared_face_dof%array()
    do i = 1, shared_face_dof%size()
       gs%shared_dof_gs(i + j) = sp(i)
    end do
    call shared_face_dof%free()
    
    ! Finalize shared gather-scatter index to dof
    allocate(gs%shared_gs_dof(gs%nshared))

    ! Add dofs on points and edges 
    sp => dof_shared%array()
    j = dof_shared%size()
    do i = 1, j
       gs%shared_gs_dof(i) = sp(i)
    end do
    call dof_shared%free()

    sp => face_dof_shared%array()
    do i = 1, face_dof_shared%size()
       gs%shared_gs_dof(i + j) = sp(i)
    end do
    call face_dof_shared%free()

    ! Allocate buffer for shared gs-ops
    allocate(gs%shared_gs(gs%nshared))

    call gs_qsort_dofmap(gs%shared_dof_gs, gs%shared_gs_dof, &
         gs%nshared, 1, gs%nshared)

  contains
    
    !> Register a unique dof
    function gs_mapping_add_dof(map_, dof, max_id) result(id)
      type(htable_i8_t), intent(inout) :: map_
      integer(kind=8), intent(inout) :: dof
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
      integer, dimension(n), intent(inout) :: dg
      integer, dimension(n), intent(inout) :: gd
      integer :: n, lo, hi
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

  end subroutine gs_init_mapping

  !> Schedule shared gather-scatter operations
  subroutine gs_schedule(gs)
    type(gs_t), intent(inout) :: gs
    integer(kind=8), allocatable :: send_buf(:), recv_buf(:)
    integer(kind=2), allocatable :: shared_flg(:), recv_flg(:)
    type(htable_iter_i8_t) :: it
    type(stack_i4_t) :: send_pe, recv_pe
    integer :: i, j, max_recv, src, dst, ierr, n_recv
    integer :: shared_id, tmp, shared_gs_id
    integer :: status(MPI_STATUS_SIZE)
    integer :: nshared_unique
    integer, pointer :: sp(:), rp(:)


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
    do i = 1, pe_size - 1
       src = modulo(pe_rank - i + pe_size, pe_size)
       dst = modulo(pe_rank + i, pe_size)

       call MPI_Sendrecv(send_buf, nshared_unique, MPI_INTEGER8, dst, 0, &
            recv_buf, max_recv, MPI_INTEGER8, src, 0, NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_INTEGER8, n_recv, ierr)

       do j = 1, n_recv
          shared_flg(j) = gs%shared_dofs%get(recv_buf(j), shared_gs_id)
          if (shared_flg(j) .eq. 0) then
             call gs%recv_dof(src)%push(shared_gs_id)
          end if
       end do

       if (gs%recv_dof(src)%size() .gt. 0) then
          call recv_pe%push(src)
       end if

       call MPI_Sendrecv(shared_flg, n_recv, MPI_INTEGER2, src, 1, &
            recv_flg, max_recv, MPI_INTEGER2, dst, 1, NEKO_COMM, status, ierr)
       call MPI_Get_count(status, MPI_INTEGER2, n_recv, ierr)

       do j = 1, n_recv
          if (recv_flg(j) .eq. 0) then
             tmp = gs%shared_dofs%get(send_buf(j), shared_gs_id) 
             call gs%send_dof(dst)%push(shared_gs_id)
          end if
       end do

       if (gs%send_dof(dst)%size() .gt. 0) then
          call send_pe%push(dst)
       end if
       
    end do

    allocate(gs%send_pe(send_pe%size()))
    allocate(gs%send_buf(send_pe%size()))
    
    sp => send_pe%array()
    do i = 1, send_pe%size()
       gs%send_pe(i) = sp(i)
       allocate(gs%send_buf(i)%data(gs%send_dof(sp(i))%size()))
    end do

    allocate(gs%recv_pe(recv_pe%size()))
    allocate(gs%recv_buf(recv_pe%size()))
    
    rp => recv_pe%array()
    do i = 1, recv_pe%size()
       gs%recv_pe(i) = rp(i)
       allocate(gs%recv_buf(i)%data(gs%recv_dof(rp(i))%size()))
    end do
    
    call send_pe%free()
    call recv_pe%free()

    deallocate(send_buf)
    deallocate(recv_flg)
    deallocate(shared_flg)

  end subroutine gs_schedule

  !> Gather-scatter operation on a field @a u with op @a op
  subroutine gs_op_fld(gs, u, op)
    type(gs_t), intent(inout) :: gs
    type(field_t), intent(inout) :: u
    integer :: n, op
    
    n = u%msh%nelv * u%Xh%lx * u%Xh%ly * u%Xh%lz
    call gs_op_vector(gs, u%x, n, op)
    
  end subroutine gs_op_fld
  
  !> Gather-scatter operation on a vector @a u with op @a op
  subroutine gs_op_vector(gs, u, n, op)
    type(gs_t), intent(inout) :: gs
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, intent(inout) :: n
    integer :: m, l, op, lo, so
    
    lo = gs%local_facet_offset
    so = -gs%shared_facet_offset
    m = gs%nlocal
    l = gs%nshared

    ! Gather shared dofs
    if (pe_size .gt. 1) then

       call gs_nbrecv(gs)

       call gs_gather(gs%shared_gs, l, so, gs%shared_dof_gs, &
            u, n, gs%shared_gs_dof, op)

       call gs_nbsend(gs, gs%shared_gs, l)
       
    end if
    
    ! Gather-scatter local dofs

    call gs_gather(gs%local_gs, m, lo, gs%local_dof_gs, &
         u, n, gs%local_gs_dof, op)
    call gs_scatter(gs%local_gs, m, gs%local_dof_gs, &
         u, n, gs%local_gs_dof)

    ! Scatter shared dofs
    if (pe_size .gt. 1) then

       call gs_nbwait(gs, gs%shared_gs, l, op)

       call gs_scatter(gs%shared_gs, l, gs%shared_dof_gs, &
            u, n, gs%shared_gs_dof)
    end if
       
  end subroutine gs_op_vector
  
  !> Gather kernel
  subroutine gs_gather(v, m, o, dg, u, n, gd, op)
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, intent(in) :: m
    integer, intent(in) :: o
    integer, intent(in) :: n
    integer :: op
    
    select case(op)
    case (GS_OP_ADD)
       call gs_gather_kernel_add(v, m, o, dg, u, n, gd)
    case (GS_OP_MUL)
       call gs_gather_kernel_mul(v, m, o, dg, u, n, gd)
    case (GS_OP_MIN)
       call gs_gather_kernel_min(v, m, o, dg, u, n, gd)
    case (GS_OP_MAX)
       call gs_gather_kernel_max(v, m, o, dg, u, n, gd)
    end select
    
  end subroutine gs_gather
 
  !> Gather kernel for addition of data
  !! \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
  subroutine gs_gather_kernel_add(v, m, o, dg, u, n, gd)
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, intent(in) :: m
    integer, intent(in) :: o
    integer, intent(in) :: n
    integer :: i, no
    real(kind=dp) :: tmp

    v = 0d0

    if (o .lt. 0) then

       no = abs(o)

       do i = 1, no - 1
          v(dg(i)) = v(dg(i)) + u(gd(i))
       end do
       
       do i = no, m
          v(dg(i)) = u(gd(i))
       end do
       
    else

       do i = 1, o - 1
          v(dg(i)) = v(dg(i)) + u(gd(i))
       end do

       do i = o, m, 2
          tmp  = u(gd(i)) + u(gd(i+1))
          v(dg(i)) = tmp
          v(dg(i+1)) = tmp
       end do
       
    end if
  end subroutine gs_gather_kernel_add

  !> Gather kernel for multiplication of data
  !! \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
  subroutine gs_gather_kernel_mul(v, m, o, dg, u, n, gd)
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, intent(in) :: m
    integer, intent(in) :: o
    integer, intent(in) :: n
    integer :: i, no
    real(kind=dp) :: tmp

    if (o .lt. 0) then

       no = abs(o)
       
       do i = 1, no - 1
          v(dg(i)) = v(dg(i)) * u(gd(i))
       end do
       
       do i = no, m
          v(dg(i)) = u(gd(i))
       end do
       
    else

       do i = 1, o - 1
          v(dg(i)) = v(dg(i)) * u(gd(i))
       end do

       do i = o, m, 2
          tmp  = u(gd(i)) * u(gd(i+1))
          v(dg(i)) = tmp
          v(dg(i+1)) = tmp          
       end do
       
    end if
  end subroutine gs_gather_kernel_mul
  
  !> Gather kernel for minimum of data
  !! \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_min(v, m, o, dg, u, n, gd)
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, intent(in) :: m
    integer, intent(in) :: o
    integer, intent(in) :: n
    integer :: i
    v = 1e24
    do i = 1, m
       v(dg(i)) = min(v(dg(i)), u(gd(i)))
    end do
  end subroutine gs_gather_kernel_min

  !> Gather kernel for maximum of data
  !! \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
  subroutine gs_gather_kernel_max(v, m, o, dg, u, n, gd)
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, intent(in) :: m
    integer, intent(in) :: o
    integer, intent(in) :: n
    integer :: i
    v = -1e24
    do i = 1, m
       v(dg(i)) = max(v(dg(i)), u(gd(i)))
    end do
  end subroutine gs_gather_kernel_max

  !> Scatter kernel  @todo Make the kernel abstract
  subroutine gs_scatter(v, m, dg, u, n, gd)
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, intent(in) :: m
    integer, intent(in) :: n
        
    call gs_scatter_kernel(v, m, dg, u, n, gd)

  end subroutine gs_scatter

  !> Scatter kernel \f$ u(gd(i) = v(dg(i)) \f$
  subroutine gs_scatter_kernel(v, m, dg, u, n, gd)
    real(kind=dp), dimension(m), intent(inout) :: v
    integer, dimension(m), intent(inout) :: dg
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, dimension(m), intent(inout) :: gd
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer :: i
    do i = 1, m
       u(gd(i)) = v(dg(i))
    end do
  end subroutine gs_scatter_kernel


  !> Post non-blocking receive operations
  subroutine gs_nbrecv(gs)
    type(gs_t), intent(inout) :: gs
    integer :: i, ierr

    do i = 1, size(gs%recv_pe)
       call MPI_IRecv(gs%recv_buf(i)%data, size(gs%recv_buf(i)%data), &
            MPI_DOUBLE_PRECISION, gs%recv_pe(i), 0, &
            NEKO_COMM, gs%recv_buf(i)%request, ierr)
       gs%recv_buf(i)%flag = .false.
    end do
    
  end subroutine gs_nbrecv

  !> Post non-blocking send operations
  subroutine gs_nbsend(gs, u, n)
    type(gs_t), intent(inout) :: gs
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, intent(in) :: n
    integer ::  i, j, ierr, dst
    integer , pointer :: sp(:)

    do i = 1, size(gs%send_pe)
       dst = gs%send_pe(i)
       sp => gs%send_dof(dst)%array()
       do j = 1, gs%send_dof(dst)%size()
          gs%send_buf(i)%data(j) = u(sp(j))
       end do

       call MPI_Isend(gs%send_buf(i)%data, size(gs%send_buf(i)%data), &
            MPI_DOUBLE_PRECISION, gs%send_pe(i), 0, &
            NEKO_COMM, gs%send_buf(i)%request, ierr)
       gs%send_buf(i)%flag = .false.
    end do
    
  end subroutine gs_nbsend

  !> Wait for non-blocking operations
  subroutine gs_nbwait(gs, u, n, op)
    type(gs_t), intent(inout) ::gs
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, intent(in) :: n
    integer :: i, j, src, ierr
    integer :: op
    integer , pointer :: sp(:)
    integer :: nreqs, ndone

    nreqs = size(gs%recv_pe)

    do while (nreqs .gt. 0) 
       do i = 1, size(gs%recv_pe)
          if (.not. gs%recv_buf(i)%flag) then
             call MPI_Test(gs%recv_buf(i)%request, gs%recv_buf(i)%flag, &
                  gs%recv_buf(i)%status, ierr)
             if (gs%recv_buf(i)%flag) then
                nreqs = nreqs - 1
                !> @todo Check size etc against status
                src = gs%recv_pe(i)
                sp => gs%recv_dof(src)%array()
                select case(op)
                case (GS_OP_ADD)
                   do j = 1, gs%send_dof(src)%size()
                      u(sp(j)) = u(sp(j)) + gs%recv_buf(i)%data(j)
                   end do
                case (GS_OP_MUL)
                   do j = 1, gs%send_dof(src)%size()
                      u(sp(j)) = u(sp(j)) * gs%recv_buf(i)%data(j)
                   end do
                case (GS_OP_MIN)
                   do j = 1, gs%send_dof(src)%size()
                      u(sp(j)) = min(u(sp(j)), gs%recv_buf(i)%data(j))
                   end do
                case (GS_OP_MAX)
                   do j = 1, gs%send_dof(src)%size()
                      u(sp(j)) = max(u(sp(j)), gs%recv_buf(i)%data(j))
                   end do
                end select
             end if
          end if
       end do
    end do

    nreqs = size(gs%send_pe)
    do while (nreqs .gt. 0) 
       do i = 1, size(gs%send_pe)
          if (.not. gs%send_buf(i)%flag) then
             call MPI_Test(gs%send_buf(i)%request, gs%send_buf(i)%flag, &
                  MPI_STATUS_IGNORE, ierr)
             if (gs%send_buf(i)%flag) nreqs = nreqs - 1
          end if
       end do
    end do
    
  end subroutine gs_nbwait
  
end module gather_scatter
