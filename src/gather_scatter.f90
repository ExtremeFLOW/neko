!> Gather-scatter
module gather_scatter
  use mesh
  use dofmap
  use field
  use num_types
  use htable
  use stack
  use utils
  implicit none

  integer, parameter :: GS_OP_ADD = 1, GS_OP_MUL = 2, &
       GS_OP_MIN = 3, GS_OP_MAX = 4
  
  type gs_t
     real(kind=dp), allocatable :: local_gs(:)        !< Buffer for local gs-ops
     integer, allocatable :: local_dof_gs(:)          !< Local dof to gs mapping
     integer, allocatable :: local_gs_dof(:)          !< Local gs to dof mapping
     type(dofmap_t), pointer ::dofmap                 !< Dofmap for gs-ops
     integer :: nlocal                                !< Local gs-ops
  end type gs_t

  private :: gs_init_mapping

contains

  !> Initialize a gather-scatter kernel
  subroutine gs_init(gs, dofmap)
    type(gs_t), intent(inout) :: gs
    type(dofmap_t), target, intent(inout) :: dofmap

    call gs_free(gs)

    gs%dofmap => dofmap

    call gs_init_mapping(gs)
    
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

    gs%nlocal = 0
    
  end subroutine gs_free

  !> Setup mapping of dofs to gather-scatter operations
  subroutine gs_init_mapping(gs)
    type(gs_t), intent(inout) :: gs
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dofmap
    type(stack_i4_t) :: local_dof, dof_local
    integer :: i, j, k, l, lx, ly, lz, max_id, id
    integer, pointer :: sp(:)
    type(htable_i8_t) :: dm

    dofmap => gs%dofmap
    msh => dofmap%msh

    call dm%init(msh%nelv, i)
    
    lx = dofmap%Xh%lx
    ly = dofmap%Xh%ly
    lz = dofmap%Xh%lz


    call local_dof%init()
    call dof_local%init()

    !
    ! Setup mapping for dofs points
    !
    
    max_id = 0
    do i = 1, msh%nelv
       id = gs_mapping_add_dof(dm, dofmap%dof(1, 1, 1, i), max_id)
       call local_dof%push(id)
       id = linear_index(1, 1, 1, i, lx, ly, lz)
       call dof_local%push(id)

       id = gs_mapping_add_dof(dm, dofmap%dof(lx, 1, 1, i), max_id)
       call local_dof%push(id)
       id = linear_index(lx, 1, 1, i, lx, ly, lz)
       call dof_local%push(id)
              
       id = gs_mapping_add_dof(dm, dofmap%dof(1, ly, 1, i), max_id)
       call local_dof%push(id)
       id = linear_index(1, ly, 1, i, lx, ly, lz)
       call dof_local%push(id)
       
       id = gs_mapping_add_dof(dm, dofmap%dof(lx, ly, 1, i), max_id)
       call local_dof%push(id)
       id = linear_index(lx, ly, 1, i, lx, ly, lz)
       call dof_local%push(id)
       
       id = gs_mapping_add_dof(dm, dofmap%dof(1, 1, lz, i), max_id)
       call local_dof%push(id)
       id = linear_index(1, 1, lz, i, lx, ly, lz)
       call dof_local%push(id)
       
       id = gs_mapping_add_dof(dm, dofmap%dof(lx, 1, lz, i), max_id)
       call local_dof%push(id)
       id = linear_index(lx, 1, lz, i, lx, ly, lz)
       call dof_local%push(id)
       
       id = gs_mapping_add_dof(dm, dofmap%dof(1, ly, lz, i), max_id)
       call local_dof%push(id)
       id = linear_index(1, ly, lz, i, lx, ly, lz)
       call dof_local%push(id)
             
       id = gs_mapping_add_dof(dm, dofmap%dof(lx, ly, lz, i), max_id)
       call local_dof%push(id)
       id = linear_index(lx, ly, lz, i, lx, ly, lz)
       call dof_local%push(id)

    end do

    !
    ! Setup mapping for dofs on edges
    !
    do i = 1, msh%nelv

       !
       ! dofs on edges in x-direction
       !
       do j = 2, lx - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(j, 1, 1, i), max_id)
          call local_dof%push(id)
          id = linear_index(j, 1, 1, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       do j = 2, lx - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(j, 1, lz, i), max_id)
          call local_dof%push(id)
          id = linear_index(j, 1, lz, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       do j = 2, lx - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, 1, i), max_id)
          call local_dof%push(id)
          id = linear_index(j, ly, 1, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       do j = 2, lx - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, lz, i), max_id)
          call local_dof%push(id)
          id = linear_index(j, ly, lz, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       !
       ! dofs on edges in y-direction
       !
       do k = 2, ly - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(1, k, 1, i), max_id)
          call local_dof%push(id)
          id = linear_index(1, k, 1, i, lx, ly, lz)
          call dof_local%push(id)
       end do
       
       do k = 2, ly - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(1, k, lz, i), max_id)
          call local_dof%push(id)
          id = linear_index(1, k, lz, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       do k = 2, ly - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, lz, i), max_id)
          call local_dof%push(id)
          id = linear_index(lx, k, 1, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       do k = 2, ly - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, lz, i), max_id)
          call local_dof%push(id)
          id = linear_index(lx, k, lz, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       !
       ! dofs on edges in z-direction
       !
       do l = 2, lz - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(1, 1, l, i), max_id)
          call local_dof%push(id)
          id = linear_index(1, 1, l, i, lx, ly, lz)
          call dof_local%push(id)
       end do
    
       do l = 2, lz - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(lx, 1, l, i), max_id)
          call local_dof%push(id)
          id = linear_index(lx, 1, l, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       do l = 2, lz - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(1, ly, l, i), max_id)
          call local_dof%push(id)
          id = linear_index(1, ly, l, i, lx, ly, lz)
          call dof_local%push(id)
       end do

       do l = 2, lz - 1
          id = gs_mapping_add_dof(dm, dofmap%dof(lx, ly, l, i), max_id)
          call local_dof%push(id)
          id = linear_index(lx, ly, l, i, lx, ly, lz)
          call dof_local%push(id)
       end do

    end do

    !
    ! Setup mapping for dofs on facets
    !
    do i = 1, msh%nelv

       ! Facets in x-direction (s, t)-plane
       do l = 2, lz - 1
          do k = 2, ly - 1
             id = gs_mapping_add_dof(dm, dofmap%dof(1, k, l,  i), max_id)
             call local_dof%push(id)
             id = linear_index(1, k, l, i, lx, ly, lz)
             call dof_local%push(id)
          end do
       end do

       do l = 2, lz - 1
          do k = 2, ly - 1
             id = gs_mapping_add_dof(dm, dofmap%dof(lx, k, l,  i), max_id)
             call local_dof%push(id)
             id = linear_index(lx, k, l, i, lx, ly, lz)
             call dof_local%push(id)
          end do
       end do
       
       ! Facets in y-direction (r, t)-plane
       do l = 2, lz - 1
          do j = 2, lx - 1
             id = gs_mapping_add_dof(dm, dofmap%dof(j, 1, l, i), max_id)
             call local_dof%push(id)
             id = linear_index(j, 1, l, i, lx, ly, lz)
             call dof_local%push(id)
          end do
       end do

       do l = 2, lz - 1
          do j = 2, lx - 1
             id = gs_mapping_add_dof(dm, dofmap%dof(j, ly, l, i), max_id)
             call local_dof%push(id)
             id = linear_index(j, ly, l, i, lx, ly, lz)
             call dof_local%push(id)
          end do
       end do

       ! Facets in z-direction (r, s)-plane
       do k = 2, ly - 1
          do j = 2, lx - 1
             id = gs_mapping_add_dof(dm, dofmap%dof(j, k, 1, i), max_id)
             call local_dof%push(id)
             id = linear_index(j, k, 1, i, lx, ly, lz)
             call dof_local%push(id)
          end do
       end do

       do k = 2, ly - 1
          do j = 2, lx - 1
             id = gs_mapping_add_dof(dm, dofmap%dof(j, k, lz, i), max_id)
             call local_dof%push(id)
             id = linear_index(j, k, lz, i, lx, ly, lz)
             call dof_local%push(id)
          end do
       end do
    end do


    call dm%free()
    
    gs%nlocal = local_dof%size()

    ! Finalize local dof to gather-scatter index
    allocate(gs%local_dof_gs(gs%nlocal))
    sp => local_dof%array()
    do i = 1, local_dof%size()
       gs%local_dof_gs(i) = sp(i)
    end do
    call local_dof%free()

    ! Finalize local gather-scatter index to dof
    allocate(gs%local_gs_dof(gs%nlocal))
    sp => dof_local%array()
    do i = 1, dof_local%size()
       gs%local_gs_dof(i) = sp(i)
    end do
    call dof_local%free()

    ! Allocate buffer for local gs-ops
    allocate(gs%local_gs(gs%nlocal))

  contains
    
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
    
  end subroutine gs_init_mapping

  !> Gather-scatter operation with op @a op
  subroutine gs_op(gs, u, op)
    type(gs_t), intent(inout) :: gs
    type(field_t), intent(inout) :: u
    integer :: op

    call gs_gather(gs, u, op)
    call gs_scatter(gs, u)
    
  end subroutine gs_op
  
  !> Gather kernel
  subroutine gs_gather(gs, u, op)
    type(gs_t), intent(inout) :: gs
    type(field_t), intent(inout) :: u
    integer, intent(in) :: op
    integer :: n
    
    n = u%msh%nelv * u%Xh%lx * u%Xh%ly * u%Xh%lz
    
    select case(op)
    case (GS_OP_ADD)
       call gs_gather_local_add(gs, u%x, n)
    case (GS_OP_MUL)
       call gs_gather_local_mul(gs, u%x, n)
    case (GS_OP_MIN)
       call gs_gather_local_min(gs, u%x, n)
    case (GS_OP_MAX)
       call gs_gather_local_max(gs, u%x, n)
    end select
    
  end subroutine gs_gather

  !> Gather kernel for addition of local data
  subroutine gs_gather_local_add(gs, u, n)
    type(gs_t), intent(inout) :: gs
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, intent(inout) :: n
    integer :: i
    gs%local_gs = 0d0
    do i = 1, gs%nlocal
       gs%local_gs(gs%local_dof_gs(i)) = &
            gs%local_gs(gs%local_dof_gs(i)) + u(gs%local_gs_dof(i))
    end do
  end subroutine gs_gather_local_add

  !> Gather kernel for multiplication of local data
  subroutine gs_gather_local_mul(gs, u, n)
    type(gs_t), intent(inout) :: gs
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, intent(inout) :: n
    integer :: i
    do i = 1, gs%nlocal
       gs%local_gs(gs%local_dof_gs(i)) = &
            gs%local_gs(gs%local_dof_gs(i)) * u(gs%local_gs_dof(i))
    end do
  end subroutine gs_gather_local_mul
  
  !> Gather kernel for minimum of local data
  subroutine gs_gather_local_min(gs, u, n)
    type(gs_t), intent(inout) :: gs
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, intent(inout) :: n
    integer :: i
    do i = 1, gs%nlocal
       gs%local_gs(gs%local_dof_gs(i)) = &
            min(gs%local_gs(gs%local_dof_gs(i)), u(gs%local_gs_dof(i)))
    end do
  end subroutine gs_gather_local_min

    !> Gather kernel for maximum of local data
  subroutine gs_gather_local_max(gs, u, n)
    type(gs_t), intent(inout) :: gs
    real(kind=dp), dimension(n), intent(inout) :: u
    integer, intent(inout) :: n
    integer :: i
    do i = 1, gs%nlocal
       gs%local_gs(gs%local_dof_gs(i)) = &
            max(gs%local_gs(gs%local_dof_gs(i)), u(gs%local_gs_dof(i)))
    end do
  end subroutine gs_gather_local_max

  !> Scatter kernel
  subroutine gs_scatter(gs, u)
    type(gs_t), intent(inout) :: gs
    type(field_t), intent(inout) :: u
    integer :: n
    
    n = u%msh%nelv * u%Xh%lx * u%Xh%ly * u%Xh%lz
    call gs_scatter_local(gs, u%x, n)

  end subroutine gs_scatter

  !> Scatter kernel for local data
  subroutine gs_scatter_local(gs, u, n)
    type(gs_t), intent(inout) :: gs
    real(kind=dp), dimension(n), intent(inout) :: u
    integer :: n
    integer :: i
    do i = 1, gs%nlocal
       u(gs%local_gs_dof(i)) = gs%local_gs(gs%local_dof_gs(i))
    end do
  end subroutine gs_scatter_local

  
end module gather_scatter
