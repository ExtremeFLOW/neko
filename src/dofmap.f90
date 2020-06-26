module dofmap
  use mesh
  use space
  use tuple
  implicit none
  private

  type, public :: dofmap_t
     integer(kind=8), allocatable :: dof(:,:,:,:)
   contains
     final :: dofmap_free
  end type dofmap_t

  interface dofmap_t
     module procedure dofmap_init
  end interface dofmap_t
  
contains

  function dofmap_init(msh, Xh) result(this)
    type(mesh_t), intent(inout) :: msh !< Mesh
    type(space_t), intent(in) :: Xh !< Function space \f$ X_h \f$
    type(dofmap_t) :: this
    type(tuple_i4_t) :: edge
    integer :: i,j,k,l
    integer :: num_dofs_edges(3) ! #dofs for each dir (r, s, t)
    integer :: local_id
    integer(kind=8) :: edge_id, edge_offset
    
    call dofmap_free(this)

    allocate(this%dof(Xh%lx, Xh%ly, Xh%lz, msh%nelv))

    this%dof = 0

    !> @todo implement for 2d elements

    ! Assign numbers to each dofs on points
    do i = 1, msh%nelv
       this%dof(1, 1, 1, i) = &
            int(msh%elements(i)%e%pts(1)%p%id(), 8)
       this%dof(Xh%lx, 1, 1, i) = &
            int(msh%elements(i)%e%pts(2)%p%id(), 8)
       this%dof(1, Xh%ly, 1, i) = &
            int(msh%elements(i)%e%pts(4)%p%id(), 8)
       this%dof(Xh%lz, Xh%ly, 1, i) = &
            int(msh%elements(i)%e%pts(3)%p%id(), 8)

       this%dof(1, 1, Xh%lz, i) = &
            int(msh%elements(i)%e%pts(5)%p%id(), 8)
       this%dof(Xh%lx, 1, Xh%lz, i) = &
            int(msh%elements(i)%e%pts(6)%p%id(), 8)
       this%dof(1, Xh%ly, Xh%lz, i) = &
            int(msh%elements(i)%e%pts(8)%p%id(), 8)
       this%dof(Xh%lz, Xh%ly, Xh%lz, i) = &
            int(msh%elements(i)%e%pts(7)%p%id(), 8)
    end do
    
    !
    ! Assign numbers to dofs on edges
    !

    ! Number of dofs on an edge excluding end-points
    num_dofs_edges(1) =  Xh%lx - 2
    num_dofs_edges(2) =  Xh%ly - 2
    num_dofs_edges(3) =  Xh%lz - 2

    edge_offset = int(msh%mpts, 8) + int(1, 8) !< @todo Only serial for now...
    
    do i = 1, msh%nelv
       
       select type(ep=>msh%elements(1)%e)
       type is (hex_t)

          !
          ! Number edges in x-direction
          ! 
          call ep%edge_id(edge, 1)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 2
             this%dof(j, 1, 1, i) = edge_id
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 3)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 2
             this%dof(j, 1, Xh%lz, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 2)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 2
             this%dof(j, Xh%ly, 1, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 4)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 2
             this%dof(j, Xh%ly, Xh%lz, i) = edge_id
             edge_id = edge_id + 1
          end do


          !
          ! Number edges in y-direction
          ! 
          call ep%edge_id(edge, 5)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 2
             this%dof(1, j, 1, i) = edge_id
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 7)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 2
             this%dof(1, j, Xh%lz, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 6)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 2
             this%dof(Xh%lx, j, 1, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 8)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 2
             this%dof(Xh%lx, j, Xh%lz, i) = edge_id
             edge_id = edge_id + 1
          end do

          !
          ! Number edges in z-direction
          ! 
          call ep%edge_id(edge, 9)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 2
             this%dof(1, 1, j, i) = edge_id
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 10)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 2
             this%dof(Xh%lx, 1, j, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 11)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 2
             this%dof(1, Xh%ly, j, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 12)
          local_id = mesh_get_local_edge(msh, edge)
          edge_id = edge_offset + int((local_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 2
             this%dof(Xh%lx, Xh%ly, j, i) = edge_id
             edge_id = edge_id + 1
          end do
          
       end select
       
    end do

    ! Assign numbers to dofs on facets 
    
  end function dofmap_init


  subroutine dofmap_free(this)
    type(dofmap_t), intent(inout) :: this

    if (allocated(this%dof)) then
       deallocate(this%dof)
    end if
    
  end subroutine dofmap_free
  
end module dofmap
