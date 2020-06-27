!> Defines a mapping of the degrees of freedom
!! @details A mapping defined based on a function space and a mesh
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
    type(tuple4_i4_t) :: face
    integer :: i,j,k,l
    integer(kind=8) :: num_dofs_edges(3) ! #dofs for each dir (r, s, t)
    integer(kind=8) :: num_dofs_faces(3) ! #dofs for each dir (r, s, t)
    integer :: global_id
    integer(kind=8) :: edge_id, edge_offset, facet_offset, facet_id
    
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
    num_dofs_edges(1) =  int(Xh%lx - 2, 8)
    num_dofs_edges(2) =  int(Xh%ly - 2, 8)
    num_dofs_edges(3) =  int(Xh%lz - 2, 8)

    edge_offset = int(msh%glb_mpts, 8) + int(1, 8)

    do i = 1, msh%nelv
       
       select type(ep=>msh%elements(i)%e)
       type is (hex_t)

          !
          ! Number edges in x-direction
          ! 
          call ep%edge_id(edge, 1)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             this%dof(j, 1, 1, i) = edge_id
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 3)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             this%dof(j, 1, Xh%lz, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 2)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             this%dof(j, Xh%ly, 1, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 4)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             this%dof(j, Xh%ly, Xh%lz, i) = edge_id
             edge_id = edge_id + 1
          end do


          !
          ! Number edges in y-direction
          ! 
          call ep%edge_id(edge, 5)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             this%dof(1, j, 1, i) = edge_id
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 7)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             this%dof(1, j, Xh%lz, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 6)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             this%dof(Xh%lx, j, 1, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 8)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             this%dof(Xh%lx, j, Xh%lz, i) = edge_id
             edge_id = edge_id + 1
          end do

          !
          ! Number edges in z-direction
          ! 
          call ep%edge_id(edge, 9)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             this%dof(1, 1, j, i) = edge_id
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 10)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             this%dof(Xh%lx, 1, j, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 11)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             this%dof(1, Xh%ly, j, i) = edge_id
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 12)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             this%dof(Xh%lx, Xh%ly, j, i) = edge_id
             edge_id = edge_id + 1
          end do
          
       end select
       
    end do

    ! Assign numbers to dofs on facets

    !> @todo don't assume lx = ly = lz
    facet_offset = int(msh%glb_mpts, 8) + &
         int(msh%glb_meds, 8) * int(Xh%lx-2, 8)

    ! Number of dofs on an edge excluding end-points
    num_dofs_faces(1) =  int((Xh%ly - 2) * (Xh%lz - 2), 8)
    num_dofs_faces(2) =  int((Xh%lx - 2) * (Xh%lz - 2), 8)
    num_dofs_faces(3) =  int((Xh%lx - 2) * (Xh%ly - 2), 8)

    do i = 1, msh%nelv
       
       !
       ! Number facets in x-direction (s, t)-plane
       !
       call msh%elements(i)%e%facet_id(face, 1)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_edges(1)
       do k = 2, Xh%lz -1
          do j = 2, Xh%ly - 1
             this%dof(1, j, k, i) = facet_id
             facet_id = facet_id + 1
          end do
       end do
       
       call msh%elements(i)%e%facet_id(face, 2)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_edges(1)
       do k = 2, Xh%lz -1
          do j = 2, Xh%ly - 1
             this%dof(Xh%lx, j, k, i) = facet_id
             facet_id = facet_id + 1
          end do
       end do


       !
       ! Number facets in y-direction (r, t)-plane
       !
       call msh%elements(i)%e%facet_id(face, 3)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_edges(2)
       do k = 2, Xh%lz - 1
          do j = 2, Xh%lx - 1
             this%dof(j, 1, k, i) = facet_id
             facet_id = facet_id + 1
          end do
       end do
       
       call msh%elements(i)%e%facet_id(face, 4)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_edges(2)
       do k = 2, Xh%lz - 1
          do j = 2, Xh%lx - 1
             this%dof(j, Xh%ly, k, i) = facet_id
             facet_id = facet_id + 1
          end do
       end do


       !
       ! Number facets in z-direction (r, s)-plane
       !
       call msh%elements(i)%e%facet_id(face, 5)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_edges(3)
       do k = 2, Xh%ly - 1
          do j = 2, Xh%lx - 1
             this%dof(j, k, 1, i) = facet_id
             facet_id = facet_id + 1
          end do
       end do
       
       call msh%elements(i)%e%facet_id(face, 6)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_edges(3)
       do k = 2, Xh%lz - 1
          do j = 2, Xh%lx - 1
             this%dof(j, k, Xh%lz, i) = facet_id
             facet_id = facet_id + 1
          end do
       end do
    end do

  end function dofmap_init


  subroutine dofmap_free(this)
    type(dofmap_t), intent(inout) :: this

    if (allocated(this%dof)) then
       deallocate(this%dof)
    end if
    
  end subroutine dofmap_free
  
end module dofmap
