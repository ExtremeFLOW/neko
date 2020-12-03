!> Defines a mapping of the degrees of freedom
!! @details A mapping defined based on a function space and a mesh
module dofmap
  use mesh
  use space
  use tuple
  use num_types
  use utils
  use fast3d
  use tensor
  use math
  implicit none
  private

  type, public :: dofmap_t
     integer(kind=8), allocatable :: dof(:,:,:,:)  !< Mapping to unique dof
     logical, allocatable :: shared_dof(:,:,:,:)   !< True if the dof is shared
     real(kind=dp), allocatable :: x(:,:,:,:)      !< Mapping to x-coordinates
     real(kind=dp), allocatable :: y(:,:,:,:)      !< Mapping to y-coordinates
     real(kind=dp), allocatable :: z(:,:,:,:)      !< Mapping to z-coordinates
     real(kind=dp), allocatable :: yinv(:,:,:,:)   !< Needed for axissymmetric case
     integer :: n_dofs                             !< Total number of dofs

     type(mesh_t), pointer :: msh
     type(space_t), pointer :: Xh

   contains
     procedure, pass(this) :: size => dofmap_size
     final :: dofmap_free
  end type dofmap_t

  interface dofmap_t
     module procedure dofmap_init
  end interface dofmap_t
  
contains

  function dofmap_init(msh, Xh) result(this)
    type(mesh_t), target, intent(inout) :: msh !< Mesh
    type(space_t), target, intent(inout) :: Xh !< Function space \f$ X_h \f$
    type(dofmap_t) :: this

    if ((msh%gdim .eq. 3 .and. Xh%lz .eq. 1) .or. &
         (msh%gdim .eq. 2 .and. Xh%lz .gt. 1)) then
       call neko_error("Invalid dimension of function space for the given mesh")
    end if


    call dofmap_free(this)

    this%msh => msh
    this%Xh => Xh

    this%n_dofs = Xh%lx* Xh%ly * Xh%lz * msh%nelv
        
    !
    ! Assign a unique id for all dofs
    ! 
    
    allocate(this%dof(Xh%lx, Xh%ly, Xh%lz, msh%nelv))
    allocate(this%shared_dof(Xh%lx, Xh%ly, Xh%lz, msh%nelv))

    this%dof = 0
    this%shared_dof = .false.

    !> @todo implement for 2d elements
    if (msh%gdim .eq. 3) then
       call dofmap_number_points(this)
       call dofmap_number_edges(this)
       call dofmap_number_facets(this)
    else
       call neko_warning('Dofmap not implemented for 2D (yet...)')
    end if

    !
    ! Generate x,y,z-coordinates for all dofs
    !
    
    allocate(this%x(Xh%lx, Xh%ly, Xh%lz, msh%nelv))
    allocate(this%y(Xh%lx, Xh%ly, Xh%lz, msh%nelv))
    allocate(this%z(Xh%lx, Xh%ly, Xh%lz, msh%nelv))
    allocate(this%yinv(Xh%lx, Xh%ly, Xh%lz, msh%nelv))
    
    this%x = 0d0
    this%y = 0d0
    this%z = 0d0
    !> @note should be intialised differently in acissymmetric case
    this%yinv = 1d0

    call dofmap_generate_xyz(this)    
    
  end function dofmap_init

  !> Deallocate the dofmap
  subroutine dofmap_free(this)
    type(dofmap_t), intent(inout) :: this

    if (allocated(this%dof)) then
       deallocate(this%dof)
    end if

    if (allocated(this%shared_dof)) then
       deallocate(this%shared_dof)
    end if

    if (allocated(this%x)) then
       deallocate(this%x)
    end if

    if (allocated(this%y)) then
       deallocate(this%y)
    end if

    if (allocated(this%z)) then
       deallocate(this%z)
    end if

    nullify(this%msh)
    nullify(this%Xh)
    
  end subroutine dofmap_free

  !> Return the total number of dofs in the dofmap
  pure function dofmap_size(this) result(res)
    class(dofmap_t), intent(in) :: this
    integer :: res
    res = this%n_dofs
  end function dofmap_size

  !> Assign numbers to each dofs on points
  subroutine dofmap_number_points(this)
    type(dofmap_t), target :: this
    integer :: i,j,k,l
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh

    msh => this%msh
    Xh => this%Xh

    do i = 1, msh%nelv
       this%dof(1, 1, 1, i) = &
            int(msh%elements(i)%e%pts(1)%p%id(), 8)
       this%dof(Xh%lx, 1, 1, i) = &
            int(msh%elements(i)%e%pts(2)%p%id(), 8)
       this%dof(1, Xh%ly, 1, i) = &
            int(msh%elements(i)%e%pts(4)%p%id(), 8)
       this%dof(Xh%lx, Xh%ly, 1, i) = &
            int(msh%elements(i)%e%pts(3)%p%id(), 8)

       this%dof(1, 1, Xh%lz, i) = &
            int(msh%elements(i)%e%pts(5)%p%id(), 8)
       this%dof(Xh%lx, 1, Xh%lz, i) = &
            int(msh%elements(i)%e%pts(6)%p%id(), 8)
       this%dof(1, Xh%ly, Xh%lz, i) = &
            int(msh%elements(i)%e%pts(8)%p%id(), 8)
       this%dof(Xh%lx, Xh%ly, Xh%lz, i) = &
            int(msh%elements(i)%e%pts(7)%p%id(), 8)

       this%shared_dof(1, 1, 1, i) = &
            mesh_is_shared(msh, msh%elements(i)%e%pts(1)%p)
       
       this%shared_dof(Xh%lx, 1, 1, i) = &
            mesh_is_shared(msh, msh%elements(i)%e%pts(2)%p)

       this%shared_dof(1, Xh%ly, 1, i) = &
            mesh_is_shared(msh, msh%elements(i)%e%pts(4)%p)
       
       this%shared_dof(Xh%lx, Xh%ly, 1, i) = &
            mesh_is_shared(msh, msh%elements(i)%e%pts(3)%p)

       this%shared_dof(1, 1, Xh%lz, i) = &
            mesh_is_shared(msh, msh%elements(i)%e%pts(5)%p)
       
       this%shared_dof(Xh%lx, 1, Xh%lz, i) = &
            mesh_is_shared(msh, msh%elements(i)%e%pts(6)%p)

       this%shared_dof(1, Xh%ly, Xh%lz, i) = &
            mesh_is_shared(msh, msh%elements(i)%e%pts(8)%p)
       
       this%shared_dof(Xh%lx, Xh%ly, Xh%lz, i) = &
            mesh_is_shared(msh, msh%elements(i)%e%pts(7)%p)
       
    end do
  end subroutine dofmap_number_points

  !> Assing numbers to dofs on edges
  subroutine dofmap_number_edges(this)
    type(dofmap_t), target :: this
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh
    integer :: i,j,k,l
    integer :: global_id        
    type(tuple_i4_t) :: edge
    integer(kind=8) :: num_dofs_edges(3) ! #dofs for each dir (r, s, t)
    integer(kind=8) :: edge_id, edge_offset
    logical :: shared_dof
        
    msh => this%msh
    Xh => this%Xh
    
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
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(1)
          !Reverse order of tranversal if edge is reversed
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(1,1,1,i)) k = Xh%lx+1-j
             this%dof(k, 1, 1, i) = edge_id
             this%shared_dof(k, 1, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 3)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(1,1,Xh%lz,i)) k = Xh%lx+1-j
             this%dof(k, 1, Xh%lz, i) = edge_id
             this%shared_dof(k, 1, Xh%lz, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 2)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(1,Xh%ly,1,i)) k = Xh%lx+1-j
             this%dof(k, Xh%ly, 1, i) = edge_id
             this%shared_dof(k, Xh%ly, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 4)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(1,Xh%ly,Xh%lz,i)) k = Xh%lx+1-j
             this%dof(k, Xh%ly, Xh%lz, i) = edge_id
             this%shared_dof(k, Xh%ly, Xh%lz, i) = shared_dof
             edge_id = edge_id + 1
          end do


          !
          ! Number edges in y-direction
          ! 
          call ep%edge_id(edge, 5)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(1,1,1,i)) k = Xh%ly+1-j
             this%dof(1, k, 1, i) = edge_id
             this%shared_dof(1, k, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 7)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(1,1,Xh%lz,i)) k = Xh%ly+1-j
             this%dof(1, k, Xh%lz, i) = edge_id
             this%shared_dof(1, k, Xh%lz, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 6)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(Xh%lx,1,1,i)) k = Xh%ly+1-j
             this%dof(Xh%lx, k, 1, i) = edge_id
             this%shared_dof(Xh%lx, k, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 8)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(Xh%lx,1,Xh%lz,i)) k = Xh%lz+1-j 
             this%dof(Xh%lx, k, Xh%lz, i) = edge_id
             this%shared_dof(Xh%lx, k, Xh%lz, i) = shared_dof
             edge_id = edge_id + 1
          end do

          !
          ! Number edges in z-direction
          ! 
          call ep%edge_id(edge, 9)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(1,1,1,i)) k = Xh%lz+1-j 
             this%dof(1, 1, k, i) = edge_id
             this%shared_dof(1, 1, k, i) = shared_dof
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 10)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(Xh%lx,1,1,i)) k = Xh%lz+1-j 
             this%dof(Xh%lx, 1, k, i) = edge_id
             this%shared_dof(Xh%lx, 1, k, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 11)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(1,Xh%ly,1,i)) k = Xh%lz+1-j 
             this%dof(1, Xh%ly, k, i) = edge_id
             this%shared_dof(1, Xh%ly, k, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 12)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), 8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             k = j
             if(int(edge%x(1),8) .ne. this%dof(Xh%lx,Xh%ly,1,i)) k = Xh%lz+1-j 
             this%dof(Xh%lx, Xh%ly, k, i) = edge_id
             this%shared_dof(Xh%lx, Xh%ly, k, i) = shared_dof
             edge_id = edge_id + 1
          end do
          
       end select
       
    end do
  end subroutine dofmap_number_edges

  !> Assign numbers to dofs on facets
  subroutine dofmap_number_facets(this)
    type(dofmap_t), target :: this
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh    
    integer :: i,j,k,l
    integer :: global_id
    type(tuple4_i4_t) :: face, face_order
    integer(kind=8) :: num_dofs_faces(3) ! #dofs for each dir (r, s, t)        
    integer(kind=8) :: facet_offset, facet_id
    logical :: shared_dof        

    msh => this%msh
    Xh => this%Xh

    !> @todo don't assume lx = ly = lz
    facet_offset = int(msh%glb_mpts, 8) + &
         int(msh%glb_meds, 8) * int(Xh%lx-2, 8) + int(1,8)

    ! Number of dofs on an edge excluding end-points
    num_dofs_faces(1) =  int((Xh%ly - 2) * (Xh%lz - 2), 8)
    num_dofs_faces(2) =  int((Xh%lx - 2) * (Xh%lz - 2), 8)
    num_dofs_faces(3) =  int((Xh%lx - 2) * (Xh%ly - 2), 8)

    do i = 1, msh%nelv
       
       !
       ! Number facets in x-direction (s, t)-plane
       !
       call msh%elements(i)%e%facet_id(face, 1)
       call msh%elements(i)%e%facet_order(face_order, 1)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_faces(1)
       do k = 2, Xh%lz -1
          do j = 2, Xh%ly - 1
             this%dof(1, j, k, i) = &
                  dofmap_facetidx(face_order,face,facet_id,j,k,Xh%lz,Xh%ly)
             this%shared_dof(1, j, k, i) = shared_dof
          end do
       end do
       
       call msh%elements(i)%e%facet_id(face, 2)
       call msh%elements(i)%e%facet_order(face_order, 2)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_faces(1)
       do k = 2, Xh%lz -1
          do j = 2, Xh%ly - 1
             this%dof(Xh%lx, j, k, i) = &
                  dofmap_facetidx(face_order,face,facet_id,j,k,Xh%lz,Xh%ly)
             this%shared_dof(Xh%lx, j, k, i) = shared_dof
          end do
       end do


       !
       ! Number facets in y-direction (r, t)-plane
       !
       call msh%elements(i)%e%facet_id(face, 3)
       call msh%elements(i)%e%facet_order(face_order, 3)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_faces(2)
       do k = 2, Xh%lz - 1
          do j = 2, Xh%lx - 1
             this%dof(j, 1, k, i) = &
                  dofmap_facetidx(face_order,face,facet_id,k,j,Xh%lz,Xh%lx)
             this%shared_dof(j, 1, k, i) = shared_dof
          end do
       end do
       
       call msh%elements(i)%e%facet_id(face, 4)
       call msh%elements(i)%e%facet_order(face_order, 4)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_faces(2)
       do k = 2, Xh%lz - 1
          do j = 2, Xh%lx - 1
             this%dof(j, Xh%ly, k, i) = &
                  dofmap_facetidx(face_order,face,facet_id,k,j,Xh%lz,Xh%lx)
             this%shared_dof(j, Xh%ly, k, i) = shared_dof
          end do
       end do


       !
       ! Number facets in z-direction (r, s)-plane
       !
       call msh%elements(i)%e%facet_id(face, 5)
       call msh%elements(i)%e%facet_order(face_order, 5)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_faces(3)
       do k = 2, Xh%ly - 1
          do j = 2, Xh%lx - 1
             this%dof(j, k, 1, i) = &
                  dofmap_facetidx(face_order,face,facet_id,k,j,Xh%ly,Xh%lx)
             this%shared_dof(j, k, 1, i) = shared_dof
          end do
       end do
       
       call msh%elements(i)%e%facet_id(face, 6)
       call msh%elements(i)%e%facet_order(face_order, 6)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), 8) * num_dofs_faces(3)
       do k = 2, Xh%ly - 1
          do j = 2, Xh%lx - 1
             this%dof(j, k, Xh%lz, i) = &
                  dofmap_facetidx(face_order,face,facet_id,k,j,Xh%lz,Xh%lx)
             this%shared_dof(j, k, Xh%lz, i) = shared_dof
          end do
       end do
    end do

  end subroutine dofmap_number_facets

 !> Get idx for GLL point on face depending on face ordering k and j
  function dofmap_facetidx(face_order, face, facet_id, k1, j1, lk1, lj1) result(facet_idx)
     type(tuple4_i4_t) :: face_order, face
     integer(8) :: facet_idx, facet_id
     integer :: k1, j1, lk1, lj1
     integer :: k,j,lk,lj
     
     k = k1 - 2
     j = j1 - 2
     lk = lk1 - 2
     lj = lj1 - 2

     ! Given the indexes k,j for a GLL point on the inner part of the
     ! face, we assign a unique number to it that depends on the
     ! corner with the lowest id and its neighbour with the lowest
     ! id. The id is assigned in this way to be consistent regardless
     ! of how the faces are rotated or mirrored.
     !
     !   4 -------- 3
     !     |      |      k
     !     |----->|      ^
     !     |----->|      |
     !     |----->|      |
     !   1 -------- 2    0--->j
     

     if(face_order%x(1) .eq. face%x(1)) then
       if(face_order%x(2) .lt. face_order%x(4)) then
         facet_idx = facet_id + j + k*lj 
       else
         facet_idx = facet_id + j*lk + k 
       endif
     else  if(face_order%x(2) .eq. face%x(1)) then
       if(face_order%x(3) .lt. face_order%x(1)) then
         facet_idx = facet_id + lk*(lj-1-j) + k 
       else
         facet_idx = facet_id + (lj-1-j) + k*lj
       endif
     else if(face_order%x(3) .eq. face%x(1)) then
       if(face_order%x(4) .lt. face_order%x(2)) then
         facet_idx = facet_id + (lj-1-j) + lj*(lk-1-k) 
       else
         facet_idx = facet_id + lk*(lj-1-j) + (lk-1-k)
       endif
     else if(face_order%x(4) .eq. face%x(1)) then
       if(face_order%x(1) .lt. face_order%x(3)) then
         facet_idx = facet_id + lk*j + (lk-1-k) 
       else
         facet_idx = facet_id + j + lj*(lk-1-k)
       endif     
     endif 

  end function dofmap_facetidx


  !> Generate x,y,z-coordinates for all dofs
  !! @note Assumes \f$ X_h_x = X_h_y = X_h_z \f$
  subroutine dofmap_generate_xyz(this)
    type(dofmap_t), target :: this
    integer :: i,j,k
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh
    real(kind=dp) :: xyzb(2,2,2,3), zgml(this%Xh%lx, 3)
    real(kind=dp) :: jx(this%Xh%lx*2), jy(this%Xh%lx*2), jz(this%Xh%lx*2)
    real(kind=dp) :: jxt(this%Xh%lx*2), jyt(this%Xh%lx*2), jzt(this%Xh%lx*2)
    real(kind=dp) :: w(4*this%Xh%lx**3), tmp(this%Xh%lx, this%Xh%lx, this%Xh%lx)
    real(kind=dp), dimension(2), parameter :: zlin = (/-1d0, 1d0/)
    

    msh => this%msh
    Xh => this%Xh
    zgml = 0d0
    xyzb = 0d0
    do i = 1, msh%nelv       

       w = 0d0
       call copy(zgml(1,1), Xh%zg(1,1), Xh%lx)                               
       call copy(zgml(1,2), Xh%zg(1,2), Xh%ly)
       if (msh%gdim .gt. 2) then
          call copy(zgml(1,3), Xh%zg(1,3), Xh%lz)
       end if

       k = 1
       do j = 1, Xh%lx
          call fd_weights_full(zgml(j,1),zlin,1,0,jxt(k))
          call fd_weights_full(zgml(j,2),zlin,1,0,jyt(k))
          if (msh%gdim .gt. 2) then
             call fd_weights_full(zgml(j,3),zlin,1,0,jzt(k))
          end if
          k = k + 2
       end do
       call trsp(jx, Xh%lx, jxt, 2)

       if (msh%gdim .eq. 2) then
          jzt = 1d0
       end if

       do j = 1, msh%gdim
          xyzb(1,1,1,j) = msh%elements(i)%e%pts(1)%p%x(j)
          xyzb(2,1,1,j) = msh%elements(i)%e%pts(2)%p%x(j)
          xyzb(1,2,1,j) = msh%elements(i)%e%pts(4)%p%x(j)
          xyzb(2,2,1,j) = msh%elements(i)%e%pts(3)%p%x(j)

          if (msh%gdim .gt. 2) then
             xyzb(1,1,2,j) = msh%elements(i)%e%pts(5)%p%x(j)
             xyzb(2,1,2,j) = msh%elements(i)%e%pts(6)%p%x(j)
             xyzb(1,2,2,j) = msh%elements(i)%e%pts(8)%p%x(j)
             xyzb(2,2,2,j) = msh%elements(i)%e%pts(7)%p%x(j)
          end if
       end do

       call tensr3(tmp, Xh%lx, xyzb(1,1,1,1), 2, jx, jyt, jzt, w)
       call copy(this%x(1,1,1,i), tmp, Xh%lx*Xh%ly*Xh%lz)
       call tensr3(tmp, Xh%ly, xyzb(1,1,1,2), 2, jx, jyt, jzt, w)
       call copy(this%y(1,1,1,i), tmp, Xh%lx*Xh%ly*Xh%lz)
       if (msh%gdim .gt. 2) then
          call tensr3(tmp, Xh%lz, xyzb(1,1,1,3), 2, jx, jyt, jzt, w)
          call copy(this%z(1,1,1,i), tmp, Xh%lx*Xh%ly*Xh%lz)
       end if
    end do
  end subroutine dofmap_generate_xyz
  
end module dofmap
