! Copyright (c) 2020-2021, The Neko Authors
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
!> Defines a mapping of the degrees of freedom
!! @details A mapping defined based on a function space and a mesh
module dofmap
  use neko_config
  use mesh
  use space
  use tuple
  use num_types
  use utils
  use fast3d
  use tensor
  use device
  use math
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public :: dofmap_t
     integer(kind=i8), allocatable :: dof(:,:,:,:)  !< Mapping to unique dof
     logical, allocatable :: shared_dof(:,:,:,:)   !< True if the dof is shared
     real(kind=rp), allocatable :: x(:,:,:,:)      !< Mapping to x-coordinates
     real(kind=rp), allocatable :: y(:,:,:,:)      !< Mapping to y-coordinates
     real(kind=rp), allocatable :: z(:,:,:,:)      !< Mapping to z-coordinates
     integer :: n_dofs                             !< Total number of dofs

     type(mesh_t), pointer :: msh
     type(space_t), pointer :: Xh

     !
     ! Device pointers (if present)
     !
     type(c_ptr) :: x_d = C_NULL_PTR
     type(c_ptr) :: y_d = C_NULL_PTR
     type(c_ptr) :: z_d = C_NULL_PTR

   contains
     procedure, pass(this) :: size => dofmap_size
!     final :: dofmap_free
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
       call dofmap_number_faces(this)
    else
       call dofmap_number_points(this)
       call dofmap_number_edges(this)
    end if

    !
    ! Generate x,y,z-coordinates for all dofs
    !
    
    allocate(this%x(Xh%lx, Xh%ly, Xh%lz, msh%nelv))
    allocate(this%y(Xh%lx, Xh%ly, Xh%lz, msh%nelv))
    allocate(this%z(Xh%lx, Xh%ly, Xh%lz, msh%nelv))
    
    this%x = 0d0
    this%y = 0d0
    this%z = 0d0
    !> @note should be intialised differently in acissymmetric case

    call dofmap_generate_xyz(this)    

    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
        (NEKO_BCKND_OPENCL .eq. 1)) then 
       call device_map(this%x, this%x_d, this%n_dofs)
       call device_map(this%y, this%y_d, this%n_dofs)
       call device_map(this%z, this%z_d, this%n_dofs)

       call device_memcpy(this%x, this%x_d, this%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(this%y, this%y_d, this%n_dofs, HOST_TO_DEVICE)
       call device_memcpy(this%z, this%z_d, this%n_dofs, HOST_TO_DEVICE)
    end if
    
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

    !
    ! Cleanup the device (if present)
    !
    if (c_associated(this%x_d)) then
       call device_free(this%x_d)
    end if

    if (c_associated(this%y_d)) then
       call device_free(this%y_d)
    end if

    if (c_associated(this%z_d)) then
       call device_free(this%z_d)
    end if
    
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
    integer :: i
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh

    msh => this%msh
    Xh => this%Xh
    if (msh%gdim .eq. 2) then
       do i = 1, msh%nelv
          this%dof(1, 1, 1, i) = &
               int(msh%elements(i)%e%pts(1)%p%id(), i8)
          this%dof(Xh%lx, 1, 1, i) = &
               int(msh%elements(i)%e%pts(2)%p%id(), i8)
          this%dof(1, Xh%ly, 1, i) = &
               int(msh%elements(i)%e%pts(4)%p%id(), i8)
          this%dof(Xh%lx, Xh%ly, 1, i) = &
               int(msh%elements(i)%e%pts(3)%p%id(), i8)

          this%shared_dof(1, 1, 1, i) = &
               mesh_is_shared(msh, msh%elements(i)%e%pts(1)%p)
          
          this%shared_dof(Xh%lx, 1, 1, i) = &
               mesh_is_shared(msh, msh%elements(i)%e%pts(2)%p)

          this%shared_dof(1, Xh%ly, 1, i) = &
               mesh_is_shared(msh, msh%elements(i)%e%pts(4)%p)
          
          this%shared_dof(Xh%lx, Xh%ly, 1, i) = &
               mesh_is_shared(msh, msh%elements(i)%e%pts(3)%p)
       end do
    else    
       do i = 1, msh%nelv
          this%dof(1, 1, 1, i) = &
               int(msh%elements(i)%e%pts(1)%p%id(), i8)
          this%dof(Xh%lx, 1, 1, i) = &
               int(msh%elements(i)%e%pts(2)%p%id(), i8)
          this%dof(1, Xh%ly, 1, i) = &
               int(msh%elements(i)%e%pts(4)%p%id(), i8)
          this%dof(Xh%lx, Xh%ly, 1, i) = &
               int(msh%elements(i)%e%pts(3)%p%id(), i8)

          this%dof(1, 1, Xh%lz, i) = &
               int(msh%elements(i)%e%pts(5)%p%id(), i8)
          this%dof(Xh%lx, 1, Xh%lz, i) = &
               int(msh%elements(i)%e%pts(6)%p%id(), i8)
          this%dof(1, Xh%ly, Xh%lz, i) = &
               int(msh%elements(i)%e%pts(8)%p%id(), i8)
          this%dof(Xh%lx, Xh%ly, Xh%lz, i) = &
               int(msh%elements(i)%e%pts(7)%p%id(), i8)

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
    end if
  end subroutine dofmap_number_points

  !> Assing numbers to dofs on edges
  subroutine dofmap_number_edges(this)
    type(dofmap_t), target :: this
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh
    integer :: i,j,k
    integer :: global_id        
    type(tuple_i4_t) :: edge
    integer(kind=i8) :: num_dofs_edges(3) ! #dofs for each dir (r, s, t)
    integer(kind=i8) :: edge_id, edge_offset
    logical :: shared_dof
        
    msh => this%msh
    Xh => this%Xh
    
    ! Number of dofs on an edge excluding end-points
    num_dofs_edges(1) =  int(Xh%lx - 2, i8)
    num_dofs_edges(2) =  int(Xh%ly - 2, i8)
    num_dofs_edges(3) =  int(Xh%lz - 2, i8)
    edge_offset = int(msh%glb_mpts, i8) + int(1, 4)

    do i = 1, msh%nelv
       
       select type(ep=>msh%elements(i)%e)
       type is (hex_t)
          !
          ! Number edges in r-direction
          ! 
          call ep%edge_id(edge, 1)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(1)
          !Reverse order of tranversal if edge is reversed
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,1,1,i)) k = Xh%lx+1-j
             this%dof(k, 1, 1, i) = edge_id
             this%shared_dof(k, 1, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 3)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,1,Xh%lz,i)) k = Xh%lx+1-j
             this%dof(k, 1, Xh%lz, i) = edge_id
             this%shared_dof(k, 1, Xh%lz, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 2)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,Xh%ly,1,i)) k = Xh%lx+1-j
             this%dof(k, Xh%ly, 1, i) = edge_id
             this%shared_dof(k, Xh%ly, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 4)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,Xh%ly,Xh%lz,i)) k = Xh%lx+1-j
             this%dof(k, Xh%ly, Xh%lz, i) = edge_id
             this%shared_dof(k, Xh%ly, Xh%lz, i) = shared_dof
             edge_id = edge_id + 1
          end do


          !
          ! Number edges in s-direction
          ! 
          call ep%edge_id(edge, 5)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,1,1,i)) k = Xh%ly+1-j
             this%dof(1, k, 1, i) = edge_id
             this%shared_dof(1, k, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 7)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,1,Xh%lz,i)) k = Xh%ly+1-j
             this%dof(1, k, Xh%lz, i) = edge_id
             this%shared_dof(1, k, Xh%lz, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 6)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(Xh%lx,1,1,i)) k = Xh%ly+1-j
             this%dof(Xh%lx, k, 1, i) = edge_id
             this%shared_dof(Xh%lx, k, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, i8)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(Xh%lx,1,Xh%lz,i)) k = Xh%lz+1-j 
             this%dof(Xh%lx, k, Xh%lz, i) = edge_id
             this%shared_dof(Xh%lx, k, Xh%lz, i) = shared_dof
             edge_id = edge_id + 1
          end do

          !
          ! Number edges in t-direction
          ! 
          call ep%edge_id(edge, 9)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,1,1,i)) k = Xh%lz+1-j 
             this%dof(1, 1, k, i) = edge_id
             this%shared_dof(1, 1, k, i) = shared_dof
             edge_id = edge_id + 1
          end do
          
          call ep%edge_id(edge, 10)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(Xh%lx,1,1,i)) k = Xh%lz+1-j 
             this%dof(Xh%lx, 1, k, i) = edge_id
             this%shared_dof(Xh%lx, 1, k, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 11)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,Xh%ly,1,i)) k = Xh%lz+1-j 
             this%dof(1, Xh%ly, k, i) = edge_id
             this%shared_dof(1, Xh%ly, k, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%edge_id(edge, 12)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(3)
          do j = 2, Xh%lz - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(Xh%lx,Xh%ly,1,i)) k = Xh%lz+1-j 
             this%dof(Xh%lx, Xh%ly, k, i) = edge_id
             this%shared_dof(Xh%lx, Xh%ly, k, i) = shared_dof
             edge_id = edge_id + 1
          end do
       type is (quad_t)
          !
          ! Number edges in r-direction
          ! 
          call ep%facet_id(edge, 3)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(1)
          !Reverse order of tranversal if edge is reversed
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,1,1,i)) k = Xh%lx+1-j
             this%dof(k, 1, 1, i) = edge_id
             this%shared_dof(k, 1, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%facet_id(edge, 4)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(1)
          do j = 2, Xh%lx - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,Xh%ly,1,i)) k = Xh%lx+1-j
             this%dof(k, Xh%ly, 1, i) = edge_id
             this%shared_dof(k, Xh%ly, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do

          !
          ! Number edges in s-direction
          ! 
          call ep%facet_id(edge, 1)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(1,1,1,i)) k = Xh%ly+1-j
             this%dof(1, k, 1, i) = edge_id
             this%shared_dof(1, k, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do

          call ep%facet_id(edge, 2)
          shared_dof = mesh_is_shared(msh, edge)
          global_id = mesh_get_global(msh, edge)
          edge_id = edge_offset + int((global_id - 1), i8) * num_dofs_edges(2)
          do j = 2, Xh%ly - 1
             k = j
             if(int(edge%x(1), i8) .ne. this%dof(Xh%lx,1,1,i)) k = Xh%ly+1-j
             this%dof(Xh%lx, k, 1, i) = edge_id
             this%shared_dof(Xh%lx, k, 1, i) = shared_dof
             edge_id = edge_id + 1
          end do
          
       end select
       
    end do
  end subroutine dofmap_number_edges

  !> Assign numbers to dofs on faces
  subroutine dofmap_number_faces(this)
    type(dofmap_t), target :: this
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh    
    integer :: i,j,k
    integer :: global_id
    type(tuple4_i4_t) :: face, face_order
    integer(kind=i8) :: num_dofs_faces(3) ! #dofs for each dir (r, s, t)        
    integer(kind=i8) :: facet_offset, facet_id
    logical :: shared_dof        

    msh => this%msh
    Xh => this%Xh

    !> @todo don't assume lx = ly = lz
    facet_offset = int(msh%glb_mpts, i8) + &
         int(msh%glb_meds, i8) * int(Xh%lx-2, i8) + int(1, i8)

    ! Number of dofs on an edge excluding end-points
    num_dofs_faces(1) =  int((Xh%ly - 2) * (Xh%lz - 2), i8)
    num_dofs_faces(2) =  int((Xh%lx - 2) * (Xh%lz - 2), i8)
    num_dofs_faces(3) =  int((Xh%lx - 2) * (Xh%ly - 2), i8)

    do i = 1, msh%nelv
       
       !
       ! Number facets in r-direction (s, t)-plane
       !
       call msh%elements(i)%e%facet_id(face, 1)
       call msh%elements(i)%e%facet_order(face_order, 1)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), i8) * num_dofs_faces(1)
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
       facet_id = facet_offset + int((global_id - 1), i8) * num_dofs_faces(1)
       do k = 2, Xh%lz -1
          do j = 2, Xh%ly - 1
             this%dof(Xh%lx, j, k, i) = &
                  dofmap_facetidx(face_order,face,facet_id,j,k,Xh%lz,Xh%ly)
             this%shared_dof(Xh%lx, j, k, i) = shared_dof
          end do
       end do


       !
       ! Number facets in s-direction (r, t)-plane
       !
       call msh%elements(i)%e%facet_id(face, 3)
       call msh%elements(i)%e%facet_order(face_order, 3)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), i8) * num_dofs_faces(2)
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
       facet_id = facet_offset + int((global_id - 1), i8) * num_dofs_faces(2)
       do k = 2, Xh%lz - 1
          do j = 2, Xh%lx - 1
             this%dof(j, Xh%ly, k, i) = &
                  dofmap_facetidx(face_order,face,facet_id,k,j,Xh%lz,Xh%lx)
             this%shared_dof(j, Xh%ly, k, i) = shared_dof
          end do
       end do


       !
       ! Number facets in t-direction (r, s)-plane
       !
       call msh%elements(i)%e%facet_id(face, 5)
       call msh%elements(i)%e%facet_order(face_order, 5)
       shared_dof = mesh_is_shared(msh, face)
       global_id = mesh_get_global(msh, face)
       facet_id = facet_offset + int((global_id - 1), i8) * num_dofs_faces(3)
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
       facet_id = facet_offset + int((global_id - 1), i8) * num_dofs_faces(3)
       do k = 2, Xh%ly - 1
          do j = 2, Xh%lx - 1
             this%dof(j, k, Xh%lz, i) = &
                  dofmap_facetidx(face_order,face,facet_id,k,j,Xh%lz,Xh%lx)
             this%shared_dof(j, k, Xh%lz, i) = shared_dof
          end do
       end do
    end do

  end subroutine dofmap_number_faces

 !> Get idx for GLL point on face depending on face ordering k and j
  function dofmap_facetidx(face_order, face, facet_id, k1, j1, lk1, lj1) result(facet_idx)
     type(tuple4_i4_t) :: face_order, face
     integer(kind=i8) :: facet_idx, facet_id
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
      end if
     else  if(face_order%x(2) .eq. face%x(1)) then
       if(face_order%x(3) .lt. face_order%x(1)) then
         facet_idx = facet_id + lk*(lj-1-j) + k 
       else
         facet_idx = facet_id + (lj-1-j) + k*lj
      end if
     else if(face_order%x(3) .eq. face%x(1)) then
       if(face_order%x(4) .lt. face_order%x(2)) then
         facet_idx = facet_id + (lj-1-j) + lj*(lk-1-k) 
       else
         facet_idx = facet_id + lk*(lj-1-j) + (lk-1-k)
      end if
     else if(face_order%x(4) .eq. face%x(1)) then
       if(face_order%x(1) .lt. face_order%x(3)) then
         facet_idx = facet_id + lk*j + (lk-1-k) 
       else
         facet_idx = facet_id + j + lj*(lk-1-k)
      end if
   end if

  end function dofmap_facetidx

  !> Generate x,y,z-coordinates for all dofs
  !! @note Assumes \f$ X_{h_x} = X_{h_y} = X_{h_z} \f$
  subroutine dofmap_generate_xyz(this)
    type(dofmap_t), target :: this
    integer :: i, j, el_idx
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh
    real(kind=rp) :: rp_curve_data(5), curve_data_tot(5,12)
    logical :: midpoint
    integer :: n_edge, curve_type(12)

    msh => this%msh
    Xh => this%Xh
    
    if (msh%gdim .eq. 3) then
       n_edge = 12
    else
       n_edge = 4
    end if

    do i = 1, msh%nelv       
       call dofmap_xyzlin(Xh,msh, msh%elements(i)%e,this%x(1,1,1,i),this%y(1,1,1,i), this%z(1,1,1,i)) 
    end do
    do i =1, msh%curve%size 
       midpoint = .false.
       el_idx = msh%curve%curve_el(i)%el_idx
       curve_type = msh%curve%curve_el(i)%curve_type
       curve_data_tot = msh%curve%curve_el(i)%curve_data
       do j = 1, n_edge
          if (curve_type(j) .eq. 4) then
             midpoint = .true.
          end if
       end do
       if (midpoint .and. Xh%lx .gt. 2) then
          call dofmap_xyzquad(Xh, msh, msh%elements(el_idx)%e, &
                              this%x(1,1,1,el_idx), this%y(1,1,1,el_idx),&
                              this%z(1,1,1,el_idx),curve_type, curve_data_tot)
       end if
    end do
    do i =1, msh%curve%size 
       el_idx = msh%curve%curve_el(i)%el_idx
       do j = 1, 8
          if (msh%curve%curve_el(i)%curve_type(j) .eq. 3) then
             rp_curve_data = msh%curve%curve_el(i)%curve_data(1:5,j)
             call arc_surface(j, rp_curve_data, &
                              this%x(1,1,1,el_idx), &
                              this%y(1,1,1,el_idx), &
                              this%z(1,1,1, el_idx), &
                              Xh, msh%elements(el_idx)%e, msh%gdim) 
          end if
       end do
    end do
    if (associated(msh%apply_deform)) then
       call msh%apply_deform(this%x, this%y, this%z, Xh%lx, Xh%ly, Xh%lz)
    end if
  end subroutine dofmap_generate_xyz

  subroutine dofmap_xyzlin(Xh, msh, element, x, y, z)
    type(mesh_t), pointer, intent(in) :: msh
    type(space_t), intent(in) :: Xh
    class(element_t), intent(in) :: element
    real(kind=rp), intent(inout) :: x(Xh%lx,Xh%ly,Xh%lz), y(Xh%lx,Xh%ly,Xh%lz), z(Xh%lx,Xh%ly,Xh%lz)
    real(kind=rp) :: xyzb(2,2,2,3), zgml(Xh%lx, 3)
    real(kind=rp) :: jx(Xh%lx*2)
    real(kind=rp) :: jxt(Xh%lx*2), jyt(Xh%lx*2), jzt(Xh%lx*2)
    real(kind=rp) :: w(4*Xh%lx**3), tmp(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp), dimension(2), parameter :: zlin = (/-1d0, 1d0/)
    
    integer :: j, k

    zgml = 0d0
    xyzb = 0d0
 
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
       xyzb(1,1,1,j) = element%pts(1)%p%x(j)
       xyzb(2,1,1,j) = element%pts(2)%p%x(j)
       xyzb(1,2,1,j) = element%pts(4)%p%x(j)
       xyzb(2,2,1,j) = element%pts(3)%p%x(j)

       if (msh%gdim .gt. 2) then
          xyzb(1,1,2,j) = element%pts(5)%p%x(j)
          xyzb(2,1,2,j) = element%pts(6)%p%x(j)
          xyzb(1,2,2,j) = element%pts(8)%p%x(j)
          xyzb(2,2,2,j) = element%pts(7)%p%x(j)
       end if
    end do
    if (msh%gdim .eq. 3) then
       call tensr3(tmp, Xh%lx, xyzb(1,1,1,1), 2, jx, jyt, jzt, w)
       call copy(x, tmp, Xh%lx*Xh%ly*Xh%lz)
       call tensr3(tmp, Xh%ly, xyzb(1,1,1,2), 2, jx, jyt, jzt, w)
       call copy(y, tmp, Xh%lx*Xh%ly*Xh%lz)
       call tensr3(tmp, Xh%lz, xyzb(1,1,1,3), 2, jx, jyt, jzt, w)
       call copy(z, tmp, Xh%lx*Xh%ly*Xh%lz)
    else
       call tnsr2d_el(tmp, Xh%lx, xyzb(1,1,1,1), 2, jx, jyt)
       call copy(x, tmp, Xh%lx*Xh%ly*Xh%lz)
       call tnsr2d_el(tmp, Xh%ly, xyzb(1,1,1,2), 2, jx, jyt)
       call copy(y, tmp, Xh%lx*Xh%ly*Xh%lz)
    end if
  end subroutine dofmap_xyzlin

  subroutine dofmap_xyzquad(Xh, msh, element, x, y, z, curve_type, curve_data)
    type(mesh_t), pointer, intent(in) :: msh
    type(space_t), intent(in) :: Xh
    class(element_t), intent(in) :: element
    real(kind=rp), dimension(Xh%lx, Xh%ly, Xh%lz), intent(inout) :: x, y, z
    integer :: curve_type(12), eindx(12)
    real(kind=rp) :: curve_data(5,12), x3(3,3,3), y3(3,3,3), z3(3,3,3)
    type(space_t), target :: Xh3
    real(kind=rp), dimension(3), parameter :: zquad = (/-1d0, 0d0,1d0/)
    real(kind=rp) :: zg(3)
    real(kind=rp), dimension(Xh%lx,Xh%lx,Xh%lx) :: tmp
    real(kind=rp) :: jx(Xh%lx*3)
    real(kind=rp) :: jxt(Xh%lx*3), jyt(Xh%lx*3), jzt(Xh%lx*3)
    real(kind=rp) :: w(4*Xh%lxyz,2)
    integer :: j, k, n_edges
    eindx = (/2 ,  6 ,  8 ,  4, &
              20 , 24 , 26 , 22, &
              10 , 12 , 18 , 16 /)

    w = 0d0
    if (msh%gdim .eq. 3) then
       n_edges = 12
       call space_init(Xh3, GLL, 3, 3, 3)
    else
       n_edges = 4
       call space_init(Xh3, GLL, 3, 3)
    end if
    call dofmap_xyzlin(Xh3, msh, element, x3, y3, z3)

    do k = 1, n_edges
       if (curve_type(k) .eq. 4) then
          x3(eindx(k),1,1) = curve_data(1,k)
          y3(eindx(k),1,1) = curve_data(2,k)
          z3(eindx(k),1,1) = curve_data(3,k)
       end if
    end do 
    zg(1) = -1
    zg(2) =  0 
    zg(3) =  1
    if (msh%gdim .eq. 3) then
       call gh_face_extend_3d(x3,zg,3,2,w(1,1),w(1,2)) ! 2 --> edge extend
       call gh_face_extend_3d(y3,zg,3,2,w(1,1),w(1,2))
       call gh_face_extend_3d(z3,zg,3,2,w(1,1),w(1,2))
    else
       call neko_warning(' m deformation not supported for 2d yet')
       call gh_face_extend_2d(x3,zg,3,2,w(1,1),w(1,2)) ! 2 --> edge extend
       call gh_face_extend_2d(y3,zg,3,2,w(1,1),w(1,2))
    end if
    k =1 
    do j = 1, Xh%lx
       call fd_weights_full(Xh%zg(j,1),zquad,2,0,jxt(k))
       call fd_weights_full(Xh%zg(j,2),zquad,2,0,jyt(k))
       if (msh%gdim .gt. 2) then
          call fd_weights_full(Xh%zg(j,3),zquad,2,0,jzt(k))
       end if
       k = k + 3
    end do
    call trsp(jx, Xh%lx, jxt, 3)
    if (msh%gdim .eq. 3) then
       call tensr3(tmp, Xh%lx, x3, 3, jx, jyt, jzt, w)
       call copy(x, tmp, Xh%lx*Xh%ly*Xh%lz)
       call tensr3(tmp, Xh%ly, y3, 3, jx, jyt, jzt, w)
       call copy(y, tmp, Xh%lx*Xh%ly*Xh%lz)
       call tensr3(tmp, Xh%lz, z3, 3, jx, jyt, jzt, w)
       call copy(z, tmp, Xh%lx*Xh%ly*Xh%lz)
    else
       call tnsr2d_el(tmp, Xh%lx, x3, 3, jx, jyt)
       call copy(x, tmp, Xh%lx*Xh%ly*Xh%lz)
       call tnsr2d_el(tmp, Xh%ly, y3, 3, jx, jyt)
       call copy(y, tmp, Xh%lx*Xh%ly*Xh%lz)
    end if

    call space_free(Xh3)
  end subroutine dofmap_xyzquad
 !> Extend faces into interior via gordon hall
 !! gh_type:  1 - vertex only                   
 !!           2 - vertex and edges
 !!           3 - vertex, edges, and faces      
 !! Original in Nek5000/core/navier5.f
    
 subroutine gh_face_extend_3d(x, zg, n, gh_type, e, v)
   integer, intent(in) :: n
   real(kind=rp), intent(inout) ::  x(n, n, n)
   real(kind=rp), intent(in) ::  zg(n)
   real(kind=rp), intent(inout) ::  e(n, n, n)
   real(kind=rp), intent(inout) ::  v(n, n, n)
   integer :: gh_type, ntot, kk, jj, ii, k, j, i
   real(kind=rp) :: si, sj, sk, hi, hj, hk
   
!
!  Build vertex interpolant
!
   ntot = n**3
   call rzero(v, ntot)
   do kk = 1, n, n-1
      do jj = 1, n, n-1
         do ii = 1, n, n-1
            do k = 1, n
               do j = 1, n
                  do i = 1, n
                     si       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
                     sj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
                     sk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
                     v(i,j,k) = v(i,j,k) + si*sj*sk*x(ii,jj,kk)
                  end do
               end do
            end do
         end do
      end do
   end do

   if (gh_type .eq. 1) then
      call copy(x, v, ntot)
      return
   end if
!
!
!  Extend 12 edges
   call rzero(e, ntot)
!
!  x-edges
!
   do kk = 1, n, n-1
      do jj = 1, n, n-1
         do k = 1, n
            do j = 1, n
               do i = 1, n
                  hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
                  hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
                  e(i,j,k) = e(i,j,k) + hj*hk*(x(i,jj,kk)-v(i,jj,kk))
               end do
            end do
         end do
      end do
   end do
!
!  y-edges
!
   do kk = 1, n, n-1
      do ii = 1, n, n-1
         do k = 1, n
            do j = 1, n
               do i = 1, n
                  hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
                  hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
                  e(i,j,k) = e(i,j,k) + hi*hk*(x(ii,j,kk)-v(ii,j,kk))
               end do
            end do
         end do
      end do
   end do
!
!  z-edges
!
   do jj = 1, n, n-1
      do ii = 1, n, n-1
         do k = 1, n
            do j = 1, n
               do i = 1, n
                  hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
                  hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
                  e(i,j,k) = e(i,j,k) + hi*hj*(x(ii,jj,k)-v(ii,jj,k))
               end do
            end do
         end do
      end do
   end do
   
   call add2(e, v, ntot)
   
   if (gh_type .eq. 2) then
      call copy(x, e, ntot)
      return
   end if
!
!  Extend faces
!
   call rzero(v, ntot)
!
!  x-edges
!
   do ii = 1, n, n-1
      do k = 1, n
         do j = 1, n
            do i = 1, n
               hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
               v(i,j,k) = v(i,j,k) + hi*(x(ii,j,k)-e(ii,j,k))
            end do
         end do
      end do
   end do
!
! y-edges
!
   do jj = 1, n, n-1
      do k = 1, n
         do j = 1, n
            do i = 1, n
               hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
               v(i,j,k) = v(i,j,k) + hj*(x(i,jj,k)-e(i,jj,k))
            end do
         end do
      end do
   end do
! 
!  z-edges
! 
   do kk= 1 , n, n-1
      do k = 1, n
         do j = 1, n
            do i = 1, n
               hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
               v(i,j,k) = v(i,j,k) + hk*(x(i,j,kk)-e(i,j,kk))
            end do
         end do
      end do
   end do
   
   call add2(v, e, ntot)
   call copy(x, v ,ntot)

  end subroutine gh_face_extend_3d

  !> Extend 2D faces into interior via gordon hall
  !! gh_type:  1 - vertex only
  !!           2 - vertex and faces
  subroutine gh_face_extend_2d(x, zg, n, gh_type, e, v)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n, n)
    real(kind=rp), intent(in) :: zg(n)
    real(kind=rp), intent(inout) :: e(n, n)
    real(kind=rp), intent(inout) :: v(n, n)
    integer, intent(in) :: gh_type
    integer :: i,j , jj, ii, ntot
    real(kind=rp) :: si, sj, hi, hj 

    !Build vertex interpolant

    ntot=n*n
    call rzero(v,ntot)
    do jj = 1, n, n-1
       do ii = 1, n, n-1
          do j = 1, n
             do i = 1, n
                si     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
                sj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
                v(i,j) = v(i,j) + si*sj*x(ii,jj)
             end do
          end do
       end do
    end do
    if (gh_type .eq. 1) then
       call copy(x,v,ntot)
       return
    end if

    !Extend 4 edges
    call rzero(e,ntot)

    !x-edges

    do jj = 1, n, n-1
       do j = 1, n
          do i = 1, n
             hj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
             e(i,j) = e(i,j) + hj*(x(i,jj)-v(i,jj))
          end do
       end do
    end do

    !y-edges

    do ii = 1, n, n-1
       do j = 1, n
          do i = 1, n
             hi     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
             e(i,j) = e(i,j) + hi*(x(ii,j)-v(ii,j))
          end do
       end do
    end do

    call add3(x, e, v, ntot)

  end subroutine gh_face_extend_2d



  subroutine arc_surface(isid, curve_data, x, y, z, Xh, element, gdim)
    integer, intent(in) :: isid, gdim
    type(space_t), intent(in) :: Xh
    class(element_t) :: element
    real(kind=rp), dimension(5), intent(in) :: curve_data
    real(kind=rp), dimension(Xh%lx, Xh%ly, Xh%lz), intent(inout) :: x, y, z 
    real(kind=rp) :: pt1x, pt1y, pt2x, pt2y, pt12x, pt12y
    real(kind=rp) :: radius, gap, dtheta, r, xys 
    real(kind=rp) :: theta0, xcenn, ycenn, h(Xh%lx, 3, 2)
    real(kind=rp) :: xcrved(Xh%lx), ycrved(Xh%lx), xs, ys
    integer :: isid1, ixt, iyt, izt, ix, isidt

    pt1x  = element%pts(isid)%p%x(1)
    pt1y  = element%pts(isid)%p%x(2)
    if(isid.eq.4) then
       pt2x  = element%pts(1)%p%x(1)
       pt2y  = element%pts(1)%p%x(2)
    else if(isid.eq.8) then
       pt2x  = element%pts(5)%p%x(1)
       pt2y  = element%pts(5)%p%x(2)
    else
       pt2x  = element%pts(isid+1)%p%x(1)
       pt2y  = element%pts(isid+1)%p%x(2)
    end if
!   find slope of perpendicular
    radius=curve_data(1)
    gap=sqrt( (pt1x-pt2x)**2 + (pt1y-pt2y)**2 )
    if (abs(2.0 * radius) .le. gap * 1.00001) then
       call neko_error('Radius to small for arced element surface')
    end if
    xs = pt2y-pt1y
    ys = pt1x-pt2x
!   make length radius
    xys=sqrt(xs**2+ys**2)
!   find center
    dtheta = abs(asin(0.5*gap/radius))
    pt12x  = (pt1x + pt2x)/2.0
    pt12y  = (pt1y + pt2y)/2.0
    xcenn  = pt12x - xs/xys * radius*cos(dtheta)
    ycenn  = pt12y - ys/xys * radius*cos(dtheta)
    theta0 = atan2((pt12y-ycenn),(pt12x-xcenn))
!   compute perturbation of geometry
    isid1 = mod(isid+4-1, 4)+1
    call compute_h(h, Xh%zg, gdim, Xh%lx) 
    do ix=1,Xh%lx
       ixt=ix
       if (isid1.gt.2) ixt=Xh%lx+1-ix
       r=Xh%zg(ix,1)
       if (radius.lt.0.0) r=-r
       xcrved(ixt) = xcenn + abs(radius) * cos(theta0 + r*dtheta) &
                           - ( h(ix,1,1)*pt1x + h(ix,1,2)*pt2x )
       ycrved(ixt) = ycenn + abs(radius) * sin(theta0 + r*dtheta) &
                           - ( h(ix,1,1)*pt1y + h(ix,1,2)*pt2y )
    end do
!   points all set, add perturbation to current mesh.
!   LEGACY WARNING
!   I dont want to dive in this again, Martin Karp 2/3 - 2021
    isidt = isid1
    select case(isidt)
    case (1)
       isid1 = 3
    case (2)
       isid1 = 2
    case (3)
       isid1 = 4
    case (4)
       isid1 = 1
    case (5)
       isid1 = 5
    case (6)
       isid1 = 6
    end select
    izt = (isid-1)/4+1
    iyt = isid1-2
    ixt = isid1
    if (isid1 .le. 2) then
       call addtnsr(x, h(1,1,ixt), xcrved, h(1,3,izt) &
                   ,Xh%lx, Xh%ly, Xh%lz)
       call addtnsr(y, h(1,1,ixt), ycrved, h(1,3,izt) &
                   ,Xh%lx, Xh%ly, Xh%lz)
    else
       call addtnsr(x, xcrved, h(1,2,iyt), h(1,3,izt) &
                   ,Xh%lx, Xh%ly, Xh%lz)
       call addtnsr(y, ycrved, h(1,2,iyt), h(1,3,izt) &
                   ,Xh%lx, Xh%ly, Xh%lz)
    end if
  end subroutine arc_surface

  subroutine compute_h(h, zgml, gdim, lx)
    integer, intent(in) :: lx, gdim
    real(kind=rp), intent(inout) ::  h(lx, 3, 2)
    real(kind=rp), intent(in) :: zgml(lx, 3)
    integer :: ix, iy, iz 

    do ix = 1, lx
       h(ix,1,1) = (1.0_rp - zgml(ix, 1)) * 0.5_rp
       h(ix,1,2) = (1.0_rp + zgml(ix, 1)) * 0.5_rp
    end do
    
    do iy = 1, lx
       h(iy,2,1) = (1.0_rp - zgml(iy, 2)) * 0.5_rp
       h(iy,2,2) = (1.0_rp + zgml(iy, 2)) * 0.5_rp
    end do
    
    if (gdim .eq. 3) then
       do iz = 1, lx
          h(iz,3,1) = (1.0_rp - zgml(iz, 3)) * 0.5_rp
          h(iz,3,2) = (1.0_rp + zgml(iz, 3)) * 0.5_rp
       end do
    else
       call rone(h(1,3,1), lx)
       call rone(h(1,3,2), lx)
    end if
    
  end subroutine compute_h
  
end module dofmap
