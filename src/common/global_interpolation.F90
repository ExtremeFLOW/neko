! Copyright (c) 2020-2023, The Neko Authors
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
!> Implements global_interpolation given a dofmap.
!! @note This modules uses functions from `gslib`, namely `findpts_setup`,
!! `findpts`, and `findpts_eval`. A full description of these subroutines can
!! be found at https://github.com/Nek5000/gslib/blob/master/src/findpts.c
!! 
module global_interpolation
  use num_types, only: rp
  use space, only: space_t
  use dofmap, only: dofmap_t
  use mesh, only: mesh_t
  use logger, only: neko_log, LOG_SIZE
  use utils, only: neko_error, neko_warning
  use comm
  use math, only: copy
  use mpi_types
  use, intrinsic :: iso_c_binding
  implicit none
  private
  !> Implements global interpolation for arbitrary points in the domain. 
  type, public :: global_interpolation_t
     !> Dofmap from which we interpolate the points
     type(dofmap_t), pointer :: dof
     !> Mesh on which we interpolate
     type(mesh_t), pointer :: mesh
     !> Space 
     type(space_t), pointer :: Xh
     !! Gslib handle
     !! @note: Remove when we remove gslib
     integer :: gs_handle
     !> Components related to the points we want to evalute
     !> Number of points we want to evaluate
     integer :: n_points
     !> x,y,z coordinates, findpts format
     !! @note: When replacing gs we can change format
     real(kind=rp), allocatable :: xyz(:,:)
     !> List of owning processes
     integer, allocatable :: proc_owner(:)
     !> List of owning elements
     integer, allocatable :: el_owner(:)
     !> r,s,t coordinates findpts format
     !! @note: When replacing gs we can change format
     real(kind=rp), allocatable :: rst(:,:)
     !> Distance squared between original and interpolated point
     !! (in xyz space) (according to gslib)
     real(kind=rp), allocatable :: dist2(:)
     !> Error code for each point, needed for gslib
     integer, allocatable :: error_code(:)
     !> Tolerance for distance squared between original and interpolated point
     real(kind=rp) :: tol = 5e-13
   contains
     !> Initialize the global interpolation object on a dofmap.
     procedure, pass(this) :: init => global_interpolation_init
     !> Destructor
     procedure, pass(this) :: free => global_interpolation_free
     !> Destructor for arrays related to evaluation points
     procedure, pass(this) :: free_points => global_interpolation_free_points
     !> Finds the process owner, global element number,
     !! and local rst coordinates for each point.
     !! Sets up correct values to be able to evalute the points
     procedure, pass(this) :: find_points_coords => global_interpolation_find_coords
     procedure, pass(this) :: find_points_xyz => global_interpolation_find_xyz
     generic :: find_points => find_points_xyz, find_points_coords
     !> Evaluate the value of the field in each point.
     procedure, pass(this) :: evaluate => global_interpolation_evaluate
  end type global_interpolation_t


contains
  !> Initialize user defined variables.
  !! @param dof Dofmap on which the interpolation is to be carried out.
  !! @param tol Tolerance for Newton iterations.
  subroutine global_interpolation_init(this, dof, tol)
    class(global_interpolation_t), intent(inout) :: this
    type(dofmap_t), target :: dof
    real(kind=rp), optional :: tol
    integer :: lx, ly, lz, nelv, max_pts_per_iter

    this%dof => dof
    this%Xh => dof%Xh
    this%mesh => dof%msh
    if(present(tol)) this%tol = tol

#ifdef HAVE_GSLIB
    
    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz
    nelv = this%mesh%nelv
    !Number of points to iterate on simultaneosuly
    max_pts_per_iter = 128
    
    call fgslib_findpts_setup(this%gs_handle, &
         NEKO_COMM, pe_size, &
         this%mesh%gdim, &
         dof%x, dof%y, dof%z, & ! Physical nodal values
         lx, ly, lz, nelv, & ! Mesh dimensions
         2*lx, 2*ly, 2*lz, & ! Mesh size for bounding box computation
         0.01, & ! relative size to expand bounding boxes by
         lx*ly*lz*nelv, lx*ly*lz*nelv, & ! local/global hash mesh sizes
         max_pts_per_iter, this%tol)
#else
    call neko_error('Neko needs to be built with GSLIB support')
#endif

  end subroutine global_interpolation_init


  !> Destructor
  subroutine global_interpolation_free(this)
    class(global_interpolation_t), intent(inout) :: this

    nullify(this%mesh)
    nullify(this%dof)
    nullify(this%Xh)

    call this%free_points()

#ifdef HAVE_GSLIB
    call fgslib_findpts_free(this%gs_handle)
#else
    call neko_error('Neko needs to be built with GSLIB support')
#endif

  end subroutine global_interpolation_free

  !> Destructor for point arrays
  subroutine global_interpolation_free_points(this)
    class(global_interpolation_t), intent(inout) :: this

    this%n_points = 0

    if (allocated(this%xyz))        deallocate(this%xyz)
    if (allocated(this%rst))        deallocate(this%rst)
    if (allocated(this%proc_owner)) deallocate(this%proc_owner)
    if (allocated(this%el_owner))   deallocate(this%el_owner)
    if (allocated(this%dist2))       deallocate(this%dist2)
    if (allocated(this%error_code)) deallocate(this%error_code)

  end subroutine global_interpolation_free_points
  
  !> Common routine for finding the points. 
  subroutine global_interpolation_find_common(this)
    class(global_interpolation_t), intent(inout) :: this
    !!Perhaps this should be kind dp
    real(kind=rp) :: xdiff, ydiff, zdiff
    character(len=8000) :: log_buf
    real(kind=rp), allocatable :: x_check(:)
    real(kind=rp), allocatable :: y_check(:)
    real(kind=rp), allocatable :: z_check(:)
    logical :: isdiff
    integer :: i

#ifdef HAVE_GSLIB

    ! gslib find points, which element they belong, to process etc.
    call fgslib_findpts(this%gs_handle, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         this%rst, this%mesh%gdim, &
         this%dist2, 1, &
         this%xyz(1,1), this%mesh%gdim, &
         this%xyz(2,1), this%mesh%gdim, &
         this%xyz(3,1), this%mesh%gdim, this%n_points)

    do i=1,this%n_points

       !
       ! Check validity of points
       !
       if (this%error_code(i) .eq. 1) then
          if (this%dist2(i) .gt. this%tol) then
             write(*,*) 'Point with coords: ',&
                this%xyz(1,i),&
                this%xyz(2,i),&
                this%xyz(3,i),&
                'Did not converge to tol. Absolute differences squared: ',&
                this%dist2(i)
          end if
       end if

       if (this%error_code(i) .eq. 2) &
             write(*,*) 'Point with coords: ',&
                this%xyz(1,i), this%xyz(2,i), this%xyz(3,i),&
                'Outside the mesh!',&
                ' Interpolation on these points will return 0.0. dist2: ', &
                this%dist2(i),&
                'el_owner, rst coords: ',&
                this%el_owner(i), this%rst(1,i), this%rst(2,i), this%rst(3,i)

    end do

    allocate(x_check(this%n_points))
    allocate(y_check(this%n_points))
    allocate(z_check(this%n_points))
    
    call fgslib_findpts_eval(this%gs_handle, x_check, &
                             1, this%error_code, 1, &
                             this%proc_owner, 1, this%el_owner, 1, &  
                             this%rst, this%mesh%gdim, &
                             this%n_points, this%dof%x)
 
    call fgslib_findpts_eval(this%gs_handle, y_check, &
                             1, this%error_code, 1, &
                             this%proc_owner, 1, this%el_owner, 1, &  
                             this%rst, this%mesh%gdim, &
                             this%n_points, this%dof%y)
 
    call fgslib_findpts_eval(this%gs_handle, z_check, &
                             1, this%error_code, 1, &
                             this%proc_owner, 1, this%el_owner, 1, &  
                             this%rst, this%mesh%gdim, &
                             this%n_points, this%dof%z)


    do i=1,this%n_points

       !
       ! Check validity of points
       !
       isdiff = .false.
       xdiff = (x_check(i)-this%xyz(1,i))**2
       if ( xdiff .gt. this%tol) isdiff = .true.
       ydiff = (y_check(i)-this%xyz(2,i))**2
       if ( ydiff .gt. this%tol) isdiff = .true.
       zdiff = (z_check(i)-this%xyz(3,i))**2
       if ( zdiff .gt. this%tol) isdiff = .true.

       if (isdiff) then
          write(*,*) 'Points with coords: ',&
                this%xyz(1,i),this%xyz(2,i),this%xyz(3,i), &
                'Differ from interpolated coords: ',&
                x_check(i), y_check(i), z_check(i),&
                'Distance squared: ',&
                xdiff, ydiff, zdiff
       end if
    end do
    
    deallocate(x_check)
    deallocate(y_check)
    deallocate(z_check)


#else
    call neko_error('Neko needs to be built with GSLIB support')
#endif
  end subroutine global_interpolation_find_common
  
  !> Finds the corresponding r,s,t coordinates 
  !! in the correct global element as well as which process that owns the point.
  !! After this the values at these points can be evaluated.
  !! If the locations of the points change this must be called again.
  !! - `error_code`: returns `0` if point found, `1` if closest point on a border
  !! (check dist2), `2` if not found
  !! - `dist2`: distance squared (used to compare the points found by each
  !! processor)
  !! @param x The x-coordinates of the points.
  !! @param y The y-coordinates of the points.
  !! @param z The z-coordinates of the points.
  !! @param n_points The number of points.
  subroutine global_interpolation_find_coords(this, x, y, z, n_points)
    class(global_interpolation_t), intent(inout) :: this
    integer :: n_points
    !!Perhaps this should be kind dp
    !!this is to get around that x,y,z often is 4 dimensional...
    !!Should maybe add interface for 1d aswell
    real(kind=rp) :: x(n_points,1,1,1)
    real(kind=rp) :: y(n_points,1,1,1)
    real(kind=rp) :: z(n_points,1,1,1)
    integer :: i

    call this%free_points()

    this%n_points = n_points
    
    allocate(this%xyz(3,n_points))
    allocate(this%rst(3,n_points))
    allocate(this%proc_owner(n_points))
    allocate(this%el_owner(n_points))
    allocate(this%dist2(n_points))
    allocate(this%error_code(n_points))

    do i = 1, n_points
       this%xyz(1,i) = x(i,1,1,1)
       this%xyz(2,i) = y(i,1,1,1)
       this%xyz(3,i) = z(i,1,1,1)
    end do

    call global_interpolation_find_common(this)

  end subroutine global_interpolation_find_coords

  !> Finds the corresponding r,s,t coordinates 
  !! in the correct global element as well as which process that owns the point.
  !! After this the values at these points can be evaluated.
  !! If the locations of the points change this must be called again.
  !! - `error_code`: returns `0` if point found, `1` if closest point on a border
  !! (check dist2), `2` if not found
  !! - `dist2`: distance squared (used to compare the points found by each
  !! processor)
  !! @param xyz The coordinates of the points.
  !! @param n_points The number of points.
  subroutine global_interpolation_find_xyz(this, xyz, n_points)
    class(global_interpolation_t), intent(inout) :: this
    integer, intent(in) :: n_points
    !!Perhaps this should be kind dp
    real(kind=rp), intent(inout) :: xyz(3,n_points)


    call this%free_points()

    this%n_points = n_points
    
    allocate(this%xyz(3,n_points))
    allocate(this%rst(3,n_points))
    allocate(this%proc_owner(n_points))
    allocate(this%el_owner(n_points))
    allocate(this%dist2(n_points))
    allocate(this%error_code(n_points))
    !> make deep copy incase xyz goes out of scope or deallocated
    call copy(this%xyz,xyz,3*n_points)

    call global_interpolation_find_common(this)

  end subroutine global_interpolation_find_xyz


  !> Evalute the interpolated value in the points given a field on the dofmap
  !! @param interp_values Array of values in the given points.
  !! @param field Array of values used for interpolation.
  subroutine global_interpolation_evaluate(this, interp_values, field)
    class(global_interpolation_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: interp_values(this%n_points)
    real(kind=rp), intent(in) :: field(this%dof%size())

#ifdef HAVE_GSLIB
       
    call fgslib_findpts_eval(this%gs_handle, interp_values, &
                             1, this%error_code, 1, &
                             this%proc_owner, 1, this%el_owner, 1, &  
                             this%rst, this%mesh%gdim, &
                             this%n_points, field)

#else
    call neko_error('Neko needs to be built with GSLIB support')
#endif

  end subroutine global_interpolation_evaluate

end module global_interpolation
