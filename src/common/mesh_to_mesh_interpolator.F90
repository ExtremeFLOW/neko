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
!> Implements mesh to mesh interpolation.
!! @note This modules uses functions from `gslib`, namely `findpts_setup`,
!! `findpts`, and `findpts_eval`. A full description of these subroutines can
!! be found at https://github.com/Nek5000/gslib/blob/master/src/findpts.c
module mesh_to_mesh_interpolator
  use num_types, only: rp
  use matrix, only: matrix_t
  use logger, only: neko_log, LOG_SIZE
  use utils, only: neko_error, neko_warning
  use comm
  use mpi_types
  use field, only: field_t
  use coefs, only: coef_t
  use field_list
  use time_based_controller, only : time_based_controller_t
  use field_registry, only : neko_field_registry
  use json_module, only : json_file
  use json_utils, only : json_get
  use device
  use file
  use math, only: rzero, copy
  use csv_file
  use tensor
  use point_interpolator
  use point, only: point_t
  use space, only: space_t
  implicit none
  private
  
  type, public :: mesh_to_mesh_interpolator_t
     !> handle to pass in each findpts call
     integer :: handle
     !> Number of probes
     integer :: n_probes
     !> Number of output fields
     integer :: n_fields = 0
     !> List of owning processes
     integer, allocatable :: proc_owner(:)
     !> List of owning elements
     integer, allocatable :: el_owner(:)
     !> Distance squared between original and interpolated point (in xyz space)
     real(kind=rp), allocatable :: dist2(:)
     !> Error code for each point
     integer, allocatable :: error_code(:)
     !> r,s,t coordinates, findpts format
     type(point_t), allocatable :: rst(:)
     !> x,y,z coordinates, findpts format
     real(kind=rp), allocatable :: xyz(:,:)
     !> interpolated fields
     real(kind=rp), allocatable :: out_fields(:,:)
     !> Time based controller for sampling
     type(time_based_controller_t) :: controller
     !> Fields to be probed
     type(field_list_t), pointer :: sampled_fields
     character(len=20), allocatable  :: which_fields(:)
     !> Interpolator instead of findpts_eval
     type(point_interpolator_t) :: interpolator
     !> Coordinated of the input mesh
     real(kind=rp), allocatable :: x_in(:,:,:,:)       !< Mapping to x-coordinates
     real(kind=rp), allocatable :: y_in(:,:,:,:)       !< Mapping to y-coordinates
     real(kind=rp), allocatable :: z_in(:,:,:,:)       !< Mapping to z-coordinates
     integer :: lx_in = 0
     integer :: ly_in = 0
     integer :: lz_in = 0
     integer :: nelv_in = 0
     integer :: lx_out = 0
     integer :: ly_out = 0
     integer :: lz_out = 0
     integer :: nelv_out = 0

     contains
       !> Initialize.
       procedure, pass(this) :: init => msh_to_msh_int_init
       !> Destructor
       procedure, pass(this) :: free => msh_to_msh_int_free
       !> Set up
       procedure, pass(this) :: setup => msh_to_msh_int_setup
       !> Maps `x,y,z` to `r,s,t` coordinates.
       procedure, pass(this) :: map => msh_to_msh_int_map
       !> Interpolate each probe from its `r,s,t` coordinates.
       procedure, pass(this) :: interpolate => msh_to_msh_int_interpolate

    end type mesh_to_mesh_interpolator_t

contains

  !> Initialize the objects.
  subroutine msh_to_msh_int_init(this, x_in,  y_in,  z_in,  lx_in,  ly_in,  lz_in,  nelv_in,  &
                                       x_out, y_out, z_out, n_out, field_list, Xh)

    class(mesh_to_mesh_interpolator_t), intent(inout) :: this
    integer, intent(in) :: lx_in, ly_in, lz_in, nelv_in, n_out
    real(kind=rp), intent(in) :: x_in(lx_in, ly_in, lz_in, nelv_in)
    real(kind=rp), intent(in) :: y_in(lx_in, ly_in, lz_in, nelv_in)
    real(kind=rp), intent(in) :: z_in(lx_in, ly_in, lz_in, nelv_in)
    real(kind=rp), intent(in) :: x_out(n_out)
    real(kind=rp), intent(in) :: y_out(n_out)
    real(kind=rp), intent(in) :: z_out(n_out)
    type(field_list_t), intent(in), target :: field_list
    type(space_t), intent(in) :: Xh
    integer :: n_probes, n_in, i, n_fields

    ! Calculate parameters
    n_in           = lx_in * ly_in * lz_in * nelv_in
    n_probes       = n_out
    n_fields       = size(field_list%fields)

    ! Assign parameters
    this%n_probes  = n_probes
    this%n_fields  = n_fields
    this%lx_in     = lx_in
    this%ly_in     = ly_in
    this%lz_in     = lz_in
    this%nelv_in   = nelv_in

    ! Assign pointers
    this%sampled_fields => field_list

    !> Allocate memory
    ! Size of xyz as used in gslib
    allocate(this%x_in(lx_in, ly_in, lz_in, nelv_in)) !< x in    
    allocate(this%y_in(lx_in, ly_in, lz_in, nelv_in)) !< y in    
    allocate(this%z_in(lx_in, ly_in, lz_in, nelv_in)) !< z in    
    allocate(this%xyz(3, n_probes))                   !< out coordinates    
    allocate(this%proc_owner(n_probes))
    allocate(this%el_owner(n_probes))
    allocate(this%dist2(n_probes))
    allocate(this%error_code(n_probes))
    allocate(this%rst(n_probes))
    allocate(this%out_fields(n_probes, n_fields))

    !> Copy the mesh data to appropiate place in object
    call copy(this%x_in, x_in, n_in)
    call copy(this%y_in, y_in, n_in)
    call copy(this%z_in, z_in, n_in)

    do i = 1, n_probes
       this%xyz(1,i) = x_out(i)
       this%xyz(2,i) = y_out(i)
       this%xyz(3,i) = z_out(i)
    end do


    !> Initialize the interpolator with a function space
    call this%interpolator%init(xh)

  end subroutine msh_to_msh_int_init

  !> Destructor
  subroutine msh_to_msh_int_free(this)
    class(mesh_to_mesh_interpolator_t), intent(inout) :: this

    if (allocated(this%x_in))        deallocate(this%x_in)
    if (allocated(this%y_in))        deallocate(this%y_in)
    if (allocated(this%z_in))        deallocate(this%z_in)
    if (allocated(this%xyz))        deallocate(this%xyz)
    if (allocated(this%rst))        deallocate(this%rst)
    if (allocated(this%proc_owner)) deallocate(this%proc_owner)
    if (allocated(this%el_owner))   deallocate(this%el_owner)
    if (allocated(this%dist2))       deallocate(this%dist2)
    if (allocated(this%error_code)) deallocate(this%error_code)
    if (allocated(this%out_fields)) deallocate(this%out_fields)
    if (allocated(this%sampled_fields%fields)) deallocate(this%sampled_fields%fields)

    call this%interpolator%free

#ifdef HAVE_GSLIB
    call fgslib_findpts_free(this%handle)
#else
    call neko_error('NEKO needs to be built with GSLIB support')
#endif
    
  end subroutine msh_to_msh_int_free

  !> Setup gslib for mapping process (with fgslib_findpts_setup).
  subroutine msh_to_msh_int_setup(this, coef)
    class(mesh_to_mesh_interpolator_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    real(kind=rp) :: tolerance
    integer :: lx, ly, lz, nelv, max_pts_per_iter

#ifdef HAVE_GSLIB
    
    ! Tolerance for Newton iterations
    tolerance = 5d-13
    lx   = this%lx_in
    ly   = this%ly_in
    lz   = this%lz_in
    nelv = this%nelv_in
    !Number of points to iterate on simultaneosuly
    max_pts_per_iter = 128

    call fgslib_findpts_setup(this%handle, &
         NEKO_COMM, pe_size, &
         coef%msh%gdim, &
         this%x_in, this%y_in, this%z_in, & ! Physical nodal values
         lx, ly, lz, nelv, & ! Mesh dimensions
         2*lx, 2*ly, 2*lz, & ! Mesh size for bounding box computation
         0.01, & ! relative size to expand bounding boxes by
         lx*ly*lz*nelv, lx*ly*lz*nelv, & ! local/global hash mesh sizes
         max_pts_per_iter, tolerance) 
#else
    call neko_error('NEKO needs to be built with GSLIB support')
#endif

  end subroutine msh_to_msh_int_setup


  !> Maps `x,y,z` to `r,s,t` coordinates. The result of the mapping for each
  !! point can be found in the following variables:
  !! - `error_code`: returns `0` if point found, `1` if closest point on a border
  !! (check dist2), `2` if not found
  !! - `dist2`: distance squared (used to compare the points found by each
  !! processor)
  !! @param coef Coefficients associated with the probe's fields.
  subroutine msh_to_msh_int_map(this, coef, tol_x, tol_y, tol_z)
    class(mesh_to_mesh_interpolator_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef  
    real(kind=rp), intent(in) :: tol_x
    real(kind=rp), intent(in) :: tol_y
    real(kind=rp), intent(in) :: tol_z

    type(logical) :: try_again = .false.
    real(kind=rp) :: tol_dist = 5d-6
    real(kind=rp) :: placeholder
    real(kind=rp) :: rst_gslib_raw(3*this%n_probes)
    integer :: i

#ifdef HAVE_GSLIB
    call fgslib_findpts(this%handle, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         rst_gslib_raw, coef%msh%gdim, &
         this%dist2, 1, &
         this%xyz(1,1), coef%msh%gdim, &
         this%xyz(2,1), coef%msh%gdim, &
         this%xyz(3,1), coef%msh%gdim, this%n_probes)

    ! Final check to see if there are any problems
    do i=1,this%n_probes
       if (this%error_code(i) .eq. 1) then
          if (this%dist2(i) .gt. tol_dist) then
             call neko_warning("Point on boundary or outside the mesh!")
          end if
       end if

       if (this%error_code(i) .eq. 2) then
          call neko_warning("Point not within the mesh!")
          write(*,*) 'point is:', this%xyz(:,i)
          try_again = .true.
          
          placeholder = this%xyz(1,i) 
          if (placeholder.gt.0) then
                this%xyz(1,i) = this%xyz(1,i) - tol_x
          else
                this%xyz(1,i) = this%xyz(1,i) + tol_x
          end if
          
          placeholder = this%xyz(2,i) 
          if (placeholder.gt.0) then
                this%xyz(2,i) = this%xyz(2,i) - tol_y
          else
                this%xyz(2,i) = this%xyz(2,i) + tol_y
          end if
          
          placeholder = this%xyz(3,i) 
          if (placeholder.gt.0) then
                this%xyz(3,i) = this%xyz(3,i) - tol_z
          else
                this%xyz(3,i) = this%xyz(3,i) + tol_z
          end if

       end if
    end do

    if (try_again) then
      
      call neko_warning("Trying mapping again substracting user tols")
    
      call fgslib_findpts(this%handle, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         rst_gslib_raw, coef%msh%gdim, &
         this%dist2, 1, &
         this%xyz(1,1), coef%msh%gdim, &
         this%xyz(2,1), coef%msh%gdim, &
         this%xyz(3,1), coef%msh%gdim, this%n_probes)

      ! Final check to see if there are any problems
      do i=1,this%n_probes
             
         ! Reformat the rst array into points
         this%rst(i)%x(1) = rst_gslib_raw(3*(i-1) + 1)
         this%rst(i)%x(2) = rst_gslib_raw(3*(i-1) + 2)
         this%rst(i)%x(3) = rst_gslib_raw(3*(i-1) + 3)

         if (this%error_code(i) .eq. 1) then
            if (this%dist2(i) .gt. tol_dist) then
               call neko_warning("Point on boundary or outside the mesh!")
            end if
         end if

         if (this%error_code(i) .eq. 2) then
            call neko_warning("Point not within the mesh!")
            write(*,*) 'point is:', this%xyz(:,i)
         end if
      end do

      call neko_warning("New mapping done")
    end if

#else
    call neko_error('NEKO needs to be built with GSLIB support')
#endif

  end subroutine msh_to_msh_int_map


  !> Interpolate each probe from its `r,s,t` coordinates.
  !! @param t Current simulation time.
  !! @param tstep Current time step.
  !! @param write_output Flag for writing output data.
  subroutine msh_to_msh_int_interpolate(this)
    class(mesh_to_mesh_interpolator_t), intent(inout) :: this

    !> Supporting variables
    integer :: il 
    integer :: n, i
    integer :: ierr
    !> Will store the local interpolated field for each rank
    real(kind=rp) :: local_out_field(this%n_probes)
    type(point_t) :: my_rst(1)
    real(kind=rp) :: my_interpolated_field(1)

#ifdef HAVE_GSLIB
    n = this%sampled_fields%fields(1)%f%dof%size()

    !! Interpolate the fields
    do il = 1, this%n_fields

       ! Copy the field to the CPU if the data is in the device
       if (NEKO_BCKND_DEVICE .eq. 1) then 
          call device_memcpy(this%sampled_fields%fields(il)%f%x, &
                          this%sampled_fields%fields(il)%f%x_d, &
                          n, DEVICE_TO_HOST)
       end if

       ! Initialize the local interpolated field with 0 since we are
       ! Reducing later with MPI_MAX
       call rzero(local_out_field, this%n_probes)

       do i = 1, this%n_probes

          ! Each rank interpolates their own points
          if (pe_rank .eq. this%proc_owner(i)) then

             my_rst(1) = this%rst(i)%x

             my_interpolated_field = this%interpolator%interpolate(my_rst, &
                     this%sampled_fields%fields(il)%f%x(:,:,:,this%el_owner(i)+1))

             local_out_field(i) = my_interpolated_field(1)
          end if
       end do

       ! Copy the data to the output - be careful with comunicating it if needed
       call copy(this%out_fields(1,il), local_out_field, this%n_probes)

    end do


#else
    call neko_error('NEKO needs to be built with GSLIB support')
#endif

  end subroutine msh_to_msh_int_interpolate

end module mesh_to_mesh_interpolator


