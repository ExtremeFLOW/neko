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
!> Implements probes.
!! @note This modules uses functions from `gslib`, namely `findpts_setup`,
!! `findpts`, and `findpts_eval`. A full description of these subroutines can
!! be found at https://github.com/Nek5000/gslib/blob/master/src/findpts.c
module probes
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
  use math, only: rzero
  use csv_file
  use tensor
  use point_interpolator
  use point, only: point_t
  use space, only: space_t
  use device
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public :: probes_t
     !> Number of probes
     integer :: n_probes
     !> Number of output fields
     integer :: n_fields = 0
     !! ===================== gslib variables ==================================
     !> handle to pass in each findpts call
     integer :: handle
     !> List of owning processes
     integer, allocatable :: proc_owner(:)
     !> List of owning elements
     integer, allocatable :: el_owner(:)
     !> Distance squared between original and interpolated point (in xyz space)
     real(kind=rp), allocatable :: dist2(:)
     !> Error code for each point
     integer, allocatable :: error_code(:)
     !> x,y,z coordinates, findpts format
     real(kind=rp), allocatable :: xyz(:,:)
     !! ========================================================================
     !> r,s,t coordinates
     type(point_t), allocatable :: rst(:)
     !> interpolated fields
     real(kind=rp), allocatable :: out_fields(:,:)
     !> Local element ownership
     integer, allocatable :: local_el_owner(:)
     !> Number of local elements per rank
     integer :: n_local_probes
     !> Local to global mapping to retrieve the fields
     integer, allocatable :: local_to_global(:)
     !> Local rst
     type(point_t), allocatable :: local_rst(:)
     !> Interpolator instead of findpts_eval
     type(point_interpolator_t) :: interpolator
     !> Weights in the r direction on the host
     real(kind=rp), allocatable :: weights_r(:,:)
     !> Weights in the s direction on the host
     real(kind=rp), allocatable :: weights_s(:,:)
     !> Weights in the s direction on the host
     real(kind=rp), allocatable :: weights_t(:,:)
     !> Local element ids on the device
     type(c_ptr) :: local_el_owner_d = C_NULL_PTR
     !> Weights in the r direction on the device
     type(c_ptr) :: weights_r_d = C_NULL_PTR
     !> Weights in the s direction on the device
     type(c_ptr) :: weights_s_d = C_NULL_PTR
     !> Weights in the t direction on the device
     type(c_ptr) :: weights_t_d = C_NULL_PTR
     !> To know if we need to recompute the weights
     logical :: update_weights = .true.
     !> Time based controller for sampling
     type(time_based_controller_t) :: controller
     !> Fields to be probed
     type(field_list_t) :: sampled_fields
     character(len=20), allocatable  :: which_fields(:)


     contains
       !> Initialize probes object.
       procedure, pass(this) :: allocate_fields => probes_allocate_fields
       !> Initialize user defined parameters.
       procedure, pass(this) :: init => probes_init
       !> Destructor
       procedure, pass(this) :: free => probes_free
       !> Print current probe status, with number of probes and coordinates
       procedure, pass(this) :: show => probes_show
       !> Show the status of processor/element owner and error code for each point
       procedure, pass(this) :: debug => probes_debug
       !> Setup the probes for mapping process (with fgslib_findpts_setup).
       procedure, pass(this) :: setup => probes_setup
       !> Maps `x,y,z` to `r,s,t` coordinates.
       procedure, pass(this) :: map => probes_map
       !> Interpolate each probe from its `r,s,t` coordinates.
       procedure, pass(this) :: interpolate => probes_interpolate

    end type probes_t

contains

  !> Initialize user defined variables.
  !! @param t Current simulation time.
  !! @param params case file.
  !! @coef Xh Function space.
  subroutine probes_init(this, t, params, xh)
    class(probes_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    type(json_file), intent(inout) :: params
    type(space_t), intent(in) :: xh

    ! Counter
    integer :: i
    ! Controller parameters 
    real(kind=rp)                  :: T_end
    real(kind=rp)                  :: output_value
    character(len=:), allocatable  :: output_control
    character(len=:), allocatable  :: points_file

    !> Read from case file
    call params%info('case.probes.fields', n_children=this%n_fields)
    call json_get(params, 'case.probes.fields', this%which_fields) 
    call json_get(params, 'case.end_time', T_end)
    call json_get(params, 'case.probes.output_control', &
         output_control)
    call json_get(params, 'case.probes.output_value', &
         output_value)
    call json_get(params, 'case.probes.points_file', points_file)

    !> Fields
    ! List with fields
    allocate(this%sampled_fields%fields(this%n_fields))
    do i = 1, this%n_fields
       this%sampled_fields%fields(i)%f => neko_field_registry%get_field(&
                                          trim(this%which_fields(i)))
    end do

    !> Controllers
    ! Calculate relevant parameters and restart                     
    call this%controller%init(T_end, output_control, &
                                  output_value)
    ! Update nexecutions in restarts
    if (this%controller%nsteps .eq. 0) then
        this%controller%nexecutions = &
               int(t / this%controller%time_interval) + 1
    end if

    !> Read file
    call read_probe_locations(this, points_file)

    !> Initialize the interpolator with a function space
    call this%interpolator%init(xh)

  end subroutine probes_init


  !> Destructor
  subroutine probes_free(this)
    class(probes_t), intent(inout) :: this

    if (allocated(this%xyz))        deallocate(this%xyz)
    if (allocated(this%rst))        deallocate(this%rst)
    if (allocated(this%proc_owner)) deallocate(this%proc_owner)
    if (allocated(this%el_owner))   deallocate(this%el_owner)
    if (allocated(this%dist2))       deallocate(this%dist2)
    if (allocated(this%error_code)) deallocate(this%error_code)
    if (allocated(this%out_fields)) deallocate(this%out_fields)
    if (allocated(this%sampled_fields%fields)) deallocate(this%sampled_fields%fields)
    if (allocated(this%local_el_owner)) deallocate(this%local_el_owner)
    if (allocated(this%local_rst)) deallocate(this%local_rst)
    if (allocated(this%local_to_global)) deallocate(this%local_to_global)

    if (allocated(this%weights_r)) deallocate(this%weights_r)
    if (allocated(this%weights_s)) deallocate(this%weights_s)
    if (allocated(this%weights_t)) deallocate(this%weights_t)

    if (c_associated(this%local_el_owner_d)) call device_free(this%local_el_owner_d)
    if (c_associated(this%weights_r_d)) call device_free(this%weights_r_d)
    if (c_associated(this%weights_s_d)) call device_free(this%weights_s_d)
    if (c_associated(this%weights_t_d)) call device_free(this%weights_t_d)

    call this%interpolator%free

#ifdef HAVE_GSLIB
    call fgslib_findpts_free(this%handle)
#else
    call neko_error('NEKO needs to be built with GSLIB support')
#endif

  end subroutine probes_free

  !> Print current probe status, with number of probes and coordinates
  subroutine probes_show(this)
    class(probes_t), intent(in) :: this
    character(len=LOG_SIZE) :: log_buf ! For logging status
    integer :: i

    !> Probes summary
    call neko_log%section('Probes')
    write(log_buf, '(A,I6)') "Number of probes: ", this%n_probes
    call neko_log%message(log_buf)
    call neko_log%message("xyz-coordinates:")
    do i=1,this%n_probes
       write(log_buf, '("(",F10.6,",",F10.6,",",F10.6,")")') this%xyz(:,i)
       call neko_log%message(log_buf)
    end do
    !> Field summary
    write(log_buf, '(A,I6)') "Number of fields: ", this%n_fields
    call neko_log%message(log_buf)
    do i=1,this%n_fields
       write(log_buf, '(A,I6, A ,A)') "Field: ", i, " ", trim(this%which_fields(i))
       call neko_log%message(log_buf)
    end do
    call neko_log%end_section()
    call neko_log%newline()
    
  end subroutine probes_show

  !> Show the status of processor/element owner and error code for each point
  subroutine probes_debug(this)
    class(probes_t) :: this

    character(len=LOG_SIZE) :: log_buf ! For logging status
    integer :: i

    do i = 1, this%n_probes
       write (log_buf, *) pe_rank, "/", this%proc_owner(i), "/" , this%el_owner(i), "/", this%error_code(i)
       call neko_log%message(log_buf)
       write(log_buf, '(A5,"(",F10.6,",",F10.6,",",F10.6,")")') "rst: ", this%rst(i)%x
       call neko_log%message(log_buf)
    end do
  end subroutine probes_debug

  !> Setup the probes for mapping process (with fgslib_findpts_setup).
  !! @param coef Coefficients associated with the probe's fields.
  subroutine probes_setup(this, coef)
    class(probes_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef

    real(kind=rp) :: tolerance
    integer :: lx, ly, lz, nelv, max_pts_per_iter

#ifdef HAVE_GSLIB

    ! Tolerance for Newton iterations
    tolerance = 5d-13
    lx = coef%xh%lx
    ly = coef%xh%ly
    lz = coef%xh%lz
    nelv = coef%msh%nelv
    !Number of points to iterate on simultaneosuly
    max_pts_per_iter = 128

    call fgslib_findpts_setup(this%handle, &
         NEKO_COMM, pe_size, &
         coef%msh%gdim, &
         coef%dof%x, coef%dof%y, coef%dof%z, & ! Physical nodal values
         lx, ly, lz, nelv, & ! Mesh dimensions
         2*lx, 2*ly, 2*lz, & ! Mesh size for bounding box computation
         0.01, & ! relative size to expand bounding boxes by
         lx*ly*lz*nelv, lx*ly*lz*nelv, & ! local/global hash mesh sizes
         max_pts_per_iter, tolerance)
#else
    call neko_error('NEKO needs to be built with GSLIB support')
#endif

  end subroutine probes_setup


  !> Maps `x,y,z` to `r,s,t` coordinates. The result of the mapping for each
  !! point can be found in the following variables:
  !! - `error_code`: returns `0` if point found, `1` if closest point on a border
  !! (check dist2), `2` if not found
  !! - `dist2`: distance squared (used to compare the points found by each
  !! processor)
  !! @param coef Coefficients associated with the probe's fields.
  subroutine probes_map(this, coef)
    class(probes_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef 

    real(kind=rp) :: rst_gslib_raw(3*this%n_probes)

    real(kind=rp) :: tol_dist = 5d-6
    integer :: i, local_index

#ifdef HAVE_GSLIB

    !
    ! ------------------------ gslib -----------------------------------------
    !
    call fgslib_findpts(this%handle, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         rst_gslib_raw, coef%msh%gdim, &
         this%dist2, 1, &
         this%xyz(1,1), coef%msh%gdim, &
         this%xyz(2,1), coef%msh%gdim, &
         this%xyz(3,1), coef%msh%gdim, this%n_probes)

    ! Number of points owned by a rank
    this%n_local_probes = 0

    do i=1,this%n_probes

       !
       ! Reformat the rst array into point_t
       !
       this%rst(i)%x(1) = rst_gslib_raw(3*(i-1) + 1)
       this%rst(i)%x(2) = rst_gslib_raw(3*(i-1) + 2)
       this%rst(i)%x(3) = rst_gslib_raw(3*(i-1) + 3)

       !
       ! Count the number of local points
       !
       if (pe_rank .eq. this%proc_owner(i)) this%n_local_probes = this%n_local_probes+1

       !
       ! Check validity of points
       !
       if (this%error_code(i) .eq. 1) then
          if (this%dist2(i) .gt. tol_dist) then
             call neko_warning("Point on boundary or outside the mesh!")
          end if
       end if

       if (this%error_code(i) .eq. 2) call neko_warning("Point not within the mesh!")
    end do

    ! -------------------------- END GSLIB ------------------------------------

    !
    ! Allocate local stuff
    !
    if (this%n_local_probes .ne. 0) then

       !
       ! Allocate and associate the element ownership locally
       !
       allocate(this%local_el_owner(this%n_local_probes))
       allocate(this%local_to_global(this%n_local_probes))
       allocate(this%local_rst(this%n_local_probes))

       local_index = 1

       do i = 1, this%n_probes
          if (pe_rank .eq. this%proc_owner(i)) then
             this%local_el_owner(local_index) = this%el_owner(i)
             this%local_to_global(local_index) = i
             this%local_rst(local_index) = this%rst(i)%x
             local_index = local_index + 1
          end if
       end do

       !
       ! This could be another way to generate the local array
       ! this%local_el_owner = pack(this%el_owner, this%proc_owner .eq. pe_rank)
       !

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(this%local_el_owner, this%local_el_owner_d, &
               this%n_local_probes)
          call device_memcpy(this%local_el_owner, this%local_el_owner_d, &
               this%n_local_probes, HOST_TO_DEVICE, sync = .true.)
       end if

    end if

#else
    call neko_error('Neko needs to be built with GSLIB support')
#endif

  end subroutine probes_map

  !> Interpolate each probe from its `r,s,t` coordinates.
  !! @note The final interpolated field is only available on rank 0.
  !! @param t Current simulation time.
  !! @param tstep Current time step.
  !! @param write_output Flag for writing output data.
  subroutine probes_interpolate(this, t, tstep, write_output)
    class(probes_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    logical, intent(inout) :: write_output
    real(kind=rp), allocatable :: tmp(:,:)
    integer :: i, ierr, lx
    integer :: size_outfields, size_weights

#ifdef HAVE_GSLIB

    lx = this%interpolator%Xh%lx
    size_outfields = this%n_probes * this%n_fields
    size_weights = this%n_local_probes * lx

    !> Check controller to determine if we must write
    if (this%controller%check(t, tstep, .false.)) then

       ! Fill the out_field array with 0s as we are reducing later
       ! with MPI_SUM
       call rzero(this%out_fields, size_outfields)

       ! Do not allocate/compute if current has no local probes
       if (this%n_local_probes .ne. 0) then

          !
          ! Only update weights if necessary
          !
          if (this%update_weights) then

             if (.not. allocated(this%weights_r)) allocate(this%weights_r(lx,this%n_local_probes))
             if (.not. allocated(this%weights_s)) allocate(this%weights_s(lx,this%n_local_probes))
             if (.not. allocated(this%weights_t)) allocate(this%weights_t(lx,this%n_local_probes))

             !
             ! Build weights
             !
             do i = 1, this%n_local_probes
                call fd_weights_full(this%rst(i)%x(1), this%interpolator%Xh%zg(:,1), &
                     lx-1, 0, this%weights_r(:,i))
                call fd_weights_full(this%rst(i)%x(2), this%interpolator%Xh%zg(:,2), &
                     lx-1, 0, this%weights_s(:,i))
                call fd_weights_full(this%rst(i)%x(3), this%interpolator%Xh%zg(:,3), &
                     lx-1, 0, this%weights_t(:,i))
             end do

             !
             ! Associate device pointers
             !
             if (NEKO_BCKND_DEVICE .eq. 1) then

                call device_map(this%weights_r, this%weights_r_d, size_weights)
                call device_map(this%weights_s, this%weights_s_d, size_weights)
                call device_map(this%weights_t, this%weights_t_d, size_weights)
                call device_memcpy(this%weights_r, this%weights_r_d, size_weights, HOST_TO_DEVICE, sync = .true.)
                call device_memcpy(this%weights_s, this%weights_s_d, size_weights, HOST_TO_DEVICE, sync = .true.)
                call device_memcpy(this%weights_t, this%weights_t_d, size_weights, HOST_TO_DEVICE, sync = .true.)

             end if

             this%update_weights = .false.

          end if

          allocate(tmp(this%n_local_probes, this%n_fields))

          !
          ! Interpolate
          !
          tmp = this%interpolator%interpolate(this%local_rst, &
               this%local_el_owner, this%sampled_fields, &
               this%weights_r, this%weights_s, this%weights_t)

          !
          ! Reconstruct the global array using global to local mapping
          !
          do i = 1, this%n_local_probes
             this%out_fields(this%local_to_global(i),:) = tmp(i,:)
          end do

          deallocate(tmp)

       end if

       ! Artificial way to gather all values to rank 0
       if (pe_rank .ne. 0) then
          call MPI_Reduce(this%out_fields(1,1), this%out_fields(1,1), size_outfields, &
               MPI_REAL_PRECISION, MPI_SUM, 0, NEKO_COMM, ierr)
       else
          call MPI_Reduce(MPI_IN_PLACE, this%out_fields(1,1), size_outfields, &
               MPI_REAL_PRECISION, MPI_SUM, 0, NEKO_COMM, ierr)
       end if

       !! Turn on flag to write output
       write_output = .true.

       !! Register the execution of the activity
       call this%controller%register_execution()
    end if

#else
    call neko_error('NEKO needs to be built with GSLIB support')
#endif

  end subroutine probes_interpolate

  !> Initialize the physical coordinates from a `csv` input file
  !! @param points_file A cvs file containing probes.
  subroutine read_probe_locations(this, points_file)
    class(probes_t), intent(inout) :: this
    character(len=:), allocatable  :: points_file
    type(matrix_t) :: mat_in
    !> Supporting variables
    type(file_t) :: file_in
    integer :: ierr, file_unit, n_lines

    file_in = file_t(trim(points_file))

    select type(ft => file_in%file_type)
      type is (csv_file_t)
          call read_from_csv(this, ft, mat_in)
      class default
          call neko_error("Invalid data. Expected csv_file_t.")
    end select

    !> After reading the number of probes is know, as well as number of fields
    call this%allocate_fields()

    ! Transpose mat%x and put it in this%xyz
    call trsp(this%xyz, 3, mat_in%x, this%n_probes)

    ! Broadcast the xyz coordinates to all processes
    call MPI_Bcast(this%xyz, 3*this%n_probes, MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr)

    !> Close the file
    call mat_in%free
    call file_free(file_in)

  end subroutine read_probe_locations

  !> Read and initialize the number of probes from a `csv` input file
  subroutine read_from_csv(this, f, mat_in)
    class(probes_t), intent(inout) :: this
    type(csv_file_t), intent(inout) :: f
    type(matrix_t), intent(inout) :: mat_in
    integer :: ierr, file_unit, n_lines

    if (pe_rank .eq. 0) n_lines = f%count_lines()
    call MPI_Bcast(n_lines, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    ! Update the number of probes
    this%n_probes = n_lines

    ! Initialize the temporal array
    call mat_in%init(this%n_probes,3)

    ! Read the data
    call f%read(mat_in)

  end subroutine read_from_csv

  !> Allocate arrays according to previously initialized # of probes and fields.
  subroutine probes_allocate_fields(this)
    class(probes_t), intent(inout) :: this
    integer :: n_probes
    integer :: n_fields

    n_probes = this%n_probes
    n_fields = this%n_fields

    ! Size of xyz as used in gslib
    allocate(this%xyz(3, n_probes))    
    allocate(this%proc_owner(n_probes))
    allocate(this%el_owner(n_probes))
    allocate(this%dist2(n_probes))
    allocate(this%error_code(n_probes))
    allocate(this%rst(n_probes))
    allocate(this%out_fields(n_probes, n_fields))

  end subroutine probes_allocate_fields

end module probes
