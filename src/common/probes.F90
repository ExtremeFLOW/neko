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
  use csv_file
  use tensor
  implicit none
  private
  
  type, public :: probes_t
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
     real(kind=rp), allocatable :: rst(:)
     !> x,y,z coordinates, findpts format
     real(kind=rp), allocatable :: xyz(:,:)
     !> interpolated fields
     real(kind=rp), allocatable :: out_fields(:,:)
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
  subroutine probes_init(this, t, params)
    class(probes_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    type(json_file), intent(inout) :: params
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

    write (*, *) pe_rank, "/", this%proc_owner, "/" , this%el_owner, "/", this%error_code

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

    real(kind=rp) :: tol_dist = 5d-6
    integer :: i

#ifdef HAVE_GSLIB
    call fgslib_findpts(this%handle, &
         this%error_code, 1, &
         this%proc_owner, 1, &
         this%el_owner, 1, &
         this%rst, coef%msh%gdim, &
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

       if (this%error_code(i) .eq. 2) call neko_warning("Point not within the mesh!")
    end do

#else
    call neko_error('NEKO needs to be built with GSLIB support')
#endif

  end subroutine probes_map


  !> Interpolate each probe from its `r,s,t` coordinates.
  !! @param t Current simulation time.
  !! @param tstep Current time step.
  !! @param write_output Flag for writing output data.
  subroutine probes_interpolate(this, t, tstep, write_output)
    class(probes_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep    
    logical, intent(inout) :: write_output

    !> Supporting variables
    integer :: il 
    integer :: n

#ifdef HAVE_GSLIB
    n = this%sampled_fields%fields(1)%f%dof%size()

    !> Check controller to determine if we must write
    if (this%controller%check(t, tstep, .false.)) then

       !! Interpolate the fields
       do il = 1, this%n_fields

          ! Copy the field to the CPU if the data is in the device
          if (NEKO_BCKND_DEVICE .eq. 1) then 
             call device_memcpy(this%sampled_fields%fields(il)%f%x, &
                             this%sampled_fields%fields(il)%f%x_d, &
                             n, DEVICE_TO_HOST)
          end if

          call fgslib_findpts_eval(this%handle, this%out_fields(il,1), this%n_fields, &
                                   this%error_code, 1, &
                                   this%proc_owner, 1, &
                                   this%el_owner, 1, &
                                   this%rst, this%sampled_fields%fields(il)%f%msh%gdim, &
                                   this%n_probes, this%sampled_fields%fields(il)%f%x)
       end do

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
    ! Size of rst as used in gslib
    allocate(this%rst(3*n_probes))     
    allocate(this%out_fields(n_fields, n_probes))

  end subroutine probes_allocate_fields


end module probes
