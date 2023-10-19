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
  use field, only: field_t
  use field_list, only: field_list_t
  use simulation_component
  use field_registry, only : neko_field_registry
  use dofmap, only: dofmap_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use global_interpolation, only: global_interpolation_t
  use math, only: rzero, copy
  use tensor, only: trsp
  use comm
  use mpi_types
  use device
  use file
  use csv_file
  use case
  use device
  use, intrinsic :: iso_c_binding
  implicit none
  private
  
  type, public, extends(simulation_component_t) :: probes_t
     !> Number of output fields
     integer :: n_fields = 0
     type(global_interpolation_t) :: global_interp
     !> Global number of probes (needed for i/o)
     integer :: n_global_probes
     !> global offset for writing
     integer :: n_probes_offset
     !> x,y,z coordinates
     real(kind=rp), allocatable :: xyz(:,:)
     !> Interpolated values
     real(kind=rp), allocatable :: out_values(:,:)
     type(c_ptr), allocatable :: out_values_d(:)
     real(kind=rp), allocatable :: out_vals_trsp(:,:)
     !> Number of local elements per rank
     integer :: n_local_probes
     !> Fields to be probed
     type(field_list_t) :: sampled_fields
     character(len=20), allocatable  :: which_fields(:)
     !> Allocated on rank 0
     integer, allocatable :: n_local_probes_tot_offset(:)
     integer, allocatable :: n_local_probes_tot(:)
     !>  For output on rank 0
     logical :: seq_io
     real(kind=rp), allocatable :: global_output_values(:,:)
     !> Output variables
     type(file_t) :: fout
     type(matrix_t) :: mat_out     
   contains
     !> Initialize from json
     procedure, pass(this) :: init => probes_init_from_json
     ! Actual constructor
     procedure, pass(this) :: init_from_attributes => &
          probes_init_from_attributes
     !> Destructor
     procedure, pass(this) :: free => probes_free
     !> Setup offset for I/O when using sequential write/read from rank 0
     procedure, pass(this) :: setup_offset => probes_setup_offset
     !> Interpolate each probe from its `r,s,t` coordinates.
     procedure, pass(this) :: compute_ => probes_evaluate_and_write
     
  end type probes_t
  
contains
  
  !> Constructor from json.
  subroutine probes_init_from_json(this, json, case)
    class(probes_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case 
    real(kind=rp), allocatable :: xyz(:,:)
    character(len=:), allocatable  :: output_file
    character(len=:), allocatable  :: points_file
    integer :: i, n_local_probes, n_global_probes
    call this%free()

    call this%init_base(json, case)

    !> Read from case file
    call json%info('fields', n_children=this%n_fields)
    call json_get(json, 'fields', this%which_fields) 
    !> Should be extended to not only csv
    !! but also be possible to define in userfile for example
    call json_get(json, 'points_file', points_file)
    call json_get(json, 'output_file', output_file) 

    allocate(this%sampled_fields%fields(this%n_fields))
    do i = 1, this%n_fields
       this%sampled_fields%fields(i)%f => neko_field_registry%get_field(&
                                          trim(this%which_fields(i)))
    end do
    !> This is distributed as to make it similar to parallel file 
    !! formats latera
    !! Reads all into rank 0
    call read_probe_locations(this, this%xyz, this%n_local_probes, &
         this%n_global_probes, points_file)
    call probes_show(this)
    call this%init_from_attributes(case%fluid%dm_Xh, output_file)
    if(allocated(xyz)) deallocate(xyz)

  end subroutine probes_init_from_json


  !> Initialize without json things
  !! @param dof Dofmap to probe
  !! @output_file Name of output file, current must be CSV
  subroutine probes_init_from_attributes(this, dof, output_file)
    class(probes_t), intent(inout) :: this
    type(dofmap_t), intent(in) :: dof
    character(len=:), allocatable, intent(inout)  :: output_file
    real(kind=rp), allocatable :: global_output_coords(:,:)
    integer :: i, ierr
    type(matrix_t) :: mat_coords



    !> Init interpolator
    call this%global_interp%init(dof)

    !> find probes and redistribute them
    call this%global_interp%find_points_and_redist(this%xyz, this%n_local_probes)

    !> Allocate output array
    allocate(this%out_values(this%n_local_probes,this%n_fields))
    allocate(this%out_values_d(this%n_fields))
    allocate(this%out_vals_trsp(this%n_fields,this%n_local_probes))
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%n_fields
          this%out_values_d(i) = c_null_ptr
          call device_map(this%out_values(:,i), this%out_values_d(i),&
                          this%n_local_probes)
       end do
    end if



    !> Initialize the output file
    this%fout = file_t(trim(output_file))

    select type(ft => this%fout%file_type)
      type is (csv_file_t)
          !> Necessary for not-parallel csv format...
          !! offsets and n points per pe
          !! Needed at root for sequential csv i/o
          allocate(this%n_local_probes_tot(pe_size))
          allocate(this%n_local_probes_tot_offset(pe_size))
          call this%setup_offset()
          if (pe_rank .eq. 0) then
             allocate(global_output_coords(3,&
                      this%n_global_probes))
             call this%mat_out%init(this%n_global_probes, this%n_fields)
             allocate(this%global_output_values(this%n_fields,&
                      this%n_global_probes))
             call mat_coords%init(this%n_global_probes,3)
          end if
          call MPI_Gatherv(this%xyz, 3*this%n_local_probes,&
                           MPI_DOUBLE_PRECISION, global_output_coords,&
                           3*this%n_local_probes_tot,&
                           3*this%n_local_probes_tot_offset,&
                           MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr) 
          if (pe_rank .eq. 0) then
             call trsp(mat_coords%x, this%n_global_probes,&
                       global_output_coords, 3)
             !! Write the data to the file
             call this%fout%write(mat_coords)
          end if
    class default
        call neko_error("Invalid data. Expected csv_file_t.")
    end select

  end subroutine probes_init_from_attributes

  !> Destructor
  subroutine probes_free(this)
    class(probes_t), intent(inout) :: this

    if (allocated(this%xyz)) then
       deallocate(this%xyz)
    end if

    if (allocated(this%out_values)) then
       deallocate(this%out_values)
    end if
    
    if (allocated(this%out_vals_trsp)) then
       deallocate(this%out_vals_trsp)
    end if
    
    if (allocated(this%sampled_fields%fields)) then
       deallocate(this%sampled_fields%fields)
    end if
    
    if (allocated(this%n_local_probes_tot)) then
       deallocate(this%n_local_probes_tot)
    end if

    if (allocated(this%n_local_probes_tot_offset)) then
       deallocate(this%n_local_probes_tot_offset)
    end if
    
    if (allocated(this%global_output_values)) then
       deallocate(this%global_output_values)
    end if

    call this%global_interp%free()

  end subroutine probes_free

  !> Print current probe status, with number of probes and coordinates
  subroutine probes_show(this)
    class(probes_t), intent(in) :: this
    character(len=LOG_SIZE) :: log_buf ! For logging status
    integer :: i

    ! Probes summary
    call neko_log%section('Probes')
    write(log_buf, '(A,I6)') "Number of probes: ", this%n_global_probes
    call neko_log%message(log_buf)
    call neko_log%message("xyz-coordinates:")
    do i=1,this%n_local_probes
       write(log_buf, '("(",F10.6,",",F10.6,",",F10.6,")")') this%xyz(:,i)
       call neko_log%message(log_buf)
    end do
    ! Field summary
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

    do i = 1, this%n_local_probes
       write (log_buf, *) pe_rank, "/", this%global_interp%proc_owner(i), "/" ,&
        this%global_interp%el_owner(i), "/",this%global_interp%error_code(i)
       call neko_log%message(log_buf)
       write(log_buf, '(A5,"(",F10.6,",",F10.6,",",F10.6,")")') "rst: ", this%global_interp%rst(:,i)
       call neko_log%message(log_buf)
    end do
  end subroutine probes_debug

  !> Setup offset for rank 0
  subroutine probes_setup_offset(this)
    class(probes_t) :: this
    integer :: ierr
    this%n_local_probes_tot = 0
    this%n_local_probes_tot_offset = 0
    this%n_probes_offset = 0
    call MPI_Gather(this%n_local_probes, 1, MPI_INTEGER,&
                    this%n_local_probes_tot, 1, MPI_INTEGER,&
                    0, NEKO_COMM, ierr)

    call MPI_Exscan(this%n_local_probes, this%n_probes_offset, 1, &
         MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
    call MPI_Gather(this%n_probes_offset, 1, MPI_INTEGER,&
                    this%n_local_probes_tot_offset, 1, MPI_INTEGER,&
                    0, NEKO_COMM, ierr)



  end subroutine probes_setup_offset

  !> Interpolate each probe from its `r,s,t` coordinates.
  !! @note The final interpolated field is only available on rank 0.
  !! @param t Current simulation time.
  !! @param tstep Current time step.
  subroutine probes_evaluate_and_write(this, t, tstep)
    class(probes_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), allocatable :: tmp(:,:)
    integer :: i, ierr, lx
    integer :: size_outfields

    !> Check controller to determine if we must write

       
    do i = 1,this%n_fields
       call this%global_interp%evaluate(this%out_values(:,i), &
                                        this%sampled_fields%fields(i)%f%x)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%n_fields
          call device_memcpy(this%out_values(:,i),this%out_values_d(i),&
               this%n_local_probes, DEVICE_TO_HOST, sync=.true.)
       end do
    end if

    if (this%output_controller%check(t, tstep)) then
       ! Gather all values to rank 0
       ! If io is only done at root
       if (this%seq_io) then
          call trsp(this%out_vals_trsp, this%n_fields, &
                    this%out_values,this%n_local_probes)
          call MPI_Gatherv(this%out_vals_trsp, this%n_fields*this%n_local_probes,&
                           MPI_DOUBLE_PRECISION, this%global_output_values,&
                           this%n_fields*this%n_local_probes_tot,&
                           this%n_fields*this%n_local_probes_tot_offset,&
                           MPI_DOUBLE_PRECISION, 0, NEKO_COMM, ierr) 
          if (pe_rank .eq. 0) then
             call trsp(this%mat_out%x, this%n_global_probes, &
                       this%global_output_values, this%n_fields)
             call this%fout%write(this%mat_out, t)
          end if
       else 
          call neko_error('probes sim comp, parallel io need implementation')
       end if

       !! Register the execution of the activity
       call this%output_controller%register_execution()
    end if

  end subroutine probes_evaluate_and_write

  !> Initialize the physical coordinates from a `csv` input file
  !! @param points_file A csv file containing probes.
  subroutine read_probe_locations(this, xyz, n_local_probes, n_global_probes, points_file)
    class(probes_t), intent(inout) :: this
    character(len=:), allocatable  :: points_file
    real(kind=rp), allocatable :: xyz(:,:)
    integer, intent(inout) :: n_local_probes, n_global_probes

    !> Supporting variables
    type(file_t) :: file_in
    integer :: ierr, file_unit, n_lines

    file_in = file_t(trim(points_file))
    !> Reads on rank 0 and distributes the probes across the different ranks
    select type(ft => file_in%file_type)
      type is (csv_file_t)
          call read_xyz_from_csv(this, xyz, n_local_probes, n_global_probes, ft)
          this%seq_io = .true.
      class default
          call neko_error("Invalid data. Expected csv_file_t.")
    end select

    !> Close the file
    call file_free(file_in)

  end subroutine read_probe_locations 

  !> Read and initialize the number of probes from a `csv` input file
  !! @param xyz xyz coordinates of the probes
  !! @param n_local_probes The number of probes local to this process
  !! @param n_global_probes The number of total probes on all processes
  !! @param f The csv file we read from
  subroutine read_xyz_from_csv(this, xyz, n_local_probes, n_global_probes, f)
    class(probes_t), intent(inout) :: this
    type(csv_file_t), intent(inout) :: f
    real(kind=rp), allocatable :: xyz(:,:)
    integer, intent(inout) :: n_local_probes, n_global_probes
    type(matrix_t) :: mat_in, mat_in2
    integer :: ierr, file_unit, n_lines

    if (pe_rank .eq. 0) n_lines = f%count_lines()
    call MPI_Bcast(n_lines, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)

    ! Update the number of probes
    n_global_probes = n_lines
    this%n_global_probes = n_lines

    ! Initialize the temporal array
    if (pe_rank .eq. 0) then
       this%n_local_probes = this%n_global_probes
       n_local_probes = n_global_probes
       allocate(xyz(3,this%n_local_probes))
       call mat_in%init(this%n_global_probes,3)
       call mat_in2%init(3,this%n_global_probes)
       call f%read(mat_in)
       call trsp(xyz, 3, mat_in%x, this%n_global_probes)
    else 
       n_local_probes = 0
       this%n_local_probes = 0
       allocate(xyz(3,this%n_local_probes))
    end if

  end subroutine read_xyz_from_csv
end module probes
