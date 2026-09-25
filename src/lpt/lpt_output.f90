! Copyright (c) 2026, The Neko Authors
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
!> Output support for Lagrangian particle tracking.
module lpt_output
  use num_types, only : rp
  use utils, only : filename_split, filename_suffix, neko_error
  use file, only : file_t
  use hdf5_file, only : hdf5_file_t
  use matrix, only : matrix_t
  use tensor, only : trsp
  use comm, only : pe_rank, pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Barrier, MPI_Gather, MPI_Gatherv, MPI_INTEGER
  implicit none
  private

  type, public :: lpt_output_t
     type(file_t) :: output_file
     character(len=1024) :: output_path = ""
     logical :: hdf5_output = .false.
     logical :: inertia = .false.
     integer :: snapshots_per_file = 0
     integer :: snapshots_in_file = 0
     integer :: output_file_index = 1
     integer :: hdf5_output_count = 0
     integer :: n_data = 0
   contains
     procedure, pass(this) :: init => lpt_output_init
     procedure, pass(this) :: free => lpt_output_free
     procedure, pass(this) :: write => lpt_output_write
     procedure, private, pass(this) :: init_file => lpt_output_init_file
     procedure, private, pass(this) :: current_path => lpt_output_current_path
     procedure, private, pass(this) :: init_hdf5 => lpt_output_init_hdf5
     procedure, private, pass(this) :: write_csv => lpt_output_write_csv
     procedure, private, pass(this) :: write_hdf5 => lpt_output_write_hdf5
  end type lpt_output_t

contains

  !> Initialise the LPT output writer.
  !! @param output_path Base output path, including csv/hdf5 suffix.
  !! @param inertia Whether particle records include diameter and density.
  !! @param snapshots_per_file Number of snapshots per file, or 0 for all.
  subroutine lpt_output_init(this, output_path, inertia, snapshots_per_file)
    class(lpt_output_t), intent(inout) :: this
    character(len=*), intent(in) :: output_path
    logical, intent(in) :: inertia
    integer, intent(in) :: snapshots_per_file

    call this%free()
    this%output_path = trim(output_path)
    this%inertia = inertia
    this%snapshots_per_file = snapshots_per_file

    if (inertia) then
       this%n_data = 11
    else
       this%n_data = 9
    end if

    call this%init_file()
  end subroutine lpt_output_init

  !> Initialise the currently active LPT output file.
  subroutine lpt_output_init_file(this)
    class(lpt_output_t), intent(inout) :: this
    character(len=1024) :: output_path
    character(len=80) :: output_suffix

    call this%output_file%free()
    call this%current_path(output_path)

    call filename_suffix(output_path, output_suffix)
    select case (trim(output_suffix))
    case ("csv")
       if (this%inertia) then
          call this%output_file%init(output_path, &
               header = "tstep,time,particle_id,x,y,z,u,v,w,d,rho", &
               overwrite = .true.)
       else
          call this%output_file%init(output_path, &
               header = "tstep,time,particle_id,x,y,z,u,v,w", &
               overwrite = .true.)
       end if
       this%hdf5_output = .false.
    case ("h5", "hdf5")
       this%hdf5_output = .true.
       call this%init_hdf5(output_path, this%inertia)
    case default
       call neko_error("lpt output_filename must end in .csv, .h5, or .hdf5")
    end select
    this%snapshots_in_file = 0
  end subroutine lpt_output_init_file

  !> Get the active output path, adding a file index when chunked output is on.
  !! @param output_path Path selected for the current output file.
  subroutine lpt_output_current_path(this, output_path)
    class(lpt_output_t), intent(in) :: this
    character(len=*), intent(out) :: output_path
    character(len=1024) :: path
    character(len=1024) :: name
    character(len=1024) :: suffix


    call filename_split(trim(this%output_path), path, name, suffix)
    write(output_path, '(A,A,A,I0,A)') trim(path), trim(name), &
         "_", this%output_file_index, trim(suffix)

  end subroutine lpt_output_current_path

  !> Initialise an HDF5 LPT trajectory file.
  !! @param output_path HDF5 file path to initialise.
  !! @param inertia Whether particle records include diameter and density.
  subroutine lpt_output_init_hdf5(this, output_path, inertia)
    class(lpt_output_t), intent(inout) :: this
    character(len=*), intent(in) :: output_path
    logical, intent(in) :: inertia
    integer :: file_unit
    integer :: ierr
    integer :: out_int

    ! The granular HDF5 API appends to existing datasets, so ensure this run
    ! starts with a fresh file before hdf5_file_t creates it collectively.
    if (pe_rank .eq. 0) then
       open(newunit = file_unit, file = trim(output_path), status = "replace", &
            action = "write", iostat = ierr)
       if (ierr .ne. 0) then
          call neko_error("Error while creating " // trim(output_path))
       end if
       close(file_unit, status = "delete")
    end if
    call MPI_Barrier(NEKO_COMM, ierr)

    call this%output_file%init(output_path)

    select type (ft => this%output_file%file_type)
    type is (hdf5_file_t)
       call ft%set_overwrite(.false.)
       call ft%open("w")
       call ft%set_active_group("lpt")
       out_int = 3
       call ft%write_attribute("FormatVersion", out_int)
       out_int = this%n_data
       call ft%write_attribute("NColumns", out_int)
       if (inertia) then
          out_int = 7
       else
          out_int = 5
       end if
       call ft%write_attribute("NCategories", out_int)
       if (inertia) then
          out_int = 1
       else
          out_int = 0
       end if
       call ft%write_attribute("Inertia", out_int)
       out_int = 0
       call ft%write_attribute("NSteps", out_int)
       call ft%close()
    class default
       call neko_error("lpt internal error: expected hdf5_file_t")
    end select

    this%hdf5_output_count = 0
  end subroutine lpt_output_init_hdf5

  !> Write one trajectory snapshot.
  !! @param local_data Particle data owned by this rank.
  !! @param n_local Number of local particles.
  subroutine lpt_output_write(this, local_data, n_local)
    class(lpt_output_t), intent(inout) :: this
    integer, intent(in) :: n_local
    real(kind=rp), intent(in) :: local_data(this%n_data, n_local)

    if (this%snapshots_per_file .gt. 0 .and. &
         this%snapshots_in_file .ge. this%snapshots_per_file) then
       this%output_file_index = this%output_file_index + 1
       call this%init_file()
    end if

    if (this%hdf5_output) then
       call this%write_hdf5(local_data, n_local)
    else
       call this%write_csv(local_data, n_local)
    end if
    this%snapshots_in_file = this%snapshots_in_file + 1
  end subroutine lpt_output_write

  !> Append one local LPT trajectory snapshot to an HDF5 file.
  !! @param local_data Particle data owned by this rank.
  !! @param n_local Number of local particles.
  subroutine lpt_output_write_hdf5(this, local_data, n_local)
    class(lpt_output_t), intent(inout) :: this
    integer, intent(in) :: n_local
    real(kind=rp), intent(in) :: local_data(this%n_data, n_local)
    type(matrix_t) :: tsteps
    type(matrix_t) :: t
    type(matrix_t) :: ids
    type(matrix_t) :: position
    type(matrix_t) :: velocity
    type(matrix_t) :: diameter
    type(matrix_t) :: density
    integer :: out_int

    call tsteps%init(1, n_local, "tsteps")
    call t%init(1, n_local, "t")
    call ids%init(1, n_local, "ids")
    call position%init(3, n_local, "position")
    call velocity%init(3, n_local, "velocity")
    if (this%inertia) then
       call diameter%init(1, n_local, "diameter")
       call density%init(1, n_local, "density")
    end if

    if (n_local .gt. 0) then
       tsteps%x(1, :) = local_data(1, :)
       t%x(1, :) = local_data(2, :)
       ids%x(1, :) = local_data(3, :)
       position%x = local_data(4:6, :)
       velocity%x = local_data(7:9, :)
       if (this%inertia) then
          diameter%x(1, :) = local_data(10, :)
          density%x(1, :) = local_data(11, :)
       end if
    end if

    select type (ft => this%output_file%file_type)
    type is (hdf5_file_t)
       call ft%open("w")
       call ft%set_active_group("lpt")
       out_int = this%hdf5_output_count + 1
       call ft%write_attribute("NSteps", out_int)
       call ft%write_dataset(tsteps)
       call ft%write_dataset(t)
       call ft%write_dataset(ids)
       call ft%write_dataset(position)
       call ft%write_dataset(velocity)
       if (this%inertia) then
          call ft%write_dataset(diameter)
          call ft%write_dataset(density)
       end if
       call ft%close()
    class default
       call neko_error("lpt internal error: expected hdf5_file_t")
    end select

    this%hdf5_output_count = this%hdf5_output_count + 1
    call tsteps%free()
    call t%free()
    call ids%free()
    call position%free()
    call velocity%free()
    if (this%inertia) then
       call diameter%free()
       call density%free()
    end if
  end subroutine lpt_output_write_hdf5

  !> Write one trajectory snapshot to CSV by gathering particle data to rank 0.
  !! @param local_data Particle data owned by this rank.
  !! @param n_local Number of local particles.
  subroutine lpt_output_write_csv(this, local_data, n_local)
    class(lpt_output_t), intent(inout) :: this
    integer, intent(in) :: n_local
    real(kind=rp), intent(in) :: local_data(this%n_data, n_local)
    type(matrix_t) :: block
    real(kind=rp), allocatable :: global_data(:, :)
    integer, allocatable :: n_local_particles_per_rank(:)
    integer, allocatable :: recvcounts(:)
    integer, allocatable :: displs(:)
    integer :: total_particles
    integer :: i
    integer :: ierr

    if (pe_rank .eq. 0) then
       allocate(n_local_particles_per_rank(pe_size))
    else
       allocate(n_local_particles_per_rank(0))
    end if
    call MPI_Gather(n_local, 1, MPI_INTEGER, n_local_particles_per_rank, 1, &
         MPI_INTEGER, 0, NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       allocate(recvcounts(pe_size))
       allocate(displs(pe_size))
       recvcounts = 0
       displs = 0

       total_particles = 0
       do i = 1, pe_size
          total_particles = total_particles + n_local_particles_per_rank(i)
          recvcounts(i) = this%n_data * n_local_particles_per_rank(i)
          displs(i) = this%n_data * &
               (total_particles - n_local_particles_per_rank(i))
       end do

       allocate(global_data(this%n_data, total_particles))
    else
       allocate(recvcounts(0))
       allocate(displs(0))
       allocate(global_data(this%n_data, 0))
    end if

    call MPI_Gatherv(local_data, this%n_data * n_local, MPI_REAL_PRECISION, &
         global_data, recvcounts, displs, MPI_REAL_PRECISION, 0, &
         NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       call block%init(total_particles, this%n_data)
       call trsp(block%x, total_particles, global_data, this%n_data)
       call this%output_file%write(block)
       call block%free()
    end if

    deallocate(global_data)
    deallocate(n_local_particles_per_rank)
    deallocate(recvcounts)
    deallocate(displs)
  end subroutine lpt_output_write_csv

  !> Free the output writer and reset file/output counters.
  subroutine lpt_output_free(this)
    class(lpt_output_t), intent(inout) :: this

    call this%output_file%free()
    this%output_path = ""
    this%hdf5_output = .false.
    this%inertia = .false.
    this%snapshots_per_file = 0
    this%snapshots_in_file = 0
    this%output_file_index = 0
    this%hdf5_output_count = 0
    this%n_data = 0
  end subroutine lpt_output_free

end module lpt_output
