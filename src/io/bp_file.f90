! Copyright (c) 2024, Gregor Weiss (HLRS)
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
!> ADIOS2 bp file format
!! @details this module defines the interface to write ADIOS2 files
!! to use compression
module bp_file
  use generic_file
  use field_list
  use structs, only: array_ptr_t
  use fld_file_data
  use utils
  use comm
  use adios2

  real(kind=dp), private, allocatable :: tmp_dp(:)
  real(kind=sp), private, allocatable :: tmp_sp(:)

  type(adios2_adios)  :: adios
  type(adios2_io)     :: ioWriter
  type(adios2_engine) :: bpWriter

  !> Interface for ADIOS2 bp files
  type, public, extends(generic_file_t) :: bp_file_t
     logical :: dp_precision = .false. !< Precision of output data
   contains
     procedure :: read => bp_file_read
     procedure :: write => bp_file_write
     procedure :: set_precision => bp_file_set_precision
  end type bp_file_t

contains

  subroutine bp_file_write(this, data, t)
    class(bp_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    type(array_ptr_t) :: u, v, w, p, tem
    real(kind=dp) :: time
    character(len=8) :: id_str
    character(len=1024) :: fname
    character(len=1024) :: start_field
    integer :: suffix_pos, tslash_pos
    logical :: write_velocity, write_pressure, write_temperature
  
    if (present(t)) then
       time = real(t,dp)
    else
       time = 0d0
    end if

    write_velocity = .false.
    write_pressure = .false.
    write_temperature = .false.

    !> @todo support for other input data types like in fld_file
    select type(data)
    type is (field_list_t)
       select case (data%size())
       case (1)
          p%ptr => data%items(1)%ptr%x(:,1,1,1)
          write_pressure = .true.
       case (2)
          p%ptr => data%items(1)%ptr%x(:,1,1,1)
          tem%ptr => data%items(2)%ptr%x(:,1,1,1)
          write_pressure = .true.
          write_temperature = .true.
       case (3)
          u%ptr => data%items(1)%ptr%x(:,1,1,1)
          v%ptr => data%items(2)%ptr%x(:,1,1,1)
          w%ptr => data%items(3)%ptr%x(:,1,1,1)
          write_velocity = .true.
       case (4)
          p%ptr => data%items(1)%ptr%x(:,1,1,1)
          u%ptr => data%items(2)%ptr%x(:,1,1,1)
          v%ptr => data%items(3)%ptr%x(:,1,1,1)
          w%ptr => data%items(4)%ptr%x(:,1,1,1)
          write_pressure = .true.
          write_velocity = .true.
       case (5)
          p%ptr => data%items(1)%ptr%x(:,1,1,1)
          u%ptr => data%items(2)%ptr%x(:,1,1,1)
          v%ptr => data%items(3)%ptr%x(:,1,1,1)
          w%ptr => data%items(4)%ptr%x(:,1,1,1)
          tem%ptr => data%items(5)%ptr%x(:,1,1,1)
          !> @todo incorporate scalar fields again (see below)
          write_pressure = .true.
          write_velocity = .true.
          write_temperature = .true.
       case default
          !> @todo incorporate scalar fields again (see above)
          !call neko_error('This many fields not supported yet, bp_file')
       end select
    class default
       call neko_error('Invalid data')
    end select

    ! Adapt filename with counter
    !> @todo write into single file with multiple steps
    suffix_pos = filename_suffix_pos(this%fname)
    write(id_str, '(i5.5,a)') this%counter, '.bp'
    fname = trim(this%fname(1:suffix_pos-1))//id_str

    if (.not. adios%valid) then
       !> @todo enable parsing XML filename
       call adios2_init(adios, 'adios2.xml', NEKO_COMM%mpi_val, ierr)
       call adios2_declare_io(ioWriter, adios, 'writer', ierr)
       call adios2_set_engine(ioWriter, 'BP5', ierr)
    end if

    call adios2_open(bpWriter, ioWriter, trim(fname), adios2_mode_write, NEKO_COMM%mpi_val, ierr)
    call adios2_begin_step(bpWriter, ierr)

    call adios2_end_step(bpWriter, ierr)
    call adios2_close(bpWriter, ierr)

    ! Write metadata file
    if (pe_rank .eq. 0) then
       tslash_pos = filename_tslash_pos(this%fname)
       write(start_field,"(I5,A7)") this%start_counter,'.adios2'
       write(*,*) "START_FIELD ", trim(start_field)
       open(unit=9, file=trim(this%fname(1:suffix_pos-1))//trim(adjustl(start_field)), &
            status='replace')
       write(9, fmt='(A,A,A)') 'filetemplate:         ', &
            this%fname(tslash_pos+1:suffix_pos-1),'%05d.bp'
       write(9, fmt='(A,i5)') 'firsttimestep: ', this%start_counter
       write(9, fmt='(A,i5)') 'numtimesteps: ', (this%counter + 1)-this%start_counter
       write(9, fmt='(A)') 'type: adios2-bp'
       close(9)
    end if

    this%counter = this%counter + 1

  end subroutine bp_file_write

  !> Load a field from a ADIOS2 bp file
  subroutine bp_file_read(this, data)
    class(bp_file_t) :: this
    class(*), target, intent(inout) :: data

  end subroutine bp_file_read

  subroutine bp_file_set_precision(this, precision)
    class(bp_file_t) :: this
    integer, intent(in) :: precision

    if (precision .eq. dp) then
       this%dp_precision = .true.
    else if (precision .eq. sp) then
       this%dp_precision = .false.
    else
       call neko_error('Invalid precision')
    end if

  end subroutine bp_file_set_precision

end module bp_file
