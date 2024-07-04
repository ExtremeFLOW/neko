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
  use dofmap, only: dofmap_t
  use vector, only: vector_t
  use space, only: space_t
  use mesh, only: mesh_t
  use structs, only: array_ptr_t
  use fld_file_data
  use utils
  use datadist
  use comm
  use adios2
  use buffer_1d
  use buffer_4d
  use buffer_4d_npar

  class(buffer_t), private, allocatable :: outbuf_points
  class(buffer_t), private, allocatable :: outbuf_npar

  type(adios2_adios)    :: adios
  type(adios2_io)       :: ioWriter
  type(adios2_io)       :: ioReader

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
    type(array_ptr_t) :: x, y, z, u, v, w, p, tem
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(space_t), pointer :: Xh
    real(kind=dp) :: time
    character(len=132) :: hdr
    character :: rdcode(10)
    character(len=8) :: id_str
    character(len=1024) :: fname
    character(len=1024) :: start_field
    integer :: i, ierr, n, suffix_pos, tslash_pos
    integer :: lx, ly, lz, lxyz, gdim, glb_nelv, nelv, offset_el
    integer :: npar
    integer, allocatable :: idx(:)
    type(array_ptr_t), allocatable :: scalar_fields(:)
    integer :: n_scalar_fields
    logical :: write_mesh, write_velocity, write_pressure, write_temperature
    integer :: adios2_type, layout
    type(adios2_engine)   :: bpWriter
    type(adios2_variable) :: variable_idx, variable_hdr, variable, variable_msh
    type(adios2_variable) :: variable_v, variable_p, variable_temp
    integer(kind=8), dimension(1) :: shape_dims, start_dims, count_dims
  
    if (present(t)) then
       time = real(t,dp)
    else
       time = 0d0
    end if

    nullify(msh)
    nullify(dof)
    nullify(Xh)
    n_scalar_fields = 0
    write_velocity = .false.
    write_pressure = .false.
    write_temperature = .false.

    !> @todo support for other input data types like in fld_file
    select type(data)
    type is (field_list_t)
       npar = data%size()
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
       case (5:99)
          p%ptr => data%items(1)%ptr%x(:,1,1,1)
          u%ptr => data%items(2)%ptr%x(:,1,1,1)
          v%ptr => data%items(3)%ptr%x(:,1,1,1)
          w%ptr => data%items(4)%ptr%x(:,1,1,1)
          tem%ptr => data%items(5)%ptr%x(:,1,1,1)
          n_scalar_fields = data%size() - 5
          allocate(scalar_fields(n_scalar_fields))
          do i = 1,n_scalar_fields
             scalar_fields(i)%ptr => data%items(i+5)%ptr%x(:,1,1,1)
          end do
          write_pressure = .true.
          write_velocity = .true.
          write_temperature = .true.
       case default
          call neko_error('This many fields not supported yet, bp_file')
       end select
       dof => data%dof(1)
    class default
       call neko_error('Invalid data')
    end select

    ! Fix things for pointers that do not exist in all data types...
    if (associated(dof)) then
       x%ptr => dof%x(:,1,1,1)
       y%ptr => dof%y(:,1,1,1)
       z%ptr => dof%z(:,1,1,1)
       msh => dof%msh
       Xh => dof%Xh
    end if

    if (associated(msh)) then
       nelv = msh%nelv
       glb_nelv = msh%glb_nelv
       offset_el = msh%offset_el
       gdim = msh%gdim
       ! Store global idx of each element
       allocate(idx(msh%nelv))
       do i = 1, msh%nelv
          idx(i) = msh%elements(i)%e%id()
       end do
    end if

    if (associated(Xh)) then
       lx = Xh%lx
       ly = Xh%ly
       lz = Xh%lz
    end if

    lxyz = lx*ly*lz
    n = nelv*lxyz

    if (this%dp_precision) then
       adios2_type = adios2_type_dp
    else
       adios2_type = adios2_type_real
    end if

    if (.not. allocated(outbuf_points)) allocate(buffer_1d_t::outbuf_points)
    call outbuf_points%init(this%dp_precision, gdim, glb_nelv, offset_el, nelv, lx, ly, lz)

    if (.not. allocated(outbuf_npar)) allocate(buffer_4d_npar_t::outbuf_npar)
    select type(outbuf_npar)
    type is (buffer_1d_t)
       call outbuf_npar%init(this%dp_precision, gdim, glb_nelv, offset_el, nelv, lx, ly, lz)
       layout = 1
    type is (buffer_4d_t)
       call outbuf_npar%init(this%dp_precision, gdim, glb_nelv, offset_el, nelv, lx, ly, lz)
       layout = 2
    type is (buffer_4d_npar_t)
       call outbuf_npar%init(this%dp_precision, npar, glb_nelv, offset_el, nelv, lx, ly, lz)
       layout = 4
    class default
       call neko_error('Invalid buffer')
    end select

    !
    ! Create fld header for NEKTON's multifile output
    !

    write_mesh = (this%counter .eq. this%start_counter)

    ! Build rdcode note that for field_t, we only support scalar
    ! fields at the moment
    rdcode = ' '
    i = 1
    if (write_mesh) then
       rdcode(i) = 'X'
       i = i + 1
    end if
    if (write_velocity) then
       rdcode(i) = 'U'
       i = i + 1
    end if
    if (write_pressure) then
       rdcode(i) = 'P'
       i = i + 1
    end if
    if (write_temperature) then
       rdcode(i) = 'T'
       i = i + 1
    end if
    if (n_scalar_fields .gt. 0 ) then
       rdcode(i) = 'S'
       i = i + 1
       write(rdcode(i), '(i1)') (n_scalar_fields)/10
       i = i + 1
       write(rdcode(i), '(i1)') (n_scalar_fields) - 10*((n_scalar_fields)/10)
       i = i + 1
    end if

    ! Create binary header information
    write(hdr, 1) adios2_type, lx, ly, lz, layout, glb_nelv,&
         time, this%counter, npar, (rdcode(i),i=1,10)
1   format('#std',1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,&
         1x,i9,1x,i6,1x,10a)

    ! Adapt filename with counter
    !> @todo write into single file with multiple steps
    !> @todo write into function because of code duplication
    suffix_pos = filename_suffix_pos(this%fname)
    write(id_str, '(i5.5,a)') this%counter, '.bp'
    fname = trim(this%fname(1:suffix_pos-1))//"0."//id_str

    if (.not. adios%valid) then
       !> @todo enable parsing XML filename
       call adios2_init(adios, 'adios2.xml', NEKO_COMM%mpi_val, ierr)
       call adios2_declare_io(ioWriter, adios, 'writer', ierr)
       call adios2_set_engine(ioWriter, 'BP5', ierr)
    end if

    call adios2_open(bpWriter, ioWriter, trim(fname), adios2_mode_write, NEKO_COMM%mpi_val, ierr)
    call adios2_begin_step(bpWriter, ierr)

    ! Write header
    call adios2_inquire_variable(variable_hdr, ioWriter, 'header', ierr)
    if (.not.variable_hdr%valid) then
       call adios2_define_variable(variable_hdr, ioWriter, 'header', adios2_type_character, ierr)
    end if
    call adios2_put(bpWriter, variable_hdr, hdr, adios2_mode_sync, ierr)

    ! Write element idxs
    shape_dims = (int(glb_nelv, i8))
    start_dims = (int(offset_el, i8))
    count_dims = (int(nelv, i8))
    call adios2_inquire_variable(variable_idx, ioWriter, 'idx', ierr)
    if (.not.variable_idx%valid) then
       call adios2_define_variable(variable_idx, ioWriter, 'idx', adios2_type_integer4,&
            size(shape_dims), shape_dims, start_dims, count_dims, .false., ierr)
    else
       call adios2_set_shape(variable_idx, size(shape_dims), shape_dims, ierr)
       call adios2_set_selection(variable_idx, size(start_dims), &
            start_dims, count_dims, ierr)
    end if
    call adios2_put(bpWriter, variable_idx, idx, adios2_mode_sync, ierr)

    deallocate(idx)

    if (write_mesh) then
       call outbuf_points%fill(x%ptr, n)
       call outbuf_points%define(variable, ioWriter, 'points-x', ierr)
       call outbuf_points%write(bpWriter, variable, ierr)
       call outbuf_points%fill(y%ptr, n)
       call outbuf_points%define(variable, ioWriter, 'points-y', ierr)
       call outbuf_points%write(bpWriter, variable, ierr)
       call outbuf_points%fill(z%ptr, n)
       call outbuf_points%define(variable, ioWriter, 'points-z', ierr)
       call outbuf_points%write(bpWriter, variable, ierr)
    end if

    if (write_velocity) then
       call outbuf_npar%fill(u%ptr, n)
       if (layout .le. 3) then
          call outbuf_npar%define(variable, ioWriter, 'velocity-u', ierr)
          call outbuf_npar%write(bpWriter, variable, ierr)
       endif
       call outbuf_npar%fill(v%ptr, n)
       if (layout .le. 3) then
          call outbuf_npar%define(variable, ioWriter, 'velocity-v', ierr)
          call outbuf_npar%write(bpWriter, variable, ierr)
       endif
       call outbuf_npar%fill(w%ptr, n)
       if (layout .le. 3) then
          call outbuf_npar%define(variable, ioWriter, 'velocity-w', ierr)
          call outbuf_npar%write(bpWriter, variable, ierr)
       endif
    end if

    if (write_pressure) then
       call outbuf_npar%fill(p%ptr, n)
       if (layout .le. 3) then
          call outbuf_npar%define(variable, ioWriter, 'pressure', ierr)
          call outbuf_npar%write(bpWriter, variable, ierr)
       endif
    end if

    if (write_temperature) then
       call outbuf_npar%fill(tem%ptr, n)
       if (layout .le. 3) then
          call outbuf_npar%define(variable, ioWriter, 'temperature', ierr)
          call outbuf_npar%write(bpWriter, variable, ierr)
       endif
    end if

    do i = 1, n_scalar_fields
       call outbuf_npar%fill(scalar_fields(i)%ptr, n)
       if (layout .le. 3) then
          write(id_str, '(a,i1,i1)') 's', i/10, i-10*(i/10)
          call outbuf_npar%define(variable, ioWriter, trim(id_str), ierr)
          call outbuf_npar%write(bpWriter, variable, ierr)
       endif
    end do

    if (layout .gt. 3) then
       call outbuf_npar%define(variable, ioWriter, 'fields', ierr)
       call outbuf_npar%write(bpWriter, variable, ierr)
    end if

    call adios2_end_step(bpWriter, ierr)
    call adios2_close(bpWriter, ierr)

    ! Write metadata file
    if (pe_rank .eq. 0) then
       tslash_pos = filename_tslash_pos(this%fname)
       write(start_field,"(I5,A7)") this%start_counter,'.adios2'
       open(unit=9, file=trim(this%fname(1:suffix_pos-1))//trim(adjustl(start_field)), &
            status='replace')
       write(9, fmt='(A,A,A)') 'filetemplate:         ', &
            this%fname(tslash_pos+1:suffix_pos-1),'%01d.%05d.bp'
       write(9, fmt='(A,i5)') 'firsttimestep: ', this%start_counter
       write(9, fmt='(A,i5)') 'numtimesteps: ', (this%counter + 1)-this%start_counter
       write(9, fmt='(A)') 'type: adios2-bp'
       close(9)
    end if

    this%counter = this%counter + 1

    if (allocated(outbuf_points)) deallocate(outbuf_points)
    if (allocated(outbuf_npar)) deallocate(outbuf_npar)
  end subroutine bp_file_write

  !> Load a field from a ADIOS2 bp file
  subroutine bp_file_read(this, data)
    class(bp_file_t) :: this
    class(*), target, intent(inout) :: data
    character(len=132) :: hdr
    integer :: ierr, suffix_pos, i, j
    character(len=1024) :: fname, meta_fname, string
    logical :: meta_file, read_mesh, read_velocity, read_pressure
    logical :: read_temp
    character(len=8) :: id_str
    integer :: lx, ly, lz, glb_nelv, counter, lxyz
    integer :: layout, npar
    integer :: adios2_type, n_scalars, n
    real(kind=rp) :: time
    type(linear_dist_t) :: dist
    character :: rdcode(10), temp_str(4)
    class(buffer_t), allocatable :: inpbuf_points, inpbuf
    type(adios2_engine)   :: bpReader
    type(adios2_variable) :: variable_hdr, variable_idx, variable
    integer(kind=8), dimension(1) :: start_dims, count_dims

    select type(data)
    type is (fld_file_data_t)
       suffix_pos = filename_suffix_pos(this%fname)
       meta_fname = trim(this%fname(1:suffix_pos-1))
       call filename_chsuffix(meta_fname, meta_fname,'adios2')

       !> @ todo debug and check if correct strings and filenames are extracted
       inquire(file=trim(meta_fname), exist=meta_file)
       if (meta_file .and. data%meta_nsamples .eq. 0) then
          if (pe_rank .eq. 0) then
             open(unit=9, file=trim(meta_fname))
             read(9, fmt='(A)') string
             read(string(14:),fmt='(A)') string
             string = trim(string)
             data%fld_series_fname = string(:scan(trim(string), '%')-1)
             data%fld_series_fname = trim(data%fld_series_fname)//'0'
             read(9, fmt='(A)') string
             read(string(scan(string,':')+1:),*) data%meta_start_counter
             read(9, fmt='(A)') string
             read(string(scan(string,':')+1:),*) data%meta_nsamples

             close(9)
             write(*,*) 'Reading meta file for bp series'
             write(*,*) 'Name: ', trim(data%fld_series_fname)
             write(*,*) 'Start counter: ', data%meta_start_counter, 'Nsamples: ', data%meta_nsamples
          end if
          call MPI_Bcast(data%fld_series_fname, 1024, MPI_CHARACTER, 0, NEKO_COMM, ierr)
          call MPI_Bcast(data%meta_start_counter, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
          call MPI_Bcast(data%meta_nsamples, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
          if(this%counter .eq. 0) this%counter = data%meta_start_counter
       end if

       if (meta_file) then
          write(id_str, '(i5.5,a)') this%counter, '.bp'
          fname = trim(data%fld_series_fname)//'.'//id_str
          if (this%counter .ge. data%meta_nsamples+data%meta_start_counter) then
             call neko_error('Trying to read more bp files than exist')
          end if
       else
          !> @todo write into function because of code duplication
          !suffix_pos = filename_suffix_pos(this%fname)
          !write(id_str, '(i5.5,a)') this%counter, '.bp'
          !fname = trim(this%fname(1:suffix_pos-1))//'.'//id_str
          fname = this%fname
       end if

       if (.not.adios%valid) then
          !> @todo enable parsing XML filename
          call adios2_init(adios, 'adios2.xml', NEKO_COMM%mpi_val, ierr)
       end if
       if (.not.ioReader%valid) then
          call adios2_declare_io(ioReader, adios, 'reader', ierr)
          call adios2_set_engine(ioReader, 'BP5', ierr)
       end if

       !> @todo check if engines and variables should be brought in to local subroutine scope
       call adios2_open(bpReader, ioReader, trim(fname), adios2_mode_read, NEKO_COMM%mpi_val, ierr)
       call adios2_begin_step(bpReader, ierr)

       ! Read header and adjust data accordingly
       if (.not.variable_hdr%valid) then
          call adios2_inquire_variable(variable_hdr, ioReader, 'header', ierr)
       end if
       call adios2_get(bpReader, variable_hdr, hdr, adios2_mode_sync, ierr)

       read(hdr, 1) temp_str, adios2_type, lx, ly, lz, layout, glb_nelv,&
          time, counter, npar, (rdcode(i),i=1,10)
1      format(4a,1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,&
         1x,i9,1x,i6,1x,10a)
       if (data%nelv .eq. 0) then
          dist = linear_dist_t(glb_nelv, pe_rank, pe_size, NEKO_COMM)
          data%nelv = dist%num_local()
          data%offset_el = dist%start_idx()
       end if
       data%lx = lx
       data%ly = ly
       data%lz = lz
       data%glb_nelv = glb_nelv
       data%t_counter = counter
       data%time = time
       lxyz = lx * ly * lz
       n = lxyz * data%nelv

       if (lz .eq. 1) then
          data%gdim = 2
       else
          data%gdim = 3
       end if

       if (adios2_type .eq. adios2_type_double) then
          this%dp_precision = .true.
       else
          this%dp_precision = .false.
       end if

       if (.not. allocated(inpbuf_points)) allocate(buffer_1d_t::inpbuf_points)
       call inpbuf_points%init(this%dp_precision, data%gdim, data%glb_nelv, data%offset_el, &
            data%nelv, lx, ly, lz)

       write(*,*) "layout ", layout
       if (layout .eq. 1) then
          if (.not. allocated(inpbuf)) allocate(buffer_1d_t::inpbuf)
       else if (layout .eq. 2) then
          if (.not. allocated(inpbuf)) allocate(buffer_4d_t::inpbuf)
       else if (layout .eq. 4) then
          if (.not. allocated(inpbuf)) allocate(buffer_4d_npar_t::inpbuf)
       end if

       select type(inpbuf)
       type is (buffer_1d_t)
          call inpbuf%init(this%dp_precision, data%gdim, data%glb_nelv, data%offset_el, &
               data%nelv, lx, ly, lz)
       type is (buffer_4d_t)
          call inpbuf%init(this%dp_precision, data%gdim, data%glb_nelv, data%offset_el, &
               data%nelv, lx, ly, lz)
       type is (buffer_4d_npar_t)
          call inpbuf%init(this%dp_precision, npar, data%glb_nelv, data%offset_el, &
               data%nelv, lx, ly, lz)
       class default
          call neko_error('Invalid buffer')
       end select

       i = 1
       read_mesh = .false.
       read_velocity = .false.
       read_pressure = .false.
       read_temp = .false.
       if (rdcode(i) .eq. 'X') then
          read_mesh = .true.
          if (data%x%n .ne. n) call data%x%init(n)
          if (data%y%n .ne. n) call data%y%init(n)
          if (data%z%n .ne. n) call data%z%init(n)
          i = i + 1
       end if
       if (rdcode(i) .eq. 'U') then
          read_velocity = .true.
          if (data%u%n .ne. n) call data%u%init(n)
          if (data%v%n .ne. n) call data%v%init(n)
          if (data%w%n .ne. n) call data%w%init(n)
          i = i + 1
       end if
       if (rdcode(i) .eq. 'P') then
          read_pressure = .true.
          if (data%p%n .ne. n) call data%p%init(n)
          i = i + 1
       end if
       if (rdcode(i) .eq. 'T') then
          read_temp = .true.
          if (data%t%n .ne. n) call data%t%init(n)
          i = i + 1
       end if
       n_scalars = 0
       if (rdcode(i) .eq. 'S') then
          i = i + 1
          read(rdcode(i),*) n_scalars
          n_scalars = n_scalars*10
          i = i + 1
          read(rdcode(i),*) j
          n_scalars = n_scalars+j
          i = i + 1
          if (allocated(data%s)) then
             if (data%n_scalars .ne. n_scalars) then
                do j = 1, data%n_scalars
                   call data%s(j)%free()
                end do
                deallocate(data%s)
                data%n_scalars = n_scalars
                allocate(data%s(n_scalars))
                do j = 1, data%n_scalars
                   call data%s(j)%init(n)
                end do
             end if
          else
             data%n_scalars = n_scalars
             allocate(data%s(data%n_scalars))
             do j = 1, data%n_scalars
                call data%s(j)%init(n)
             end do
          end if
          i = i + 1
       end if

       if (allocated(data%idx)) then
          if (size(data%idx) .ne. data%nelv) then
             deallocate(data%idx)
             allocate(data%idx(data%nelv))
          end if
       else
          allocate(data%idx(data%nelv))
       end if

       ! Read element idxs
       start_dims = (int(data%offset_el, i8))
       count_dims = (int(data%nelv, i8))
       call adios2_inquire_variable(variable_idx, ioReader, 'idx', ierr)
       if (variable_idx%valid) then
          call adios2_set_selection(variable_idx, size(start_dims), &
               start_dims, count_dims, ierr)
       end if
       call adios2_get(bpReader, variable_idx, data%idx, adios2_mode_sync, ierr)

       if (read_mesh) then
          call inpbuf_points%inquire(variable, ioReader, 'points-x', ierr)
          call inpbuf_points%read(bpReader, variable, ierr)
          call inpbuf_points%copy(data%x)
          call inpbuf_points%inquire(variable, ioReader, 'points-y', ierr)
          call inpbuf_points%read(bpReader, variable, ierr)
          call inpbuf_points%copy(data%y)
          call inpbuf_points%inquire(variable, ioReader, 'points-z', ierr)
          call inpbuf_points%read(bpReader, variable, ierr)
          call inpbuf_points%copy(data%z)
       end if

       if (layout .gt. 3) then
          call inpbuf%inquire(variable, ioReader, 'fields', ierr)
          call inpbuf%read(bpReader, variable, ierr)
       end if

       if (read_velocity) then
          if (layout .le. 3) then
             call inpbuf%inquire(variable, ioReader, 'velocity-u', ierr)
             call inpbuf%read(bpReader, variable, ierr)
          end if
          call inpbuf%copy(data%u)
          if (layout .le. 3) then
             call inpbuf%inquire(variable, ioReader, 'velocity-v', ierr)
             call inpbuf%read(bpReader, variable, ierr)
          end if
          call inpbuf%copy(data%v)
          if (layout .le. 3) then
             call inpbuf%inquire(variable, ioReader, 'velocity-w', ierr)
             call inpbuf%read(bpReader, variable, ierr)
          end if
          call inpbuf%copy(data%w)
       end if

       if (read_pressure) then
          if (layout .le. 3) then
             call inpbuf%inquire(variable, ioReader, 'pressure', ierr)
             call inpbuf%read(bpReader, variable, ierr)
          endif
          call inpbuf%copy(data%p)
       end if

       if (read_temp) then
          if (layout .le. 3) then
             call inpbuf%inquire(variable, ioReader, 'temperature', ierr)
             call inpbuf%read(bpReader, variable, ierr)
          endif
          call inpbuf%copy(data%t)
       end if

       do i = 1, n_scalars
          if (layout .le. 3) then
             write(id_str, '(a,i1,i1)') 's', i/10, i-10*(i/10)
             call inpbuf%inquire(variable, ioReader, trim(id_str), ierr)
             call inpbuf%read(bpReader, variable, ierr)
          endif
          call inpbuf%copy(data%s(i))
       end do

       call adios2_end_step(bpReader, ierr)
       call adios2_close(bpReader, ierr)

       this%counter = this%counter + 1

       if (allocated(inpbuf_points)) deallocate(inpbuf_points)
       if (allocated(inpbuf)) deallocate(inpbuf)
    class default
       call neko_error('Currently we only read into fld_file_data_t,&
                        please use that data structure instead.&
                        (output_format.adios2)')
    end select

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
