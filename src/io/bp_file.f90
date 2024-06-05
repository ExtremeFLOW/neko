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
  use space, only: space_t
  use mesh, only: mesh_t
  use structs, only: array_ptr_t
  use fld_file_data
  use utils
  use comm
  use adios2

  real(kind=dp), private, allocatable :: tmp_dp(:)
  real(kind=sp), private, allocatable :: tmp_sp(:)

  type(adios2_adios)    :: adios
  type(adios2_io)       :: ioWriter
  type(adios2_engine)   :: bpWriter
  type(adios2_engine)   :: bpReader

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
    integer, allocatable :: idx(:)
    logical :: write_mesh, write_velocity, write_pressure, write_temperature
    integer :: adios2_type
    type(adios2_variable) :: variable_idx, variable_hdr, variable_msh
    type(adios2_variable) :: variable_v, variable_p, variable_temp
    integer(kind=8), dimension(1) :: shape_dims_idx, start_dims_idx, count_dims_idx
    integer(kind=8), dimension(1) :: shape_dims_msh, start_dims_msh, count_dims_msh
    integer(kind=8), dimension(1) :: shape_dims_v, start_dims_v, count_dims_v
    integer(kind=8), dimension(1) :: shape_dims_p, start_dims_p, count_dims_p
    integer(kind=8), dimension(1) :: shape_dims_temp, start_dims_temp, count_dims_temp
  
    if (present(t)) then
       time = real(t,dp)
    else
       time = 0d0
    end if

    nullify(msh)
    nullify(dof)
    nullify(Xh)
    !n_scalar_fields = 0
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
       case (5:99)
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

    if (this%dp_precision) then
       if (.not. allocated(tmp_dp)) allocate(tmp_dp(gdim*n))
    else
       if (.not. allocated(tmp_sp)) allocate(tmp_sp(gdim*n))
    end if

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
    !if (n_scalar_fields .gt. 0 ) then
       !rdcode(i) = 'S'
       !i = i + 1
       !write(rdcode(i), '(i1)') (n_scalar_fields)/10
       !i = i + 1
       !write(rdcode(i), '(i1)') (n_scalar_fields) - 10*((n_scalar_fields)/10)
       !i = i + 1
    !end if

    ! Create binary header information
    write(hdr, 1) adios2_type, lx, ly, lz, glb_nelv, glb_nelv,&
         time, this%counter, 1, 1, (rdcode(i),i=1,10)
1   format('#std',1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,&
         1x,i9,1x,i6,1x,i6,1x,10a)

    ! Adapt filename with counter
    !> @todo write into single file with multiple steps
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
    shape_dims_idx = (int(glb_nelv, i8))
    start_dims_idx = (int(offset_el, i8))
    count_dims_idx = (int(nelv, i8))
    call adios2_inquire_variable(variable_idx, ioWriter, 'idx', ierr)
    if (.not.variable_idx%valid) then
       call adios2_define_variable(variable_idx, ioWriter, 'idx', adios2_type_integer4,&
            size(shape_dims_idx), shape_dims_idx, start_dims_idx, count_dims_idx, .false., ierr)
    else
       call adios2_set_shape(variable_idx, size(shape_dims_idx), shape_dims_idx, ierr)
       call adios2_set_selection(variable_idx, size(start_dims_idx), &
            start_dims_idx, count_dims_idx, ierr)
    end if
    call adios2_put(bpWriter, variable_idx, idx, adios2_mode_sync, ierr)

    deallocate(idx)

    if (write_mesh) then
       shape_dims_msh = (int(glb_nelv, i8) * int(lxyz, i8) * int(gdim, i8))
       start_dims_msh = (int(offset_el, i8) * int(lxyz, i8) * int(gdim, i8))
       count_dims_msh = (int(nelv, i8) * int(lxyz, i8) * int(gdim, i8))
       call adios2_inquire_variable(variable_msh, ioWriter, 'mesh', ierr)
       if (.not.variable_msh%valid) then
          call adios2_define_variable(variable_msh, ioWriter, 'mesh', adios2_type,&
               size(shape_dims_msh), shape_dims_msh, start_dims_msh, count_dims_msh, .false., ierr)
       else
          call adios2_set_shape(variable_msh, size(shape_dims_msh), shape_dims_msh, ierr)
          call adios2_set_selection(variable_msh, size(start_dims_msh), &
               start_dims_msh, count_dims_msh, ierr)
       end if
       call bp_file_buffer_vector_field(this, x%ptr, y%ptr, z%ptr, n, gdim, lxyz, nelv)
       if (this%dp_precision) then
          call adios2_put(bpWriter, variable_msh, tmp_dp, adios2_mode_sync, ierr)
       else
          call adios2_put(bpWriter, variable_msh, tmp_sp, adios2_mode_sync, ierr)
       end if
    end if

    if (write_velocity) then
       shape_dims_v = (int(glb_nelv, i8) * int(lxyz, i8) * int(gdim, i8))
       start_dims_v = (int(offset_el, i8) * int(lxyz, i8) * int(gdim, i8))
       count_dims_v = (int(nelv, i8) * int(lxyz, i8) * int(gdim, i8))
       call adios2_inquire_variable(variable_v, ioWriter, 'velocity', ierr)
       if (.not.variable_v%valid) then
          call adios2_define_variable(variable_v, ioWriter, 'velocity', adios2_type,&
               size(shape_dims_v), shape_dims_v, start_dims_v, count_dims_v, .false., ierr)
       else
          call adios2_set_shape(variable_v, size(shape_dims_v), shape_dims_v, ierr)
          call adios2_set_selection(variable_v, size(start_dims_v), &
               start_dims_v, count_dims_v, ierr)
       end if
       call bp_file_buffer_vector_field(this, u%ptr, v%ptr, w%ptr, n, gdim, lxyz, nelv)
       if (this%dp_precision) then
          call adios2_put(bpWriter, variable_v, tmp_dp, adios2_mode_sync, ierr)
       else
          call adios2_put(bpWriter, variable_v, tmp_sp, adios2_mode_sync, ierr)
       end if
    end if

    if (write_pressure) then
       shape_dims_p = (int(glb_nelv, i8) * int(lxyz, i8))
       start_dims_p = (int(offset_el, i8) * int(lxyz, i8))
       count_dims_p = (int(nelv, i8) * int(lxyz, i8))
       call adios2_inquire_variable(variable_p, ioWriter, 'pressure', ierr)
       if (.not.variable_p%valid) then
          call adios2_define_variable(variable_p, ioWriter, 'pressure', adios2_type,&
               size(shape_dims_p), shape_dims_p, start_dims_p, count_dims_p, .false., ierr)
       else
          call adios2_set_shape(variable_p, size(shape_dims_p), shape_dims_p, ierr)
          call adios2_set_selection(variable_p, size(start_dims_p), &
               start_dims_p, count_dims_p, ierr)
       end if
       call bp_file_buffer_field(this, p%ptr, n)
       if (this%dp_precision) then
          call adios2_put(bpWriter, variable_p, tmp_dp, adios2_mode_sync, ierr)
       else
          call adios2_put(bpWriter, variable_p, tmp_sp, adios2_mode_sync, ierr)
       end if
    end if

    if (write_temperature) then
       shape_dims_temp = (int(glb_nelv, i8) * int(lxyz, i8))
       start_dims_temp = (int(offset_el, i8) * int(lxyz, i8))
       count_dims_temp = (int(nelv, i8) * int(lxyz, i8))
       call adios2_inquire_variable(variable_temp, ioWriter, 'temperature', ierr)
       if (.not.variable_temp%valid) then
          call adios2_define_variable(variable_temp, ioWriter, 'temperature', adios2_type,&
               size(shape_dims_temp), shape_dims_temp, start_dims_temp, count_dims_temp, .false., ierr)
       else
          call adios2_set_shape(variable_temp, size(shape_dims_temp), shape_dims_temp, ierr)
          call adios2_set_selection(variable_temp, size(start_dims_temp), &
               start_dims_temp, count_dims_temp, ierr)
       end if
       call bp_file_buffer_field(this, tem%ptr, n)
       if (this%dp_precision) then
          call adios2_put(bpWriter, variable_temp, tmp_dp, adios2_mode_sync, ierr)
       else
          call adios2_put(bpWriter, variable_temp, tmp_sp, adios2_mode_sync, ierr)
       end if
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

    if (allocated(tmp_dp)) deallocate(tmp_dp)
    if (allocated(tmp_sp)) deallocate(tmp_sp)
  end subroutine bp_file_write

  subroutine bp_file_buffer_field(this, p, n)
    class(bp_file_t), intent(inout) :: this
    integer, intent(inout) :: n
    real(kind=rp), intent(inout) :: p(n)
    integer :: i

    if (this%dp_precision) then
       do i = 1, n
          tmp_dp(i) = real(p(i),dp)
       end do
    else
       do i = 1, n
          tmp_sp(i) = real(p(i),sp)
       end do
    end if

  end subroutine bp_file_buffer_field

  subroutine bp_file_buffer_vector_field(this, x, y, z, n, gdim, lxyz, nelv)
    class(bp_file_t), intent(inout) :: this
    integer, intent(in) :: n, gdim, lxyz, nelv
    real(kind=rp), intent(in) :: x(lxyz,nelv), y(lxyz,nelv), z(lxyz,nelv)
    integer :: i, el, j

    if (this%dp_precision) then
       i = 1
       do el = 1, nelv
          do j = 1, lxyz
             tmp_dp(i) = real(x(j,el),dp)
             i = i +1
          end do
          do j = 1, lxyz
             tmp_dp(i) = real(y(j,el),dp)
             i = i +1
          end do
          if (gdim .eq. 3) then
             do j = 1, lxyz
                tmp_dp(i) = real(z(j,el),dp)
                i = i +1
             end do
          end if
       end do
    else
       i = 1
       do el = 1, nelv
          do j = 1, lxyz
             tmp_sp(i) = real(x(j,el),sp)
             i = i +1
          end do
          do j = 1, lxyz
             tmp_sp(i) = real(y(j,el),sp)
             i = i +1
          end do
          if (gdim .eq. 3) then
             do j = 1, lxyz
                tmp_sp(i) = real(z(j,el),sp)
                i = i +1
             end do
          end if
       end do
    end if

  end subroutine bp_file_buffer_vector_field

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
