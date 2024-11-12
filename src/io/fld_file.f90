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
!> NEKTON fld file format
!! @details this module defines interface to write NEKTON's fld fields
module fld_file
  use num_types, only: rp, dp, sp, i8
  use generic_file, only: generic_file_t
  use field, only: field_t
  use field_list, only: field_list_t
  use dofmap, only: dofmap_t
  use space, only: space_t
  use structs, only: array_ptr_t
  use vector, only: vector_t
  use fld_file_data, only: fld_file_data_t
  use mean_flow, only: mean_flow_t
  use mean_sqr_flow, only: mean_sqr_flow_t
  use vector, only : vector_t
  use space, only : space_t
  use mesh, only : mesh_t
  use utils, only: filename_suffix_pos, filename_tslash_pos, filename_chsuffix
  use utils, only: neko_error
  use comm
  use datadist, only: linear_dist_t
  use math, only: vlmin, vlmax
  use neko_mpi_types, only: MPI_CHARACTER_SIZE, MPI_DOUBLE_PRECISION_SIZE, &
       MPI_REAL_SIZE, MPI_INTEGER_SIZE
  implicit none
  private

  real(kind=dp), private, allocatable :: tmp_dp(:)
  real(kind=sp), private, allocatable :: tmp_sp(:)

  !> Interface for NEKTON fld files
  type, public, extends(generic_file_t) :: fld_file_t
     logical :: dp_precision = .false. !< Precision of output data
   contains
     procedure :: read => fld_file_read
     procedure :: write => fld_file_write
     procedure :: set_precision => fld_file_set_precision
  end type fld_file_t


contains

  !> Write fields to a NEKTON fld file
  !! @note currently limited to double precision data
  subroutine fld_file_write(this, data, t)
    class(fld_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    type(array_ptr_t) :: x, y, z, u, v, w, p, tem
    real(kind=rp), allocatable, target :: tempo(:)
    type(mesh_t), pointer :: msh
    type(dofmap_t), pointer :: dof
    type(space_t), pointer :: Xh
    real(kind=dp) :: time
    character(len= 132) :: hdr
    character :: rdcode(10)
    character(len=6) :: id_str
    character(len= 1024) :: fname
    character(len= 1024) :: start_field
    integer :: i, ierr, n, suffix_pos, tslash_pos
    integer :: lx, ly, lz, lxyz, gdim, glb_nelv, nelv, offset_el
    integer, allocatable :: idx(:)
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset, temp_offset
    real(kind=sp), parameter :: test_pattern = 6.54321
    type(array_ptr_t), allocatable :: scalar_fields(:)
    logical :: write_mesh, write_velocity, write_pressure, write_temperature
    integer :: FLD_DATA_SIZE, n_scalar_fields

    if (present(t)) then
       time = real(t, dp)
    else
       time = 0d0
    end if

    nullify(msh)
    nullify(dof)
    nullify(Xh)
    n_scalar_fields = 0
    write_pressure = .false.
    write_velocity = .false.
    write_temperature = .false.

    select type (data)
      type is (fld_file_data_t)
       nelv = data%nelv
       lx = data%lx
       ly = data%ly
       lz = data%lz
       gdim = data%gdim
       glb_nelv = data%glb_nelv
       offset_el = data%offset_el
 
       if (data%x%n .gt. 0) x%ptr => data%x%x
       if (data%y%n .gt. 0) y%ptr => data%y%x
       if (data%z%n .gt. 0) z%ptr => data%z%x
       if (gdim .eq. 2) z%ptr => data%y%x
       if (data%u%n .gt. 0) then
          u%ptr => data%u%x
          write_velocity = .true.
       end if
       if (data%v%n .gt. 0) v%ptr => data%v%x
       if (data%w%n .gt. 0) w%ptr => data%w%x
       if (data%p%n .gt. 0) then
          p%ptr => data%p%x
          write_pressure = .true.
       end if
       if (data%t%n .gt. 0) then
          write_temperature = .true.
          tem%ptr => data%t%x
       end if
       ! If gdim = 2 and Z-velocity component exists,
       ! it is stored in last scalar field 
       if (gdim .eq. 2 .and. data%w%n .gt. 0) then
          n_scalar_fields = data%n_scalars + 1
          allocate(scalar_fields(n_scalar_fields))
          do i = 1, n_scalar_fields -1
             scalar_fields(i)%ptr => data%s(i)%x
          end do
          scalar_fields(n_scalar_fields)%ptr => data%w%x
       else 
          n_scalar_fields = data%n_scalars 
          allocate(scalar_fields(n_scalar_fields+1))
          do i = 1, n_scalar_fields 
             scalar_fields(i)%ptr => data%s(i)%x
          end do
          scalar_fields(n_scalar_fields+1)%ptr => data%w%x
       end if
       ! This is very stupid...
       ! Some compilers cannot handle that these pointers dont point to anything
       ! (although they are not used) this fixes a segfault due to this.
       if (nelv .eq. 0) then
          allocate(tempo(1))
          x%ptr => tempo
          y%ptr => tempo
          z%ptr => tempo
          u%ptr => tempo
          v%ptr => tempo
          w%ptr => tempo
          p%ptr => tempo
          tem%ptr => tempo
       end if

       allocate(idx(nelv))
       do i = 1, nelv
          idx(i) = data%idx(i)
       end do
      type is (field_t)
       p%ptr => data%x(:,1,1,1)
       dof => data%dof
       write_pressure = .true.
       write_velocity = .false.
      type is (field_list_t)
       select case (data%size())
         case (1)
          p%ptr => data%items(1)%ptr%x(:,1,1,1)
          write_pressure = .true.
          write_velocity = .false.
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
          do i = 1, n_scalar_fields
             scalar_fields(i)%ptr => data%items(i+5)%ptr%x(:,1,1,1)
          end do
          write_pressure = .true.
          write_velocity = .true.
          write_temperature = .true.
         case default
          call neko_error('This many fields not supported yet, fld_file')
       end select
       dof => data%dof(1)

      type is (mean_flow_t)
       u%ptr => data%u%mf%x(:,1,1,1)
       v%ptr => data%v%mf%x(:,1,1,1)
       w%ptr => data%w%mf%x(:,1,1,1)
       p%ptr => data%p%mf%x(:,1,1,1)
       dof => data%u%mf%dof
       write_pressure = .true.
       write_velocity = .true.
      type is (mean_sqr_flow_t)
       u%ptr => data%uu%mf%x(:,1,1,1)
       v%ptr => data%vv%mf%x(:,1,1,1)
       w%ptr => data%ww%mf%x(:,1,1,1)
       p%ptr => data%pp%mf%x(:,1,1,1)
       dof => data%pp%mf%dof
       write_pressure = .true.
       write_velocity = .true.
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
       FLD_DATA_SIZE = MPI_DOUBLE_PRECISION_SIZE
    else
       FLD_DATA_SIZE = MPI_REAL_SIZE
    end if
    if (this%dp_precision) then
       allocate(tmp_dp(gdim*n))
    else
       allocate(tmp_sp(gdim*n))
    end if


    !
    ! Create fld header for NEKTON's multifile output
    !

    write_mesh = (this%counter .eq. this%start_counter)
    call MPI_Allreduce(MPI_IN_PLACE, write_mesh, 1, &
         MPI_LOGICAL, MPI_LOR, NEKO_COMM) 
    call MPI_Allreduce(MPI_IN_PLACE, write_velocity, 1, &
         MPI_LOGICAL, MPI_LOR, NEKO_COMM) 
    call MPI_Allreduce(MPI_IN_PLACE, write_pressure, 1, &
         MPI_LOGICAL, MPI_LOR, NEKO_COMM) 
    call MPI_Allreduce(MPI_IN_PLACE, write_temperature, 1, & 
         MPI_LOGICAL, MPI_LOR, NEKO_COMM) 
    call MPI_Allreduce(MPI_IN_PLACE, n_scalar_fields, 1, &
         MPI_INTEGER, MPI_MAX, NEKO_COMM) 
 
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

    !> @todo fix support for single precision output?
    write(hdr, 1) FLD_DATA_SIZE, lx, ly, lz, glb_nelv, glb_nelv,&
         time, this%counter, 1, 1, (rdcode(i),i = 1, 10)
1   format('#std', 1x, i1, 1x, i2, 1x, i2, 1x, i2, 1x, i10, 1x, i10, &
         1x, e20.13, 1x, i9, 1x, i6, 1x, i6, 1x, 10a)

    ! Change to NEKTON's fld file format
    suffix_pos = filename_suffix_pos(this%fname)
    write(id_str, '(a,i5.5)') 'f', this%counter
    fname = trim(this%fname(1:suffix_pos-1))//'0.'//id_str

    call MPI_File_open(NEKO_COMM, trim(fname), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, &
         ierr)

    call MPI_File_write_all(fh, hdr, 132, MPI_CHARACTER, status, ierr)
    mpi_offset = 132 * MPI_CHARACTER_SIZE

    call MPI_File_write_all(fh, test_pattern, 1, MPI_REAL, status, ierr)
    mpi_offset = mpi_offset + MPI_REAL_SIZE

    byte_offset = mpi_offset + &
         int(offset_el, i8) * int(MPI_INTEGER_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, idx, nelv, &
         MPI_INTEGER, status, ierr)
    mpi_offset = mpi_offset + int(glb_nelv, i8) * int(MPI_INTEGER_SIZE, i8)
    deallocate(idx)
    if (write_mesh) then

       byte_offset = mpi_offset + int(offset_el, i8) * &
            (int(gdim*lxyz, i8) * &
            int(FLD_DATA_SIZE, i8))
       call fld_file_write_vector_field(this, fh, byte_offset, &
            x%ptr, y%ptr, z%ptr, &
            n, gdim, lxyz, nelv)
       mpi_offset = mpi_offset + int(glb_nelv, i8) * &
            (int(gdim *lxyz, i8) * &
            int(FLD_DATA_SIZE, i8))
    end if
    if (write_velocity) then
       byte_offset = mpi_offset + int(offset_el, i8) * &
            (int(gdim * (lxyz), i8) * int(FLD_DATA_SIZE, i8))
       call fld_file_write_vector_field(this, fh, byte_offset, &
            u%ptr, v%ptr, w%ptr, n, gdim, lxyz, nelv)

       mpi_offset = mpi_offset + int(glb_nelv, i8) * &
            (int(gdim * (lxyz), i8) * &
            int(FLD_DATA_SIZE, i8))

    end if

    if (write_pressure) then
       byte_offset = mpi_offset + int(offset_el, i8) * &
            (int((lxyz), i8) * int(FLD_DATA_SIZE, i8))
       call fld_file_write_field(this, fh, byte_offset, p%ptr, n)
       mpi_offset = mpi_offset + int(glb_nelv, i8) * &
            (int((lxyz), i8) * int(FLD_DATA_SIZE, i8))
    end if

    if (write_temperature) then
       byte_offset = mpi_offset + int(offset_el, i8) * &
            (int((lxyz), i8) * &
            int(FLD_DATA_SIZE, i8))
       call fld_file_write_field(this, fh, byte_offset, tem%ptr, n)
       mpi_offset = mpi_offset + int(glb_nelv, i8) * &
            (int((lxyz), i8) * &
            int(FLD_DATA_SIZE, i8))
    end if

    temp_offset = mpi_offset

    do i = 1, n_scalar_fields
       ! Without this redundant if statement Cray optimizes this loop to
       ! Oblivion
       if (i .eq. 2) then
          mpi_offset = int(temp_offset, i8) + int(1_i8*glb_nelv, i8) * &
               (int(lxyz, i8) * int(FLD_DATA_SIZE, i8))
       end if
       byte_offset = int(mpi_offset, i8) + int(offset_el, i8) * &
            (int((lxyz), i8) * &
            int(FLD_DATA_SIZE, i8))
       call fld_file_write_field(this, fh, byte_offset, scalar_fields(i)%ptr, n)
       mpi_offset = int(mpi_offset, i8) + int(glb_nelv, i8) * &
            (int(lxyz, i8) * &
            int(FLD_DATA_SIZE, i8))
    end do
    
    if (gdim .eq. 3) then

       !> Include metadata with bounding boxes (Just copying from nek5000)
       if (write_mesh) then
          !The offset is:
          ! mpioff + element_off * 2(min max value)
          !           * 4(single precision) * gdim(dimensions)
          byte_offset = int(mpi_offset, i8) + &
               int(offset_el, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8) * &
               int(gdim, i8)
          call fld_file_write_metadata_vector(this, fh, byte_offset, &
               x%ptr, y%ptr, z%ptr, gdim, lxyz, nelv)
          mpi_offset = int(mpi_offset, i8) + &
               int(glb_nelv, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8) * &
               int(gdim, i8)
       end if

       if (write_velocity) then
          byte_offset = int(mpi_offset, i8) + &
               int(offset_el, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8) * &
               int(gdim, i8)
          call fld_file_write_metadata_vector(this, fh, byte_offset, &
               u%ptr, v%ptr, w%ptr, gdim, lxyz, nelv)
          mpi_offset = int(mpi_offset, i8) + &
               int(glb_nelv, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8) * &
               int(gdim, i8)

       end if

       if (write_pressure) then
          byte_offset = int(mpi_offset, i8) + &
               int(offset_el, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8)
          call fld_file_write_metadata_scalar(this, fh, byte_offset, &
               p%ptr, lxyz, nelv)
          mpi_offset = int(mpi_offset, i8) + &
               int(glb_nelv, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8)

       end if

       if (write_temperature) then
          byte_offset = int(mpi_offset, i8) + &
               int(offset_el, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8)
          call fld_file_write_metadata_scalar(this, fh, byte_offset, &
               tem%ptr, lxyz, nelv)
          mpi_offset = int(mpi_offset, i8) + &
               int(glb_nelv, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8)

       end if



       temp_offset = mpi_offset

       do i = 1, n_scalar_fields
          ! Without this redundant if statement, Cray optimizes this loop to
          ! Oblivion
          if (i .eq. 2) then
             mpi_offset = int(temp_offset, i8) + &
                  int(1_i8*glb_nelv, i8) * &
                  int(2, i8) * &
                  int(MPI_REAL_SIZE, i8)
          end if

          byte_offset = int(mpi_offset, i8) + &
               int(offset_el, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8)
          call fld_file_write_metadata_scalar(this, fh, byte_offset, &
               scalar_fields(i)%ptr, lxyz, nelv)
          mpi_offset = int(mpi_offset, i8) + &
               int(glb_nelv, i8) * &
               int(2, i8) * &
               int(MPI_REAL_SIZE, i8)
       end do
    end if


    call MPI_File_sync(fh, ierr)
    call MPI_File_close(fh, ierr)
    ! Write metadata file
    if (pe_rank .eq. 0) then
       tslash_pos = filename_tslash_pos(this%fname)
       write(start_field, "(I5,A8)") this%start_counter, '.nek5000'
       open(unit = 9, &
            file = trim(this%fname(1:suffix_pos-1)) // &
            trim(adjustl(start_field)), status = 'replace')
       write(9, fmt = '(A,A,A)') 'filetemplate:         ', &
            this%fname(tslash_pos+1:suffix_pos-1), '%01d.f%05d'
       write(9, fmt = '(A,i5)') 'firsttimestep: ', this%start_counter
       write(9, fmt = '(A,i5)') 'numtimesteps: ', &
            (this%counter + 1) - this%start_counter
       write(9, fmt = '(A)') 'type: binary'
       close(9)
    end if

    this%counter = this%counter + 1

    if (allocated(tmp_dp)) deallocate(tmp_dp)
    if (allocated(tmp_sp)) deallocate(tmp_sp)
  end subroutine fld_file_write

  subroutine fld_file_write_metadata_vector(this, fh, byte_offset, x, y, z, &
       gdim, lxyz, nelv)
    class(fld_file_t), intent(inout) :: this
    type(MPI_File), intent(inout) :: fh
    integer, intent(in) :: gdim, lxyz, nelv
    real(kind=rp), intent(in) :: x(lxyz, nelv), y(lxyz, nelv), z(lxyz, nelv)
    integer (kind=MPI_OFFSET_KIND), intent(in) :: byte_offset
    integer :: el, j, ierr, nout
    type(MPI_Status) :: status
    real(kind=sp) :: buffer(2*gdim*nelv)

    j = 1
    do el = 1, nelv
       buffer(j+0) = real(vlmin(x(1, el), lxyz), sp)
       buffer(j+1) = real(vlmax(x(1, el), lxyz), sp)
       buffer(j+2) = real(vlmin(y(1, el), lxyz), sp)
       buffer(j+3) = real(vlmax(y(1, el), lxyz), sp)
       j = j + 4
       if (gdim .eq. 3) then
          buffer(j+0) = real(vlmin(z(1, el), lxyz), sp)
          buffer(j+1) = real(vlmax(z(1, el), lxyz), sp)
          j = j + 2 
       end if
    end do

    ! write out data
    nout = 2*gdim*nelv

    call MPI_File_write_at_all(fh, byte_offset, buffer, nout, &
         MPI_REAL, status, ierr)

  end subroutine fld_file_write_metadata_vector

  subroutine fld_file_write_metadata_scalar(this, fh, byte_offset, x, lxyz, &
       nelv)
    class(fld_file_t), intent(inout) :: this
    type(MPI_File), intent(inout) :: fh
    integer, intent(in) :: lxyz, nelv
    real(kind=rp), intent(in) :: x(lxyz, nelv)
    integer (kind=MPI_OFFSET_KIND), intent(in) :: byte_offset
    integer :: el, j, ierr, nout
    type(MPI_Status) :: status
    real(kind=sp) :: buffer(2*nelv)

    j = 1
    do el = 1, nelv
       buffer(j+0) = real(vlmin(x(1, el), lxyz), sp)
       buffer(j+1) = real(vlmax(x(1, el), lxyz), sp)
       j = j + 2
    end do

    ! write out data
    nout = 2*nelv

    call MPI_File_write_at_all(fh, byte_offset, buffer, nout, &
         MPI_REAL, status, ierr)

  end subroutine fld_file_write_metadata_scalar

  subroutine fld_file_write_field(this, fh, byte_offset, p, n)
    class(fld_file_t), intent(inout) :: this
    type(MPI_File), intent(inout) :: fh
    integer, intent(inout) :: n
    real(kind=rp), intent(inout) :: p(n)
    integer (kind=MPI_OFFSET_KIND), intent(in) :: byte_offset
    integer :: i, ierr
    type(MPI_Status) :: status

    if ( this%dp_precision) then
       do i = 1, n
          tmp_dp(i) = real(p(i), dp)
       end do

       call MPI_File_write_at_all(fh, byte_offset, tmp_dp, n, &
            MPI_DOUBLE_PRECISION, status, ierr)
    else
       do i = 1, n
          tmp_sp(i) = real(p(i), sp)
       end do
       call MPI_File_write_at_all(fh, byte_offset, tmp_sp, n, &
            MPI_REAL, status, ierr)
    end if

  end subroutine fld_file_write_field

  subroutine fld_file_write_vector_field(this, fh, byte_offset, x, y, z, n, &
       gdim, lxyz, nelv)
    class(fld_file_t), intent(inout) :: this
    type(MPI_File), intent(inout) :: fh
    integer, intent(in) :: n, gdim, lxyz, nelv
    real(kind=rp), intent(in) :: x(lxyz, nelv), y(lxyz, nelv), z(lxyz, nelv)
    integer (kind=MPI_OFFSET_KIND), intent(in) :: byte_offset
    integer :: i, el, j, ierr
    type(MPI_Status) :: status

    if (this%dp_precision) then
       i = 1
       do el = 1, nelv
          do j = 1, lxyz
             tmp_dp(i) = real(x(j, el), dp)
             i = i +1
          end do
          do j = 1, lxyz
             tmp_dp(i) = real(y(j, el), dp)
             i = i +1
          end do
          if (gdim .eq. 3) then
             do j = 1, lxyz
                tmp_dp(i) = real(z(j, el), dp)
                i = i +1
             end do
          end if
       end do
       call MPI_File_write_at_all(fh, byte_offset, tmp_dp, gdim*n, &
            MPI_DOUBLE_PRECISION, status, ierr)
    else
       i = 1
       do el = 1, nelv
          do j = 1, lxyz
             tmp_sp(i) = real(x(j, el), sp)
             i = i +1
          end do
          do j = 1, lxyz
             tmp_sp(i) = real(y(j, el), sp)
             i = i +1
          end do
          if (gdim .eq. 3) then
             do j = 1, lxyz
                tmp_sp(i) = real(z(j, el), sp)
                i = i +1
             end do
          end if
       end do
       call MPI_File_write_at_all(fh, byte_offset, tmp_sp, gdim*n, &
            MPI_REAL, status, ierr)
    end if


  end subroutine fld_file_write_vector_field

  !> Load a field from a NEKTON fld file
  subroutine fld_file_read(this, data)
    class(fld_file_t) :: this
    class(*), target, intent(inout) :: data
    character(len= 132) :: hdr
    integer :: ierr, suffix_pos, i, j
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    character(len= 1024) :: fname, meta_fname, string, path
    logical :: meta_file, read_mesh, read_velocity, read_pressure
    logical :: read_temp
    character(len=6) :: id_str
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer :: lx, ly, lz, glb_nelv, counter, lxyz
    integer :: FLD_DATA_SIZE, n_scalars, n
    real(kind=rp) :: time
    real(kind=sp) :: temp
    type(linear_dist_t) :: dist
    real(kind=sp), parameter :: test_pattern = 6.54321
    character :: rdcode(10), temp_str(4)

    select type (data)
      type is (fld_file_data_t)
       call filename_chsuffix(this%fname, meta_fname, 'nek5000')

       inquire(file = trim(meta_fname), exist = meta_file)
       if (meta_file .and. data%meta_nsamples .eq. 0) then
          if (pe_rank .eq. 0) then
             open(unit = 9, file = trim(meta_fname))
             read(9, fmt = '(A)') string
             read(string(14:), fmt = '(A)') string
             string = trim(string)
             data%fld_series_fname = string(:scan(trim(string), '%')-1)
             data%fld_series_fname = adjustl(data%fld_series_fname)
             data%fld_series_fname = trim(data%fld_series_fname)//'0'
             read(9, fmt = '(A)') string
             read(string(scan(string, ':')+1:), *) data%meta_start_counter
             read(9, fmt = '(A)') string
             read(string(scan(string, ':')+1:), *) data%meta_nsamples

             close(9)
             write(*,*) 'Reading meta file for fld series'
             write(*,*) 'Name: ', trim(data%fld_series_fname)
             write(*,*) 'Start counter: ', data%meta_start_counter, &
                  'Nsamples: ', data%meta_nsamples
          end if
          call MPI_Bcast(data%fld_series_fname, 1024, MPI_CHARACTER, 0, &
               NEKO_COMM, ierr)
          call MPI_Bcast(data%meta_start_counter, 1, MPI_INTEGER, 0, &
               NEKO_COMM, ierr)
          call MPI_Bcast(data%meta_nsamples, 1, MPI_INTEGER, 0, &
               NEKO_COMM, ierr)
          if (this%counter .eq. 0) this%counter = data%meta_start_counter
       end if

       if (meta_file) then
          write(id_str, '(a,i5.5)') 'f', this%counter
          path = trim(meta_fname(1:scan(meta_fname, '/', .true. )))
          fname = trim(path)//trim(data%fld_series_fname)//'.'//id_str
          if (this%counter .ge. data%meta_nsamples+data%meta_start_counter) then
             call neko_error('Trying to read more fld files than exist')
          end if
       else
          suffix_pos = filename_suffix_pos(this%fname)
          write(id_str, '(a,i5.5)') 'f', this%counter
          fname = trim(this%fname(1:suffix_pos-1))//'.'//id_str
       end if
       call MPI_File_open(NEKO_COMM, trim(fname), &
            MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

       if (ierr .ne. 0) call neko_error("Could not read "//trim(fname))

       call MPI_File_read_all(fh, hdr, 132, MPI_CHARACTER, status, ierr)
       ! This read can prorbably be done wihtout the temp variables,
       !      temp_str, i, j

       read(hdr, 1) temp_str, FLD_DATA_SIZE, lx, ly, lz, glb_nelv, glb_nelv, &
            time, counter, i, j, (rdcode(i), i = 1, 10)
1      format(4a, 1x, i1, 1x, i2, 1x, i2, 1x, i2, 1x, i10, 1x, i10, &
            1x, e20.13, 1x, i9, 1x, i6, 1x, i6, 1x, 10a)
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


       if (FLD_DATA_SIZE .eq. MPI_DOUBLE_PRECISION_SIZE) then
          this%dp_precision = .true.
       else
          this%dp_precision = .false.
       end if
       if (this%dp_precision) then
          allocate(tmp_dp(data%gdim*n))
       else
          allocate(tmp_sp(data%gdim*n))
       end if


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

       mpi_offset = 132 * MPI_CHARACTER_SIZE
       call MPI_File_read_at_all(fh, mpi_offset, temp, 1, &
            MPI_REAL, status, ierr)
       if (temp .ne. test_pattern) then
          call neko_error('Incorrect format for fld file, &
               &test pattern does not match.')
       end if
       mpi_offset = mpi_offset + MPI_REAL_SIZE


       if (allocated(data%idx)) then
          if (size(data%idx) .ne. data%nelv) then
             deallocate(data%idx)
             allocate(data%idx(data%nelv))
          end if
       else
          allocate(data%idx(data%nelv))
       end if

       byte_offset = mpi_offset + &
            int(data%offset_el, i8) * int(MPI_INTEGER_SIZE, i8)

       call MPI_File_read_at_all(fh, byte_offset, data%idx, data%nelv, &
            MPI_INTEGER, status, ierr)

       mpi_offset = mpi_offset + &
            int(data%glb_nelv, i8) * int(MPI_INTEGER_SIZE, i8)

       if (read_mesh) then
          byte_offset = mpi_offset + int(data%offset_el, i8) * &
               (int(data%gdim*lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
          call fld_file_read_vector_field(this, fh, byte_offset, &
               data%x, data%y, data%z, data)
          mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
               (int(data%gdim *lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
       end if

       if (read_velocity) then
          byte_offset = mpi_offset + int(data%offset_el, i8) * &
               (int(data%gdim*lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
          call fld_file_read_vector_field(this, fh, byte_offset, &
               data%u, data%v, data%w, data)
          mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
               (int(data%gdim *lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
       end if

       if (read_pressure) then
          byte_offset = mpi_offset + int(data%offset_el, i8) * &
               (int(lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
          call fld_file_read_field(this, fh, byte_offset, data%p, data)
          mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
               (int(lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
       end if

       if (read_temp) then
          byte_offset = mpi_offset + int(data%offset_el, i8) * &
               (int(lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
          call fld_file_read_field(this, fh, byte_offset, data%t, data)
          mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
               (int(lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
       end if

       do i = 1, n_scalars
          byte_offset = mpi_offset + int(data%offset_el, i8) * &
               (int(lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
          call fld_file_read_field(this, fh, byte_offset, data%s(i), data)
          mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
               (int(lxyz, i8) * &
               int(FLD_DATA_SIZE, i8))
       end do

       this%counter = this%counter + 1

       if (allocated(tmp_dp)) deallocate(tmp_dp)
       if (allocated(tmp_sp)) deallocate(tmp_sp)
      class default
       call neko_error('Currently we only read into fld_file_data_t, &
            &please use that data structure instead.')
    end select

  end subroutine fld_file_read

  subroutine fld_file_read_field(this, fh, byte_offset, x, fld_data)
    class(fld_file_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x
    type(fld_file_data_t) :: fld_data
    integer(kind=MPI_OFFSET_KIND) :: byte_offset
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    integer :: n, ierr, lxyz, i

    n = x%n
    lxyz = fld_data%lx*fld_data%ly*fld_data%lz

    if (this%dp_precision) then
       call MPI_File_read_at_all(fh, byte_offset, tmp_dp, n, &
            MPI_DOUBLE_PRECISION, status, ierr)
    else
       call MPI_File_read_at_all(fh, byte_offset, tmp_sp, n, &
            MPI_REAL, status, ierr)
    end if

    if (this%dp_precision) then
       do i = 1, n
          x%x(i) = tmp_dp(i)
       end do
    else
       do i = 1, n
          x%x(i) = tmp_sp(i)
       end do
    end if


  end subroutine fld_file_read_field


  subroutine fld_file_read_vector_field(this, fh, byte_offset, &
       x, y, z, fld_data)
    class(fld_file_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x, y, z
    type(fld_file_data_t) :: fld_data
    integer(kind=MPI_OFFSET_KIND) :: byte_offset
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    integer :: n, ierr, lxyz, i, j, e, nd

    n = x%n
    nd = n*fld_data%gdim
    lxyz = fld_data%lx*fld_data%ly*fld_data%lz

    if (this%dp_precision) then
       call MPI_File_read_at_all(fh, byte_offset, tmp_dp, nd, &
            MPI_DOUBLE_PRECISION, status, ierr)
    else
       call MPI_File_read_at_all(fh, byte_offset, tmp_sp, nd, &
            MPI_REAL, status, ierr)
    end if


    if (this%dp_precision) then
       i = 1
       do e = 1, fld_data%nelv
          do j = 1, lxyz
             x%x((e-1)*lxyz+j) = tmp_dp(i)
             i = i +1
          end do
          do j = 1, lxyz
             y%x((e-1)*lxyz+j) = tmp_dp(i)
             i = i +1
          end do
          if (fld_data%gdim .eq. 3) then
             do j = 1, lxyz
                z%x((e-1)*lxyz+j) = tmp_dp(i)
                i = i +1
             end do
          end if
       end do
    else
       i = 1
       do e = 1, fld_data%nelv
          do j = 1, lxyz
             x%x((e-1)*lxyz+j) = tmp_sp(i)
             i = i +1
          end do
          do j = 1, lxyz
             y%x((e-1)*lxyz+j) = tmp_sp(i)
             i = i +1
          end do
          if (fld_data%gdim .eq. 3) then
             do j = 1, lxyz
                z%x((e-1)*lxyz+j) = tmp_sp(i)
                i = i +1
             end do
          end if
       end do
    end if

  end subroutine fld_file_read_vector_field

  subroutine fld_file_set_precision(this, precision)
    class(fld_file_t) :: this
    integer, intent(in) :: precision

    if (precision .eq. dp) then
       this%dp_precision = .true.
    else if (precision .eq. sp) then
       this%dp_precision = .false.
    else
       call neko_error('Invalid precision')
    end if

  end subroutine fld_file_set_precision


end module fld_file
