! Copyright (c) 2020-2022, The Neko Authors
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
  use generic_file
  use field
  use dofmap
  use fluid_method
  use scalar
  use fld_file_data
  use mean_flow
  use mean_sqr_flow
  use mesh
  use utils
  use comm
  use mpi_types
  use mpi_f08    
  implicit none
  


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
    type(field_t), pointer :: u, v, w, p
    type(mesh_t), pointer :: msh
    type(space_t), pointer :: Xh
    type(dofmap_t), pointer :: dof
    real(kind=dp) :: time
    character(len=132) :: hdr
    character :: rdcode(10)
    character(len=6) :: id_str
    character(len=1024) :: fname
    character(len=1024) :: start_field
    integer :: i, ierr, n, j,k,l,el, suffix_pos,tslash_pos
    integer, allocatable :: idx(:)
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    real(kind=dp), allocatable :: tmp_dp(:)
    real(kind=sp), allocatable :: tmp_sp(:)
    real(kind=sp), parameter :: test_pattern = 6.54321
    logical :: write_mesh, write_velocity, write_pressure
    integer :: FLD_DATA_SIZE

    if (present(t)) then
       time = real(t,dp)
    else
       time = 0d0
    end if
    
    select type(data)
    type is (field_t)
       p => data
       msh => p%msh
       Xh => p%Xh
       dof => p%dof
       write_pressure = .true.
       write_velocity = .false.
    class is (fluid_scheme_t)
       u => data%u
       v => data%v
       w => data%w
       p => data%p
       msh => p%msh              
       Xh => p%Xh
       dof => p%dof
       write_pressure = .true.
       write_velocity = .true.
    type is (mean_flow_t)
       u => data%u%mf
       v => data%v%mf
       w => data%w%mf
       p => data%p%mf
       msh => p%msh              
       Xh => p%Xh
       dof => p%dof
       write_pressure = .true.
       write_velocity = .true.
    type is (mean_sqr_flow_t)
       u => data%uu%mf
       v => data%vv%mf
       w => data%ww%mf
       p => data%pp%mf
       msh => p%msh              
       Xh => p%Xh
       dof => p%dof
       write_pressure = .true.
       write_velocity = .true.
    class is (scalar_scheme_t)
       u => data%u
       v => data%v
       w => data%w
       p => data%s
       msh => p%msh
       Xh => p%Xh
       dof => p%dof
       write_pressure = .true.
       write_velocity = .false.
    class default
       call neko_error('Invalid data')
    end select

    if (this%dp_precision) then
       FLD_DATA_SIZE = MPI_DOUBLE_PRECISION_SIZE
    else
       FLD_DATA_SIZE = MPI_REAL_SIZE
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
    end if


    !> @todo fix support for single precision output?
    write(hdr, 1) FLD_DATA_SIZE, Xh%lx, Xh%ly, Xh%lz,msh%glb_nelv,msh%glb_nelv,&
         time, this%counter, 1, 1, (rdcode(i),i=1,10)
1   format('#std',1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,&
         1x,i9,1x,i6,1x,i6,1x,10a)

    ! Store global idx of each element
    allocate(idx(msh%nelv))
    do i = 1, msh%nelv
       idx(i) = msh%elements(i)%e%id()
    end do

    ! Change to NEKTON's fld file format
    suffix_pos = filename_suffix_pos(this%fname)
    write(id_str, '(a,i5.5)') 'f', this%counter
    fname = trim(this%fname(1:suffix_pos-1))//'0.'//id_str
    
    call MPI_File_open(NEKO_COMM, trim(fname), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)

    call MPI_File_write_all(fh, hdr, 132, MPI_CHARACTER,  status, ierr)
    mpi_offset = 132 * MPI_CHARACTER_SIZE

    call MPI_File_write_all(fh, test_pattern, 1, MPI_REAL, status, ierr)
    mpi_offset = mpi_offset  + MPI_REAL_SIZE

    byte_offset = mpi_offset + &
         int(msh%offset_el, i8) * int(MPI_INTEGER_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, idx, msh%nelv, &
         MPI_INTEGER, status, ierr)
    mpi_offset = mpi_offset + int(msh%glb_nelv, i8) * int(MPI_INTEGER_SIZE, i8)

    deallocate(idx)
    
    n = msh%gdim*(Xh%lx * Xh%ly * Xh%lz * msh%nelv)

    if (this%dp_precision) then
       allocate(tmp_dp(n))
    else
       allocate(tmp_sp(n))
    end if
       
    if (write_mesh) then

       byte_offset = mpi_offset + int(msh%offset_el, i8) * &
            (int(msh%gdim*Xh%lxyz, i8) * &
            int(FLD_DATA_SIZE, i8))
       
       if (this%dp_precision) then
          i = 1
          do el = 1, msh%nelv
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_dp(i) = real(dof%x(j,k,l,el),dp)
                      i = i +1
                   end do
                end do
             end do
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_dp(i) = real(dof%y(j,k,l,el),dp)
                      i = i +1
                   end do
                end do
             end do
             if (msh%gdim .eq. 3) then
                do l = 1, Xh%lz
                   do k = 1, Xh%ly
                      do j = 1, Xh%lx
                         tmp_dp(i) = real(dof%z(j,k,l,el),dp)
                         i = i +1
                      end do
                   end do
                end do
             end if
          end do
          
          call MPI_File_write_at_all(fh, byte_offset, tmp_dp, n, &
               MPI_DOUBLE_PRECISION, status, ierr)
          
       else
          i = 1
          do el = 1, msh%nelv
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_sp(i) = real(dof%x(j,k,l,el),sp)
                      i = i +1
                   end do
                end do
             end do
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_sp(i) = real(dof%y(j,k,l,el),sp)
                      i = i +1
                   end do
                end do
             end do
             if (msh%gdim .eq. 3) then
                do l = 1, Xh%lz
                   do k = 1, Xh%ly
                      do j = 1, Xh%lx
                         tmp_sp(i) = real(dof%z(j,k,l,el),sp)
                         i = i +1
                      end do
                   end do
                end do
             end if
          end do
          
          call MPI_File_write_at_all(fh, byte_offset, tmp_sp, n, &
               MPI_REAL, status, ierr)
          
       end if

       mpi_offset = mpi_offset + int(msh%glb_nelv, i8) * &
            (int(msh%gdim *Xh%lxyz, i8) * & 
            int(FLD_DATA_SIZE, i8))
    end if

    if (write_velocity) then
       byte_offset = mpi_offset + int(msh%offset_el, i8) * &
            (int(msh%gdim * (Xh%lx * Xh%ly * Xh%lz), i8) * &
            int(FLD_DATA_SIZE, i8))

       if (this%dp_precision) then
          i = 1
          do el = 1, msh%nelv
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_dp(i) = real(u%x(j,k,l,el),dp)
                      i = i +1
                   end do
                end do
             end do
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_dp(i) = real(v%x(j,k,l,el),dp)
                      i = i +1
                   end do
                end do
             end do
             if (msh%gdim .eq. 3) then
                do l = 1, Xh%lz
                   do k = 1, Xh%ly
                      do j = 1, Xh%lx
                         tmp_dp(i) = real(w%x(j,k,l,el),dp)
                         i = i +1
                      end do
                   end do
                end do
             end if
          end do
          
          call MPI_File_write_at_all(fh, byte_offset, tmp_dp, n, &
               MPI_DOUBLE_PRECISION, status, ierr)
          
       else
          i = 1
          do el = 1, msh%nelv
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_sp(i) = real(u%x(j,k,l,el),sp)
                      i = i +1
                   end do
                end do
             end do
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_sp(i) = real(v%x(j,k,l,el),sp)
                      i = i +1
                   end do
                end do
             end do
             if (msh%gdim .eq. 3) then
                do l = 1, Xh%lz
                   do k = 1, Xh%ly
                      do j = 1, Xh%lx
                         tmp_sp(i) = real(w%x(j,k,l,el),sp)
                         i = i +1
                      end do
                   end do
                end do
             end if
          end do

          call MPI_File_write_at_all(fh, byte_offset, tmp_sp, n, &
               MPI_REAL, status, ierr)

       end if
       
       mpi_offset = mpi_offset + int(msh%glb_nelv, i8) * &
            (int(msh%gdim * (Xh%lx * Xh%ly * Xh%lz), i8) * &
            int(FLD_DATA_SIZE, i8))
       
    end if
 
    if (write_pressure) then
       byte_offset = mpi_offset + int(msh%offset_el, i8) * &
            (int((Xh%lx * Xh%ly * Xh%lz), i8) * &
            int(FLD_DATA_SIZE, i8))
      
       if (.not. this%dp_precision) then

          i = 1
          do el = 1, msh%nelv
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_sp(i) = real(p%x(j,k,l,el),sp)
                      i = i + 1
                   end do
                end do
             end do
          end do
          
          call MPI_File_write_at_all(fh, byte_offset, tmp_sp, n/msh%gdim, &
               MPI_REAL, status, ierr)
       else
          i = 1
          do el = 1, msh%nelv
             do l = 1, Xh%lz
                do k = 1, Xh%ly
                   do j = 1, Xh%lx
                      tmp_dp(i) = real(p%x(j,k,l,el),dp)
                      i = i + 1
                   end do
                end do
             end do
          end do
          call MPI_File_write_at_all(fh, byte_offset, tmp_dp, n/msh%gdim, &
               MPI_DOUBLE_PRECISION, status, ierr)
       end if
    end if
    
    if (allocated(tmp_dp)) then
       deallocate(tmp_dp)
    end if

    if (allocated(tmp_sp)) then
       deallocate(tmp_sp)
    end if
    
    call MPI_File_sync(fh, ierr)
    call MPI_File_close(fh, ierr)

    ! Write metadata file 
    if (pe_rank .eq. 0) then
       tslash_pos = filename_tslash_pos(this%fname)
       write(start_field,"(I5,A8)") this%start_counter,'.nek5000'
       open(unit=9, file=trim(this%fname(1:suffix_pos-1))//trim(adjustl(start_field)), &
            status='replace')
       write(9, fmt='(A,A,A)') 'filetemplate:         ', &
            this%fname(tslash_pos+1:suffix_pos-1),'%01d.f%05d'
       write(9, fmt='(A,i5)') 'firsttimestep: ', this%start_counter
       write(9, fmt='(A,i5)') 'numtimesteps: ', (this%counter + 1)-this%start_counter
       write(9, fmt='(A)') 'type: binary'
       close(9)
    end if

    this%counter = this%counter + 1
    
  end subroutine fld_file_write
  
  !> Load a field from a NEKTON fld file
  subroutine fld_file_read(this, data)
    class(fld_file_t) :: this
    class(*), target, intent(inout) :: data
    character(len=132) :: hdr
    integer :: ierr, suffix_pos, i, j
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    character(len=1024) :: fname, meta_fname, string
    logical :: meta_file, read_mesh, read_velocity, read_pressure
    logical :: read_temp
    character(len=6) :: id_str
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer :: lx, ly, lz, glb_nelv, counter, lxyz
    integer :: FLD_DATA_SIZE, n_scalars, n, nd
    real(kind=rp) ::  time 
    real(kind=sp) :: temp
    real(kind=sp), parameter :: test_pattern = 6.54321
    character :: rdcode(10),temp_str(4)
    select type(data)
    type is (fld_file_data_t)
      call filename_chsuffix(this%fname, meta_fname,'nek5000')

      inquire(file=meta_fname, exist=meta_file)
      if (meta_file .and. data%meta_nsamples .eq. 0) then
         if (pe_rank .eq. 0) then
            open(unit=9, file=meta_fname)
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
         end if
         call MPI_Bcast(data%fld_series_fname, 1024, MPI_CHARACTER, 0, NEKO_COMM, ierr)
         call MPI_Bcast(data%meta_start_counter, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
         call MPI_Bcast(data%meta_nsamples, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
         this%counter = data%meta_start_counter
      end if

      if (meta_file) then
         write(id_str, '(a,i5.5)') 'f', this%counter
         fname = trim(data%fld_series_fname)//'.'//id_str
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
      call MPI_File_read_all(fh, hdr, 132, MPI_CHARACTER,  status, ierr)
      !This read can prorbably be done wihtout the temp variables, temp_str, i, j       
      
      read(hdr, 1) temp_str,FLD_DATA_SIZE, lx, ly, lz, glb_nelv, glb_nelv,&
          time, counter, i, j, (rdcode(i),i=1,10)
1     format(4a,1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,&
         1x,i9,1x,i6,1x,i6,1x,10a)
       
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
      

      if (this%dp_precision) then
         FLD_DATA_SIZE = MPI_DOUBLE_PRECISION_SIZE
      else
         FLD_DATA_SIZE = MPI_REAL_SIZE
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
         read(rdcode(i:),*) n_scalars
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
            allocate(data%idx(data%nelv))
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
         call neko_error('Incorrect format for fld file, test pattern does not match.')
      end if
      mpi_offset = mpi_offset  + MPI_REAL_SIZE


      if (allocated(data%idx)) then
         if (size(data%idx) .ne. data%nelv) then
            deallocate(data%idx)
            allocate(data%idx(data%nelv))
         end if
      else
         allocate(data%idx(data%nelv))
      end if
      
      byte_offset = mpi_offset + &
           int(data%offset_nelv, i8) * int(MPI_INTEGER_SIZE, i8)

      call MPI_File_read_at_all(fh, byte_offset, data%idx, data%nelv, &
           MPI_INTEGER, status, ierr)
      
      mpi_offset = mpi_offset + int(data%glb_nelv, i8) * int(MPI_INTEGER_SIZE, i8)
      
      if (read_mesh) then
         byte_offset = mpi_offset + int(data%offset_nelv, i8) * &
             (int(data%gdim*lxyz, i8) * &
              int(FLD_DATA_SIZE, i8))
         call fld_file_read_3dfields(this, fh, byte_offset, data%x, data%y, data%z, data)
         mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
             (int(data%gdim *lxyz, i8) * & 
              int(FLD_DATA_SIZE, i8))
      end if

      if (read_velocity) then
         byte_offset = mpi_offset + int(data%offset_nelv, i8) * &
             (int(data%gdim*lxyz, i8) * &
              int(FLD_DATA_SIZE, i8))
         call fld_file_read_3dfields(this, fh, byte_offset, data%u, data%v, data%w, data)
         mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
             (int(data%gdim *lxyz, i8) * & 
              int(FLD_DATA_SIZE, i8))
      end if

      if (read_pressure) then
         byte_offset = mpi_offset + int(data%offset_nelv, i8) * &
             (int(lxyz, i8) * &
              int(FLD_DATA_SIZE, i8))
         call fld_file_read_field(this, fh, byte_offset, data%p, data)
         mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
             (int(lxyz, i8) * & 
              int(FLD_DATA_SIZE, i8))
      end if

      if (read_temp) then
         byte_offset = mpi_offset + int(data%offset_nelv, i8) * &
             (int(lxyz, i8) * &
              int(FLD_DATA_SIZE, i8))
         call fld_file_read_field(this, fh, byte_offset, data%t, data)
         mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
             (int(lxyz, i8) * & 
              int(FLD_DATA_SIZE, i8))
      end if

      do i = 1, n_scalars
         byte_offset = mpi_offset + int(data%offset_nelv, i8) * &
             (int(lxyz, i8) * &
              int(FLD_DATA_SIZE, i8))
         call fld_file_read_field(this, fh, byte_offset, data%s(i), data)
         mpi_offset = mpi_offset + int(data%glb_nelv, i8) * &
             (int(lxyz, i8) * & 
              int(FLD_DATA_SIZE, i8))
      end do

      this%counter = this%counter + 1

    class default 
       call neko_error('Currently we only read into fld_file_data_t, please use that data structure instead.')
    end select

  end subroutine fld_file_read

  subroutine fld_file_read_field(this, fh, byte_offset, x, fld_data)
    class(fld_file_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x
    type(fld_file_data_t) :: fld_data
    integer(kind=MPI_OFFSET_KIND) :: byte_offset
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    real(kind=dp), allocatable :: tmp_dp(:)
    real(kind=sp), allocatable :: tmp_sp(:)
    integer :: n, ierr, lxyz, i, j, e

    n = x%n
    lxyz = fld_data%lx*fld_data%ly*fld_data%lz

    if (this%dp_precision) then
       allocate(tmp_dp(n))
       call MPI_File_read_at_all(fh, byte_offset, tmp_dp, n, &
            MPI_DOUBLE_PRECISION, status, ierr)
    else
       allocate(tmp_sp(n))
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


  subroutine fld_file_read_3dfields(this, fh, byte_offset, x, y, z, fld_data)
    class(fld_file_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x, y, z
    type(fld_file_data_t) :: fld_data
    integer(kind=MPI_OFFSET_KIND) :: byte_offset
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    real(kind=dp), allocatable :: tmp_dp(:)
    real(kind=sp), allocatable :: tmp_sp(:)
    integer :: n, ierr, lxyz, i, j, e, nd

    n = x%n
    nd = n*fld_data%gdim
    lxyz = fld_data%lx*fld_data%ly*fld_data%lz

    if (this%dp_precision) then
       allocate(tmp_dp(fld_data%gdim*n))
       call MPI_File_read_at_all(fh, byte_offset, tmp_dp, nd, &
            MPI_DOUBLE_PRECISION, status, ierr)
    else
       allocate(tmp_sp(fld_data%gdim*n))
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
  end subroutine fld_file_read_3dfields

  subroutine fld_file_set_precision(this, precision)
    class(fld_file_t) :: this
    integer, intent(inout) :: precision

    if (precision .eq. dp) then
       this%dp_precision = .true.
    else if (precision .eq. sp) then
       this%dp_precision = .false.
    else
       call neko_error('Invalid precision')
    end if
    
  end subroutine fld_file_set_precision

  
end module fld_file
