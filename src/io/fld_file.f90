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
  use mean_flow
  use mean_sqr_flow
  use mesh
  use utils
  use comm
  use mpi_types
  use mpi_f08    
  implicit none
  private


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

    call neko_error('Not implemented yet!')
    
  end subroutine fld_file_read


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
