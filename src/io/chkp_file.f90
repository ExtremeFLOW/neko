! Copyright (c) 2021-2022, The Neko Authors
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
!> Neko checkpoint file format
!! @details this module defines interface to read/write Neko's ceckpoint files
module chkp_file
  use generic_file
  use field_series
  use checkpoint    
  use num_types
  use field
  use utils
  use mesh
  use interpolation
  use mpi_types
  use mpi_f08
  implicit none
  private

  !> Interface for Neko checkpoint files
  type, public, extends(generic_file_t) :: chkp_file_t
     type(space_t) :: chkp_Xh
     type(space_t), pointer :: sim_Xh
     type(interpolator_t) :: interp
   contains
     procedure :: read => chkp_file_read
     procedure :: read_field => chkp_read_field
     procedure :: write => chkp_file_write
  end type chkp_file_t

contains
  
  !> Write a Neko checkpoint
  subroutine chkp_file_write(this, data, t)
    class(chkp_file_t), intent(inout) :: this
    class(*), target, intent(in) :: data
    real(kind=rp), intent(in), optional :: t
    real(kind=dp) :: time
    character(len=5) :: id_str
    character(len=1024) :: fname
    integer :: ierr, suffix_pos, have_lag
    type(field_t), pointer :: u, v, w, p
    type(field_series_t), pointer :: ulag => null()
    type(field_series_t), pointer :: vlag => null()
    type(field_series_t), pointer :: wlag => null()
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer(kind=i8) :: n_glb_dofs, dof_offset
    logical write_lag
    integer :: i

    if (present(t)) then
       time = real(t,dp)
    else
       time = 0d0
    end if
    
    select type(data)
    type is (chkp_t)

       if ( .not. associated(data%u) .or. &
            .not. associated(data%v) .or. &
            .not. associated(data%w) .or. &
            .not. associated(data%p) ) then
          call neko_error('Checkpoint not initialized')
       end if
    
       u => data%u
       v => data%v
       w => data%w
       p => data%p
       msh => u%msh

       if (associated(data%ulag)) then       
          ulag => data%ulag
          vlag => data%vlag
          wlag => data%wlag
          write_lag = .true.
          have_lag = 1
       else
          write_lag = .false.
          have_lag = 0
       end if
       
    class default
       call neko_error('Invalid data')
    end select


    
    suffix_pos = filename_suffix_pos(this%fname)
    write(id_str, '(i5.5)') this%counter
    fname = trim(this%fname(1:suffix_pos-1))//id_str//'.chkp'


    dof_offset = int(msh%offset_el, i8) * int(u%Xh%lx * u%Xh%ly * u%Xh%lz, i8)
    n_glb_dofs = int(u%Xh%lx * u%Xh%ly * u%Xh%lz, i8) * int(msh%glb_nelv, i8)
    
    call MPI_File_open(NEKO_COMM, trim(fname), &
         MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)
    call MPI_File_write_all(fh, msh%glb_nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, msh%gdim, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, u%Xh%lx, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, have_lag, 1, MPI_INTEGER, status, ierr)
    call MPI_File_write_all(fh, time, 1, MPI_DOUBLE_PRECISION, status, ierr)
    
    
    !
    ! Dump mandatory checkpoint data
    !
    
    byte_offset = 4 * MPI_INTEGER_SIZE + MPI_DOUBLE_PRECISION_SIZE + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, u%x, u%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = 4 * MPI_INTEGER_SIZE + MPI_DOUBLE_PRECISION_SIZE + &
         n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    
    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, v%x, v%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, w%x, w%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call MPI_File_write_at_all(fh, byte_offset, p%x, p%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    !
    ! Dump optional payload
    !

    if (write_lag) then

       do i = 1, ulag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call MPI_File_write_at_all(fh, byte_offset, ulag%lf(i)%x, &
               ulag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, vlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call MPI_File_write_at_all(fh, byte_offset, vlag%lf(i)%x, &
               vlag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, wlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call MPI_File_write_at_all(fh, byte_offset, wlag%lf(i)%x, &
               wlag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
              
    end if
    
    call MPI_File_close(fh, ierr)

    this%counter = this%counter + 1
    
  end subroutine chkp_file_write
  
  !> Load a checkpoint from file
  subroutine chkp_file_read(this, data)
    class(chkp_file_t) :: this
    class(*), target, intent(inout) :: data
    type(chkp_t), pointer :: chkp
    character(len=5) :: id_str
    character(len=1024) :: fname
    integer :: ierr, suffix_pos
    type(field_t), pointer :: u, v, w, p
    type(field_series_t), pointer :: ulag => null()
    type(field_series_t), pointer :: vlag => null()
    type(field_series_t), pointer :: wlag => null()
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer(kind=i8) :: n_glb_dofs, dof_offset
    integer :: glb_nelv, gdim, lx, have_lag, nel
    logical read_lag
    integer :: i
    
    select type(data)
    type is (chkp_t)       

       if ( .not. associated(data%u) .or. &
            .not. associated(data%v) .or. &
            .not. associated(data%w) .or. &
            .not. associated(data%p) ) then
          call neko_error('Checkpoint not initialized')
       end if
    
       u => data%u
       v => data%v
       w => data%w
       p => data%p
       msh => u%msh

       if (associated(data%ulag)) then       
          ulag => data%ulag
          vlag => data%vlag
          wlag => data%wlag
          read_lag = .true.
       else
          read_lag = .false.
       end if

       chkp => data
       
    class default
       call neko_error('Invalid data')
    end select


    
    call MPI_File_open(NEKO_COMM, trim(this%fname), &
         MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_File_read_all(fh, glb_nelv, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, gdim, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, lx, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, have_lag, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, chkp%t, 1, MPI_DOUBLE_PRECISION, status, ierr)

    if ( ( glb_nelv .ne. msh%glb_nelv ) .or. &
         ( gdim .ne. msh%gdim) .or. &
         ( (have_lag .eq. 1) .and. (.not. read_lag) ) ) then
       call neko_error('Checkpoint does not match case')
    end if
    nel = msh%nelv
    this%sim_Xh => u%Xh
    if (gdim .eq. 3) then
       call space_init(this%chkp_Xh, GLL, lx, lx, lx)
    else
       call space_init(this%chkp_Xh, GLL, lx, lx)
    end if
    call this%interp%init(this%sim_Xh, this%chkp_Xh) 
    dof_offset = int(msh%offset_el, i8) * int(this%chkp_Xh%lxyz, i8)
    n_glb_dofs = int(this%chkp_Xh%lxyz, i8) * int(msh%glb_nelv, i8)    
    
    !
    ! Read mandatory checkpoint data
    !
    
    byte_offset = 4 * MPI_INTEGER_SIZE + MPI_DOUBLE_PRECISION_SIZE + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call this%read_field(fh, byte_offset, u%x, nel)
    mpi_offset = 4 * MPI_INTEGER_SIZE + MPI_DOUBLE_PRECISION_SIZE + &
         n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    
    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call this%read_field(fh, byte_offset, v%x, nel)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call this%read_field(fh, byte_offset, w%x, nel)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, i8)
    call this%read_field(fh, byte_offset, p%x, nel)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)

    !
    ! Read optional payload
    !

    if (read_lag) then

       do i = 1, ulag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, ulag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, vlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, vlag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
       
       do i = 1, wlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          call this%read_field(fh, byte_offset, wlag%lf(i)%x, nel)
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
       
    end if
    
    call MPI_File_close(fh, ierr)      
    
  end subroutine chkp_file_read

  subroutine chkp_read_field(this, fh, byte_offset, x, nel)
    class(chkp_file_t) :: this
    type(MPI_File) :: fh
    integer(kind=MPI_OFFSET_KIND) :: byte_offset
    integer :: nel
    real(kind=rp) :: x(this%sim_Xh%lxyz, nel)
    real(kind=rp), allocatable :: read_array(:)
    integer :: nel_stride, frac_space
    type(MPI_Status) :: status
    integer :: ierr

    allocate(read_array(this%chkp_Xh%lxyz*nel)) 
    call MPI_File_read_at_all(fh, byte_offset, read_array, &
               nel*this%chkp_Xh%lxyz, MPI_REAL_PRECISION, status, ierr)
    call this%interp%map_host(x, read_array, nel, this%sim_Xh)
    deallocate(read_array)
  end subroutine chkp_read_field
  
end module chkp_file
