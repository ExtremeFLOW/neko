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
  use dofmap, only: dofmap_t
  use utils
  use space
  use mesh
  use math
  use interpolation
  use mpi_types
  use comm
  use global_interpolation
  implicit none
  private

  !> Interface for Neko checkpoint files
  type, public, extends(generic_file_t) :: chkp_file_t
     type(space_t) :: chkp_Xh !< Function space in the loaded checkpoint file
     type(space_t), pointer :: sim_Xh !< Function space used in the simulation
     type(interpolator_t) :: space_interp !< Interpolation when only changing lx
     type(global_interpolation_t) :: global_interp !< Interpolation for different meshes
     logical :: mesh2mesh !< Flag if previous mesh difers from current.
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
    integer :: ierr, suffix_pos, optional_fields
    type(field_t), pointer :: u, v, w, p, s
    type(field_series_t), pointer :: ulag => null()
    type(field_series_t), pointer :: vlag => null()
    type(field_series_t), pointer :: wlag => null()
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer(kind=i8) :: n_glb_dofs, dof_offset
    logical :: write_lag, write_scalar
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
       
       optional_fields = 0

       if (associated(data%ulag)) then       
          ulag => data%ulag
          vlag => data%vlag
          wlag => data%wlag
          write_lag = .true.
          optional_fields = optional_fields + 1
       else
          write_lag = .false.
       end if
 
       if (associated(data%s)) then       
          s => data%s
          write_scalar = .true.
          optional_fields = optional_fields + 2
       else
          write_scalar = .false.
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
    call MPI_File_write_all(fh, optional_fields, 1, MPI_INTEGER, status, ierr)
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
          ! We should not need this extra associate block, ant it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => ulag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 ulag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, vlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, ant it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => vlag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 vlag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do

       do i = 1, wlag%size()
          byte_offset = mpi_offset + &
               dof_offset * int(MPI_REAL_PREC_SIZE, i8)
          ! We should not need this extra associate block, ant it works
          ! great without it for GNU, Intel, NEC and Cray, but throws an
          ! ICE with NAG.
          associate (x => wlag%lf(i)%x)
            call MPI_File_write_at_all(fh, byte_offset, x, &
                 wlag%lf(i)%dof%size(), MPI_REAL_PRECISION, status, ierr)
          end associate
          mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
       end do
              
    end if

    if (write_scalar) then 
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call MPI_File_write_at_all(fh, byte_offset, s%x, p%dof%size(), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
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
    type(field_t), pointer :: u, v, w, p, s
    type(field_series_t), pointer :: ulag => null()
    type(field_series_t), pointer :: vlag => null()
    type(field_series_t), pointer :: wlag => null()
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    real(kind=rp), allocatable :: x_coord(:,:,:,:)
    real(kind=rp), allocatable :: y_coord(:,:,:,:)
    real(kind=rp), allocatable :: z_coord(:,:,:,:)
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer(kind=i8) :: n_glb_dofs, dof_offset
    integer :: glb_nelv, gdim, lx, have_lag, have_scalar, nel, optional_fields
    logical :: read_lag, read_scalar
    real(kind=rp) :: tol
    real(kind=rp) :: center_x, center_y, center_z
    integer :: i, e
    type(dofmap_t) :: dof

    
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
       !> If checkpoint used another mesh
       if (allocated(data%previous_mesh%elements)) then
          msh => data%previous_mesh
          this%mesh2mesh = .true.
          tol = data%mesh2mesh_tol
       else !< The checkpoint was written on the same mesh
          msh => u%msh
          this%mesh2mesh = .false.
       end if 

       if (associated(data%ulag)) then       
          ulag => data%ulag
          vlag => data%vlag
          wlag => data%wlag
          read_lag = .true.
       else
          read_lag = .false.
       end if

       if (associated(data%s)) then       
          s => data%s
          read_scalar = .true.
       else
          read_scalar = .false.
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
    call MPI_File_read_all(fh, optional_fields, 1, MPI_INTEGER, status, ierr)
    call MPI_File_read_all(fh, chkp%t, 1, MPI_DOUBLE_PRECISION, status, ierr)

    have_lag = mod(optional_fields,2)/1
    have_scalar = mod(optional_fields,4)/2

    if ( ( glb_nelv .ne. msh%glb_nelv ) .or. &
         ( gdim .ne. msh%gdim) .or. &
         ( (have_lag .eq. 0) .and. (read_lag) ) .or. &
        ( (have_scalar .eq. 0) .and. (read_scalar) ) ) then
       call neko_error('Checkpoint does not match case')
    end if
    nel = msh%nelv
    this%sim_Xh => u%Xh
    if (gdim .eq. 3) then
       call this%chkp_Xh%init(GLL, lx, lx, lx)
    else
       call this%chkp_Xh%init(GLL, lx, lx)
    end if
    if (this%mesh2mesh) then
       dof = dofmap_t(msh, this%chkp_Xh)
       allocate(x_coord(u%Xh%lx,u%Xh%ly,u%Xh%lz,u%msh%nelv))
       allocate(y_coord(u%Xh%lx,u%Xh%ly,u%Xh%lz,u%msh%nelv))
       allocate(z_coord(u%Xh%lx,u%Xh%ly,u%Xh%lz,u%msh%nelv))
       !> To ensure that each point is within an element
       !! Remedies issue with points on the boundary
       !! Technically gives each point a slightly different value
       !! but still within the specified tolerance
       do e = 1, u%dof%msh%nelv
          center_x = 0d0
          center_y = 0d0
          center_z = 0d0
          do i = 1,u%dof%Xh%lxyz
             center_x = center_x + u%dof%x(i,1,1,e)
             center_y = center_y + u%dof%y(i,1,1,e)
             center_z = center_z + u%dof%z(i,1,1,e)
          end do
          center_x = center_x/u%Xh%lxyz
          center_y = center_y/u%Xh%lxyz
          center_z = center_z/u%Xh%lxyz
          do i = 1,u%dof%Xh%lxyz
             x_coord(i,1,1,e) = u%dof%x(i,1,1,e) - tol*(u%dof%x(i,1,1,e)-center_x)
             y_coord(i,1,1,e) = u%dof%y(i,1,1,e) - tol*(u%dof%y(i,1,1,e)-center_y)
             z_coord(i,1,1,e) = u%dof%z(i,1,1,e) - tol*(u%dof%z(i,1,1,e)-center_z)
          end do
       end do
       call this%global_interp%init(dof,tol=tol)
       call this%global_interp%find_points(x_coord,y_coord,z_coord,u%dof%size())
       deallocate(x_coord)
       deallocate(y_coord)
       deallocate(z_coord)
    else
       call this%space_interp%init(this%sim_Xh, this%chkp_Xh) 
    end if
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

    if (read_scalar) then 
       byte_offset = mpi_offset + &
            dof_offset * int(MPI_REAL_PREC_SIZE, i8)
       call this%read_field(fh, byte_offset, s%x, nel)
       mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, i8)
    end if
    
    call MPI_File_close(fh, ierr)      

    call this%global_interp%free()
    call this%space_interp%free()
    
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
    if (this%mesh2mesh) then
       x = 0.0_rp
       call this%global_interp%evaluate(x,read_array)

    else
       call this%space_interp%map_host(x, read_array, nel, this%sim_Xh)
    end if
    deallocate(read_array)
  end subroutine chkp_read_field
  
end module chkp_file
