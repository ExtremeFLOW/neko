!> Neko checkpoint file format
!! @details this module defines interface to read/write Neko's ceckpoint files
module chkp_file
  use generic_file
  use checkpoint    
  use num_types
  use field
  use utils
  use mesh
  use mpi_types
  use mpi_f08
  implicit none
  private

  !> Interface for Neko checkpoint files
  type, public, extends(generic_file_t) :: chkp_file_t
   contains
     procedure :: read => chkp_file_read
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
    character(len=80) :: fname
    integer :: ierr, suffix_pos, have_lag
    type(field_t), pointer :: u, v, w, p
    real(kind=rp), pointer :: ulag(:,:,:,:,:) => null()
    real(kind=rp), pointer :: vlag(:,:,:,:,:) => null()
    real(kind=rp), pointer :: wlag(:,:,:,:,:) => null()
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer(kind=8) :: n_glb_dofs, dof_offset
    integer(kind=8) :: n_lag, n_glb_lag, lag_offset
    logical write_lag

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


    dof_offset = int(msh%offset_el, 8) * int(u%Xh%lx * u%Xh%ly * u%Xh%lz, 8)
    n_glb_dofs = int(u%Xh%lx * u%Xh%ly * u%Xh%lz, 8) * int(msh%glb_nelv, 8)
    
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
         dof_offset * int(MPI_REAL_PREC_SIZE, 8)
    call MPI_File_write_at_all(fh, byte_offset, u%x, u%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = 4 * MPI_INTEGER_SIZE + MPI_DOUBLE_PRECISION_SIZE + &
         n_glb_dofs * int(MPI_REAL_PREC_SIZE, 8)
    
    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, 8)
    call MPI_File_write_at_all(fh, byte_offset, v%x, v%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, 8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, 8)
    call MPI_File_write_at_all(fh, byte_offset, w%x, w%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, 8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, 8)
    call MPI_File_write_at_all(fh, byte_offset, p%x, p%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, 8)

    !
    ! Dump optional payload
    !

    if (write_lag) then
       
       n_lag = int(size(ulag) / u%dof%size(), 8)
       n_glb_lag = int(n_lag * (u%Xh%lx * u%Xh%ly * u%Xh%lz), 8) * &
            int(msh%glb_nelv, 8)
       lag_offset = int(msh%offset_el, 8) * &
            int(n_lag * (u%Xh%lx * u%Xh%ly * u%Xh%lz), 8)

       byte_offset = mpi_offset + &            
                lag_offset * int(MPI_REAL_PREC_SIZE, 8)
       call MPI_File_write_at_all(fh, byte_offset, ulag, size(ulag), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_lag * int(MPI_REAL_PREC_SIZE, 8)

       byte_offset = mpi_offset + &            
                lag_offset * int(MPI_REAL_PREC_SIZE, 8)
       call MPI_File_write_at_all(fh, byte_offset, vlag, size(vlag), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_lag * int(MPI_REAL_PREC_SIZE, 8)

       byte_offset = mpi_offset + &            
            lag_offset * int(MPI_REAL_PREC_SIZE, 8)
       call MPI_File_write_at_all(fh, byte_offset, wlag, size(wlag), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_lag * int(MPI_REAL_PREC_SIZE, 8)
       
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
    character(len=80) :: fname
    integer :: ierr, suffix_pos
    type(field_t), pointer :: u, v, w, p
    real(kind=rp), pointer :: ulag(:,:,:,:,:) => null()
    real(kind=rp), pointer :: vlag(:,:,:,:,:) => null()
    real(kind=rp), pointer :: wlag(:,:,:,:,:) => null()
    type(mesh_t), pointer :: msh
    type(MPI_Status) :: status
    type(MPI_File) :: fh
    integer (kind=MPI_OFFSET_KIND) :: mpi_offset, byte_offset
    integer(kind=8) :: n_glb_dofs, dof_offset
    integer(kind=8) :: n_lag, n_glb_lag, lag_offset
    integer :: glb_nelv, gdim, lx, have_lag
    logical read_lag
    
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
         ( lx .ne. u%Xh%lx .or. lx .ne. u%Xh%ly .or. lx .ne. u%Xh%lz) .or. &
         ( (have_lag .eq. 1) .and. (.not. read_lag) ) ) then
       call neko_error('Checkpoint does not match case')
    end if
    

    dof_offset = int(msh%offset_el, 8) * int(u%Xh%lx * u%Xh%ly * u%Xh%lz, 8)
    n_glb_dofs = int(u%Xh%lx * u%Xh%ly * u%Xh%lz, 8) * int(msh%glb_nelv, 8)    
    
    !
    ! Read mandatory checkpoint data
    !
    
    byte_offset = 4 * MPI_INTEGER_SIZE + MPI_DOUBLE_PRECISION_SIZE + &
         dof_offset * int(MPI_REAL_PREC_SIZE, 8)
    call MPI_File_read_at_all(fh, byte_offset, u%x, u%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = 4 * MPI_INTEGER_SIZE + MPI_DOUBLE_PRECISION_SIZE + &
         n_glb_dofs * int(MPI_REAL_PREC_SIZE, 8)
    
    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, 8)
    call MPI_File_read_at_all(fh, byte_offset, v%x, v%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, 8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, 8)
    call MPI_File_read_at_all(fh, byte_offset, w%x, w%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, 8)

    byte_offset = mpi_offset + &
         dof_offset * int(MPI_REAL_PREC_SIZE, 8)
    call MPI_File_read_at_all(fh, byte_offset, p%x, p%dof%size(), &
         MPI_REAL_PRECISION, status, ierr)
    mpi_offset = mpi_offset + n_glb_dofs * int(MPI_REAL_PREC_SIZE, 8)

    !
    ! Read optional payload
    !

    if (read_lag) then
       
       n_lag = int(size(ulag) / u%dof%size(), 8)
       n_glb_lag = int(n_lag * (u%Xh%lx * u%Xh%ly * u%Xh%lz), 8) * &
            int(msh%glb_nelv, 8)
       lag_offset = int(msh%offset_el, 8) * &
            int(n_lag * (u%Xh%lx * u%Xh%ly * u%Xh%lz), 8)

       byte_offset = mpi_offset + &            
                lag_offset * int(MPI_REAL_PREC_SIZE, 8)
       call MPI_File_read_at_all(fh, byte_offset, ulag, size(ulag), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_lag * int(MPI_REAL_PREC_SIZE, 8)

       byte_offset = mpi_offset + &            
                lag_offset * int(MPI_REAL_PREC_SIZE, 8)
       call MPI_File_read_at_all(fh, byte_offset, vlag, size(vlag), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_lag * int(MPI_REAL_PREC_SIZE, 8)

       byte_offset = mpi_offset + &            
            lag_offset * int(MPI_REAL_PREC_SIZE, 8)
       call MPI_File_read_at_all(fh, byte_offset, wlag, size(wlag), &
            MPI_REAL_PRECISION, status, ierr)
       mpi_offset = mpi_offset + n_glb_lag * int(MPI_REAL_PREC_SIZE, 8)
       
    end if
    
    call MPI_File_close(fh, ierr)      
    
  end subroutine chkp_file_read
  
  
end module chkp_file
