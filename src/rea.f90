!
!> NEKTON session data
!! @details This module is used to represent NEKTON session data
module rea
  use generic_file
  use num_types
  use utils
  use mesh
  use point
  use re2_file
  implicit none
  private

  !> Interface for NEKTON ascii files
  type, public, extends(generic_file_t) :: rea_file_t
   contains
     procedure :: read => rea_file_read
     procedure :: write => rea_file_write
  end type rea_file_t

  !> NEKTON session data struct.
  !! @todo add missing data fields
  type, public :: rea_t
     type(mesh_t) :: msh                     !< Mesh (rep. as a Neko mesh)
     real(kind=dp), allocatable :: params(:) !< Parameters
  end type rea_t
  
contains
  
  !> Free a NEKTON session data
  subroutine rea_free(r)
    type(rea_t), intent(inout) :: r

    call mesh_free(r%msh)
    
    if (allocated(r%params)) then
       deallocate(r%params)
    end if

  end subroutine rea_free

  !> Load NEKTON session data from an ascii file
  subroutine rea_file_read(this, data)
    class(rea_file_t) :: this
    class(*), target, intent(inout) :: data
    type(mesh_t), pointer :: msh
    real(kind=dp), pointer :: params(:)
    integer :: ndim, nparam, nskip, nlogic
    integer :: nelgs, nelgv, i, j, ierr
    integer :: el_idx, pt_idx
    logical :: read_param
    real(kind=dp) :: xc(8), yc(8), zc(8)
    type(point_t) :: p(8)
    type(re2_file_t) :: re2_file
    character(len=80) :: re2_fname

    select type(data)
    type is (rea_t)
       call rea_free(data)       
       msh => data%msh
       params => data%params
       read_param = .true.
    type is (mesh_t)    
       msh => data
       read_param = .false.
    class default
       call neko_error('Invalid output data')
    end select

    open(unit=9,file=trim(this%fname), status='old', iostat=ierr)
    write(*, '(A,A)') " Reading NEKTON file ", this%fname
    
    read(9, *)
    read(9, *)
    read(9, *) ndim
    read(9, *) nparam
    
    if (.not. read_param) then
       ! Skip parameters
       do i = 1, nparam
          read(9, *)
       end do
    else       
       allocate(params(nparam))
       do i = 1, nparam
          read(9, *) params(i)
       end do
    end if
    
    ! Skip passive scalars
    read(9, *) nskip
    do i = 1, nskip
       read(9, *)
    end do
    
    ! Skip logic switches
    read(9, *) nlogic
    do i = 1, nlogic
       read(9, *)
    end do
    
    ! Read mesh info
    read(9, *)
    read(9, *)
    read(9, *) nelgs,ndim, nelgv
    if (nelgs .lt. 0) then
       re2_fname = trim(this%fname(1:scan(trim(this%fname), '.', back=.true.)))//'re2' 
       call re2_file%init(re2_fname)
       call re2_file%read(msh)
    else       
       write(*,1) ndim, nelgv
1      format(1x,'ndim = ', i1, ', nelements =', i7)
       call mesh_init(msh, ndim, nelgv)
       
       el_idx = 1
       pt_idx = 1
       do i = 1, nelgv
          read(9, *)
          if (ndim .eq. 2) then
             read(9, *) (xc(j),j=1,4)
             read(9, *) (yc(j),j=1,4)
             do j = 1, 4
                p(j) = point_t(xc(j), yc(j), 0d0, pt_idx)
                pt_idx = pt_idx
             end do
             call mesh_add_element(msh, el_idx, p(1), p(2), p(3), p(4))
          else if (ndim .eq. 3) then
             read(9, *) (xc(j),j=1,4)
             read(9, *) (yc(j),j=1,4)
             read(9, *) (zc(j),j=1,4)
             read(9, *) (xc(j),j=5,8)
             read(9, *) (yc(j),j=5,8)
             read(9, *) (zc(j),j=5,8)
             do j = 1, 8
                p(j) = point_t(xc(j), yc(j), zc(j), pt_idx)
                pt_idx = pt_idx + 1
             end do
             call mesh_add_element(msh, el_idx, &
                  p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))
          end if
          el_idx = el_idx + 1
       end do
       write(*,*) 'Done'
       
       
       !> @todo Add support for curved side data
       
       close(9)
    endif
    
  end subroutine rea_file_read

  subroutine rea_file_write(this, data)
    class(rea_file_t), intent(in) :: this
    class(*), target, intent(in) :: data
  end subroutine rea_file_write

end module rea
