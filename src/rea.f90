module rea
  use generic_file
  use num_types
  use mesh
  implicit none
  private

  type, extends(generic_file_t) :: rea_t
   contains
     procedure :: read => rea_read
     procedure :: write => rea_write
  end type rea_t

  public :: rea_t

contains
  
  subroutine rea_read(this, data)
    class(rea_t) :: this
    class(*), intent(inout) :: data
    integer :: ndim, nparam, nskip, nlogic
    integer :: nelgs, nelgv, i, j, ierr

    open(unit=9,file=trim(this%fname), status='old', iostat=ierr)
    write(*, '(A,A)') ' Reading ', this%fname
    
    read(9, *)
    read(9, *)
    read(9, *) ndim
    read(9, *) nparam
    
    ! Skip parameters
    do i = 1, nparam
       read(9, *)
    end do
    
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
    write(*,*) nelgs, ndim, nelgv
    
    
    select type(data)
    type is (mesh_t)    
       call mesh_init_coordinates(data, ndim, nelgv)       
       do i = 1, nelgv
          read(9, *)
          if (ndim .eq. 2) then
             read(9, *) (data%xc(j, i),j=1,4)
             read(9, *) (data%yc(j, i),j=1,4)
          else if (ndim .eq. 3) then
             read(9, *) (data%xc(j, i),j=1,4)
             read(9, *) (data%yc(j, i),j=1,4)
             read(9, *) (data%zc(j, i),j=1,4)
             read(9, *) (data%xc(j, i),j=5,8)
             read(9, *) (data%yc(j, i),j=5,8)
             read(9, *) (data%zc(j, i),j=5,8)
          end if
       end do
    class default
       write(*,*) 'Fail!'
    end select
    
  end subroutine rea_read


  subroutine rea_write(this, data)
    class(rea_t) :: this
    class(*), intent(in) :: data
  end subroutine rea_write
end module rea
