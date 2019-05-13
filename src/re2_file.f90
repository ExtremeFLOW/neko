!> NEKTON mesh data in re2 format
!! @details This module is used to read/write binary NEKTION mesh data
module re2_file
  use generic_file
  use num_types
  use utils
  use mesh
  use point
  implicit none
  private
  

  !> Interface for NEKTON re2 files
  type, public, extends(generic_file_t) :: re2_file_t
   contains
     procedure :: read => re2_file_read
     procedure :: write => re2_file_write
  end type re2_file_t

contains

  !> Load a binary NEKTON mesh from a re2 file
  subroutine re2_file_read(this, data)
    class(re2_file_t) :: this
    class(*), target, intent(inout) :: data
    type(mesh_t), pointer :: msh
    character(len=80) :: hdr_ver, hdr_str
    integer :: i, j, k, nel, ndim, nelv, ierr, el_idx, pt_idx
    real(kind=sp), allocatable :: xyz(:)
    real(kind=sp) :: test, x(8), y(8), z(8)
    type(point_t) :: p(8)

    select type(data)
    type is (mesh_t)
       msh => data
    end select

    open(unit=9,file=trim(this%fname), status='old', iostat=ierr)
    write(*, '(A,A)') " Reading binary NEKTON file ", this%fname
    read(9, '(a5,i9,i3,i9,a54)') hdr_ver, nel, ndim, nelv, hdr_str
    write(*,1) ndim, nelv
1   format(1x,'ndim = ', i1, ', nelements =', i7)
    close(9)


    call mesh_init(msh, ndim, nelv)

    allocate(xyz(nelv * 26))

    open(unit=9,file=trim(this%fname), status='old', access='stream', form='unformatted')
    read(9, pos=81) test
    read(9, pos=85) xyz

    pt_idx = 1
    el_idx = 1
    k = 2
    if (ndim .eq. 2) then
       do i = 1, nelv
          do j = 1, 8
             x(j) = xyz(k)
             k = k + 1
          end do

          do j = 1, 8
             y(j) = xyz(k)
             k = k + 1
          end do

          do j = 1, 8             
             p(j) = point_t(dble(x(j)), dble(y(j)), 0d0, pt_idx)
             pt_idx = pt_idx + 1
          end do

          call mesh_add_element(msh, el_idx, p(1), p(2), p(3), p(4))
          el_idx = el_idx + 1
       end do       
    else if (ndim .eq. 3) then
       do i = 1, nelv
          do j = 1, 8
             x(j) = xyz(k)
             k = k + 1
          end do

          do j = 1, 8
             y(j) = xyz(k)
             k = k + 1
          end do

          do j = 1, 8
             z(j) = xyz(k)
             k = k + 1
          end do
          k = k + 1

          do j = 1, 8             
             p(j) = point_t(dble(x(j)), dble(y(j)), dble(z(j)), pt_idx)
             pt_idx = pt_idx + 1
          end do

          call mesh_add_element(msh, el_idx, &
               p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8))          
          el_idx = el_idx + 1
       end do
    end if
    write(*,*) 'Done'

    !> @todo Add support for curved side data

    close(9)

    deallocate(xyz)
    
  end subroutine re2_file_read

  subroutine re2_file_write(this, data)
    class(re2_file_t), intent(in) :: this
    class(*), target, intent(in) :: data
  end subroutine re2_file_write

end module re2_file
