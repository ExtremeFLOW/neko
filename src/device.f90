module device
  use num_types
  use utils
  !$ use omp_lib
  use, intrinsic :: iso_c_binding
  implicit none

contains

  !> Allocate a Fortran array on the device
  subroutine device_alloc(x, x_d, n, d) 
    integer, intent(in), value :: n
    real(kind=dp), intent(in), target :: x(n)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), optional :: d
    integer :: dev
    integer(c_size_t) :: nb    

!$  if (present(d)) then
!$     dev = d
!$  else
!$     dev = omp_get_default_device()
!$  end if
    x_d = C_NULL_PTR
!$  nb = n * 8
!$  x_d = omp_target_alloc(nb, dev)
!$  if (x_d .eq. C_NULL_PTR) then
!$     call neko_error('Device allocation failed')
!$  end if
   
  end subroutine device_alloc
    
  !> Associate a Fortran array to an array on the device
  subroutine device_associate(x, x_d, n, d)
    integer, intent(in), value :: n
    real(kind=dp), intent(in), target :: x(n)
    type(c_ptr), intent(in) :: x_d
    integer, intent(in), optional :: d
    integer :: dev
    integer(c_size_t) :: nb    

!$  if (present(d)) then
!$     dev = d
!$  else
!$     dev = omp_get_default_device()
!$  end if
    
!$  if (x_d .eq. C_NULL_PTR) then
!$     call neko_error('Array not allocated on device')
!$  end if
    
!$  nb = size(x) * 8
!$  if (omp_target_associate_ptr(c_loc(x), x_d, &
!$       nb, 0_c_size_t, dev) .gt. 0) then
!$     call neko_error('Device association failed')
!$  end if
              
  end subroutine device_associate
  
end module device
