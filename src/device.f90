!> Simple interface to device memory
module device
  use num_types
  use utils
  !$ use omp_lib
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: HOST_TO_DEVICE = 1, DEVICE_TO_HOST = 2
  
contains

  !> Allocate a Fortran array on the device
  subroutine device_alloc(x, x_d, n, d) 
    integer, intent(in), value :: n
    real(kind=dp), intent(in), target :: x(n)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), optional :: d
    integer(c_int) :: dev
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
    integer(c_int) :: dev
    integer(c_size_t) :: nb    

!$  if (present(d)) then
!$     dev = d
!$  else
!$     dev = omp_get_default_device()
!$  end if
    
!$  if (x_d .eq. C_NULL_PTR) then
!$     call neko_error('Array not allocated on device')
!$  end if
    
!$  nb = n * 8
!$  if (omp_target_associate_ptr(c_loc(x), x_d, &
!$       nb, 0_c_size_t, dev) .gt. 0) then
!$     call neko_error('Device association failed')
!$  end if
    
  end subroutine device_associate

  subroutine device_memcpy(x, x_d, n, dir, d)
    integer, intent(in), value :: n
    real(kind=dp), intent(in), target :: x(n)
    type(c_ptr), intent(in) :: x_d
    integer, intent(in) :: dir
    integer, intent(in), optional :: d
    integer(c_int) :: dev, hst
    integer(c_size_t) :: nb

!$  if (present(d)) then
!$     dev = d
!$  else
!$     dev = omp_get_default_device()
!$  end if

!$  hst = omp_get_initial_device()

!$  if ((x_d .eq. C_NULL_PTR) .or. &
!$       (omp_target_is_present(c_loc(x), dev) .eq. C_NULL_PTR)) then
!$     call neko_error('Array not allocated on device')
!$  end if

!$  nb = n * 8
!$  if (dir .eq. HOST_TO_DEVICE) then
!$     if (omp_target_memcpy(x_d, c_loc(x), nb, &
!$          0_c_size_t, 0_c_size_t, dev, hst) .gt. 0) then
!$        call neko_error('Device memcpy (host-to-device) failed')
!$     end if
!$  else if (dir .eq. DEVICE_TO_HOST) then
!$     if (omp_target_memcpy(c_loc(x), x_d, nb, &
!$          0_c_size_t, 0_c_size_t, hst, dev) .gt. 0) then
!$        call neko_error('Device memcpy (device-to-host) failed')
!$     end if
!$  else
!$     call neko_error('Invalid direction')
!$  end if
    
  end subroutine device_memcpy

end module device
