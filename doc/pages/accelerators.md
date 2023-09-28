# Accelerators {#accelerators}

## Device abstraction layer

Neko has a device abstraction layer (device.F90) to manage device memory, data transfer and kernel invocations directly from a familiar Fortran interface, targeting all supported accelerator backends.

### Memory management
#### Allocation/deallocation
Allocating device memory can be done via the low-level device::device_alloc function, which takes a C pointer and requested size (in bytes) as arguments (see below). The pointer points to the allocated device memory if the allocation is successful.  
~~~~~~~~~~~~~~~{.f90}
  type(c_ptr) :: x_d
  integer(c_size_t) :: size
  ...
  call device_alloc(x_d, size)
~~~~~~~~~~~~~~~

To deallocate device memory, call device::device_free as shown below.
~~~~~~~~~~~~~~~{.f90}
  type(c_ptr) :: x_d
  ...
  call device_free(x_d)
~~~~~~~~~~~~~~~
#### Associate data on host and device
It is often helpful to associate a Fortran array with an array allocated on the device. This can be done via the device::device_associate routine, which takes a Fortran array `x` of type `integer`, `integer(i8)`, `real(kind=sp)` or `real(kind=dp)` (of maximum rank four) and a device pointer `x_d` as arguments.
~~~~~~~~~~~~~~~{.f90}
  integer, allocatable :: x(:)
  type(c_ptr) :: x_d
  ...
  call device_associate(x, x_d)
~~~~~~~~~~~~~~~
Once associated, the device pointer can be retrieved by calling the device::device_get_ptr function, which returns the pointer associated with an array `x`.
~~~~~~~~~~~~~~~{.f90}
  integer, allocatable :: x(:)
  type(c_ptr) :: x_d
  ...
  x_d = device_get_ptr(x)
~~~~~~~~~~~~~~~
@attention `device_get_ptr` returns a fatal error unless `x` has been associated to a device pointer. If in doubt, use the routine device::device_associated to check the status of an array.
~~~~~~~~~~~~~~~{.f90}
  integer, allocatable :: x(:)
  type(c_ptr) :: x_d
  ...
  if (device_associated(x)) then
     x_d = device_get_ptr(x)
  end if
~~~~~~~~~~~~~~~

#### Map a host array to a device
Since allocation and association is such an ordinary operation, Neko provides a combined routine device::device_map which takes a Fortran array `x` of type `integer`, `integer(i8)`, `real(kind=sp)` or `real(kind=dp)`, its size `n` (**Note** number of entries, not size in bytes as for `device_allocate`) and a C pointer `x_d` to the (to be allocated) device memory.
~~~~~~~~~~~~~~~{.f90}
  integer, allocatable :: x(:)
  type(c_ptr) :: x_d
  integer :: n
  ...
  allocate(x(n))
  ...
  call device_map(x, x_d, n)
~~~~~~~~~~~~~~~

#### Data transfer
To copy data between host and device (and device to use) use the routine device::device_memcpy which takes a Fortran array `x` of type `integer`, `integer(i8)`, `real(kind=sp)` or `real(kind=dp)`, its size `n` and the direction as the third argument which can either be `HOST_TO_DEVICE` to copy data to the device, or `DEVICE_TO_HOST` to retrieve data from the device.
~~~~~~~~~~~~~~~{.f90}
  integer, allocatable :: x(:)
  type(c_ptr) :: x_d
  integer :: n
  ...
  allocate(x(n))
  ...
  ! Fill the device with data
  call device_memcpy(x, x_d, n, HOST_TO_DEVICE)
  ...
  ! Retrieve data from device
  call device_memcpy(x, x_d, n, DEVICE_TO_HOST)
~~~~~~~~~~~~~~~

@attention device::device_memcpy defaults to asynchronous data transfer. The optional boolean argument `sync` must be true if synchronous transfers are needed.

@note that there's a special device::device_memcpy_cptr routine which only works with pointers, either on the host or the device, and the size needs to be given in bytes. This routine can be used to copy data between two arrays on a device with the direction `DEVICE_TO_DEVICE`.
~~~~~~~~~~~~~~~{.f90}
  type(c_ptr) :: x_d, y_d
  integer(c_size_t) :: s
  ...
  ! Copy the content of x into y on the device
  call device_memcpy_cptr(y_d, x_d, s, DEVICE_TO_DEVICE)
~~~~~~~~~~~~~~~

@attention It is the programmers' responsibility to make sure that device arrays are kept in sync with the associated host array. Neko does not perform any implicit data movement.

### Offload work
To offload work to a device, most routines in Neko have a device version prefixed with `device_<name of routine>`. These routines have the same arguments as the host equivalent, but one must pass device pointers instead of Fortran arrays.

For example, we call the `math::add2` routine to add two arrays together on the host.
~~~~~~~~~~~~~~~{.f90}
  integer, allocatable :: x(:), y(:)
  integer :: n
  ...
  allocate(x(n), y(n))
  ...
  call add2(x, y, n)
~~~~~~~~~~~~~~~
To offload the computation to the device, one must obtain the device pointers of `x` and `y`, and instead call device_math::device_add2
~~~~~~~~~~~~~~~{.f90}
  integer, allocatable :: x(:), y(:)
  type(c_ptr) :: x_d, y_d
  integer :: n
  ...
  allocate(x(n), y(n))
  ...
  x_d = device_get_ptr(x)
  y_d = device_get_ptr(y)
  call device_add2(x_d, y_d, n)
~~~~~~~~~~~~~~~
@note Most derived types in Neko already contain one or several device pointers associated with its internal data. Thus the `device_get_ptr` call can often be omitted.

However, for type bound procedures, such as computing the matrix-vector product derived from `ax_t`, one should **always** call the same type bound procedure (in this case `compute`) as on the host. This is because derived types contain all the logic such that the fastest backend is always selected or is instantiated as a backend-specific type during initialisation (see for example ax_helm_fctry.f90)

