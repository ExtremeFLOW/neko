!> Fortran CUDA interface
module cuda_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_CUDA

  !> Enum @a cudaError
  enum, bind(c)
     enumerator :: cudaSuccess = 0
     enumerator :: cudaErrorInvalidValue = 1
     enumerator :: cudaErrorMemoryAllocation = 2
     enumerator :: cudaErrorInitializationError = 3
  end enum

  !> Enum @a cudaMemcpyKind
  enum, bind(c)
     enumerator :: cudaMemcpyHostToHost = 0
     enumerator :: cudaMemcpyHostToDevice = 1
     enumerator :: cudaMemcpyDeviceToHost = 2
     enumerator :: cudaMemcpyDevicetoDevice = 3
     enumerator :: cudaMemcpyDefault = 4
  end enum

  interface
     integer (c_int) function cudaMalloc(ptr_d, s) &
          bind(c, name='cudaMalloc')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       integer(c_size_t), value :: s
     end function cudaMalloc
  end interface

  interface
     integer (c_int) function cudaFree(ptr_d) &
          bind(c, name='cudaFree')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function cudaFree
  end interface
  
  interface
     integer (c_int) function cudaMemcpy(ptr_dst, ptr_src, s, dir) &
          bind(c, name='cudaMemcpy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_dst, ptr_src
       integer(c_size_t), value :: s
       integer(c_int), value :: dir
     end function cudaMemcpy
  end interface

  interface
     integer (c_int) function cudaMemcpyAsync(ptr_dst, ptr_src, s, dir) &
          bind(c, name='cudaMemcpyAsync')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_dst, ptr_src
       integer(c_size_t), value :: s
       integer(c_int), value :: dir
     end function cudaMemcpyAsync
  end interface
  
  interface
     integer (c_int) function cudaDeviceSynchronize() &
          bind(c, name='cudaDeviceSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
     end function cudaDeviceSynchronize
  end interface

  interface
     integer (c_int) function cudaGetDeviceProperties(prop, device) &
          bind(c, name='cudaGetDeviceProperties')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: prop
       integer(c_int), value :: device
     end function cudaGetDeviceProperties
  end interface

contains

  subroutine cuda_device_name(name)
    character(len=*), intent(inout) :: name
    character(kind=c_char, len=8192), target :: prop
    integer :: end_pos

    !
    ! Yes this is an ugly hack!
    ! Since we're only interested in the device name (first 256 bytes)
    ! we pass down a large enough chunk of memory to the cuda runtime
    ! and extract what we need later on
    !
    ! This will of course break if sizeof(cudaDeviceProp) > 8192
    !
    
    if (cudaGetDeviceProperties(c_loc(prop), 0) .ne. cudaSuccess) then
       call neko_error('Failed to query device')
    end if
    
    end_pos = scan(prop(1:256), C_NULL_CHAR)
    name(1:end_pos) = prop(1:end_pos)
    
  end subroutine cuda_device_name
  
#endif
 
end module cuda_intf
