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

  interface
     integer (c_int) function cudaStreamCreate(stream) &
          bind(c, name='cudaStreamCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: stream
     end function cudaStreamCreate
  end interface

  interface
     integer (c_int) function cudaStreamCreateWithFlags(stream, flags) &
          bind(c, name='cudaStreamCreateWithFlags')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: stream
       integer(c_int), value :: flags
     end function cudaStreamCreateWithFlags
  end interface

  interface
     integer (c_int) function cudaStreamDestroy(steam) &
          bind(c, name='cudaStreamDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: steam
     end function cudaStreamDestroy
  end interface

  interface 
     integer (c_int) function cudaStreamSynchronize(stream) &
          bind(c, name='cudaStreamSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream
     end function cudaStreamSynchronize
  end interface

  interface 
     integer (c_int) function cudaStreamWaitEvent(stream, event, flags) &
          bind(c, name='cudaStreamWaitEvent')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream, event
       integer(c_int), value :: flags
     end function cudaStreamWaitEvent
  end interface
  
  interface
     integer (c_int) function cudaProfilerStart() &
          bind(c, name='cudaProfilerStart')
       use, intrinsic :: iso_c_binding
       implicit none
     end function cudaProfilerStart
  end interface

  interface
     integer (c_int) function cudaProfilerStop() &
          bind(c, name='cudaProfilerStop')
       use, intrinsic :: iso_c_binding
       implicit none
     end function cudaProfilerStop
  end interface

  interface
     integer (c_int) function cudaEventCreate(event) &
          bind(c, name='cudaEventCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: event
     end function cudaEventCreate
  end interface

  interface
     integer (c_int) function cudaEventDestroy(event) &
          bind(c, name='cudaEventDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function cudaEventDestroy
  end interface
  
  interface
     integer (c_int) function cudaEventCreateWithFlags(event, flags) &
          bind(c, name='cudaEventCreateWithFlags')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: event
       integer(c_int), value :: flags
     end function cudaEventCreateWithFlags
  end interface

  interface
     integer (c_int) function cudaEventRecord(event, stream) &
          bind(c, name='cudaEventRecord')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event, stream
     end function cudaEventRecord
  end interface

  interface 
     integer (c_int) function cudaEventSynchronize(event) &
          bind(c, name='cudaEventSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function cudaEventSynchronize
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
    if(end_pos .ge. 2) then
       name(1:end_pos-1) = prop(1:end_pos-1)
    endif
  end subroutine cuda_device_name
  
#endif
 
end module cuda_intf
