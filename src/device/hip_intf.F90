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
!> Fortran HIP interface
module hip_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  
  !> Enum @a hipError_t
  enum, bind(c)
     enumerator :: hipSuccess = 0
     enumerator :: hipErrorInvalidContext = 1
     enumerator :: hipErrorInvalidKernelFile = 2
     enumerator :: hipErrorMemoryAllocation = 3
     enumerator :: hipErrorInitializationError = 4
     enumerator :: hipErrorLaunchFailure = 5
     enumerator :: hipErrorLaunchOutOfResources = 6
     enumerator :: hipErrorInvalidDevice = 7
     enumerator :: hipErrorInvalidValue = 8
     enumerator :: hipErrorInvalidDevicePointer = 9
     enumerator :: hipErrorInvalidMemcpyDirection = 10
     enumerator :: hipErrorUnknown = 11
     enumerator :: hipErrorInvalidResourceHandle = 12
     enumerator :: hipErrorNotReady = 13
     enumerator :: hipErrorNoDevice = 14
     enumerator :: hipErrorPeerAccessAlreadyEnabled = 15
     enumerator :: hipErrorPeerAccessNotEnabled = 16
     enumerator :: hipErrorRuntimeMemory = 17
     enumerator :: hipErrorRuntimeOther = 18
     enumerator :: hipErrorHostMemoryAlreadyRegistered = 19
     enumerator :: hipErrorHostMemoryNotRegistered = 20
     enumerator :: hipErrorMapBufferObjectFailed = 21
     enumerator :: hipErrorTbd = 22
  end enum
  
  !> Enum @a hipMemcpyKind
  enum, bind(c)
     enumerator :: hipMemcpyHostToHost = 0
     enumerator :: hipMemcpyHostToDevice = 1
     enumerator :: hipMemcpyDeviceToHost = 2
     enumerator :: hipMemcpyDevicetoDevice = 3
     enumerator :: hipMemcpyDefault = 4
  end enum
  
  interface
     integer (c_int) function hipMalloc(ptr_d, s) &
          bind(c, name='hipMalloc')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       integer(c_size_t), value :: s
     end function hipMalloc
  end interface

  interface
     integer (c_int) function hipFree(ptr_d) &
          bind(c, name='hipFree')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function hipFree
  end interface
  
  interface
     integer (c_int) function hipMemcpy(ptr_dst, ptr_src, s, dir) &
          bind(c, name='hipMemcpy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_dst, ptr_src
       integer(c_size_t), value :: s
       integer(c_int), value :: dir
     end function hipMemcpy
  end interface

  interface
     integer (c_int) function hipMemcpyAsync(ptr_dst, ptr_src, s, dir) &
          bind(c, name='hipMemcpyAsync')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_dst, ptr_src
       integer(c_size_t), value :: s
       integer(c_int), value :: dir
     end function hipMemcpyAsync
  end interface
  
  interface
     integer (c_int) function hipDeviceSynchronize() &
          bind(c, name='hipDeviceSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
     end function hipDeviceSynchronize
  end interface

  interface
     integer (c_int) function hipDeviceGetName(name, len, device) &
          bind(c, name='hipDeviceGetName')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: name
       integer(c_int), value :: len
       integer(c_int), value :: device
     end function hipDeviceGetName
  end interface

  interface
     integer (c_int) function hipStreamCreate(stream) &
          bind(c, name='hipStreamCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: stream
     end function hipStreamCreate
  end interface

  interface
     integer (c_int) function hipStreamCreateWithFlags(stream, flags) &
          bind(c, name='hipStreamCreateWithFlags')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: stream
       integer(c_int), value :: flags
     end function hipStreamCreateWithFlags
  end interface

  interface
     integer (c_int) function hipStreamDestroy(steam) &
          bind(c, name='hipStreamDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: steam
     end function hipStreamDestroy
  end interface

  interface 
     integer (c_int) function hipStreamSynchronize(stream) &
          bind(c, name='hipStreamSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream
     end function hipStreamSynchronize
  end interface
  
  interface 
     integer (c_int) function hipStreamWaitEvent(stream, event, flags) &
          bind(c, name='hipStreamWaitEvent')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: stream, event
       integer(c_int), value :: flags
     end function hipStreamWaitEvent
  end interface

  interface
     integer (c_int) function hipEventCreate(event) &
          bind(c, name='hipEventCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: event
     end function hipEventCreate
  end interface

  interface
     integer (c_int) function hipEventDestroy(event) &
          bind(c, name='hipEventDestroy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function hipEventDestroy
  end interface

  interface
     integer (c_int) function hipEventCreateWithFlags(event, flags) &
          bind(c, name='hipEventCreateWithFlags')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: event
       integer(c_int), value :: flags
     end function hipEventCreateWithFlags
  end interface

  interface
     integer (c_int) function hipEventRecord(event, stream) &
          bind(c, name='hipEventRecord')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event, stream
     end function hipEventRecord
  end interface

  interface 
     integer (c_int) function hipEventSynchronize(event) &
          bind(c, name='hipEventSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function hipEventSynchronize
  end interface
  
contains

  subroutine hip_device_name(name)
    character(len=*), intent(inout) :: name
    character(kind=c_char, len=1024), target :: c_name
    integer :: end_pos
    
    if (hipDeviceGetName(c_loc(c_name), 1024, 0) .ne. hipSuccess) then
       call neko_error('Failed to query device')
    end if

    end_pos = scan(c_name, C_NULL_CHAR)
    if(end_pos .ge. 2) then
       name(1:end_pos-1) = c_name(1:end_pos-1)
    endif

  end subroutine hip_device_name

#endif
 
end module hip_intf
